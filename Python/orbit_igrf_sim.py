"""
Orbit + IGRF simulator for attitude filter testing.

Propagates a circular LEO orbit analytically and returns the IGRF magnetic
field unit vector in the ECI frame at each time step.  This vector serves as
the reference (known ECI direction) for the UKF measurement model.

Typical usage
-------------
    from Python.orbit_igrf_sim import OrbitIGRFSimulator
    from datetime import datetime

    sim = OrbitIGRFSimulator(epoch=datetime(2024, 6, 1))
    for _ in range(N_steps):
        b_eci, pos_eci = sim.step(dt)   # b_eci is the unit reference vector
"""

import math
from datetime import datetime, timedelta

import numpy as np
import ppigrf
import pymap3d
import pyorb.kepler as kepler_lib
from pyquaternion import Quaternion

# ── constants ──────────────────────────────────────────────────────────────────
MU_EARTH = 3.986004418e14   # m³ s⁻²


# ── IGRF helper ────────────────────────────────────────────────────────────────

def igrf_eci_vector(x_m: float, y_m: float, z_m: float,
                    time: datetime) -> np.ndarray:
    """Return the IGRF magnetic field vector in ECI at a given position.

    Parameters
    ----------
    x_m, y_m, z_m : ECI position in metres (J2000).
    time : UTC datetime for the IGRF coefficients and Earth rotation angle.

    Returns
    -------
    (3,) vector in the ECI frame [uT].  Not normalised.
    """
    # 1. ECI position → geodetic (lat/lon in degrees, alt in metres)
    lat, lon, alt_m = pymap3d.eci2geodetic(x_m, y_m, z_m, time)

    # 2. IGRF in local ENU frame [nT]  (ppigrf returns shape-(1,) arrays)
    b_e, b_n, b_u = ppigrf.igrf(float(lon), float(lat), float(alt_m) / 1e3, time)
    b_e, b_n, b_u = float(b_e), float(b_n), float(b_u)

    # 3. ENU *vector* → ECEF vector
    #    pymap3d.enu2ecef gives the ECEF *position* of the ENU point, so we
    #    subtract the observer's ECEF origin to isolate the vector part.
    x0, y0, z0 = pymap3d.geodetic2ecef(float(lat), float(lon), float(alt_m))
    bx_ecef, by_ecef, bz_ecef = pymap3d.enu2ecef(
        b_e, b_n, b_u, float(lat), float(lon), float(alt_m)
    )
    b_ecef = np.array([bx_ecef - x0, by_ecef - y0, bz_ecef - z0])

    # 4. ECEF vector → ECI vector (pure rotation — valid for vectors)
    bx_eci, by_eci, bz_eci = pymap3d.ecef2eci(b_ecef[0], b_ecef[1], b_ecef[2], time)
    b_eci = np.array([bx_eci, by_eci, bz_eci])

    norm = np.linalg.norm(b_eci)
    if norm < 1e-30:
        raise ValueError("IGRF returned a near-zero magnetic field vector.")
    return b_eci * 1e-3


# ── orbit simulator ────────────────────────────────────────────────────────────

class OrbitIGRFSimulator:
    """Circular-orbit propagator that yields IGRF reference vectors each step.

    The true anomaly is advanced analytically (ν += n·Δt), which is exact for
    a circular orbit (e = 0) and a good approximation for small eccentricities.

    Parameters
    ----------
    epoch : datetime
        UTC start time of the simulation.
    a : float
        Semi-major axis [m].  Defaults to a 400 km altitude LEO orbit.
    e : float
        Eccentricity (keep small for the analytic propagator).
    inc_deg : float
        Inclination [degrees].
    omega_deg : float
        Argument of perigee [degrees].
    raan_deg : float
        Right ascension of the ascending node [degrees].
    nu0_deg : float
        Initial true anomaly [degrees].
    """

    def __init__(
        self,
        epoch: datetime = datetime(2024, 6, 1, 0, 0, 0),
        a: float = 6_771_000.0,   # 400 km altitude
        e: float = 0.0,
        inc_deg: float = 51.6,    # ISS-like inclination
        omega_deg: float = 0.0,
        raan_deg: float = 20.0,
        nu0_deg: float = 0.0,
    ):
        self.epoch = epoch
        self.a = a
        self.e = e
        self.inc = math.radians(inc_deg)
        self.omega = math.radians(omega_deg)
        self.raan = math.radians(raan_deg)
        self.nu = math.radians(nu0_deg)
        self.n = math.sqrt(MU_EARTH / a ** 3)   # mean motion [rad/s]
        self._elapsed = 0.0                      # total elapsed time [s]

    @property
    def current_time(self) -> datetime:
        return self.epoch + timedelta(seconds=self._elapsed)

    def step(self, dt: float):
        """Advance the orbit by *dt* seconds.

        Returns
        -------
        b_eci : (3,) ndarray
            IGRF magnetic field vector in the ECI frame [uT].  Not normalised.
        pos_eci : (3,) ndarray
            Satellite ECI position [m] at the new position.
        vel_eci : (3,) ndarray
            Satellite ECI velocity [m/s] at the new position.
        """
        self._elapsed += dt

        # Advance true anomaly (circular-orbit analytic step)
        self.nu = (self.nu + self.n * dt) % (2 * math.pi)

        # Keplerian elements in pyorb convention: [a, e, i, ω, Ω, ν] (radians)
        kep = np.array([self.a, self.e, self.inc, self.omega, self.raan, self.nu])
        cart = kepler_lib.kep_to_cart(kep=kep, mu=MU_EARTH)   # [x,y,z,vx,vy,vz]
        pos_eci = cart[:3]
        vel_eci = cart[3:6]

        b_eci = igrf_eci_vector(pos_eci[0], pos_eci[1], pos_eci[2],
                                self.current_time)
        return b_eci, pos_eci, vel_eci

def _rotvec_to_quat(vec: np.ndarray) -> Quaternion:
    """Convert a rotation vector to a unit quaternion."""
    angle = np.linalg.norm(vec)
    if angle < 1e-12:
        return Quaternion()
    return Quaternion(axis=vec / angle, angle=angle)


def nadir_quaternion(pos_eci: np.ndarray, vel_eci: np.ndarray,
                     q_guess: Quaternion = None) -> Quaternion:
    """Return the body-to-ECI quaternion for perfect nadir-pointing attitude.

    Body frame convention (matches the MATLAB reference):
      - body z-axis → nadir (toward Earth centre)
      - body y-axis → flipped orbit normal (−h̃ = −(r × v)/|r × v|)
      - body x-axis → completes the right-hand system (≈ along-track)

    Parameters
    ----------
    pos_eci, vel_eci : ECI position [m] and velocity [m/s].
    q_guess : optional warm-start quaternion from the previous step.
        When supplied, a 2-vector Wahba fixed-point iteration is used so
        the output quaternion stays sign-continuous with the previous step.
        When omitted (cold start), the result is built directly from a DCM.
    """
    z_target = -pos_eci / np.linalg.norm(pos_eci)          # nadir  → body z
    h = np.cross(pos_eci, vel_eci)
    y_target = -h / np.linalg.norm(h)                       # −orbit normal → body y

    if q_guess is None:
        # Cold start: construct directly from the two target axes
        x_hat = np.cross(y_target, z_target)
        x_hat = x_hat / np.linalg.norm(x_hat)
        R = np.column_stack([x_hat, y_target, z_target])    # R_body_to_ECI
        return Quaternion(matrix=R)

    # Warm-start: iterative 2-vector Wahba alignment
    # err = 0.5 * (z_body_in_ECI × z_target + y_body_in_ECI × y_target)
    # Correction dq = rotvec(err) is left-multiplied onto the current estimate.
    z_body = np.array([0.0, 0.0, 1.0])
    y_body = np.array([0.0, 1.0, 0.0])
    q = q_guess
    for _ in range(100):
        z_now = q.rotate(z_body)
        y_now = q.rotate(y_body)
        err = 0.5 * (np.cross(z_now, z_target) + np.cross(y_now, y_target))
        if np.linalg.norm(err) < 1e-10:
            break
        q = _rotvec_to_quat(err) * q
    return q.normalised

