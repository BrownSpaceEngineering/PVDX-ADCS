import numpy as np
import pyquaternion
from pyquaternion import Quaternion
from numpy.linalg import cholesky
import math
from typing import List
from datetime import datetime

from orbit_igrf_sim import OrbitIGRFSimulator, nadir_quaternion

'''
ALL CREDITS TO https://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf FOR THE MUKF EQUATIONS

Took inspiration from Kail P Laughlin's thesis on "Single-Vector Aiding of an IMU for CubeSat Attitude
Determination" for some of the methods

This file estimates a Nadir-pointing satellite across both sun and eclipse, achieving good performance of sub-10 norm degree accuracies in both circumstances (with sigma_gyro = 0.001 and sigma_mag = 0.05 and sigma_unit_vector = 0.01)

'''

# ── simulation truth parameters ────────────────────────────────────────────────
sigma_gyro          = 0.0005   # gyro white-noise std dev [rad/s]
sigma_magnetometer  = 0.05    # magnetometer vector-sensor noise std dev in uT
sigma_unit_vector   = 0.01    # second reference vector sensor noise std dev (unit vector)
true_gyro_bias   = np.array([ 0.002, -0.001,  0.0015])  # constant gyro bias [rad/s]
dt = 0.1 

# 6-state error vector: [dq_x, dq_y, dq_z, b_x, b_y, b_z]
global_Q = np.zeros((6, 6))
global_Q[0:3, 0:3] = np.eye(3) * (sigma_gyro * dt)**2    # Attitude noise
global_Q[3:6, 3:6] = np.eye(3) * 1e-10  # Bias random walk

update_every = int(10 / dt)  # perform a measurement update every N propagation steps → 0.1 Hz
switch_every = int(60 * 45 / dt) # switch between 2 vector and 1 vector 
vector_input = 1 # 1 or 2 vector. the code will start on this mode until switch_every
simulation_time = int(60 * 180 / dt) # seconds of simulation. ~90 minutes per orbit

gyro_measurement = np.array([0.0, 0.0, 0.0])



def rotate(vec, q):
    v = Quaternion(scalar = 0, vector = vec)
    return (q * v * q.inverse).vector

def quaternion_to_rotation(x : Quaternion):
    quat = x.elements
    if quat[0] < 0:
        quat = -quat
    if(quat[0] >= 1):
        return np.zeros(3)
    theta = 2 * math.acos(quat[0])
    if(theta == 0):
        return np.zeros(3)
    vec = theta * quat[1:] / np.sqrt(1-quat[0]**2)
    return vec

def rotation_to_quat(vec):
    if(np.linalg.norm(vec) == 0):
        return Quaternion()
    axis = vec / np.linalg.norm(vec)
    angle = np.linalg.norm(vec)
    return Quaternion(angle = angle, axis = axis)

def get_error_vectors(Y: List[Quaternion], x : Quaternion):
    error_vectors = []
    for quat in Y:
        error_quaternion = quat * x.inverse
        error_vectors.append(quaternion_to_rotation(error_quaternion))
    return error_vectors

def gradient_descent(Y : List[Quaternion], x : Quaternion, weights=None):
    n_pts = len(Y)
    if weights is None:
        weights = np.ones(n_pts) / n_pts
    # Clamp negative weights and renormalise (safe for negative lambda values)
    safe_w = np.maximum(weights, 0)
    safe_w = safe_w / safe_w.sum()
    for i in range(100):
        error_vectors = get_error_vectors(Y, x)
        ave = np.average(error_vectors, weights=safe_w, axis=0)
        if(np.linalg.norm(ave) < 1e-6):
            break
        average_error_quat = rotation_to_quat(ave)
        x = average_error_quat * x
    return x, error_vectors

def get_weights(lam, x, alpha, beta):
    n = x.shape[0]
    c = .5 / (n + lam)
    covariance_weights = np.full(2*n + 1, c)
    mean_weights = np.full(2*n + 1, c)
    covariance_weights[0] = lam / (n + lam) + (1 - alpha**2 + beta)
    mean_weights[0] = lam/ (n + lam)
    return covariance_weights, mean_weights

def get_sigma_points(lam, x, P):
    n = x.shape[0]
    P = (P + P.T) / 2
    P += np.eye(n) * 1e-10
    U = cholesky((lam + n)*P).T
    sigmas = [x]
    for k in range(1, n + 1, 1):
        idx = k - 1
        sigmas.append(np.add(x, U[idx]))
    for k in range(n + 1, 2*n + 1, 1):
        idx = k - (n + 1)
        sigmas.append(np.subtract(x, U[idx]))
    return np.array(sigmas)

def calculate_lambda(alpha, x):
    n = x.shape[0]
    kappa = 9 - n
    return alpha**2 * (n + kappa) - n

def error_sigmas_to_quat_sigmas(error_sigmas, rotation : Quaternion):
    quat_sigmas = []
    for sigma in error_sigmas:
        quat_rot = rotation_to_quat(sigma[:3])
        new_quat_rot = quat_rot * rotation
        new_quat_rot = new_quat_rot.normalised
        quat_sigmas.append(np.concatenate([new_quat_rot.elements, sigma[3:6]]))
    return np.array(quat_sigmas)

def propagate_quat_sigmas(quat_sigmas):
    global gyro_measurement
    new_sigmas = []
    for sigma in quat_sigmas:
        quaternion = sigma[:4]
        bias = sigma[4:7]

        w = gyro_measurement - bias

        new_quaternion = Quaternion(quaternion)*rotation_to_quat(w * dt).inverse
        new_quaternion = new_quaternion.normalised

        new_sigmas.append(np.concatenate([new_quaternion.elements, sigma[4:7]]))
    return np.array(new_sigmas)

def get_measurements(quat_sigmas, ref_vecs):
    """Compute expected body-frame measurements for each sigma point.

    Parameters
    ----------
    ref_vecs : (3,) ndarray  OR  list of (3,) ndarrays
        One or more ECI reference unit vectors.  When multiple vectors are
        supplied the per-sigma measurements are concatenated, producing a
        row width of 3 * len(ref_vecs).
    """
    if isinstance(ref_vecs, np.ndarray) and ref_vecs.ndim == 1:
        ref_vecs = [ref_vecs]
    measurements = []
    for sigma in quat_sigmas:
        quat = Quaternion(sigma[:4])
        ECI_to_body = quat.inverse
        row = np.concatenate([rotate(rv, ECI_to_body) for rv in ref_vecs])
        measurements.append(row)
    return np.array(measurements)

def ensure_positive_definite(P, epsilon=1e-9):
    P = (P + P.T) / 2 
    eigvals = np.linalg.eigvalsh(P)
    if np.min(eigvals) < epsilon:
        P += np.eye(P.shape[0]) * (epsilon - np.min(eigvals) + epsilon)
    return P

def iterate(error_state, rotation: Quaternion, P, obs=None, ref_vecs=None, R=None):
    """Propagate the UKF by one step and, optionally, apply a measurement update.

    Parameters
    ----------
    error_state, rotation, P : current filter state.
    obs      : (3*k,) stacked body-frame observations for k reference vectors.
               Pass ``None`` (default) to perform a propagation-only step.
    ref_vecs : (3,) ndarray or list of k (3,) ndarrays.  Required when *obs*
               is given.  Each vector corresponds to a 3-element block of *obs*.
    R        : measurement noise covariance (m×m).  Defaults to
               ``global_R[0,0] * I`` sized to match the observation.

    Returns
    -------
    (x_hat, quaternion, P_new)
    """
    global global_Q
    n = error_state.shape[0]
    alpha = 0.1
    beta = 2
    P = ensure_positive_definite(P)

    lam = calculate_lambda(alpha, error_state)
    sigmas = get_sigma_points(lam, error_state, P)
    quat_sigmas = error_sigmas_to_quat_sigmas(sigmas, rotation)
    propagated_quat_sigmas = propagate_quat_sigmas(quat_sigmas)

    covariance_weights, mean_weights = get_weights(lam, error_state, alpha, beta)

    quaternions_of_propagated_sigmas = [Quaternion(y[:4]) for y in propagated_quat_sigmas]
    average_quaternion, propagated_error_vectors = gradient_descent(
        quaternions_of_propagated_sigmas, rotation, weights=mean_weights
    )
    propagated_errors = np.hstack([propagated_error_vectors, propagated_quat_sigmas[:, 4:7]])

    mean_error = np.zeros(n)
    for column in range(n):
        var_mean = 0
        for row in range(2*n + 1):
            var_mean += propagated_errors[row][column] * mean_weights[row]
        mean_error[column] = var_mean

    P_hat = np.zeros((n, n))
    for i, error in enumerate(propagated_errors):
        err = error - mean_error
        P_hat += covariance_weights[i] * np.outer(err, err)
    P_hat += global_Q

    # ── propagate-only: skip measurement update ────────────────────────────────
    if obs is None or ref_vecs is None:
        return mean_error, average_quaternion, ensure_positive_definite(P_hat)

    # ── measurement update ─────────────────────────────────────────────────────
    measurements = get_measurements(propagated_quat_sigmas, ref_vecs)

    mean_measurement = np.zeros(measurements.shape[1])
    for column in range(measurements.shape[1]):
        var_mean = 0
        for row in range(2*n + 1):
            var_mean += measurements[row][column] * mean_weights[row]
        mean_measurement[column] = var_mean

    P_xz = np.zeros((n, measurements.shape[1]))
    for i, error in enumerate(propagated_errors):
        err = error - mean_error
        msmt_err = measurements[i] - mean_measurement
        P_xz += covariance_weights[i] * np.outer(err, msmt_err)

    P_zz = np.zeros((measurements.shape[1], measurements.shape[1]))
    for i, msmt in enumerate(measurements):
        msmt_err = msmt - mean_measurement
        P_zz += covariance_weights[i] * np.outer(msmt_err, msmt_err)

    # Build R to match the observation dimension (3 per reference vector)
    if R is None:
        R = np.eye(measurements.shape[1]) * sigma_magnetometer**2
    P_vv = P_zz + R
    k = P_xz @ np.linalg.inv(P_vv)
    x_hat = mean_error + k @ (obs - mean_measurement)

    P = P_hat - k @ P_vv @ k.T
    P = ensure_positive_definite(P)

    x_hat_rot = rotation_to_quat(x_hat[:3])
    return x_hat, x_hat_rot * average_quaternion, P

def quat_diff(q1, q2):
    q_err = q1 * q2.inverse
    return np.linalg.norm(quaternion_to_rotation(q_err))

if __name__ == '__main__':
    # 6x6 Covariance and 6-element state
    P = np.eye(6) * 0.01
    state = np.zeros(6)

    orbit_sim = OrbitIGRFSimulator(epoch=datetime(2024, 6, 1, 0, 0, 0))

    # Initialise true attitude from perfect nadir pointing at t=0
    _b0, pos0, vel0 = orbit_sim.step(0.0)
    true_rot = nadir_quaternion(pos0, vel0)            # cold start (no q_guess)
    true_rot_prev = true_rot
    rot = Quaternion(true_rot.elements + np.random.normal(loc=0, scale=0.05, size=4)).normalised #sets the first guess to decently close but not too close

    np.set_printoptions(linewidth=300, threshold=np.inf, suppress=True, precision = 6)
    
    ref_vec_2 = np.array([1.0, 0.0, 0.0])   # second ECI reference vector (fixed)

    for i in range(simulation_time):
        # Toggle 1-vector ↔ 2-vector mode every switch_every steps
        if i > 0 and i % switch_every == 0:
            vector_input = 3 - vector_input   # 1 → 2 or 2 → 1

        # ── simulation environment (10 Hz) ────────────────────────────────────
        b_eci_uT, pos_eci, vel_eci = orbit_sim.step(dt)
        b_uT_mag = np.linalg.norm(b_eci_uT)
        reference_vector = b_eci_uT / b_uT_mag   # unit vector for filter

        true_rot = nadir_quaternion(pos_eci, vel_eci, q_guess=true_rot_prev)

        # True angular velocity from quaternion delta (captures actual nadir-tracking motion)
        delta_q = true_rot.inverse * true_rot_prev
        true_angular_velocity = quaternion_to_rotation(delta_q) / dt
        true_rot_prev = true_rot

        gyro_noise = np.random.normal(loc=0, scale=sigma_gyro, size=3)
        gyro_measurement = true_angular_velocity + true_gyro_bias + gyro_noise

        # ── filter step ───────────────────────────────────────────────────────
        try:
            if i % update_every == 0:
                # Propagate + measurement update (0.1 Hz)
                ref_to_body = true_rot.inverse
                if vector_input == 1:
                    b_body_uT = rotate(b_eci_uT, ref_to_body)
                    r1 = b_body_uT + np.random.normal(loc=0, scale=sigma_magnetometer, size=3)
                    measurement = r1 / np.linalg.norm(r1)
                    ref_vecs = reference_vector
                    R_update = np.eye(3) * 2e-5
                else:
                    b_body_uT = rotate(b_eci_uT, ref_to_body)
                    r1 = b_body_uT + np.random.normal(loc=0, scale=sigma_magnetometer, size=3)
                    r2 = rotate(ref_vec_2, ref_to_body) + np.random.normal(loc=0, scale=sigma_unit_vector, size=3)
                    measurement = np.concatenate([r1 / np.linalg.norm(r1), r2 / np.linalg.norm(r2)])
                    ref_vecs = [reference_vector, ref_vec_2]
                    R_update = np.block([
                        [np.eye(3) * 2e-5, np.zeros((3, 3))],
                        [np.zeros((3, 3)), np.eye(3) * sigma_unit_vector**2],
                    ])
                state, rot, P = iterate(state, rot, P, measurement, ref_vecs, R_update)
            else:
                # Propagate only (10 Hz)
                state, rot, P = iterate(state, rot, P)

            if i % 100 == 0:
                error_deg = math.degrees(quat_diff(rot, true_rot))
                print(f"{int(i * dt):4d} seconds in | Mode: {vector_input}v | Magnitude of attitude err: {error_deg:6.3f} deg | Bias Est: {state[3:6]}")

            state[:3] = np.zeros(3)

        except Exception as e:
            print(f"Filter failed at step {i}: {e}")
            break