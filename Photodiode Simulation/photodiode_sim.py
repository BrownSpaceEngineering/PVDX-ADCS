import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

n_sims           = 10000
cone_deg         = 65.0
max_voltage_mean = 1.0
assumed_voltages = [0.5, 1.0, 1.5]
std_levels       = [0.0, 0.05, 0.1, 0.25]

def noise_params(theta_deg, noise_fraction, true_max_voltage):
    """
    Flat Gaussian noise throughout. Past 65°, adds a positive voltage bias
    of 0.025 * (theta - 65) to model the cosine model breakdown seen in real data.
    Returns (mean, std) of the noise distribution.
    """
    std = noise_fraction * true_max_voltage
    if theta_deg <= 65.0:
        mean = 0.0
    else:
        mean = .005 * (theta_deg - 65.0)
    return mean, std


class PhotodiodePair:
    def __init__(self, vec_a, vec_b, noise_fraction, atol=1e-5):
        vec_a = np.array(vec_a, dtype=float)
        vec_b = np.array(vec_b, dtype=float)
        if not np.isclose(np.linalg.norm(vec_a), 1.0, atol=atol):
            raise ValueError("vec_a is not a unit vector")
        if not np.isclose(np.linalg.norm(vec_b), 1.0, atol=atol):
            raise ValueError("vec_b is not a unit vector")
        if not np.allclose(vec_a, -vec_b, atol=atol):
            raise ValueError("Vectors are not opposite.")
        self.vec_a          = vec_a
        self.vec_b          = vec_b
        self.noise_fraction = noise_fraction

    def reading(self, sun_vector, true_max_voltage):
        sun_vector = np.array(sun_vector, dtype=float)
        # Raw dot products (can be negative for back-facing diodes)
        raw_dot_a = np.dot(self.vec_a, sun_vector)
        raw_dot_b = np.dot(self.vec_b, sun_vector)
        # True angle from raw dot (not clamped) — so back-facing diodes get >90°
        theta_a = np.rad2deg(np.arccos(np.clip(raw_dot_a, -1.0, 1.0)))
        theta_b = np.rad2deg(np.arccos(np.clip(raw_dot_b, -1.0, 1.0)))
        # Voltage is clamped to 0 for back-facing diodes
        v_a = true_max_voltage * max(0.0, raw_dot_a)
        v_b = true_max_voltage * max(0.0, raw_dot_b)
        # Bias only applies if diode is within its forward hemisphere (theta <= 90°)
        # Back-facing diodes (theta > 90°) get flat noise only
        mean_a, std_a = noise_params(theta_a if theta_a <= 90.0 else 0.0,
                                     self.noise_fraction, true_max_voltage)
        mean_b, std_b = noise_params(theta_b if theta_b <= 90.0 else 0.0,
                                     self.noise_fraction, true_max_voltage)
        v_a = max(0.0, v_a + np.random.normal(mean_a, std_a))
        v_b = max(0.0, v_b + np.random.normal(mean_b, std_b))
        if v_a >= v_b:
            return v_a, self.vec_a
        else:
            return v_b, self.vec_b


def least_squares_solve(vectors, voltages, assumed_max_voltage):
    A = np.array(vectors)
    b = np.array(voltages) / assumed_max_voltage
    sun_vec = (np.linalg.inv(A.T @ A)) @ A.T @ b
    return sun_vec / np.linalg.norm(sun_vec)


def estimate_cone(vecs, volts, assumed_max, cone_deg=cone_deg):
    min_v = assumed_max * np.cos(np.deg2rad(cone_deg))
    lit_vecs = [v for v, volt in zip(vecs, volts) if volt > min_v]
    lit_volts = [volt for volt in volts if volt > min_v]
    if len(lit_vecs) < 3:
        return None
    return least_squares_solve(lit_vecs, lit_volts, assumed_max)


def estimate_top5(vecs, volts, assumed_max):
    ranked = sorted(zip(volts, vecs), key=lambda x: x[0], reverse=True)[:5]
    return least_squares_solve([v for _, v in ranked], [v for v, _ in ranked], assumed_max)


def estimate_all9(vecs, volts, assumed_max):
    return least_squares_solve(vecs, volts, assumed_max)


pairs = [
    PhotodiodePair([ 0.8660254, -0.27968387, -0.41445981], [-0.8660254,  0.27968387,  0.41445981], noise_fraction=0.025),
    PhotodiodePair([ 0.8660254, -0.16748457,  0.47111455], [-0.8660254,  0.16748457, -0.47111455], noise_fraction=0.025),
    PhotodiodePair([ 0.8660254,  0.4999767,   0.00482655], [-0.8660254, -0.4999767,  -0.00482655], noise_fraction=0.025),
    PhotodiodePair([-0.07387608,  0.8660254,  0.49451221], [ 0.07387608, -0.8660254, -0.49451221], noise_fraction=0.025),
    PhotodiodePair([-0.48736734,  0.8660254, -0.11168291], [ 0.48736734, -0.8660254,  0.11168291], noise_fraction=0.025),
    PhotodiodePair([ 0.32616791,  0.8660254, -0.37896503], [-0.32616791, -0.8660254,  0.37896503], noise_fraction=0.025),
    PhotodiodePair([-0.36819159, -0.33828235, -0.8660254], [ 0.36819159,  0.33828235,  0.8660254], noise_fraction=0.025),
    PhotodiodePair([ 0.48601695,  0.1174203,  -0.8660254], [-0.48601695, -0.1174203,   0.8660254], noise_fraction=0.025),
    PhotodiodePair([-0.23581003,  0.44090093, -0.8660254], [ 0.23581003, -0.44090093,  0.8660254], noise_fraction=0.025),
]

# results[assumed_idx][std_idx] = {'cone', 'top5', 'all9', 'fails'}
results = []
for assumed_max in assumed_voltages:
    row = []
    for mv_std in std_levels:
        errs_cone, errs_top5, errs_all9 = [], [], []
        fails = 0
        for _ in range(n_sims):
            v = np.random.randn(3)
            true_vector = v / np.linalg.norm(v)
            if mv_std == 0.0:
                true_maxes = np.full(len(pairs), max_voltage_mean)
            else:
                true_maxes = np.clip(
                    np.random.normal(max_voltage_mean, mv_std, len(pairs)), 0.05, None)
            active_vecs, active_volts = [], []
            for pair, tm in zip(pairs, true_maxes):
                voltage, active_vec = pair.reading(true_vector, tm)
                active_vecs.append(active_vec)
                active_volts.append(voltage)

            result = estimate_cone(active_vecs, active_volts, assumed_max, cone_deg)
            if result is None:
                fails += 1
            else:
                dot = np.clip(true_vector @ result, -1.0, 1.0)
                errs_cone.append(np.rad2deg(np.arccos(dot)))

            result = estimate_top5(active_vecs, active_volts, assumed_max)
            dot = np.clip(true_vector @ result, -1.0, 1.0)
            errs_top5.append(np.rad2deg(np.arccos(dot)))

            result = estimate_all9(active_vecs, active_volts, assumed_max)
            dot = np.clip(true_vector @ result, -1.0, 1.0)
            errs_all9.append(np.rad2deg(np.arccos(dot)))

        row.append({'cone': np.array(errs_cone),
                    'top5': np.array(errs_top5),
                    'all9': np.array(errs_all9),
                    'fails': fails})
    results.append(row)

# ── One figure per assumed max ────────────────────────────────────────────────
method_keys   = ['cone',     'top5',     'all9']
method_labels = ['Cone filter', 'Top 5', 'All 9']
method_colors = ['#3a9edb',  '#e8834a',  '#a855f7']
std_labels    = [f'σ = {s} V' for s in std_levels]

figs = []
for ci, assumed_max in enumerate(assumed_voltages):

    # Shared x-axis limit across all panels in this figure
    all_errs = np.concatenate([results[ci][si][k]
                                for si in range(len(std_levels))
                                for k in method_keys
                                if len(results[ci][si][k]) > 0])
    x_max = np.percentile(all_errs, 99) * 1.15 if len(all_errs) > 0 else 15.0
    bins  = np.linspace(0, x_max, 40)

    # Layout: rows = std levels, cols = methods
    n_rows = len(std_levels)
    n_cols = len(method_keys)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(5.5 * n_cols, 3.0 * n_rows),
                             facecolor='#0d0d0d',
                             sharex=True, sharey=False)
    title_color = '#88ff88' if assumed_max == max_voltage_mean else 'white'
    fig.suptitle(
        f'Assumed Max = {assumed_max} V   '
        f'(true max ~ N({max_voltage_mean}V, σ²) per diode per run  |  '
        f'noise = 2.5% flat, bias +0.005×(θ-65)V past 65°  |  {n_sims} sims)',
        color=title_color, fontsize=13, y=0.98)

    # Column headers (method names) — top row only
    for col, (mlabel, mcolor) in enumerate(zip(method_labels, method_colors)):
        axes[0][col].set_title(mlabel, color=mcolor, fontsize=12, fontweight='bold')

    for si, (mv_std, slabel) in enumerate(zip(std_levels, std_labels)):
        for col, (mkey, mcolor) in enumerate(zip(method_keys, method_colors)):
            ax   = axes[si][col]
            errs = results[ci][si][mkey]
            fails = results[ci][si]['fails']

            ax.set_facecolor('#1a1a1a')
            ax.tick_params(colors='#aaaaaa', labelsize=8)
            for spine in ax.spines.values():
                spine.set_edgecolor('#444444')

            if len(errs) > 0:
                ax.hist(errs, bins=bins, color=mcolor,
                        edgecolor='#0d0d0d', linewidth=0.3, alpha=0.88)
                ax.axvline(np.mean(errs),   color='#ffcc00', linewidth=1.4, linestyle='--')
                ax.axvline(np.median(errs), color='#ff6655', linewidth=1.2, linestyle=':')

            # Row label on left
            if col == 0:
                ax.set_ylabel(f'{slabel}\nCount', color='#cccccc',
                              fontsize=9, labelpad=4)

            # Stats box
            if len(errs) > 0:
                fail_str = f'\n{fails} fails' if (mkey == 'cone' and fails > 0) else ''
                ax.text(0.97, 0.95,
                        f'μ = {np.mean(errs):.2f}°\n'
                        f'med = {np.median(errs):.2f}°{fail_str}',
                        transform=ax.transAxes, color='#dddddd', fontsize=8,
                        ha='right', va='top',
                        bbox=dict(facecolor='#2a2a2a', alpha=0.7, edgecolor='none'))
            elif mkey == 'cone':
                ax.text(0.5, 0.5, f'ALL {fails} FAILED',
                        transform=ax.transAxes, color='#ff6655', fontsize=11,
                        ha='center', va='center', fontweight='bold')

            if si == n_rows - 1:
                ax.set_xlabel('Error (°)', color='#888888', fontsize=8)

    legend_elems = [
        Line2D([0], [0], color='#ffcc00', linestyle='--', linewidth=1.4, label='Mean'),
        Line2D([0], [0], color='#ff6655', linestyle=':',  linewidth=1.2, label='Median'),
    ]
    fig.legend(handles=legend_elems, loc='lower center', ncol=1,
               fontsize=9, framealpha=0.3, labelcolor='white',
               facecolor='#2a2a2a', edgecolor='#555555',
               bbox_to_anchor=(0.5, -0.02))

    plt.tight_layout(pad=1.5)
    plt.subplots_adjust(top=0.93, left=0.08)
    figs.append(fig)

plt.show()