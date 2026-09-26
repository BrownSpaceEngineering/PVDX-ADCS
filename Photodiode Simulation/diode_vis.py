import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

all_vectors = np.array([
    [ 0.8660254, -0.27968387, -0.41445981],
    [ 0.8660254, -0.16748457,  0.47111455],
    [ 0.8660254,  0.4999767,   0.00482655],
    [-0.8660254,  0.27968387,  0.41445981],
    [-0.8660254,  0.16748457, -0.47111455],
    [-0.8660254, -0.4999767,  -0.00482655],
    [-0.07387608,  0.8660254,  0.49451221],
    [-0.48736734,  0.8660254, -0.11168291],
    [ 0.32616791,  0.8660254, -0.37896503],
    [ 0.07387608, -0.8660254, -0.49451221],
    [ 0.48736734, -0.8660254,  0.11168291],
    [-0.32616791, -0.8660254,  0.37896503],
    [-0.36819159, -0.33828235, -0.8660254],
    [ 0.48601695,  0.1174203,  -0.8660254],
    [-0.23581003,  0.44090093, -0.8660254],
    [ 0.36819159,  0.33828235,  0.8660254],
    [-0.48601695, -0.1174203,   0.8660254],
    [ 0.23581003, -0.44090093,  0.8660254],
])

cos_cone = np.cos(np.deg2rad(70.0))

n_theta, n_phi = 60, 120
theta = np.linspace(0, np.pi,   n_theta)
phi   = np.linspace(0, 2*np.pi, n_phi)
T, P  = np.meshgrid(theta, phi, indexing='ij')

X = np.sin(T) * np.cos(P)
Y = np.sin(T) * np.sin(P)
Z = np.cos(T)

sample_dirs = np.stack([X.ravel(), Y.ravel(), Z.ravel()], axis=1)
counts = (sample_dirs @ all_vectors.T > cos_cone).sum(axis=1).reshape(n_theta, n_phi)

fig = plt.figure(figsize=(11, 9), facecolor='#0d0d0d')
ax  = fig.add_subplot(111, projection='3d', facecolor='#0d0d0d')

cmap = plt.colormaps['RdYlGn']
norm = mcolors.Normalize(vmin=counts.min(), vmax=counts.max())

ax.plot_surface(X, Y, Z,
                facecolors=cmap(norm(counts)),
                rstride=1, cstride=1,
                linewidth=0, antialiased=False,
                shade=False, alpha=0.95)

for v in all_vectors:
    ax.quiver(0, 0, 0, v[0], v[1], v[2],
              length=1.08, color='white',
              arrow_length_ratio=0.12, linewidth=1.2, alpha=0.7)

for axis, col in zip([(1,0,0),(0,1,0),(0,0,1)], ['#ff5555','#55ff55','#5599ff']):
    for sign in [1, -1]:
        a = np.array(axis, dtype=float) * sign * 1.35
        ax.quiver(0, 0, 0, a[0], a[1], a[2],
                  length=1.0, color=col, linewidth=1.0, alpha=0.4,
                  arrow_length_ratio=0.07, linestyle='dashed')

ax.set_xlim(-1.3, 1.3); ax.set_ylim(-1.3, 1.3); ax.set_zlim(-1.3, 1.3)
ax.set_xlabel('X', color='white', labelpad=8)
ax.set_ylabel('Y', color='white', labelpad=8)
ax.set_zlabel('Z', color='white', labelpad=8)
ax.tick_params(colors='#555555', labelsize=8)
for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
    pane.fill = False
    pane.set_edgecolor('#1a1a1a')
ax.grid(False)
ax.set_title('Diodes within 70° cone — coverage across full sphere',
             color='white', fontsize=12, pad=16)

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.08)
cbar.set_label('# diodes in cone', color='white', fontsize=10)
cbar.ax.yaxis.set_tick_params(color='white', labelcolor='white')
cbar.set_ticks(range(counts.min(), counts.max()+1))

plt.tight_layout()
plt.show()
print(f"Samples: {n_theta}x{n_phi} = {n_theta*n_phi}")