import numpy as np
def cartesian_to_spherical(pt):
    pt = np.array(pt)
    r = np.linalg.norm(pt)
    theta = np.arctan2(pt[1], pt[0])
    phi = np.arccos(pt[2]/r)
    return np.array([r, theta, phi])

normal_vectors = {
    "Plane 1": [1, 0, 0],
    "Plane 2": [0, 1, 0],
    "Plane 3": [0, -1, 0],
    "Plane 4": [0, 0, -1],
    "Plane 5": [0 ,0, 1],
    "Plane 6": [-1, 0, 0]
}
def angle_between(u, v, degrees=False):
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)

    # Compute norms
    nu = np.linalg.norm(u)
    nv = np.linalg.norm(v)

    if nu == 0 or nv == 0:
        raise ValueError("Angle is undefined for zero-length vectors")

    # Cosine of angle (clipped for numerical stability)
    cos_theta = np.dot(u, v) / (nu * nv)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    theta = np.arccos(cos_theta)

    return np.degrees(theta) if degrees else theta

def dowork(pts):
    original_pts = []
    for pt in range(9):
        pts.append(-1 * np.array(pts[pt]))
        original_pts.append(pts[pt])
    pts = np.array(pts)
    max_angle = 45#degrees
    mapped = {
        "Plane 1": [],
        "Plane 2": [],
        "Plane 3": [],
        "Plane 4": [],
        "Plane 5": [],
        "Plane 6": [],
    }
    original_mapped = {
        "Plane 1": [],
        "Plane 2": [],
        "Plane 3": [],
        "Plane 4": [],
        "Plane 5": [],
        "Plane 6": [],
    }
    missing = False
    worst = 0
    for pt in pts:
        pt = np.array(pt)/np.linalg.norm(np.array(pt))
        pt_mapped = False
        minimum_angle = 1000000
        best_plane = None
        for _, (plane, nv) in enumerate(normal_vectors.items()):
            nv = np.array(nv)
            angle = angle_between(nv, pt)
            if(angle < np.deg2rad(max_angle)):
                best_plane = plane
                pt_mapped = True
                minimum_angle = min(minimum_angle, angle)
        if not pt_mapped:
            missing = True
        else:
            mapped[best_plane].append((pt, minimum_angle))
            worst = max(worst, minimum_angle)
    for pt in original_pts:
        pt = np.array(pt)/np.linalg.norm(np.array(pt))
        pt_mapped = False
        minimum_angle = 1000000
        best_plane = None
        for _, (plane, nv) in enumerate(normal_vectors.items()):
            nv = np.array(nv)
            angle = angle_between(nv, pt)
            if(angle < np.deg2rad(max_angle)):
                best_plane = plane
                pt_mapped = True
                minimum_angle = min(minimum_angle, angle)
        if not pt_mapped:
            missing = True
        else:
            original_mapped[best_plane].append((pt, minimum_angle))
    return missing, mapped, worst, original_mapped

