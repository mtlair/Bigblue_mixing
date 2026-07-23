#!/usr/bin/env python3
# =============================================================================
# SEM morphology metrics for the dried product (samples 114641-114644 = up3_1..4)
# =============================================================================
# Defines and extracts three image metrics that map to the model's morphology
# outputs (up2_spray_dryer_module.R), so the SEM data can calibrate the chain:
#
#   sphericity  Omega_SEM   <-> model Omega_struct_z   (low-mag whole particles)
#   porosity    phi_SEM     <-> model phi_porosity_z   (high-mag surface pores)
#   skin        theta_SEM   <-> model theta_skin_z     (high-mag surface fusion)
#
# METRIC DEFINITIONS
# ------------------
# Sphericity (from 500x): segment bright particles on the dark stub, split
#   touching particles by watershed, and per particle take
#     solidity   = area / convex-hull area          (1 = convex/round)
#     aspect     = minor axis / major axis           (1 = circular)
#     circ       = 4*pi*area / perimeter^2           (1 = perfect circle)
#   Report the median over all accepted particles. Omega_SEM := median solidity
#   (robust to perimeter pixelation; aspect/circ reported alongside).
#
# Surface porosity (from 5kx/10kx): mask the in-focus particle body, flatten the
#   curvature shading (subtract a large-sigma Gaussian), and classify the dark
#   interstitial depressions between primary beads as pores.
#     phi_SEM := pore-pixel area / particle-body area
#
# Surface skin / fusion (from 5kx/10kx): an unfused surface is a granular packing
#   of discrete primary beads (high local texture); a fused skin is smooth
#   (low local texture). Within the particle mask, compute the local texture
#   (std in a small window, as a fraction of local mean) and
#     theta_SEM := fraction of surface with local texture below TAU_SMOOTH
#   i.e. the smooth/fused area fraction. Also report rms_rough (mean texture CV)
#   and bump_density (primary-bead peaks per 100 um^2 when scale is known).
#
# Run:  python3 sem_morphology.py   ->  data/sem_morphology.csv
# =============================================================================
import os, glob, csv, warnings, numpy as np
warnings.filterwarnings("ignore")
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation, feature

SEM_DIR = "sem"
SAMPLES = {"114641": "up3_1", "114642": "up3_2", "114643": "up3_3", "114644": "up3_4"}

# inner crop avoiding the burned-in text banners (top ~30px, bottom ~40px)
R0, R1, C0, C1 = 70, 660, 45, 920
TAU_SMOOTH = 0.06   # local texture CV below this = "smooth/fused" surface pixel

def load(path):
    a = np.asarray(Image.open(path).convert("L"), float)
    return a[R0:R1, C0:C1]

# ---- sphericity (low-mag whole particles) -----------------------------------
def sphericity(path):
    img = load(path)
    thr = filters.threshold_otsu(img)
    bw = img > thr
    bw = morphology.remove_small_objects(bw, 60)
    bw = ndi.binary_fill_holes(bw)
    # split touching particles by watershed on the distance transform
    dist = ndi.distance_transform_edt(bw)
    peaks = feature.peak_local_max(dist, min_distance=7, labels=bw)
    mrk = np.zeros(dist.shape, int)
    for i, (r, c) in enumerate(peaks, 1):
        mrk[r, c] = i
    lab = segmentation.watershed(-dist, mrk, mask=bw)
    lab = segmentation.clear_border(lab)
    sol, asp, cir, npart = [], [], [], 0
    for p in measure.regionprops(lab):
        if p.area < 80 or p.area > 0.15 * img.size:
            continue
        if p.perimeter <= 0 or p.axis_major_length <= 0:
            continue
        sol.append(p.solidity)
        asp.append(p.axis_minor_length / p.axis_major_length)
        cir.append(min(1.0, 4 * np.pi * p.area / p.perimeter ** 2))
        npart += 1
    if npart < 5:
        return None
    return dict(n=npart, solidity=np.median(sol),
                aspect=np.median(asp), circ=np.median(cir))

# ---- surface porosity + skin (high-mag) -------------------------------------
def surface(path):
    img = load(path)
    # particle body mask: bright in-focus region (drop dark background/shadow)
    body = img > filters.threshold_otsu(img) * 0.55
    body = morphology.remove_small_objects(body, 2000)
    body = ndi.binary_fill_holes(body)
    if body.sum() < 0.2 * img.size:
        body = np.ones_like(img, bool)          # frame is all particle
    # flatten curvature shading
    flat = img - ndi.gaussian_filter(img, 25)
    m = body
    # pores: dark depressions below the flattened surface
    neg = flat[m]
    pore = (flat < (neg.mean() - 1.0 * neg.std())) & m
    phi = pore.sum() / m.sum()
    # local texture CV (std/mean in a 5px window)
    loc_mean = ndi.uniform_filter(img, 5)
    loc_sq = ndi.uniform_filter(img ** 2, 5)
    loc_std = np.sqrt(np.maximum(loc_sq - loc_mean ** 2, 0))
    cv = loc_std / np.maximum(loc_mean, 1)
    rms_rough = cv[m].mean()
    theta = (cv[m] < TAU_SMOOTH).mean()          # smooth/fused area fraction
    # primary-bead peak density (per megapixel of body)
    pk = feature.peak_local_max(ndi.gaussian_filter(img, 1.5),
                                min_distance=3, labels=m)
    bump_density = len(pk) / (m.sum() / 1e6)
    return dict(phi=phi, theta=theta, rms_rough=rms_rough,
                bump_density=bump_density)

def mean_over(paths, fn, keys):
    vals = [fn(p) for p in paths]
    vals = [v for v in vals if v]
    if not vals:
        return {k: float("nan") for k in keys}
    return {k: float(np.mean([v[k] for v in vals])) for k in keys}

rows = []
for sid, cond in SAMPLES.items():
    lo = sorted(glob.glob(f"{SEM_DIR}/SEM_X{sid}*500x0.jpg"))      # 500x a & b
    hi = sorted(glob.glob(f"{SEM_DIR}/SEM_X{sid}_5kx0.jpg") +
                glob.glob(f"{SEM_DIR}/SEM_X{sid}_10kx0.jpg"))       # 5k + 10k
    sph = mean_over(lo, sphericity, ["n", "solidity", "aspect", "circ"])
    srf = mean_over(hi, surface, ["phi", "theta", "rms_rough", "bump_density"])
    rows.append(dict(cond=cond, sample=sid,
        n_particles=round(sph["n"]) if sph["n"] == sph["n"] else "",
        Omega_SEM=round(sph["solidity"], 3), aspect_SEM=round(sph["aspect"], 3),
        circ_SEM=round(sph["circ"], 3),
        phi_SEM=round(srf["phi"], 4), theta_SEM=round(srf["theta"], 3),
        rms_rough=round(srf["rms_rough"], 4),
        bump_density=round(srf["bump_density"], 1)))
    print(f"{cond} ({sid}): Omega_SEM={rows[-1]['Omega_SEM']} aspect={rows[-1]['aspect_SEM']} "
          f"phi_SEM={rows[-1]['phi_SEM']} theta_SEM={rows[-1]['theta_SEM']} "
          f"rough={rows[-1]['rms_rough']} bumps/Mpx={rows[-1]['bump_density']}  (n={rows[-1]['n_particles']})")

os.makedirs("data", exist_ok=True)
with open("data/sem_morphology.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("\nWrote data/sem_morphology.csv")
