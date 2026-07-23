#!/usr/bin/env python3
# =============================================================================
# SEM morphology batch: sphericity + surface porosity/skin for ALL samples
# =============================================================================
# Runs the same metric extraction as sem_morphology.py + sem_porosity.py over
# every distinct sample in sem/ (auto-discovered), except 114645 (excluded).
# Handles both export sizes (947x725 and 710x543) and both magnification-naming
# conventions (_5kx0 / _a5kx0 / _5000x0 for the surface; _a500x0 / _b500x0 for
# whole particles) via auto-cropping and scale-aware thresholds.
#
#   sphericity Omega_SEM  = median particle solidity        (500x)
#   porosity   = % open space between spheres (granular)    (5k/5000x)
#   skin       = smooth gap-free merged fraction            (5k/5000x)
# See review/sem_morphology_metrics.md for the definitions and the model map.
#
# Run:  python3 sem_all.py  ->  data/sem_all_morphology.csv + sem/overlay_X*.png
# =============================================================================
import glob, re, csv, warnings, numpy as np
warnings.filterwarnings("ignore")
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology, measure, segmentation, feature

EXCLUDE = {"114645"}
REF_AREA = 516000.0     # reference cropped area (947x725 export) for scale-aware params
# Surface classification. SMOOTH is defined by an ABSOLUTE local-range threshold
# (normalised by image contrast) so a fused surface reads as more skin than a
# granular one -- a percentile would self-normalise and erase that difference.
#   TAU_NORM : local range / body-sigma below this = smooth (merged skin)
#   DARK_K   : within the granular packing, flat < mean - DARK_K*sigma = open gap
TAU_NORM, DARK_K, BRIGHT_K = 1.5, 0.30, 0.40

def auto_crop(path):
    a = np.asarray(Image.open(path).convert("L"), float)
    rm, cm = a.mean(1), a.mean(0)
    rows, cols = np.where(rm < 245)[0], np.where(cm < 245)[0]
    if len(rows) < 10 or len(cols) < 10:
        return a
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    H, W = r1 - r0, c1 - c0
    # inset to drop the burned-in title (top) and scale-bar (bottom) text banners
    return a[r0 + int(0.06 * H): r1 - int(0.07 * H),
             c0 + int(0.03 * W): c1 - int(0.02 * W)]

def scales(img):
    ar = img.size / REF_AREA          # area ratio
    lin = max(ar ** 0.5, 0.4)         # linear ratio (for window/distance px)
    return ar, lin

def isolate_particles(img):
    """Particle vs non-particle (stub/background). Particle = raised (bright,
    large-scale) AND/OR textured; stub = darker + flatter. Keeps fused particles
    (bright even when smooth) and excludes the smooth substrate between grains."""
    ar, lin = scales(img)
    s = max(8, int(12 * lin)); w = max(3, int(5 * lin))
    B = ndi.gaussian_filter(img, s)
    rng = ndi.maximum_filter(img, w) - ndi.minimum_filter(img, w)
    T = ndi.gaussian_filter(rng, s)
    def nrm(a):
        lo, hi = np.percentile(a, 2), np.percentile(a, 98)
        return np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1)
    m = (0.5 * nrm(B) + 0.5 * nrm(T)) > filters.threshold_otsu(0.5 * nrm(B) + 0.5 * nrm(T))
    m = ndi.binary_fill_holes(ndi.binary_closing(m, iterations=max(1, int(2 * lin))))
    m = morphology.remove_small_objects(m, int(4000 * ar))
    m = morphology.remove_small_objects(ndi.binary_opening(m, iterations=max(1, int(lin))),
                                        int(4000 * ar))
    return m if m.sum() > 0.1 * img.size else np.ones_like(img, bool)

# ---- sphericity (whole particles, 500x) -------------------------------------
def sphericity(img):
    ar, lin = scales(img)
    thr = filters.threshold_otsu(img)
    bw = ndi.binary_fill_holes(morphology.remove_small_objects(img > thr, int(60 * ar)))
    dist = ndi.distance_transform_edt(bw)
    peaks = feature.peak_local_max(dist, min_distance=max(4, int(7 * lin)), labels=bw)
    mrk = np.zeros(dist.shape, int)
    for i, (r, c) in enumerate(peaks, 1):
        mrk[r, c] = i
    lab = segmentation.clear_border(segmentation.watershed(-dist, mrk, mask=bw))
    sol, asp = [], []
    for p in measure.regionprops(lab):
        if p.area < 80 * ar or p.area > 0.15 * img.size or p.axis_major_length <= 0:
            continue
        sol.append(p.solidity); asp.append(p.axis_minor_length / p.axis_major_length)
    if len(sol) < 5:
        return None
    return dict(n=len(sol), solidity=float(np.median(sol)), aspect=float(np.median(asp)))

# ---- surface porosity + skin (5k/5000x) -------------------------------------
def surface(img, overlay=None):
    ar, lin = scales(img)
    sig = max(8, 18 * lin); win = max(3, int(round(5 * lin)))
    flat = img - ndi.gaussian_filter(img, sig)
    rng = ndi.maximum_filter(img, win) - ndi.minimum_filter(img, win)
    # PARTICLE-FIRST: isolate the particle body (excludes the stub/background),
    # then measure skin/porosity only inside it.
    body = isolate_particles(img)
    fb = flat[body]; med = np.median(fb)
    body_sigma = max(img[body].std(), 1)
    bright_cut = fb.mean() + BRIGHT_K * fb.std()
    # SMOOTH: absolute normalised local range below TAU_NORM (fused -> more smooth)
    smooth = (rng / body_sigma) < TAU_NORM
    cap = body & (flat >= bright_cut)
    skin = morphology.remove_small_objects(body & smooth & (flat < med), int(120 * ar))
    granular = body & (~smooth)
    # OPEN GAP inside the granular packing: darker than the local average there
    gf = flat[granular]
    dark_cut = (gf.mean() - DARK_K * gf.std()) if granular.sum() > 100 else -1e9
    pore = granular & (flat < dark_cut) & (~skin)
    A, G = body.sum(), max(granular.sum(), 1)
    if overlay is not None:
        base = np.dstack([img, img, img]).astype(float)
        A = 0.5   # overlay alpha -> underlying texture stays visible
        for mask, col in ((pore, (220, 40, 40)), (skin, (40, 90, 230))):
            base[mask] = (1 - A) * base[mask] + A * np.array(col)
        Image.fromarray(np.clip(base, 0, 255).astype(np.uint8)).save(overlay)
    return dict(porosity_gran=pore.sum() / G, porosity_surf=pore.sum() / A,
                skin=skin.sum() / A, granular=G / A)

# ---- discover samples -------------------------------------------------------
sids = sorted({m.group(1) for f in glob.glob("sem/SEM_X*.jpg")
               for m in [re.search(r"SEM_X(\d+)", f)]} - EXCLUDE)
CONDMAP = {"114641": "up3_1", "114642": "up3_2", "114643": "up3_3", "114644": "up3_4"}

def pick(sid, patterns):
    out = []
    for pat in patterns:
        out += glob.glob(f"sem/SEM_X{sid}{pat}")
    return sorted(set(out))

rows = []
for sid in sids:
    lo = pick(sid, ["a_500x0.jpg", "b_500x0.jpg", "_a500x0.jpg", "_b500x0.jpg"])
    hi = pick(sid, ["_5kx0.jpg", "_a5kx0.jpg", "_b5kx0.jpg", "_5000x0.jpg", "_10kx0.jpg"])
    sph = [sphericity(auto_crop(p)) for p in lo]; sph = [s for s in sph if s]
    srf = []
    for i, p in enumerate(hi):
        ov = f"sem/overlay_X{sid}.png" if i == 0 else None
        srf.append(surface(auto_crop(p), overlay=ov))
    def avg(lst, k):
        v = [d[k] for d in lst if d]
        return float(np.mean(v)) if v else float("nan")
    row = dict(sample=sid, cond=CONDMAP.get(sid, ""),
        n_part=round(avg(sph, "n")) if sph else 0,
        Omega_SEM=round(avg(sph, "solidity"), 3), aspect_SEM=round(avg(sph, "aspect"), 3),
        porosity_pct=round(100 * avg(srf, "porosity_gran"), 1),
        porosity_surf_pct=round(100 * avg(srf, "porosity_surf"), 1),
        skin_pct=round(100 * avg(srf, "skin"), 1),
        granular_pct=round(100 * avg(srf, "granular"), 1))
    rows.append(row)
    print(f"X{sid} {row['cond']:6}: Omega={row['Omega_SEM']} aspect={row['aspect_SEM']} "
          f"| porosity={row['porosity_pct']}% skin={row['skin_pct']}% "
          f"granular={row['granular_pct']}%  (n={row['n_part']})")

with open("data/sem_all_morphology.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"\nWrote data/sem_all_morphology.csv ({len(rows)} samples) + sem/overlay_X*.png")
