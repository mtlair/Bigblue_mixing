#!/usr/bin/env python3
# =============================================================================
# SEM surface porosity: % open space between primary spheres  (skin vs gap)
# =============================================================================
# Per the reviewer's definition (114641_5kx0 annotation):
#   * OPEN PORE  = interstitial gap between packing spheres (dark thin network,
#                  bounded by bright sphere caps -> high local texture nearby)
#   * SKIN       = area where spheres have merged and leave NO visible gap:
#                  smooth AND darker than the surrounding granular packing
#                  (low local texture)
#   * SPHERE CAP = the bright top of a primary bead (solid material)
#
# Pixels inside the particle body are classified from the curvature-flattened
# image (flat = img - gaussian(img, sigma)) and the local texture CV:
#   cap  : flat >= bright cut                              (bead tops, solid)
#   pore : flat <  dark cut  AND  texture >= tau           (open interstitial gap)
#   skin : (flat <  dark cut OR mid) AND texture <  tau     (smooth merged patch)
#
#   porosity = pore_area / body_area           <-- "% open space between spheres"
#   skin_frac = skin_area / body_area
#
# Saves an overlay (pore=red, skin=blue) per image for visual verification.
# Run:  python3 sem_porosity.py  ->  data/sem_porosity.csv + sem/overlay_*.png
# =============================================================================
import glob, csv, warnings, numpy as np
warnings.filterwarnings("ignore")
from PIL import Image
from scipy import ndimage as ndi
from skimage import filters, morphology

SAMPLES = {"114641": "up3_1", "114642": "up3_2", "114643": "up3_3", "114644": "up3_4"}
R0, R1, C0, C1 = 70, 660, 45, 920
SIGMA_FLAT = 18      # curvature-shading removal scale [px]
TEX_WIN    = 5       # window for the local-range (bead-scale) measure [px]
SMOOTH_PCT = 40      # local range below this percentile = smooth (candidate skin)
DARK_PCT   = 42      # flattened intensity below this percentile = open gap
SKIN_MINPX = 120     # min connected area for a skin patch (~ the circled size)
BRIGHT_K   = 0.40    # bright cut = mean + BRIGHT_K*std (sphere caps)

def analyze(path, save_overlay=None):
    img = np.asarray(Image.open(path).convert("L"), float)[R0:R1, C0:C1]
    # particle-body mask: drop the off-particle dark background/corners
    body = ndi.gaussian_filter(img, 12) > filters.threshold_otsu(img) * 0.55
    body = morphology.remove_small_objects(body, 3000)
    body = ndi.binary_fill_holes(body)
    if body.sum() < 0.25 * img.size:
        body = np.ones_like(img, bool)
    flat = img - ndi.gaussian_filter(img, SIGMA_FLAT)      # remove curvature
    # LOCAL RANGE (max-min over a bead-scale window): high where distinct sphere
    # caps sit next to dark gaps (granular/porous); low where beads have merged
    # into a smooth skin. This is the discriminator the reviewer described.
    rng = ndi.maximum_filter(img, TEX_WIN) - ndi.minimum_filter(img, TEX_WIN)
    fb = flat[body]
    med = np.median(fb)
    dark_cut   = np.percentile(fb, DARK_PCT)      # open-gap darkness cut
    bright_cut = fb.mean() + BRIGHT_K * fb.std()
    rng_cut    = np.percentile(rng[body], SMOOTH_PCT)   # below = smooth (merged)
    smooth = rng < rng_cut
    cap  = body & (flat >= bright_cut)
    # SKIN: smooth AND darker-than-median, broad patches (the merged, gap-free areas)
    skin_raw = body & smooth & (flat < med)
    skin = morphology.remove_small_objects(skin_raw, SKIN_MINPX)
    # GRANULAR (bead-packed) region: rough, non-smooth surface = where spheres and
    # their gaps are resolved. Excludes smooth skin AND smooth off-particle substrate.
    granular = body & (~smooth)
    # OPEN PORE: dark interstitial gaps within the granular packing
    pore = granular & (flat < dark_cut) & (~skin)
    A = body.sum(); G = max(granular.sum(), 1)
    res = dict(porosity_granular=pore.sum() / G,   # <-- % open space BETWEEN spheres
               porosity_surface=pore.sum() / A,    # (of whole framed surface)
               skin=skin.sum() / A, granular_frac=G / A,
               cap=cap.sum() / A, tex_mean=rng[body].mean())
    if save_overlay:
        rgb = np.dstack([img, img, img]).astype(np.uint8)
        rgb[pore] = [220, 40, 40]      # open pores -> red
        rgb[skin] = [40, 90, 230]      # skin -> blue
        Image.fromarray(rgb).save(save_overlay)
    return res

rows = []
for sid, cond in SAMPLES.items():
    imgs = sorted(glob.glob(f"sem/SEM_X{sid}_5kx0.jpg") +
                  glob.glob(f"sem/SEM_X{sid}_10kx0.jpg"))
    ov = f"sem/overlay_{cond}_5kx.png" if sid else None
    vals = []
    for i, p in enumerate(imgs):
        vals.append(analyze(p, save_overlay=ov if p.endswith("5kx0.jpg") else None))
    agg = {k: float(np.mean([v[k] for v in vals])) for k in vals[0]}
    rows.append(dict(cond=cond, sample=sid,
        porosity_pct=round(100 * agg["porosity_granular"], 1),   # % open space between spheres
        porosity_surface_pct=round(100 * agg["porosity_surface"], 1),
        skin_pct=round(100 * agg["skin"], 1),
        granular_frac_pct=round(100 * agg["granular_frac"], 1)))
    print(f"{cond} ({sid}): POROSITY (open space between spheres)={rows[-1]['porosity_pct']}%"
          f"   [surface {rows[-1]['porosity_surface_pct']}%, skin {rows[-1]['skin_pct']}%, "
          f"granular {rows[-1]['granular_frac_pct']}%]")

with open("data/sem_porosity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("\nWrote data/sem_porosity.csv + sem/overlay_*_5kx.png")
