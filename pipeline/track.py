import glob, sys, numpy as np, warnings, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
from astropy.io import fits; from astropy.wcs import WCS; from astropy.time import Time
pat = sys.argv[1]; t1, r1, d1, t2, r2, d2 = map(float, sys.argv[2:8]); out = sys.argv[8]
vr, vd = (r2 - r1) / (t2 - t1), (d2 - d1) / (t2 - t1)
fs = sorted(glob.glob(pat), key=lambda f: fits.getheader(f)["EXPSTART"])
cols = 6; rows = int(np.ceil(len(fs) / cols))
fig, ax = plt.subplots(rows, cols, figsize=(2.6 * cols, 2.9 * rows), facecolor="k"); ax = ax.ravel()
for a in ax: a.axis("off")
for a, f in zip(ax, fs):
    h = fits.open(f); p = h[0].header; t = (p["EXPSTART"] + p["EXPEND"]) / 2
    ra, de = r1 + vr * (t - t1), d1 + vd * (t - t1)
    x, y = WCS(h["SCI"].header).all_world2pix(ra, de, 0); x, y = int(round(float(x))), int(round(float(y))); d = h["SCI"].data; r = 25
    ok = r < x < d.shape[1] - r and r < y < d.shape[0] - r
    a.set_title(f"{p['FILTER']} {p['DATE-BEG'][11:19]}" + ("" if ok else " off"), color="w", fontsize=8)
    if not ok: continue
    c = d[y - r:y + r, x - r:x + r]; v = np.nanpercentile(c, [5, 99.7])
    a.imshow(c, origin="lower", cmap="magma", vmin=v[0], vmax=v[1]); a.plot(r, r, "c+", ms=16, mew=0.8)
    print(f.split("/")[-1][:30], p["FILTER"], f"{t:.5f}", "pred", f"{ra:.6f} {de:.6f}", "pix", x, y)
plt.tight_layout(); plt.savefig(out, dpi=80, facecolor="k")
