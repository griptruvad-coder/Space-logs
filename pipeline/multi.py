import sys, warnings, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); import lightkurve as lk
from scipy.ndimage import median_filter
tic = int(sys.argv[1]); thr = float(sys.argv[2])
sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC"); print(sr)
seen = set(); rows = []
for r in sr:
    s = r.mission[0]
    if s in seen: continue
    seen.add(s); lc = r.download(download_dir="vetdl").remove_nans().normalize(); rows.append((s, lc))
fig, axs = plt.subplots(len(rows), 1, figsize=(13, 2.4 * len(rows)))
for a, (s, lc) in zip(np.atleast_1d(axs), rows):
    t, f = lc.time.value, lc.flux.value; a.plot(t, f, ".", ms=0.8, color="k"); a.set_title(f"TIC {tic} {s}", fontsize=9)
    dt = np.median(np.diff(t)); rr = f / median_filter(f, size=int(1.5 / dt) | 1, mode="nearest")
    low = np.where(rr < thr)[0]
    for g in (np.split(low, np.where(np.diff(low) > 30)[0] + 1) if len(low) else []):
        if len(g) < 3: continue
        k = g[np.argmin(rr[g])]; print(s, f"eclipse BTJD {t[k]:.3f} depth {1-rr[k]:.3f} dur {(t[g[-1]]-t[g[0]])*24:.1f} h")
plt.tight_layout(); plt.savefig(f"multi_{tic}.png", dpi=70)
