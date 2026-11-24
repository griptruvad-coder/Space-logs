"""Download all SPOC 2-min sectors for a TIC, detrend, fold, and plot: full LC, folded, odd/even, secondary."""
import sys, glob, os, subprocess, warnings
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.ndimage import median_filter
from astropy.timeseries import BoxLeastSquares
warnings.filterwarnings("ignore")
tic = int(sys.argv[1]); P = float(sys.argv[2]) if len(sys.argv) > 2 else None
import lightkurve as lk
sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
print(sr)
T, R = [], []
for row in sr:
    lc = row.download(download_dir="vetdl", flux_column="pdcsap_flux", quality_bitmask="hard")
    t = lc.time.value; f = lc.flux.value; m = np.isfinite(t) & np.isfinite(f); t, f = t[m], f[m] / np.nanmedian(f[m])
    w = int(0.75 / np.median(np.diff(t))) | 1
    r = f / median_filter(f, size=w, mode="nearest"); T.append(t); R.append(r)
t, r = np.concatenate(T), np.concatenate(R)
mad = 1.4826 * np.median(np.abs(r - 1)); k = r < 1 + 4 * mad; t, r = t[k], r[k]
bls = BoxLeastSquares(t, r)
if P is None:
    periods = np.exp(np.linspace(np.log(0.4), np.log(min(30, (t.max() - t.min()) / 2)), 40000))
else:
    periods = np.linspace(P * 0.995, P * 1.005, 4000)
res = bls.power(periods, [0.03, 0.05, 0.08, 0.12, 0.16, 0.2, 0.25])
i = np.argmax(res.power); P, T0, D, dep = res.period[i], res.transit_time[i], res.duration[i], res.depth[i]
st = bls.compute_stats(P, D, T0)
print(f"P={P:.5f} T0={T0:.4f} dur={D*24:.2f}h depth={dep*1e6:.0f}ppm odd={st['depth_odd'][0]*1e6:.0f}±{st['depth_odd'][1]*1e6:.0f} even={st['depth_even'][0]*1e6:.0f}±{st['depth_even'][1]*1e6:.0f} sec={st['depth_half'][0]*1e6:.0f}±{st['depth_half'][1]*1e6:.0f} ntr={np.sum(st['per_transit_count']>0)}")
ph = ((t - T0 + 0.5 * P) % P) / P - 0.5
def binned(x, y, n=200):
    b = np.linspace(x.min(), x.max(), n + 1); c = np.digitize(x, b)
    return 0.5 * (b[1:] + b[:-1]), np.array([np.median(y[c == j]) if np.any(c == j) else np.nan for j in range(1, n + 1)])
fig, ax = plt.subplots(2, 2, figsize=(13, 8))
ax[0, 0].plot(t, r, ",k"); ax[0, 0].set_title(f"TIC {tic} detrended ({len(sr)} sectors)")
for tt in st["transit_times"]: ax[0, 0].axvline(tt, color="r", alpha=.15)
z = np.abs(ph * P * 24) < 4 * D * 24
ax[0, 1].plot(ph[z] * P * 24, r[z], ".", ms=1, color="0.6"); bx, by = binned(ph[z] * P * 24, r[z], 80); ax[0, 1].plot(bx, by, "r-")
ax[0, 1].set_title(f"Folded P={P:.4f} d, depth {dep*1e6:.0f} ppm, dur {D*24:.2f} h"); ax[0, 1].set_xlabel("hours from mid-transit")
n = np.round((t - T0) / P).astype(int)
for par, col in [(0, "b"), (1, "orange")]:
    s = z & (n % 2 == par); bx, by = binned(ph[s] * P * 24, r[s], 50); ax[1, 0].plot(bx, by, "-", color=col, label="even" if par == 0 else "odd")
ax[1, 0].legend(); ax[1, 0].set_title("Odd vs even transits")
ph2 = ((t - T0) % P) / P - 0.5; z2 = np.abs(ph2 * P * 24) < 4 * D * 24
bx, by = binned(ph2[z2] * P * 24, r[z2], 60); ax[1, 1].plot(ph2[z2] * P * 24, r[z2], ".", ms=1, color="0.6"); ax[1, 1].plot(bx, by, "r-")
ax[1, 1].set_title("Around phase 0.5 (secondary eclipse?)")
plt.tight_layout(); plt.savefig(f"vet_{tic}.png", dpi=90); print("saved")
