import warnings, numpy as np; warnings.filterwarnings("ignore"); import lightkurve as lk
lc = lk.search_lightcurve("TIC 97628340", author="SPOC", exptime=120)[-1].download(download_dir="vetdl").remove_nans().normalize()
t, f = lc.time.value, lc.flux.value
from scipy.ndimage import median_filter
r = f / median_filter(f, size=1081, mode="nearest")
low = np.where(r < 0.975)[0]; grp = np.split(low, np.where(np.diff(low) > 30)[0] + 1)
for g in grp:
    w = (t > t[g].min() - 0.15) & (t < t[g].max() + 0.15)
    tt, rr = t[w], r[w]; m = rr < np.percentile(rr, 30)
    c = np.polyfit(tt[m] - tt[m].mean(), rr[m], 2); tmin = -c[1] / (2 * c[0]) + tt[m].mean()
    print(f"eclipse at BTJD {tmin:.4f}  depth {1-rr.min():.3f}  dur {(t[g].max()-t[g].min())*24:.2f} h")
