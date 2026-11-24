"""Prepare VSX-style parameters for eclipsing-binary candidates: true period, epoch of primary minimum, TESS-band range, phase plot."""
import json, sys, warnings, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
import lightkurve as lk
from astropy.timeseries import BoxLeastSquares
C = json.load(open("coords.json"))
res = []
tics = [int(t) for t in sys.argv[1:]]
fig, axs = plt.subplots(len(tics), 1, figsize=(9, 3.1 * len(tics)))
for ax, tic in zip(np.atleast_1d(axs), tics):
    sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    lcs = [r.download(download_dir="vetdl", quality_bitmask="hard").remove_nans().normalize() for r in sr]
    lc = lk.LightCurveCollection(lcs).stitch().remove_outliers(sigma_upper=4, sigma_lower=50)
    t, f = lc.time.value, lc.flux.value
    bls = BoxLeastSquares(t, f)
    pg = bls.power(np.exp(np.linspace(np.log(0.3), np.log(5), 60000)), [0.02, 0.04, 0.06, 0.1])
    P1 = pg.period[np.argmax(pg.power)]
    # test P1 vs 2*P1: compare depths of alternating eclipses
    best = None
    for P in [P1, 2 * P1]:
        r = bls.power(np.linspace(P * 0.999, P * 1.001, 3000), [0.02, 0.04, 0.06, 0.1])
        i = np.argmax(r.power); Pf, T0 = r.period[i], r.transit_time[i]
        ph = ((t - T0) / Pf + 0.25) % 1 - 0.25
        def mn(lo, hi):
            m = (ph > lo) & (ph < hi); return np.nanpercentile(f[m], 5) if m.sum() > 20 else np.nan
        pri, sec = mn(-0.02, 0.02), mn(0.48, 0.52)
        best = (Pf, T0, pri, sec) if P != P1 or best is None else best
        if P == 2 * P1: P2 = (Pf, T0, pri, sec)
    # choose 2P if secondary at 2P is a real eclipse (deeper than 0.25 of primary) -> equal-ish eclipses = W UMa / EB at 2P
    Pf, T0, pri, sec = P2
    mag = -2.5 * np.log10(np.clip(f, 1e-3, None)) + C[str(tic)]["Tmag"]
    ph = ((t - T0) / Pf + 0.25) % 1 - 0.25
    bins = np.linspace(-0.25, 0.75, 201); idx = np.digitize(ph, bins)
    bm = np.array([np.nanmedian(mag[idx == k]) for k in range(1, 201)])
    maxm, minm = np.nanmin(bm), np.nanmax(bm)
    sec_depth = np.nanmax(bm[(bins[:-1] > 0.4) & (bins[:-1] < 0.6)]) - maxm
    amp = minm - maxm
    # light-curve shape: continuous variation (EW) vs flat out-of-eclipse (EA)
    out = bm[((bins[:-1] > 0.12) & (bins[:-1] < 0.38)) | ((bins[:-1] > 0.62) & (bins[:-1] < 0.73))]
    flat = np.nanstd(out) < 0.1 * amp
    vtype = "EA" if flat else ("EW" if sec_depth > 0.6 * amp else "EB")
    d = dict(tic=tic, gaia=C[str(tic)]["gaia"], ra=C[str(tic)]["ra"], dec=C[str(tic)]["dec"], period=round(float(Pf), 6),
             epoch_bjd=round(float(T0) + 2457000, 5), max_T=round(float(maxm), 3), min_T=round(float(minm), 3),
             sec_min_T=round(float(maxm + sec_depth), 3), type=vtype, sectors=[str(x) for x in sr.mission])
    res.append(d); print(json.dumps(d))
    ax.plot(ph, mag, ".", ms=1, color="0.7"); ax.plot((bins[1:] + bins[:-1]) / 2, bm, "r-")
    ax.invert_yaxis(); ax.set_title(f"TIC {tic}  P={Pf:.5f} d  type {vtype}  TESS mag {maxm:.2f}-{minm:.2f}", fontsize=10)
    ax.set_xlabel("phase"); ax.set_ylabel("TESS mag")
plt.tight_layout(); plt.savefig("vsx_candidates.png", dpi=90)
json.dump(res, open("vsx_candidates.json", "w"), indent=1)
