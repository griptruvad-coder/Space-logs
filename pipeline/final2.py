import json, warnings, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); import lightkurve as lk
from astropy.timeseries import LombScargle
C = json.load(open("coords.json")); out = []
fig, axs = plt.subplots(2, 1, figsize=(9, 6.4), facecolor="#070b1a")
for ax, (tic, Pg) in zip(axs, [(168228752, 0.745457 / 2), (1159914, 0.628091 / 2)]):
    sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    lc = lk.LightCurveCollection([r.download(download_dir="vetdl", quality_bitmask="hard").remove_nans().normalize() for r in sr]).stitch()
    t, f = lc.time.value, lc.flux.value
    fr = np.linspace(2 / Pg * 0.995, 2 / Pg * 1.005, 200000)  # double-wave: LS peak at 2/P
    pw = LombScargle(t, f).power(fr); P = 2 / fr[np.argmax(pw)]
    # epoch of primary (deeper) minimum
    ph = (t / P) % 1; b = np.linspace(0, 1, 101); idx = np.digitize(ph, b)
    bf = np.array([np.nanmedian(f[idx == k]) for k in range(1, 101)]); kmin = np.nanargmin(bf)
    T0 = t.min() - (t.min() / P % 1 - (b[kmin] + 0.005)) * P
    T0 = T0 + np.round((np.median(t) - T0) / P) * P
    mag = -2.5 * np.log10(f) + C[str(tic)]["Tmag"]; ph = ((t - T0) / P + 0.25) % 1 - 0.25
    bm = np.array([np.nanmedian(mag[(ph >= a) & (ph < a + 0.01)]) for a in np.arange(-0.25, 0.75, 0.01)])
    d = dict(tic=tic, gaia_dr3=C[str(tic)]["gaia"], ra=round(C[str(tic)]["ra"], 6), dec=round(C[str(tic)]["dec"], 6), type="EW",
             period_d=round(float(P), 6), epoch_min_bjd=round(float(T0) + 2457000, 5),
             max_T=round(float(np.nanmin(bm)), 3), min_T=round(float(np.nanmax(bm)), 3),
             sec_min_T=round(float(np.nanmax(bm[(np.arange(-0.25, 0.75, 0.01) > 0.35) & (np.arange(-0.25, 0.75, 0.01) < 0.65)])), 3),
             sectors=[s.split()[-1] for s in sr.mission])
    out.append(d); print(json.dumps(d))
    ax.set_facecolor("#070b1a"); ax.plot(ph, mag, ".", ms=1, color="#5b6aa8", alpha=.5)
    ax.plot(np.arange(-0.25, 0.75, 0.01) + 0.005, bm, "-", color="#ffb347", lw=2.5); ax.invert_yaxis()
    ax.set_title(f"TIC {tic} · W UMa-type binary · P = {P:.5f} d", color="w", fontsize=11)
    ax.tick_params(colors="#c9d2f0"); ax.set_ylabel("TESS mag", color="#c9d2f0"); [s.set_color("#2a3360") for s in ax.spines.values()]
axs[-1].set_xlabel("Orbital phase", color="#c9d2f0")
plt.tight_layout(); plt.savefig("vsx_final.png", dpi=100, facecolor="#070b1a"); json.dump(out, open("vsx_final.json", "w"), indent=1)
