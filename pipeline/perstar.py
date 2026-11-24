import json, warnings, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); import lightkurve as lk
F = {d["tic"]: d for d in json.load(open("vsx_final.json"))}; C = json.load(open("coords.json"))
for tic, name in [(168228752, "2MASS J19195569-2609427"), (1159914, "2MASS J19062908-2902093")]:
    d = F[tic]; P, T0 = d["period_d"], d["epoch_min_bjd"] - 2457000
    sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    lc = lk.LightCurveCollection([r.download(download_dir="vetdl", quality_bitmask="hard").remove_nans().normalize() for r in sr]).stitch()
    t, f = lc.time.value, lc.flux.value; mag = -2.5 * np.log10(f) + C[str(tic)]["Tmag"]
    ph = ((t - T0) / P) % 1
    fig, ax = plt.subplots(2, 1, figsize=(10, 8))
    m = t < t.min() + 3
    ax[0].plot(t[m], mag[m], ".", ms=2, color="k"); ax[0].invert_yaxis(); ax[0].set_xlabel("BTJD (BJD - 2457000)"); ax[0].set_ylabel("TESS mag")
    ax[0].set_title(f"{name} = TIC {tic}: first 3 days of TESS 2-min data")
    for k in [0, 1]: ax[1].plot(ph + k, mag, ".", ms=1, color="0.4")
    ax[1].invert_yaxis(); ax[1].set_xlabel("Phase"); ax[1].set_ylabel("TESS mag")
    ax[1].set_title(f"Phased: P = {P} d, Epoch(min) = BJD {d['epoch_min_bjd']}  | EW, range {d['max_T']}-{d['min_T']} T")
    fig.text(0.01, 0.005, f"Data: TESS SPOC PDCSAP, sectors {','.join(d['sectors'])}. Mag zero-point from TIC v8 Tmag. Analysis: Anurag Kumar.", fontsize=8)
    plt.tight_layout(rect=(0, 0.02, 1, 1)); fn = f"VSX_{tic}_lightcurve.png"; plt.savefig(fn, dpi=85); print(fn)
