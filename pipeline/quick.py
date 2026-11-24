import sys, warnings, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore"); import lightkurve as lk
tics = [int(x) for x in sys.argv[1:]]
fig, axs = plt.subplots(len(tics), 2, figsize=(15, 2.6 * len(tics)), gridspec_kw=dict(width_ratios=[2, 1]))
for (a1, a2), tic in zip(axs, tics):
    sr = lk.search_lightcurve(f"TIC {tic}", author="SPOC", exptime=120)
    i = [k for k, m in enumerate(sr.mission) if m.endswith("92")]; lc = sr[i[0] if i else -1].download(download_dir="vetdl", quality_bitmask="default").remove_nans().normalize()
    t, f = lc.time.value, lc.flux.value
    a1.plot(t, f, ".", ms=0.8, color="k"); a1.set_title(f"TIC {tic} ({len(sr)} sectors avail: {','.join(s.split()[-1] for s in sr.mission)})", fontsize=9)
    k = np.argmin(f); m = np.abs(t - t[k]) < 1.2
    a2.plot(t[m], f[m], ".", ms=1.5, color="k"); a2.set_title("zoom around deepest point", fontsize=9)
plt.tight_layout(); plt.savefig("quick.png", dpi=70)
