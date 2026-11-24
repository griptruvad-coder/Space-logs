"""Rank stars by how unusual their light curves are (Isolation Forest + per-feature extremes) and plot the strangest."""
import json, sys, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import RobustScaler
S = int(sys.argv[1]); TOP = int(sys.argv[2])
df = pd.DataFrame([r for r in map(json.loads, open(f"feat_s{S}.jsonl")) if "error" not in r]).drop_duplicates("tic")
df = df[df.tmag <= 13.5]
F = ["lmad", "smooth", "skew", "kurt", "p001", "p999", "ndip", "dipmax", "dipcv", "nflare", "lspow", "lsamp", "sde", "blsdep", "rawrange"]
X = df[F].copy(); X["ndip"] = np.log1p(X.ndip); X["nflare"] = np.log1p(X.nflare); X["dipmax"] = np.log1p(X.dipmax)
X = X.replace([np.inf, -np.inf], np.nan).fillna(X.median())
Xs = np.clip(RobustScaler().fit_transform(X), -30, 30)
iso = IsolationForest(n_estimators=500, contamination="auto", random_state=0).fit(Xs)
df["score"] = -iso.score_samples(Xs)
# known-object flags
known = set()
import re, csv
for fn in ["toi.csv", "ctoi.csv"]:
    for r in csv.DictReader(open(fn)):
        try: known.add(int(r["TIC ID"]))
        except: pass
for fn in ["ps.csv", "k2.csv"]:
    for l in open(fn):
        m = re.search(r"(\d+)", l)
        if m: known.add(int(m.group(1)))
df["known"] = df.tic.isin(known)
df = df.sort_values("score", ascending=False)
pd.set_option("display.width", 250)
print(len(df), "stars")
print(df.head(TOP)[["tic", "tmag", "teff", "rad", "score", "skew", "ndip", "dipmax", "dipcv", "nflare", "lsper", "sde", "blsper", "known"]].round(2).to_string(index=False))
df.to_csv(f"anomaly_s{S}.csv", index=False)
top = df.head(TOP)
cols = 4; rows = int(np.ceil(TOP / cols))
fig, axs = plt.subplots(rows, cols, figsize=(5 * cols, 2.3 * rows)); axs = axs.ravel()
for a in axs: a.axis("off")
for a, (_, r) in zip(axs, top.iterrows()):
    t, f = np.load(f"bins{S}/{int(r.tic)}.npy"); a.axis("on")
    a.plot(t, f, ".", ms=1.5, color="k"); a.set_title(f"TIC {int(r.tic)} T={r.tmag:.1f} R={r.rad if r.rad==r.rad else 0:.1f}{' KNOWN' if r.known else ''}", fontsize=8)
    a.tick_params(labelsize=6)
plt.tight_layout(); plt.savefig(f"anomaly_top_s{S}.png", dpi=70)
