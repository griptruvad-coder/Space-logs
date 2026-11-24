import json, pandas as pd
dv = open("dv92.sh").read()
df = pd.DataFrame([r for r in map(json.loads, open("results_s92.jsonl")) if "error" not in r]).drop_duplicates("tic")
df["tce"] = df.tic.map(lambda t: f"{t:016d}" in dv)
c = df[(df.tmag <= 13.5) & (df.sde >= 7.5) & (df.snr >= 8) & (df.ntr >= 3) & (df.maxfrac < 0.5)
       & (df.depth_ppm < 30000) & (df.oe_sigma < 3) & (df.sec_snr < 3) & (df.dur_h < 0.2 * df.P * 24) & (df.P > 0.55)]
c = c[(c.rad.isna()) | (c.rad < 3)]
pd.set_option("display.width", 250)
print(len(df), "stars;", len(c), "planet-like;", (~c.tce).sum(), "of them NOT caught by NASA pipeline")
print(c.sort_values("sde", ascending=False)[["tic","tmag","teff","rad","P","dur_h","depth_ppm","snr","sde","ntr","maxfrac","tce"]].round(2).to_string(index=False))
