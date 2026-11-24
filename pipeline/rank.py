import json, sys
import pandas as pd
rows = [json.loads(l) for l in open(sys.argv[1])]
df = pd.DataFrame([r for r in rows if "error" not in r])
print("analysed", len(df), "errors", sum("error" in r for r in rows))
df = df[(df.tmag <= 13.5) & (df.sde >= 9) & (df.snr >= 10) & (df.ntr >= 2) & (df.maxfrac < 0.6)]
df["kind"] = "planet-like"
df.loc[(df.depth_ppm > 30000) | (df.oe_sigma > 3) | (df.sec_snr > 4), "kind"] = "binary-like"
df.loc[(df.dur_h >= 0.3 * df.P * 24), "kind"] = "variable (not dips)"
df = df.sort_values("sde", ascending=False)
pd.set_option("display.width", 250)
print(df[["tic","tmag","teff","rad","P","dur_h","depth_ppm","snr","sde","oe_sigma","sec_snr","ntr","maxfrac","kind"]].round(2).to_string(index=False))
df.to_csv(sys.argv[1].replace(".jsonl", "_cands.csv"), index=False)
