"""Scan TESS SPOC 2-min light curves for periodic dips (BLS) on stars with no known planet/TOI/CTOI."""
import sys, os, re, random, json, warnings, subprocess
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
warnings.filterwarnings("ignore")

SECTOR = int(sys.argv[1]); N = int(sys.argv[2]); SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 0
os.makedirs("lc", exist_ok=True)

def known_tics():
    k = set()
    import csv
    for fn, col in [("toi.csv", "TIC ID"), ("ctoi.csv", "TIC ID")]:
        with open(fn) as f:
            for r in csv.DictReader(f):
                try: k.add(int(r[col]))
                except: pass
    for fn2 in ["ps.csv","k2.csv"]:
     for line in open(fn2):
        m = re.search(r"(\d+)", line)
        if m: k.add(int(m.group(1)))
    return k

def analyse(url):
    from astropy.io import fits
    from astropy.timeseries import BoxLeastSquares
    fn = "lc/" + url.split("product/")[1]
    tic = int(fn.split("-")[2])
    try:
        if not os.path.exists(fn):
            subprocess.run(["curl", "-s", "-L", "-o", fn, url], timeout=120, check=True)
        with fits.open(fn) as h:
            d = h[1].data; hdr = h[0].header
            t, f, q = d["TIME"], d["PDCSAP_FLUX"], d["QUALITY"]
            tmag = hdr.get("TESSMAG"); teff = hdr.get("TEFF"); rad = hdr.get("RADIUS")
        m = np.isfinite(t) & np.isfinite(f) & (q == 0)
        t, f = t[m], f[m] / np.nanmedian(f[m])
        if len(t) < 2000: return None
        # detrend: running median over 0.75 d
        from scipy.ndimage import median_filter
        dt = np.median(np.diff(t)); w = int(0.75 / dt) | 1
        trend = median_filter(f, size=w, mode="nearest")
        r = f / trend
        # sigma clip upward outliers only
        mad = 1.4826 * np.median(np.abs(r - 1)); keep = r < 1 + 4 * mad
        t, r = t[keep], r[keep]
        span = t.max() - t.min()
        periods = np.exp(np.linspace(np.log(0.4), np.log(span / 2), 6000))
        bls = BoxLeastSquares(t, r)
        res = bls.power(periods, [0.04, 0.08, 0.12, 0.2], objective="snr")
        p = res.power; sde = (p.max() - p.mean()) / p.std()
        i = np.argmax(p); P, T0, D, dep = res.period[i], res.transit_time[i], res.duration[i], res.depth[i]
        st = bls.compute_stats(P, D, T0)
        odd, even = st["depth_odd"], st["depth_even"]
        oe = abs(odd[0] - even[0]) / np.hypot(odd[1], even[1])
        sec = st["depth_half"]  # depth at phase 0.5 (secondary eclipse check)
        ntr = int(np.sum(st["per_transit_count"] > 0))
        intr = bls.transit_mask(t, P, D, T0); nin = int(intr.sum())
        snr = float(dep / mad * np.sqrt(max(nin, 1)))
        # single-event check: biggest single transit share of total SNR
        ev = [np.mean(1 - r[intr & (np.abs(((t - T0 + P/2) % P) - P/2) < D) & (np.round((t - T0)/P) == k)]) for k in np.unique(np.round((t[intr]-T0)/P))]
        ev = [e for e in ev if np.isfinite(e)]
        maxfrac = float(max(ev) / max(sum(ev), 1e-9)) if ev else 1.0
        os.remove(fn)
        return dict(tic=tic, tmag=tmag, teff=teff, rad=rad, P=float(P), T0=float(T0), dur_h=float(D * 24),
                    depth_ppm=float(dep * 1e6), snr=snr, maxfrac=maxfrac, nin=nin, sde=float(sde),
                    oe_sigma=float(oe), sec_ppm=float(sec[0] * 1e6), sec_snr=float(sec[0] / sec[1]) if sec[1] > 0 else 0,
                    ntr=ntr, noise_ppm=float(mad * 1e6))
    except Exception as e:
        if os.path.exists(fn): os.remove(fn)
        return dict(tic=tic, error=str(e)[:100])

if __name__ == "__main__":
    known = known_tics()
    urls = [l.split()[-1] for l in open(f"lc{SECTOR}.sh") if l.startswith("curl")]
    done = set()
    if os.path.exists(f"results_s{SECTOR}.jsonl"):
        done = {json.loads(l)["tic"] for l in open(f"results_s{SECTOR}.jsonl")}
    urls = [u for u in urls if int(u.split("-")[2]) not in known and int(u.split("-")[2]) not in done]
    random.Random(SEED).shuffle(urls); urls = urls[:N]
    print(f"known excluded; scanning {len(urls)} stars", flush=True)
    out = open(f"results_s{SECTOR}.jsonl", "a")
    with ProcessPoolExecutor(6) as ex:
        for k, fut in enumerate(as_completed([ex.submit(analyse, u) for u in urls])):
            r = fut.result()
            if r: out.write(json.dumps(r) + "\n"); out.flush()
            if k % 50 == 0: print(k, flush=True)
