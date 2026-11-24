"""Extract anomaly-detection features + a compact binned light curve for every star in a TESS sector."""
import sys, os, json, random, subprocess, warnings
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
warnings.filterwarnings("ignore")
SECTOR, N, SEED = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
os.makedirs("lc", exist_ok=True); os.makedirs(f"bins{SECTOR}", exist_ok=True)

def feats(url):
    from astropy.io import fits
    from astropy.timeseries import LombScargle, BoxLeastSquares
    from scipy.ndimage import median_filter
    from scipy.stats import skew, kurtosis
    fn = "lc/" + url.split("product/")[1]; tic = int(fn.split("-")[2])
    try:
        if not os.path.exists(fn):
            subprocess.run(["curl", "-s", "-L", "-o", fn, url], timeout=120, check=True)
        with fits.open(fn) as h:
            d = h[1].data; hd = h[0].header
            t, f, q = d["TIME"], d["PDCSAP_FLUX"], d["QUALITY"]
            tmag, teff, rad = hd.get("TESSMAG"), hd.get("TEFF"), hd.get("RADIUS")
        os.remove(fn)
        m = np.isfinite(t) & np.isfinite(f) & (q == 0)
        t, f = t[m], f[m]
        if len(t) < 3000: return None
        f = f / np.nanmedian(f)
        dt = np.median(np.diff(t))
        # gentle detrend (3 d) keeps intrinsic shapes, removes slow systematics
        r = f / median_filter(f, size=int(3.0 / dt) | 1, mode="nearest")
        med = np.median(r); mad = 1.4826 * np.median(np.abs(r - med)) + 1e-9
        z = (r - med) / mad
        # point-to-point noise vs total scatter (smooth variability has small p2p)
        p2p = 1.4826 * np.median(np.abs(np.diff(r))) / np.sqrt(2) + 1e-9
        # dip events: consecutive runs below -5 sigma
        low = z < -5; ev = np.split(np.where(low)[0], np.where(np.diff(np.where(low)[0]) > 3)[0] + 1) if low.any() else []
        ev = [e for e in ev if len(e) >= 3]
        evd = [float(-z[e].min()) for e in ev]
        high = z > 5; fl = np.split(np.where(high)[0], np.where(np.diff(np.where(high)[0]) > 3)[0] + 1) if high.any() else []
        fl = [e for e in fl if len(e) >= 2]
        ls = LombScargle(t, r); fr, pw = ls.autopower(minimum_frequency=1 / 13, maximum_frequency=24, samples_per_peak=5)
        k = np.argmax(pw)
        # BLS
        bls = BoxLeastSquares(t, r); per = np.exp(np.linspace(np.log(0.3), np.log(13), 3000))
        res = bls.power(per, [0.04, 0.08, 0.16], objective="snr"); sde = (res.power.max() - res.power.mean()) / res.power.std()
        # binned LC for later viewing (30-min bins)
        b = np.arange(t.min(), t.max() + 0.0209, 0.0208); idx = np.digitize(t, b)
        cnt = np.bincount(idx, minlength=len(b) + 1)[1:len(b)]; s = np.bincount(idx, weights=f, minlength=len(b) + 1)[1:len(b)]
        np.save(f"bins{SECTOR}/{tic}.npy", np.vstack([b[:-1], np.where(cnt > 0, s / np.maximum(cnt, 1), np.nan)]).astype(np.float32))
        return dict(tic=tic, tmag=tmag, teff=teff, rad=rad,
                    lmad=float(np.log10(mad)), lp2p=float(np.log10(p2p)), smooth=float(np.log10(mad / p2p)),
                    skew=float(skew(z)), kurt=float(np.log10(kurtosis(z) + 3.5)),
                    p001=float(np.percentile(z, 0.1)), p999=float(np.percentile(z, 99.9)),
                    ndip=len(ev), dipmax=float(max(evd)) if evd else 0.0, dipcv=float(np.std(evd) / np.mean(evd)) if len(evd) > 1 else 0.0,
                    nflare=len(fl), lspow=float(pw[k]), lsper=float(1 / fr[k]), lsamp=float(np.log10(np.std(ls.model(t, fr[k])) + 1e-7)),
                    sde=float(sde), blsper=float(res.period[np.argmax(res.power)]), blsdep=float(res.depth[np.argmax(res.power)] / mad),
                    rawrange=float(np.log10(np.percentile(f, 99) - np.percentile(f, 1) + 1e-6)))
    except Exception as e:
        if os.path.exists(fn): os.remove(fn)
        return dict(tic=tic, error=str(e)[:80])

if __name__ == "__main__":
    urls = [l.split()[-1] for l in open(f"lc{SECTOR}.sh") if l.startswith("curl")]
    outf = f"feat_s{SECTOR}.jsonl"
    done = {json.loads(l)["tic"] for l in open(outf)} if os.path.exists(outf) else set()
    random.Random(SEED).shuffle(urls)
    urls = [u for u in urls if int(u.split("-")[2]) not in done][:N]
    print("scanning", len(urls), flush=True)
    with open(outf, "a") as out, ProcessPoolExecutor(6) as ex:
        for k, fu in enumerate(as_completed([ex.submit(feats, u) for u in urls])):
            r = fu.result()
            if r: out.write(json.dumps(r) + "\n"); out.flush()
            if k % 100 == 0: print(k, flush=True)
