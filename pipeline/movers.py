"""Find sources in individual JWST MIRI exposures that are NOT static on the sky, then link them into linear tracks."""
import glob, sys, json, warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources, SourceCatalog
from photutils.background import Background2D, MedianBackground
from astropy.convolution import Gaussian2DKernel, convolve
warnings.filterwarnings("ignore")

files = sorted(glob.glob(sys.argv[1]))
dets = []
for f in files:
    h = fits.open(f); p = h[0].header; sci = h["SCI"].data.astype(float); dq = h["DQ"].data
    bad = (dq & 1).astype(bool) | ~np.isfinite(sci)
    sci[bad] = np.nan
    if np.isfinite(sci).sum() < 1e5: continue
    bkg = Background2D(np.nan_to_num(sci, nan=np.nanmedian(sci)), 64, filter_size=3, bkg_estimator=MedianBackground(), mask=bad)
    d = sci - bkg.background; d[bad] = 0
    _, med, std = sigma_clipped_stats(d[~bad])
    dc = convolve(d, Gaussian2DKernel(1.5))
    seg = detect_sources(dc, 4 * std, npixels=10, mask=bad)
    if seg is None: continue
    cat = SourceCatalog(d, seg, convolved_data=dc)
    w = WCS(h["SCI"].header)
    t0 = Time(p["DATE-BEG"]).mjd; t1 = Time(p["DATE-END"]).mjd
    for s in cat:
        x, y = s.xcentroid, s.ycentroid
        if not np.isfinite(x): continue
        # reject sources touching bad pixels / edges
        yy, xx = int(round(y)), int(round(x))
        if bad[max(0, yy-6):yy+7, max(0, xx-6):xx+7].any(): continue
        ra, de = w.all_pix2world(x, y, 0)
        dets.append(dict(file=f.split("/")[-1], filt=p["FILTER"], tmid=(t0 + t1) / 2, texp=p["EFFEXPTM"], x=float(x), y=float(y),
                         ra=float(ra), dec=float(de), flux=float(s.segment_flux), snr=float(s.segment_flux / (std * np.sqrt(s.area.value))),
                         elong=float(s.elongation), area=float(s.area.value), pa=float(s.orientation.value)))
    print(f.split("/")[-1], p["FILTER"], len(cat), flush=True)

D = {k: np.array([d[k] for d in dets]) for k in dets[0]}
ra, de = np.radians(D["ra"]), np.radians(D["dec"])
# static test: how many OTHER exposures have a detection within 0.5"
from astropy.coordinates import SkyCoord
import astropy.units as u
c = SkyCoord(D["ra"] * u.deg, D["dec"] * u.deg)
idx1, idx2, sep, _ = c.search_around_sky(c, 0.5 * u.arcsec)
files_arr = D["file"]
nother = np.zeros(len(c), int)
for i in range(len(c)):
    m = (idx1 == i) & (files_arr[idx2] != files_arr[i])
    nother[i] = len(set(files_arr[idx2[m]]))
# how many exposures cover this position at all (same filter or any)
nexp = len(set(files_arr))
# reject detector-fixed artifacts: same pixel (within 2px) in another exposure
pixrep = np.zeros(len(c), bool)
for i in range(len(c)):
    m = (files_arr != files_arr[i]) & (np.hypot(D["x"] - D["x"][i], D["y"] - D["y"][i]) < 2.5)
    pixrep[i] = m.any()
lonely = (nother == 0) & (D["snr"] > 6) & ~pixrep
print(f"\n{len(c)} detections in {nexp} exposures; {lonely.sum()} appear in only one exposure")
out = [{k: (D[k][i].item() if hasattr(D[k][i], 'item') else D[k][i]) for k in D} for i in np.where(lonely)[0]]
json.dump(out, open(sys.argv[2], "w"))
# try to link lonely detections into linear tracks (>=3 points, consistent rate 0.5-300 "/hr)
L = [o for o in out]
L.sort(key=lambda o: o["tmid"])
tracks = []
for i in range(len(L)):
    for j in range(i + 1, len(L)):
        dt = (L[j]["tmid"] - L[i]["tmid"]) * 24
        if dt <= 0.003: continue
        dra = (L[j]["ra"] - L[i]["ra"]) * np.cos(np.radians(L[i]["dec"])) * 3600; dde = (L[j]["dec"] - L[i]["dec"]) * 3600
        rate = np.hypot(dra, dde) / dt
        if not (0.5 < rate < 300): continue
        members = [i, j]
        for k in range(len(L)):
            if k in members: continue
            dtk = (L[k]["tmid"] - L[i]["tmid"]) * 24
            pra = L[i]["ra"] + dra / dt * dtk / 3600 / np.cos(np.radians(L[i]["dec"])); pde = L[i]["dec"] + dde / dt * dtk / 3600
            off = np.hypot((L[k]["ra"] - pra) * np.cos(np.radians(L[i]["dec"])), L[k]["dec"] - pde) * 3600
            if off < 1.0 and L[k]["file"] not in [L[m]["file"] for m in members]: members.append(k)
        if len(members) >= 3:
            key = tuple(sorted(members))
            if key not in [t[0] for t in tracks]: tracks.append((key, rate, np.degrees(np.arctan2(dra, dde)) % 360))
tracks.sort(key=lambda t: -len(t[0]))
seen = set(); final = []
for key, rate, pa in tracks:
    if seen & set(key): continue
    seen |= set(key); final.append((key, rate, pa))
for key, rate, pa in final:
    print(f"TRACK n={len(key)} rate={rate:.1f}\"/hr PA={pa:.0f}deg")
    for k in key: print("   ", L[k]["file"], L[k]["filt"], f"{L[k]['tmid']:.5f}", f"{L[k]['ra']:.6f} {L[k]['dec']:.6f}", f"snr={L[k]['snr']:.1f} elong={L[k]['elong']:.2f}")
json.dump([[L[k] for k in key] for key, _, _ in final], open(sys.argv[2].replace(".json", "_tracks.json"), "w"))
