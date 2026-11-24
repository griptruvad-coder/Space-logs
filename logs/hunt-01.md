# Hunt 1: First transit search (24 Nov 2026)

- **Data:** TESS Sector 92, 2-min SPOC PDCSAP light curves, 2,823 random stars with no confirmed planet, TOI or CTOI.
- **Method:** running-median detrend, Box Least Squares over 0.4 d to half the baseline; kept SDE ≥ 9, SNR ≥ 10, ≥ 2 transits.
- **Sanity check:** WASP-18 b recovered at P = 0.9415 d.
- **Result:** 31 strong periodic signals, mostly eclipsing binaries and pulsators. The best "planet-like" survivor, TIC 377241703 (P = 1.1577 d, ~1,700 ppm, grazing), passed odd/even, secondary and centroid checks, and Gaia RVs looked flat.
- **Lesson:** it turned out to be EPIC 218711655.01, a K2 candidate already published by Petigura et al. (2018). My exclusion list did not include K2 candidates. Fixed in Hunt 2.
- **Still useful:** an independent TESS recovery 9 years after K2, with a matching period (1.15740 d vs 1.157563 d).
