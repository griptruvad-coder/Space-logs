# Hunt 4: AI anomaly hunt, TESS Sector 92 (24 Nov 2026)

- **Method:** 15 shape features per light curve (scatter, smoothness, skew, kurtosis, dip/flare counts and depth irregularity, Lomb-Scargle and BLS strength), Isolation Forest on 6,250 stars with Tmag ≤ 13.5 (12,840 processed).
- **Highlights**
  - **TIC 111716343** (TYC 5819-17-1, Tmag 8.5): one 27% eclipse per sector in Sectors 42, 70 and 92. Combining them gives P = 21.7686 d with zero timing mismatch. The 1% secondary falls at phase 0.409, so the orbit is eccentric. Gaia DR3 RUWE = 19 hints at a third body or a resolved orbit. Not a TOI, not in the Gaia DR3 EB table. VSX check pending.
  - **TIC 97628340** (HIP 97448): 5.05 d eclipsing binary, Gaia RUWE 19.8, possible triple.
  - **TIC 10686217:** M-dwarf eclipsing binary (Teff ~3160 K), Gaia RV swing 137 km/s.
  - **Asteroid crossings:** TIC 130290692 and TIC 169388003 show a dip-bump-dip shape from an asteroid crossing the photometric and background apertures.
- **Next:** Gaia DR4 (2 Dec 2026) epoch astrometry for TIC 111716343 and TIC 97628340; run the anomaly hunt on more sectors.
