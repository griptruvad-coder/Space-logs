# VSX new-variable submission packet — Anurag (prepared 24 Nov 2026)

Before submitting: search each position in VSX (https://www.aavso.org/vsx → Search → Coordinates, radius 60") and only submit if nothing is found.

Data source for both: NASA TESS 2-minute SPOC PDCSAP light curves (MAST). Magnitudes are in the TESS band (T), anchored on the TIC v8 Tmag. Amplitudes are probably diluted by nearby stars in TESS's 21" pixels, so the true amplitudes are likely larger. Gaia DR3 radial velocities independently show large velocity swings, which confirms each target star itself is a close binary.

---

## 1. TIC 168228752

| Field | Value |
|---|---|
| Name / cross-IDs | TIC 168228752; Gaia DR3 6769125881659800192 |
| RA, Dec (J2000, deg) | 289.982055, −26.161874 |
| Variability type | EW (W UMa-type contact/near-contact eclipsing binary) |
| Period | 0.372687 d |
| Epoch of primary minimum | BJD 2460817.56497 |
| Range (TESS band) | 11.527 – 11.571 T; secondary minimum 11.557 T |
| Data | TESS Sector 92 (2025), 2-min cadence |
| Supporting evidence | Gaia DR3 rv_amplitude_robust = 106.5 km/s over 11 epochs; not in Gaia DR3 vari_eclipsing_binary; not a TOI/CTOI; no ExoFOP history; not flagged by the SPOC pipeline as a TCE |
| Remarks | Unequal maxima (O'Connell effect) |

## 2. TIC 1159914

| Field | Value |
|---|---|
| Name / cross-IDs | TIC 1159914; Gaia DR3 6762282422513576576; TYC 6885-1397-1; EPIC 213550243 |
| RA, Dec (J2000, deg) | 286.621112, −29.035943 |
| Variability type | EW (W UMa-type) |
| Period | 0.313998 d |
| Epoch of primary minimum | BJD 2458680.78346 |
| Range (TESS band) | 10.770 – 10.790 T; secondary minimum 10.785 T |
| Data | TESS Sectors 13 (2019) and 92 (2025), 2-min cadence (6-year baseline) |
| Supporting evidence | Gaia DR3 rv_amplitude_robust = 47.8 km/s over 12 epochs; not in Gaia DR3 vari_eclipsing_binary; not a TOI/CTOI; not flagged by SPOC as a TCE |
| Remarks | Was a K2 target (EPIC 213550243): check whether any K2 variability catalogue already lists it |

---

Discoverer credit: Anurag (independent analysis of public TESS data with a custom Python/BLS pipeline).
