# Space Logs

Hunting for things nobody has noticed yet in public space-telescope data, using Python and AI. Every hunt, false alarm and finding is logged here.

**Blog:** https://griptruvad-coder.github.io/Space-logs/

By **Anurag Kumar** (Delhi, India).

## Findings so far

| Object | What it is | Status |
|---|---|---|
| 2MASS J19195569-2609427 (TIC 168228752) | New W UMa-type eclipsing binary, P = 0.372687 d | Submitted to AAVSO VSX, 24 Nov 2026 (pending) |
| 2MASS J19062908-2902093 (TIC 1159914) | New W UMa-type eclipsing binary, P = 0.313998 d | Not in VSX; submission queued |
| TIC 111716343 (TYC 5819-17-1) | Long-period eccentric eclipsing binary, P = 21.7686 d, secondary at phase 0.41, Gaia RUWE 19 | Candidate; VSX check pending |
| TIC 377241703 (EPIC 218711655) | Recovered the K2 hot-Jupiter candidate EPIC 218711655.01 (Petigura et al. 2018) independently in TESS | Known candidate, not new |

## Hunts

| Log | Data | Stars | Result |
|---|---|---|---|
| [Hunt 1](logs/hunt-01.md) | TESS Sector 92, BLS transit search | 2,823 | Re-found a known K2 candidate |
| [Hunt 2](logs/hunt-02.md) | TESS Sector 92, full sector, stricter filters | 12,666 | No new planets; new binaries |
| [Hunt 3](logs/hunt-03.md) | JWST MIRI archive, moving-object search | 49 exposures | Only detector artifacts so far |
| [Hunt 4](logs/hunt-04.md) | TESS Sector 92, AI anomaly detection (Isolation Forest) | 12,840 | Long-period eccentric binary, possible triples, asteroid crossings |

## Pipeline

- `pipeline/hunt.py`: downloads TESS SPOC 2-min light curves for a sector, removes known planets/TOIs/CTOIs/K2 candidates, runs Box Least Squares.
- `pipeline/features.py` + `pipeline/anomaly.py`: shape features for every light curve, Isolation Forest ranking, gallery of the strangest stars.
- `pipeline/vet.py`, `pipeline/multi.py`, `pipeline/novelty.py`: vetting (odd/even, secondary eclipse, multi-sector timing), plus "is it already known?" checks against ExoFOP, SPOC TCEs and Gaia DR3.
- `pipeline/movers.py`: finds sources that move between individual JWST exposures.

Data: NASA TESS (MAST), ESA Gaia DR3, NASA Exoplanet Archive, ExoFOP, JWST (MAST).
