# Hunt 3: JWST archive moving-object search (24 Nov 2026)

- **Idea:** moving objects (asteroids) get rejected when JWST exposures are combined, so look at individual `_cal` exposures instead.
- **Data:** 7,111 public MIRI imaging observations; 702 within 5° of the ecliptic; 49 exposures from 3 fields downloaded (IR14348-NE, vMa-2, HIGHZODI-DOY117).
- **Method:** source detection per exposure, then keep sources seen in only one exposure and not repeated at the same detector pixel, then link into linear tracks.
- **Result:** 19 single-exposure detections, and all were detector artifacts (a glowing detector row; the coronagraph-edge glow).
- **Next:** mask non-imager regions, scale to 100+ fields, train an artifact-vs-asteroid classifier, check survivors against the Minor Planet Center.
