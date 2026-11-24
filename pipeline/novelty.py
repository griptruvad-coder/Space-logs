"""For each TIC: is this signal already known anywhere? ExoFOP history, SPOC TCE, Gaia DR3 variability/EB/NSS."""
import sys, json, subprocess, warnings
warnings.filterwarnings("ignore")
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia

dv = open("dv92.sh").read()
for tic in map(int, sys.argv[1:]):
    out = {"tic": tic, "spoc_tce": f"{tic:016d}" in dv}
    try:
        j = json.loads(subprocess.run(["curl", "-s", f"https://exofop.ipac.caltech.edu/tess/target.php?id={tic}&json"],
                                      capture_output=True, text=True, timeout=60).stdout)
        out["names"] = j["basic_info"]["star_names"]
        out["exofop"] = {k: len(j.get(k) or []) for k in ["tois", "ctois", "imaging", "spectroscopy", "time_series", "files"]}
        out["k2_campaign"] = j["basic_info"].get("k2_campaign")
    except Exception as e:
        out["exofop_err"] = str(e)[:80]
    try:
        g = Catalogs.query_criteria(catalog="TIC", ID=tic)[0]["GAIA"]
        out["gaia"] = str(g)
        if g and str(g) != "--":
            r = Gaia.launch_job(f"""select g.ruwe, g.non_single_star, g.rv_amplitude_robust, g.rv_nb_transits, g.phot_variable_flag,
                v.best_class_name from gaiadr3.gaia_source g left join gaiadr3.vari_classifier_result v on v.source_id=g.source_id
                where g.source_id={g}""").get_results()
            out.update({c: (None if str(r[c][0]) == "--" else str(r[c][0])) for c in r.colnames})
            out["gaia_eb"] = len(Gaia.launch_job(f"select source_id from gaiadr3.vari_eclipsing_binary where source_id={g}").get_results())
    except Exception as e:
        out["gaia_err"] = str(e)[:80]
    print(json.dumps(out))
