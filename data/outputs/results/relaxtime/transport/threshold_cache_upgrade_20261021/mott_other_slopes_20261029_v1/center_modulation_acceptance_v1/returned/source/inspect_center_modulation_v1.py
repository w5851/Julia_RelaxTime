import csv, json, hashlib, zipfile
from pathlib import Path
d=Path(r"D:/CodexCompute/trmott-other-slopes-20261029-v1/center_modulation_acceptance_v1")
def read(name):
    with (d/(name+".csv")).open(encoding="utf-8-sig",newline="") as f: return list(csv.DictReader(f))
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
centers=read("centered_modulation")
proc=read("process_modulation")
final=read("final_species_modulation")
precision=read("centered_precision")
chosen=["tau_ubar","tau_sbar","eta_over_s","zeta_over_s","sigma_over_T"]
v={"definition":"J=right log-secant minus left log-secant, at identical half-width h; final groups sum exact species modulation contributions","half_width":.0002,"window_table":[],"K_process_table":[],"final_groups":[],"sigma_scales":[],"pi_scales":[]}
for cid in sorted({r["curve_id"] for r in centers}):
    sel=[r for r in centers if r["curve_id"]==cid and r["profile"]=="tighter" and float(r["half_width"])==.0002]
    row={k:sel[0][k] for k in ["curve_id","meson","window","T_MeV","muB_MeV","xi_mott"]}
    row.update({r["observable"]:float(r["right_minus_left"]) for r in sel})
    row["precision_max_abs_J_change"]=max(float(r["absolute_difference"]) for r in precision if r["curve_id"]==cid and float(r["half_width"])==.0002)
    v["window_table"].append(row)
for r in proc:
    if r["profile"]=="tighter" and float(r["half_width"])==.0002 and r["species"] in ("sbar","ubar") and r["process_group"]=="K_s_processes":
        v["K_process_table"].append({k:r[k] for k in ("curve_id","species","delta_log_slope_from_rate","delta_log_slope_from_density","delta_log_slope_contribution")})
for cid in sorted({r["curve_id"] for r in centers}):
    for ob in chosen[2:]:
        sel=[r for r in final if r["profile"]=="tighter" and r["curve_id"]==cid and r["observable"]==ob and float(r["half_width"])==.0002]
        out=dict(curve_id=cid,observable=ob)
        for group,species in (("light",("u","d","ubar","dbar")),("strange",("s","sbar"))):
            for field in ("delta_log_slope_from_tau","delta_log_slope_from_coefficient","delta_log_slope_contribution"):
                out[group+"_"+field]=sum(float(r[field]) for r in sel if r["species"] in species)
        v["final_groups"].append(out)
for r in centers:
    if r["profile"]=="tighter" and r["observable"] in chosen:
        out={k:r[k] for k in ("curve_id","half_width","observable","left_log_slope","right_log_slope","right_minus_left")}
        if r["curve_id"]=="other07": v["sigma_scales"].append(out)
        if r["curve_id"] in ("other05","other08"): v["pi_scales"].append(out)
v["input_sha256"]={name:sha(d/name) for name in ("centered_modulation.csv","process_modulation.csv","final_species_modulation.csv","centered_precision.csv")}
(d/"interpretation_views.json").write_text(json.dumps(v,indent=2)+"\n",encoding="utf-8")
for r in v["window_table"]:
    print(json.dumps(r))
print(json.dumps({"K_process_examples":[r for r in v["K_process_table"] if r["curve_id"] in ("other01","other03")],"final_examples":[r for r in v["final_groups"] if r["curve_id"] in ("other01","other03","other07")]}))
files=[dict(path=p.relative_to(d).as_posix(),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(d.rglob("*")) if p.is_file()]
(d/"package_manifest.json").write_text(json.dumps(dict(schema="derived_diagnostic_package_v1",files=files),indent=2)+"\n",encoding="utf-8")
pack=Path(r"D:/CodexCompute/trmott-other-slopes-20261029-v1/packages/center_modulation_acceptance_v1")
pack.mkdir(parents=True,exist_ok=False)
z=pack/"center_modulation_acceptance_v1.zip"
with zipfile.ZipFile(z,"w",zipfile.ZIP_DEFLATED) as f:
    for p in sorted(d.rglob("*")):
        if p.is_file(): f.write(p,p.relative_to(d).as_posix())
receipt=dict(archive=z.name,bytes=z.stat().st_size,sha256=sha(z),manifest_files=len(files),package_manifest_sha256=sha(d/"package_manifest.json"))
(pack/"archive_receipt.json").write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
print(json.dumps(receipt))
