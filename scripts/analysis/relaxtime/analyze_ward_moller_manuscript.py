"""Audit the revised manuscript observables on the complete frozen transport grid."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

SPECIES = ("u", "d", "s", "ubar", "dbar", "sbar")
OBSERVABLES = tuple("tau_" + item for item in SPECIES) + (
    "eta_over_s", "zeta_over_s", "sigma_over_T")


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def write(path, rows):
    if not rows:
        raise ValueError(f"empty result table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def number(row, key):
    return float(row[key])


def difference(left, right, field):
    return 100 * (number(right, field) / number(left, field) - 1)


def audit(results):
    grids = {mode: read(results / f"{mode}_grid.csv") for mode in ("mode_a", "mode_b")}
    assert len(grids["mode_a"]) == 942 and len(grids["mode_b"]) == 925
    all_rows = [row for rows in grids.values() for row in rows]
    ids = {row["point_id"]: row for row in all_rows}
    assert len(ids) == 1867
    assert all(row["a_prescription"] == "ward_matched" and
               row["relative_velocity"] == "medium_moller" and
               row["b0_direction"] == "sphere_average" for row in all_rows)
    assert all(math.isfinite(number(row, field)) and number(row, field) > 0
               for row in all_rows for field in OBSERVABLES)
    sigma_closure = max(abs(sum(number(row, "tau_"+s)*number(row, "J_"+s) for s in SPECIES)
                            /number(row, "sigma_over_T")-1) for row in all_rows)
    assert sigma_closure < 1e-10
    rates=defaultdict(dict)
    for row in read(results/'channel_rates.csv'):
        rates[row['point_id']][row['channel']]=number(row,'rate')
    inverse_error=0.
    for row in all_rows:
        w=rates[row['point_id']]
        assert len(w)==17
        u,s,ub,sb=(number(row,'n_'+x) for x in ('u','s','ubar','sbar'))
        expected={
            'u':u*(w['uu_to_uu']+w['ud_to_ud'])+
                ub*(w['uubar_to_uubar']+w['uubar_to_ddbar']+w['uubar_to_ssbar']+w['udbar_to_udbar'])+
                s*w['us_to_us']+sb*w['usbar_to_usbar'],
            's':2*u*w['us_to_us']+2*ub*w['usbar_to_usbar']+s*w['ss_to_ss']+
                sb*(w['ssbar_to_ssbar']+2*w['ssbar_to_uubar']),
            'ubar':u*(w['uubar_to_uubar']+w['uubar_to_ddbar']+w['uubar_to_ssbar']+w['dubar_to_dubar'])+
                ub*(w['ubardbar_to_ubardbar']+w['ubarubar_to_ubarubar'])+s*w['subar_to_subar']+sb*w['ubarsbar_to_ubarsbar'],
            'sbar':2*u*w['usbar_to_usbar']+2*ub*w['ubarsbar_to_ubarsbar']+sb*w['sbarsbar_to_sbarsbar']+
                s*(w['ssbar_to_ssbar']+2*w['ssbar_to_uubar'])}
        inverse_error=max(inverse_error,*(abs(v/number(row,'tauinv_'+sp)-1) for sp,v in expected.items()))
    assert inverse_error<1e-10

    groups = {}
    for mode, rows in grids.items():
        grouped = defaultdict(list)
        for row in rows:
            grouped[(row["plot_panel"], row["plot_series"])].append(row)
        for key, subset in grouped.items():
            subset.sort(key=lambda row: number(row, "xi"))
            assert len({number(row, "xi") for row in subset}) == len(subset)
            groups[(mode, *key)] = subset
    endpoints = []
    references = []
    for (mode, panel, series), subset in groups.items():
        left, right = subset[0], subset[-1]
        assert number(left, "xi") == -.5 and number(right, "xi") == .5
        for field in OBSERVABLES:
            minimum, maximum = min(subset, key=lambda row:number(row, field)), max(subset, key=lambda row:number(row, field))
            endpoints.append(dict(mode=mode, panel=panel, series=series, observable=field,
                T_MeV=number(left, "T_MeV"), muB_MeV=number(left, "muB_MeV"),
                alpha_T=number(left, "alpha_T"), left=number(left, field), right=number(right, field),
                net_change_percent=difference(left,right,field),
                minimum=number(minimum,field), minimum_xi=number(minimum,"xi"),
                maximum=number(maximum,field), maximum_xi=number(maximum,"xi")))
        if mode == "mode_a" and number(left,"alpha_T") == 1.2:
            zero = next(row for row in subset if number(row,"xi") == 0)
            for field in ("eta_over_s", "sigma_over_T"):
                fn = ((lambda x:(3+2*x)/(3*(1+x))) if field=="eta_over_s" else
                      (lambda x:(1+x)**-.5))
                relative = [(row,number(row,field)/number(zero,field)/fn(number(row,"xi"))-1)
                            for row in subset]
                worst, drift = max(relative,key=lambda item:abs(item[1]))
                references.append(dict(muB_MeV=number(left,"muB_MeV"),T_MeV=number(left,"T_MeV"),
                    observable=field,max_relative_deviation_percent=100*abs(drift),
                    signed_deviation_percent=100*drift,xi_at_max=number(worst,"xi")))
    write(results/"endpoint_summary.csv",endpoints)
    write(results/"rs_reference_comparison.csv",references)

    ordering = []
    for mode in grids:
        for fixed in ((120.,160.,200.) if mode=="mode_b" else (1.,1.1,1.2)):
            field_fixed = "T_MeV" if mode=="mode_b" else "alpha_T"
            densities = {mu:{number(r,"xi"):r for r in grids[mode]
                            if number(r,"muB_MeV")==mu and number(r,field_fixed)==fixed}
                         for mu in (0.,450.,900.)}
            common = sorted(set.intersection(*(set(r) for r in densities.values())))
            for field in ("eta_over_s","zeta_over_s","sigma_over_T"):
                ratios = [number(densities[900.][xi],field)/number(densities[0.][xi],field) for xi in common]
                ordering.append(dict(mode=mode,fixed_variable=field_fixed,fixed_value=fixed,observable=field,
                    shared_points=len(common),increasing_count=sum(number(densities[0.][x],field)<
                        number(densities[450.][x],field)<number(densities[900.][x],field) for x in common),
                    decreasing_count=sum(number(densities[0.][x],field)>
                        number(densities[450.][x],field)>number(densities[900.][x],field) for x in common),
                    high_to_zero_min=min(ratios),high_to_zero_max=max(ratios)))
    write(results/"density_ordering.csv",ordering)

    hierarchy = []
    for (mode,panel,series), subset in groups.items():
        target_rate_ratios=[2*number(r,'n_u')*rates[r['point_id']]['us_to_us']
                            /number(r,'tauinv_u') for r in subset]
        hierarchy.append(dict(mode=mode,panel=panel,series=series,
            strange_below_light_count=sum(number(r,"tau_s")<number(r,"tau_u") for r in subset),
            strange_anti_above_light_count=sum(number(r,"tau_sbar")>number(r,"tau_ubar") for r in subset),
            anti_below_quark_count=sum(number(r,"tau_ubar")<number(r,"tau_u") for r in subset),
            points=len(subset),su_sd_fraction_min=min(number(r,"su_sd_fraction") for r in subset),
            su_sd_fraction_max=max(number(r,"su_sd_fraction") for r in subset),
            su_sd_over_total_light_min=min(target_rate_ratios),
            su_sd_over_total_light_max=max(target_rate_ratios),
            su_sd_above_total_light_count=sum(ratio>1 for ratio in target_rate_ratios)))
    write(results/"flavor_hierarchy.csv",hierarchy)

    windows = []
    species = []
    for left_id,right_id in (("a0653","a0654"),("a0776","a0777")):
        l,r=ids[left_id],ids[right_id]
        dtau = sum((number(r,"tau_"+s)-number(l,"tau_"+s))*number(l,"J_"+s) for s in SPECIES)
        dJ = sum(number(r,"tau_"+s)*(number(r,"J_"+s)-number(l,"J_"+s)) for s in SPECIES)
        dsigma=number(r,"sigma_over_T")-number(l,"sigma_over_T")
        assert math.isclose(dtau+dJ,dsigma,rel_tol=1e-10,abs_tol=1e-14)
        inverse_drop=number(l,"tauinv_ubar")-number(r,"tauinv_ubar")
        contribution_drop=number(l,"dubar_inverse_contribution")-number(r,"dubar_inverse_contribution")
        wl=number(l,"dubar_inverse_contribution")/number(l,"n_u")
        wr=number(r,"dubar_inverse_contribution")/number(r,"n_u")
        windows.append(dict(left_point_id=left_id,right_point_id=right_id,
            T_MeV=number(l,"T_MeV"),alpha_T=number(l,"alpha_T"),
            xi_left=number(l,"xi"),xi_right=number(r,"xi"),
            w_left=wl,w_right=wr,w_change_percent=100*(wr/wl-1),
            tauinv_ubar_left=number(l,"tauinv_ubar"),tauinv_ubar_right=number(r,"tauinv_ubar"),
            tauinv_change_percent=difference(l,r,"tauinv_ubar"),
            dubar_share_of_inverse_drop_percent=100*contribution_drop/inverse_drop,
            tau_ubar_change_percent=difference(l,r,"tau_ubar"),
            sigma_change_percent=difference(l,r,"sigma_over_T"),
            sigma_left=number(l,"sigma_over_T"),sigma_right=number(r,"sigma_over_T"),
            delta_tau_percent=100*dtau/number(l,"sigma_over_T"),
            delta_J_percent=100*dJ/number(l,"sigma_over_T"),
            J_ubar_over_J_u_left_percent=100*number(l,"J_ubar")/number(l,"J_u")))
        for s in SPECIES:
            j0,j1=number(l,"J_"+s),number(r,"J_"+s)
            t0,t1=number(l,"tau_"+s),number(r,"tau_"+s)
            species.append(dict(T_MeV=number(l,"T_MeV"),xi_left=number(l,"xi"),xi_right=number(r,"xi"),
                species=s,tau_left=t0,tau_right=t1,J_left=j0,J_right=j1,
                tau_change_percent=100*(t1/t0-1),J_change_percent=100*(j1/j0-1),
                left_current_fraction_percent=100*t0*j0/number(l,"sigma_over_T"),
                right_current_fraction_percent=100*t1*j1/number(r,"sigma_over_T"),
                delta_tau_percent=100*(t1-t0)*j0/number(l,"sigma_over_T"),
                delta_J_percent=100*t1*(j1-j0)/number(l,"sigma_over_T")))
    write(results/"mott_window_summary.csv",windows)
    write(results/"mott_species_decomposition.csv",species)

    roots=read(results/"mott_temperatures.csv")
    slopes={}
    for meson in ("pi","K"):
        pairs=[(number(r,"xi"),number(r,"T_Mott_MeV")) for r in roots if r["meson"]==meson]
        xm=sum(x for x,y in pairs)/len(pairs);ym=sum(y for x,y in pairs)/len(pairs)
        slopes[meson]=sum((x-xm)*(y-ym) for x,y in pairs)/sum((x-xm)**2 for x,y in pairs)
    extrema={}
    for mode,rows in grids.items():
        for field in ("eta_over_s","zeta_over_s","sigma_over_T"):
            r=min(rows,key=lambda r:number(r,field))
            extrema[f"{mode}_{field}"]={key:r[key] for key in
                ("point_id","T_MeV","muB_MeV","xi","alpha_T",field)}
    soft=max(all_rows,key=lambda r:abs(number(r,"soft_pion_residual")))
    threshold_poles=read(results/'mechanism_poles.csv')
    threshold_areas=read(results/'mechanism_thermal_areas.csv')
    near_threshold=next(r for r in threshold_poles if r['point_id']=='a0777' and r['meson']=='pi')
    low_energy={r['point_id']:r for r in threshold_areas if number(r,'lower_excess_MeV')==0}
    threshold_note=dict(point_id='a0777',xi=.18,
        pion_binding_energy_MeV=number(near_threshold,'threshold_MeV')-number(near_threshold,'mass_MeV'),
        excess_energy_interval_MeV=[0,25],
        left_collision_contribution_inv_fm=number(low_energy['a0776'],'inverse_ubar_contribution'),
        right_collision_contribution_inv_fm=number(low_energy['a0777'],'inverse_ubar_contribution'),
        contribution_change_percent=difference(low_energy['a0776'],low_energy['a0777'],'inverse_ubar_contribution'))
    report=dict(status="complete_grid_analyzed",points=len(all_rows),sigma_decomposition_max_relative_error=sigma_closure,
        channel_rate_reassembly_max_relative_error=inverse_error,
        terminal_replay_max_relative_error=max(number(r,"terminal_replay_error") for r in all_rows),
        bulk_replay_max_relative_error=max(abs(number(r,"bulk_replay_error")) for r in all_rows),
        soft_pion_residual_max=abs(number(soft,"soft_pion_residual")),soft_pion_residual_point=soft["point_id"],
        isospin_tau_max_relative_error=max(abs(number(r,"tau_d")/number(r,"tau_u")-1) for r in all_rows),
        zero_density_particle_anti_max_relative_error=max(abs(number(r,"tau_ubar")/number(r,"tau_u")-1)
            for r in all_rows if number(r,"muB_MeV")==0),
        spectral_mott_fit_slopes_MeV=slopes,extrema=extrema,windows=windows,
        near_threshold_pion_thermal_area=threshold_note,
        inputs={f"{mode}_grid.csv":hashlib.sha256((results/f"{mode}_grid.csv").read_bytes()).hexdigest() for mode in grids})
    (results/"physics_audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,required=True)
    audit(parser.parse_args().results)
