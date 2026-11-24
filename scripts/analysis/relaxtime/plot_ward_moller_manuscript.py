"""Render the revised PNJL transport figures from frozen numerical results."""
from __future__ import annotations
import argparse, csv, json, math, sys, tempfile, shutil
from pathlib import Path
from dataclasses import replace

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import FixedLocator, MaxNLocator, ScalarFormatter, FormatStrFormatter, LogLocator
from scripts.plotting.plot_series import SeriesStyle, SeriesStyleRegistry
from scripts.plotting.plot_style import load_profile, configure_axis_ticks
from scripts.plotting.plot_manifest import (build_manifest,input_record,output_record,
    generator_record,runtime_record,write_manifest)
from scripts.plotting.plot_quality import export_figure,inspect_export,placement_limits
from scripts.plotting.plot_accessibility import accessibility_record
from scripts.plotting.plot_bundle import build_bundle
from scripts.plotting.validate_plot_artifact import validate_manifest,validate_manifest_record
from scripts.plotting import transport_panels as panels
from scripts.plotting.publication_review import (
    REVIEW_COLORS, figure_context, grayscale_review_record, register_shared_series_styles, write_sop_snapshot,
)

COLORS=REVIEW_COLORS
KINDS=('spectra','relaxation','transport','fixed_T','mechanism')
SPECTRAL_XI = (-.3, 0., .3)
SPECTRAL_STYLES = SeriesStyleRegistry('ward_moller_spectra_xi_v1', tuple(
    style
    for xi, color, line in zip(SPECTRAL_XI, COLORS, ('-', '--', '-.'))
    for style in (SeriesStyle(f'spectrum.{xi:g}', color, line),
                  SeriesStyle(f'mott_root.{xi:g}', color, 'None', 'o'))
))
MECHANISM_STYLES = SeriesStyleRegistry('ward_moller_mechanism_lines_v1', (
    SeriesStyle('pion_mass', COLORS[0], '-', 'o'),
    SeriesStyle('pion_threshold', '#333333', '--'),
    SeriesStyle('sigma_left', COLORS[0], '-'),
    SeriesStyle('sigma_right', COLORS[1], '--'),
    SeriesStyle('total_rate', '#222222', '-', 'o'),
    SeriesStyle('net_species', '#222222', 'None', 'o'),
))
MECHANISM_LINES = {style.series_id: style for style in MECHANISM_STYLES.styles}

def rows(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))

def n(row,field):return float(row[field])

def declared_legends(quality,scope):
    return [dict(legend_index=item['legend_index'],host_axes_index=item['host_axes_index'],
        title=item['title'],labels=item['labels'],placement='in_axes' if item['host_axes_index'] is not None else 'outside_axes',
        location='declared clear panel corner',scope=scope) for item in quality['legend_layout']]

def export_record(fig,stem,meta,profile,input_paths,calculation_sha,*,description,axes,series,registry,placements=None,panel_specs=None):
    # Legends added after a bound style context inherit Matplotlib defaults.
    # Bind their ordinary and mathematical faces to the same manuscript fonts.
    for artist in fig.findobj(match=Text):
        artist.set_fontfamily('Times New Roman')
        artist.set_math_fontfamily('stix')
    fig.set_dpi(600)
    with tempfile.TemporaryDirectory() as tmp:
        temporary_stem=Path(tmp)/stem.name
        _,quality=export_figure(fig,temporary_stem,profile,formats=['png'],intended_width_inches=6.75)
        shutil.copyfile(temporary_stem.with_suffix('.png'),stem.with_suffix('.png'))
    outputs=[output_record(stem.with_suffix('.png'),fmt='png',dpi=600,vector=False)]
    color=stem.with_suffix('.png')
    outputs[0].update(inspection=inspect_export(color),role='color_review')
    outputs.append(grayscale_review_record(color))
    series_styles=register_shared_series_styles(series, registry)
    rendering=dict(column='double_column',figure_size_inches=quality['figure_size_inches'],
        size_override_reason=description,quality=quality,placement_limits=placement_limits(quality,profile,outputs),
        legend_policy='declared_geometry_checked',legend_placements=placements or declared_legends(quality,description),
        legend_outside=all(not item['overlapped_axes'] for item in quality['legend_layout']),
        color_route='color_print_and_online',
        delivery_stage='png_review',vector_delivery_pending=True,review_contract='prd_review_v1',
        accessibility=accessibility_record([dict(value=str(i),color=c,noncolor_encoding=style)
             for i,(c,style) in enumerate(zip(COLORS,('-','--','-.')))],
            backgrounds=['#ffffff'],background_scope='opaque white panel backgrounds'),
        path_simplification=False,panel_specs=panel_specs or [],
        series_styles=series_styles,
        series_identity_scope='shared style identities refer to the panel series and selected frozen rows')
    record=build_manifest(asset_id='ward_moller_'+stem.name,figure_family='pnjl_transport_manuscript',
        case_slug=stem.parent.name,figure_mode='audit',semantic_status='new_prescription_computed_samples',
        style_profile=profile.profile_id,publication_scope='internal_review',
        generator=generator_record(__file__,command=' '.join(sys.argv),
            runtime=runtime_record({'matplotlib':matplotlib.__version__})),
        inputs=[input_record(path,role=role) for path,role in input_paths],outputs=outputs,axes=axes,series=series,
        selection_rule=description,interpolation_policy='none',connector_policy='forbidden',
        missing_value_policy='reject missing/nonfinite/duplicate rows',
        validation=dict(finite=True,duplicate_keys=True,support=True,strict_gate=False),
        rendering=rendering,calculation_sha=calculation_sha)
    record.update(manuscript_eligible=False,current_publication_layer=False,delivery_stage='png_review',
        vector_delivery_pending=True,solver_called=False,canonical_data_modified=False)
    errors=validate_manifest_record(record)
    if errors:raise ValueError(stem.name+': '+'; '.join(errors))
    print(json.dumps({'figure':stem.name,'glyph_mm':quality['minimum_capital_numeral_height_mm'],
                      'text_overlaps':quality['text_overlap_pairs'],'legend_curves':quality['legend_curve_overlap_count']}),flush=True)
    plt.close(fig)
    return record

def render_spectra(source,mott,profile):
    data=rows(source);roots=rows(mott)
    with figure_context(profile,compact=True) as context:
        fig,axs=plt.subplots(2,2,figsize=(6.75,4.2),sharex=True,sharey='row')
        context.bind(fig);fig.subplots_adjust(left=.132,right=.978,bottom=.144,top=.923,wspace=.11,hspace=.24)
        for col,meson in enumerate(('pi','K')):
            for xi, style, mott_style in zip(SPECTRAL_XI, SPECTRAL_STYLES.styles[::2], SPECTRAL_STYLES.styles[1::2]):
                subset=sorted([r for r in data if r['meson']==meson and abs(n(r,'xi')-xi)<1e-9],key=lambda r:n(r,'T_MeV'))
                assert len(subset)==201
                T=[n(r,'T_MeV') for r in subset]
                for ax,field in ((axs[0,col],'mass_MeV'),(axs[0,col],'threshold_MeV'),(axs[1,col],'width_MeV')):
                    ax.plot(T,[n(r,field) for r in subset],lw=1.2,**style.plot_kwargs())
                root=next(r for r in roots if r['meson']==meson and abs(n(r,'xi')-xi)<1e-9)
                # The marker uses a directly solved Mott point, including its own frozen threshold mass.
                axs[0,col].plot([n(root,'T_Mott_MeV')],[n(root,'threshold_MeV')],
                    mec='white',mew=.45,ms=5.3,zorder=5,**mott_style.plot_kwargs())
            for row in range(2):
                ax=axs[row,col]
                ax.text(.97,.97,f'({chr(97+row*2+col)})'+(' Pion' if row==0 and col==0 else ' Kaon' if row==0 else ''),
                    transform=ax.transAxes,ha='right',va='top',fontsize=11)
                ax.set_xlim(100,280);ax.xaxis.set_major_locator(FixedLocator([100,150,200,250]))
                ax.set_ylim((0,1000) if row==0 else (-15,650))
                ax.yaxis.set_major_locator(FixedLocator([0,400,800] if row==0 else [0,300,600]))
                configure_axis_ticks(ax,profile)
        axs[0,0].set_ylabel('Energy (MeV)',labelpad=7);axs[1,0].set_ylabel('Width (MeV)',labelpad=7)
        fig.supxlabel(r'$T$ (MeV)',x=.56,y=.012,fontsize=13)
        handles=[Line2D([],[],color=COLORS[i],ls=style,lw=1.4) for i,style in enumerate(('-','--','-.'))]
        axs[1,0].legend(handles,[r'$-0.3$',r'$0$',r'$+0.3$'],title=r'$\xi$',loc='upper left',frameon=False,
            fontsize=11,title_fontsize=13,handlelength=1.55,handletextpad=.35,labelspacing=.10,borderpad=.05,borderaxespad=.35)
        return fig

def transport_inputs(a,b):
    groups={**panels.validate_rows(rows(a),'mode_a'),**panels.validate_rows(rows(b),'mode_b')}
    gaps={('mode_a','muB900.0','alpha1.0'):[(-.003,.003)],
          ('mode_b','T120.0','muB900.0'):[(-.13,-.12)]}
    return panels.TransportInputs(groups,gaps,{'mode_a':a,'mode_b':b})

def render_transport(kind,data,profile):
    if kind=='fixed_T':
        fig,specs,placements=panels.render_fixed_T(data,profile)
        fig.axes[5].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        panels.capture_ticks(fig.axes[5],specs[5])
        return fig,specs,placements
    name='figure1_relaxation_times_comparison' if kind=='relaxation' else 'figure2_transport_coefficients_comparison'
    fig,specs,placements=panels.render_composite(name,data,profile)
    if kind=='transport':
        # Explicitly preserve the linear axes of the reviewed viscosity/conductivity composite.
        for ax,spec in zip(fig.axes,specs):
            ax.set_yscale('linear');ax.relim();ax.autoscale_view(scalex=False)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4,min_n_ticks=3,steps=[1,2,5,10]))
            ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
            spec.update(axis_scale='linear',axis_scale_reason='reviewed linear transport-ratio panels')
            configure_axis_ticks(ax,profile);panels.capture_ticks(ax,spec)
        first=fig.axes[0];lo,hi=first.get_ylim();first.set_ylim(min(lo,.07),hi)
        first.axhline(1/(4*math.pi),color='#666666',ls=':',lw=.8)
        first.text(.51,1/(4*math.pi),r'$1/(4\pi)$',transform=first.get_yaxis_transform(),fontsize=11,va='bottom',
                   bbox=dict(facecolor='white',edgecolor='none',pad=.3))
        # Recreate shared keys after the data-derived axes changed.
        for legend in panels.visible_legends(fig):legend.remove()
        fig.subplots_adjust(top=.94)
        placements=panels.place_legends(fig,specs,profile,hosts=(6,2))
        legends=panels.visible_legends(fig)
        for placement in placements:
            anchor=(.06,1.) if placement['role']=='endpoints' else (.96,1.)
            legends[placement['legend_index']].set_bbox_to_anchor(anchor)
            placement.update(anchor_axes=list(anchor),
                selection_reason='inset the shared key to separate it from adjacent tick labels')
        panels.capture_ticks(first,specs[0])
    return fig,specs,placements


def series_from_panels(specs):
    output=[]
    for spec in specs:
        for curve in spec['series']:
            output.append(dict(curve,series_id=curve['curve_id'],style_ids=[curve['series_id']],state='computed_native_series',
                support_rule='computed source rows, split only at the declared first-order pair',
                mask_rule='none'))
    return output


def spectral_metadata(source,mott):
    data,roots=rows(source),rows(mott)
    series=[];specs=[]
    for panel,meson,fields in (('a','pi',('mass_MeV','threshold_MeV')),
                             ('b','K',('mass_MeV','threshold_MeV')),
                             ('c','pi',('width_MeV',)),('d','K',('width_MeV',))):
        for xi in SPECTRAL_XI:
            selected=[i+2 for i,r in enumerate(data) if r['meson']==meson and abs(n(r,'xi')-xi)<1e-9]
            for field in fields:
                series.append(dict(series_id=f'{panel}.{meson}.{xi}.{field}',style_ids=[f'spectrum.{xi:g}'],
                    state='computed_native_series',
                    source_csv=str(source),source_column=field,source_rows=selected,point_count=len(selected),
                    panel=panel,meson=meson,xi=xi,support_rule='all 201 frozen temperatures in the input scan',
                    mask_rule='viewport displays 100 to 280 MeV without changing source vertices'))
            if panel in ('a','b'):
                selected_root=[i+2 for i,r in enumerate(roots) if r['meson']==meson and abs(n(r,'xi')-xi)<1e-9]
                series.append(dict(series_id=f'{panel}.{meson}.{xi}.direct_mott',style_ids=[f'mott_root.{xi:g}'],
                    state='confirmed_root',
                    source_csv=str(mott),source_rows=selected_root,point_count=1,panel=panel,
                    support_rule='directly solved Mott root and its own threshold mass',mask_rule='none'))
        specs.append(dict(panel=panel,meson=meson,observables=list(fields),x_limits=[100,280],
            y_limits=[0,1000] if panel in ('a','b') else [-15,650],axis_scale='linear'))
    return series,specs


def render_mechanism(meta,profile):
    spectra=rows(meta/'mechanism_spectrum.csv');poles=rows(meta/'mechanism_poles.csv')
    areas=rows(meta/'mechanism_thermal_areas.csv');grid=rows(meta/'mechanism_transport.csv')
    windows=rows(meta/'mott_window_summary.csv');split=rows(meta/'mott_species_decomposition.csv')
    specs=[]
    with figure_context(profile,compact=True) as context:
        fig,axs=plt.subplots(4,2,figsize=(6.75,7.95))
        context.bind(fig)
        fig.subplots_adjust(left=.14,right=.98,bottom=.065,top=.962,wspace=.27,hspace=.55)
        for col,win in enumerate(windows):
            T=n(win,'T_MeV');left_id,right_id=win['left_point_id'],win['right_point_id']
            local=sorted([r for r in grid if abs(n(r,'T_MeV')-T)<1e-7],key=lambda r:n(r,'xi'))
            ids=[r['point_id'] for r in local];xs=[n(r,'xi') for r in local]
            pirows=[next(r for r in poles if r['point_id']==identifier and r['meson']=='pi') for identifier in ids]
            top=axs[0,col]
            top.plot(xs,[n(r,'mass_MeV') for r in pirows],lw=1.25,ms=4.5,label=r'$M_\pi$',
                     **MECHANISM_LINES['pion_mass'].plot_kwargs())
            top.plot(xs,[n(r,'threshold_MeV') for r in pirows],lw=1.15,label=r'$M_u+M_d$',
                     **MECHANISM_LINES['pion_threshold'].plot_kwargs())
            top.set_xticks(xs,[f'{x:.2f}' for x in xs]);top.set_xlabel(r'$\xi$')
            top.set_title(rf'$T={T:.0f}$ MeV',pad=6,fontsize=13)
            top.yaxis.set_major_locator(MaxNLocator(4,min_n_ticks=3))
            top.set_ylabel('Energy (MeV)' if col==0 else '')
            if col==0:
                top.legend(loc='upper left',fontsize=11,frameon=False,handlelength=1.4,handletextpad=.35,labelspacing=.10)

            ax=axs[1,col]
            peak=0.
            for i,identifier in enumerate((left_id,right_id)):
                selected=sorted([r for r in spectra if r['point_id']==identifier],key=lambda r:n(r,'excess_energy_MeV'))
                xi=next(n(r,'xi') for r in grid if r['point_id']==identifier)
                ax.plot([n(r,'excess_energy_MeV') for r in selected],[n(r,'sigma_fm2') for r in selected],
                    lw=1.25,label=rf'$\xi={xi:.2f}$',
                    **MECHANISM_LINES[('sigma_left','sigma_right')[i]].plot_kwargs())
                peak=max(peak,*(n(r,'sigma_fm2') for r in selected if n(r,'excess_energy_MeV')<=80))
            ax.set_yscale('log');ax.set_xlim(0,80);ax.set_xlabel(r'$\sqrt{\hat s}-M_u-M_d$ (MeV)')
            low,high=ax.get_ylim();ax.set_ylim(low,max(high,2.5*peak))
            ax.tick_params(axis='y',labelsize=13)
            ax.set_ylabel(r'$\sigma_{d\bar u}$ (fm$^2$)' if col==0 else '')
            ax.legend(loc='upper right',fontsize=11,frameon=False,handlelength=1.4,handletextpad=.35,labelspacing=.10)
            # The energy-resolved thermal area is shown next to its total-rate consequence.
            ax=axs[2,col]
            positions=[0,1];bottom=[0.,0.]
            interval_names=['0–25','25–75','75–150','>150']
            fills=[COLORS[0],COLORS[1],COLORS[2],'#555555'];hatches=['','///','xx','..']
            for j,lo in enumerate((0.,25.,75.,150.)):
                vals=[sum(n(r,'inverse_ubar_contribution') for r in areas if
                    r['point_id']==identifier and n(r,'lower_excess_MeV')==lo) for identifier in (left_id,right_id)]
                ax.bar(positions,vals,bottom=bottom,width=.48,color=fills[j],hatch=hatches[j],
                    edgecolor='white',linewidth=.5,label=interval_names[j])
                bottom=[b+v for b,v in zip(bottom,vals)]
            ax.plot(positions,[n(win,'tauinv_ubar_left'),n(win,'tauinv_ubar_right')],
                ms=4.8,lw=1.1,**MECHANISM_LINES['total_rate'].plot_kwargs())
            ax.set_xticks(positions,[rf'$\xi={n(win,"xi_left"):.2f}$',rf'$\xi={n(win,"xi_right"):.2f}$'])
            ax.set_ylim(0,1.05)
            ax.yaxis.set_major_locator(FixedLocator([0,.25,.50,.75,1.]))
            ax.set_ylabel(r'Collision rate (fm$^{-1}$)' if col==0 else '')
            if col==0:
                ax.legend(title='Excess energy (MeV)',loc='upper left',bbox_to_anchor=(.02,.68),
                    fontsize=11,title_fontsize=11,frameon=False,ncol=2,
                    handlelength=1.1,handletextpad=.3,labelspacing=.10,columnspacing=.55,borderaxespad=.15)
            else:
                ax.text(.05,.62,r'Total $\tau_{\bar u}^{-1}$',transform=ax.transAxes,fontsize=13)
                ax.text(.05,.40,r'Bars: $d\bar u\to d\bar u$',transform=ax.transAxes,fontsize=11)

            ax=axs[3,col]
            grouped=[(r'$u+d$',('u','d')),(r'$s$',('s',)),(r'$\bar{u}+\bar{d}$',('ubar','dbar')),(r'$\bar{s}$',('sbar',))]
            all_entries=[r for r in split if abs(n(r,'T_MeV')-T)<1e-7]
            tau=[sum(n(r,'delta_tau_percent') for r in all_entries if r['species'] in ss) for _,ss in grouped]
            weight=[sum(n(r,'delta_J_percent') for r in all_entries if r['species'] in ss) for _,ss in grouped]
            ax.axhline(0,color='#888888',lw=.65)
            ax.bar([i-.14 for i in range(4)],tau,width=.27,color=COLORS[0],label=r'$\Delta_\tau$')
            ax.bar([i+.14 for i in range(4)],weight,width=.27,color=COLORS[1],hatch='///',label=r'$\Delta_J$')
            ax.plot(range(4),[a+b for a,b in zip(tau,weight)],ms=4.5,**MECHANISM_LINES['net_species'].plot_kwargs())
            ax.set_xticks(range(4),[label for label,_ in grouped],fontsize=13)
            ax.set_ylabel(r'$\Delta(\sigma/T)/(\sigma/T)_L$ (%)' if col==0 else '')
            ax.yaxis.set_major_locator(FixedLocator([-2,-1,0,1,2]));ax.set_ylim(-2.8,2.)
            ax.text(.98,.05,rf'Net: {n(win,"sigma_change_percent"):+.2f}%',transform=ax.transAxes,
                ha='right',va='bottom',fontsize=11)
            ax.legend(loc='upper right',fontsize=13,frameon=False,ncol=2,handlelength=.95,
                handletextpad=.3,labelspacing=.05,columnspacing=.55,borderpad=.05,borderaxespad=.25)
            for row in range(4):
                ax=axs[row,col];configure_axis_ticks(ax,profile)
                if row==1:
                    ax.yaxis.set_minor_locator(LogLocator(base=10,subs=(2,3,4,5,6,7,8,9),numticks=100))
                ax.text(0,1.015,f'({chr(97+2*row+col)})',transform=ax.transAxes,ha='left',va='bottom',fontsize=11)
                specs.append(dict(panel=chr(97+2*row+col),T_MeV=T,point_ids=ids if row==0 else [left_id,right_id],
                    observable=('mass and threshold','elastic cross section','energy-resolved collision rate','species conductivity changes')[row],
                    x_limits=list(ax.get_xlim()),y_limits=list(ax.get_ylim()),axis_scale=ax.get_yscale()))
        return fig,specs


def mechanism_series(meta):
    series=[]
    style_ids = {'mass_MeV': ['pion_mass'], 'threshold_MeV': ['pion_threshold'],
                 'sigma_fm2': ['sigma_left','sigma_right'], 'tauinv_ubar': ['total_rate'],
                 'delta_tau_percent': ['net_species'], 'delta_J_percent': ['net_species']}
    windows=rows(meta/'mott_window_summary.csv')
    for col,window in enumerate(windows):
        T=n(window,'T_MeV')
        for offset,name,fields in ((0,'mechanism_poles.csv',('mass_MeV','threshold_MeV')),
                                  (2,'mechanism_spectrum.csv',('sigma_fm2',)),
                                  (4,'mechanism_thermal_areas.csv',('inverse_ubar_contribution',)),
                                  (4,'mechanism_transport.csv',('tauinv_ubar',)),
                                  (6,'mott_species_decomposition.csv',('delta_tau_percent','delta_J_percent'))):
            data=rows(meta/name)
            source_rows=[i+2 for i,r in enumerate(data) if abs(n(r,'T_MeV')-T)<1e-7
                         and (r.get('meson','pi')=='pi')]
            for field in fields:
                series.append(dict(series_id=f'{chr(97+offset+col)}.{name}.{field}',
                    style_ids=style_ids.get(field, []),state='computed_native_series',source_csv=str(meta/name),source_column=field,
                    source_rows=source_rows,point_count=len(source_rows),panel=chr(97+offset+col),T_MeV=T,
                    support_rule='computed fixed-temperature local window and exact finite-difference decomposition',
                    mask_rule='selected pion records or fixed window; energy viewport0to80MeV',
                    aggregation_rule='species grouped as u+d,s,ubar+dbar,sbar; areas stacked by recorded energy interval'))
    return series

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--kind',choices=(*KINDS,'all'),required=True)
    parser.add_argument('--review-results-dir',type=Path,help='Case results directory for the frozen SOP snapshot')
    args=parser.parse_args()
    profile=load_profile('candidate_aps_v2')
    palette={**profile.data['palette'],'colors':[*COLORS,*profile.colors[3:]]}
    profile=replace(profile,data={**profile.data,'palette':palette})
    out=args.out.resolve();meta=args.results.resolve();out.mkdir(parents=True,exist_ok=True)
    numerical=json.loads((meta/'numerical_manifest.json').read_text(encoding='utf-8'))
    inputs=[(meta/'numerical_manifest.json','frozen numerical provenance'),(profile.path,'plot profile')]
    inputs.extend((ROOT/f'scripts/plotting/{helper}.py','plotting helper') for helper in (
        'publication_review','transport_panels','plot_style','plot_series','plot_manifest','plot_quality',
        'plot_provenance','plot_bundle','plot_accessibility','validate_plot_artifact'))
    if args.review_results_dir is not None:
        inputs.append((write_sop_snapshot(args.review_results_dir.resolve()),'production figure SOP snapshot, LF as in Git'))
    records=[];charts=[]
    wanted=KINDS if args.kind=='all' else (args.kind,)
    for kind in wanted:
      sources=list(inputs)
      if kind=='spectra':
        source,mott=meta/'spectra.csv',meta/'mott_temperatures.csv'
        fig=render_spectra(source,mott,profile);sources.extend([(source,'computed spectral series'),(mott,'direct Mott roots')])
        axes=[dict(field='T',source_unit='MeV',display_unit='MeV',label=r'$T$ (MeV)'),
              dict(field='mass_width_threshold',source_unit='MeV',display_unit='MeV',label='Energy (MeV)'),
              dict(field='width',source_unit='MeV',display_unit='MeV',label='Width (MeV)')]
        series,specs=spectral_metadata(source,mott)
        record=export_record(fig,out/'pion_kaon_spectra',meta,profile,sources,numerical['calculation_sha'],
            description='pi/K at zero density; xi=-0.3,0,0.3; 201 computed temperatures; direct Mott roots',
            axes=axes,series=series,registry=SPECTRAL_STYLES,panel_specs=specs)
      elif kind in ('relaxation','transport','fixed_T'):
        a,b=meta/'mode_a_grid.csv',meta/'mode_b_grid.csv'
        data=transport_inputs(a,b);fig,specs,placements=render_transport(kind,data,profile)
        sources.extend([(a,'recomputed phase-referenced transport'),(b,'recomputed fixed-temperature transport')])
        stems={'relaxation':'relaxation_times','transport':'transport_coefficients','fixed_T':'transport_fixed_T'}
        axes=[dict(field='xi',source_unit='dimensionless',display_unit='dimensionless',label=r'$\xi$')]
        for observable in sorted({spec['observable'] for spec in specs}):
            axes.append(dict(field=observable,source_unit='fm' if observable.startswith('tau_') else 'dimensionless',
                display_unit='fm' if observable.startswith('tau_') else 'dimensionless',label=panels.LABELS[observable]))
        record=export_record(fig,out/stems[kind],meta,profile,sources,numerical['calculation_sha'],
            description='original published scan points, recomputed with Ward A, Moller flux and sphere-averaged linear B0',
            axes=axes,placements=placements,series=series_from_panels(specs),
            registry=panels.make_registry(specs[0]['mode_key'],profile),panel_specs=specs)
      else:
        fig,specs=render_mechanism(meta,profile)
        names=['mechanism_spectrum.csv','mechanism_poles.csv','mechanism_thermal_areas.csv',
               'mechanism_transport.csv','mott_window_summary.csv','mott_species_decomposition.csv']
        sources.extend((meta/name,'computed local Mott mechanism') for name in names)
        axes=[dict(field=field,source_unit=unit,display_unit=unit,label=label) for field,unit,label in (
            ('xi','dimensionless',r'$\xi$'),('mass_threshold','MeV','Energy (MeV)'),
            ('excess_energy','MeV',r'$\sqrt{\hat s}-M_u-M_d$ (MeV)'),('sigma','fm2',r'$\sigma_{d\bar u}$ (fm$^2$)'),
            ('collision_rate','fm-1','Collision rate (fm$^{-1}$)'),('relative_conductivity_change','percent','Relative conductivity change (%)'))]
        record=export_record(fig,out/'mott_transport_mechanism',meta,profile,sources,numerical['calculation_sha'],
            description='two complete fixed-temperature Mott windows at muB=900 MeV; exact conductivity decomposition',
            axes=axes,panel_specs=specs,series=mechanism_series(meta),registry=MECHANISM_STYLES)
      records.append(record);charts.append(dict(name=kind,asset_id=record['asset_id']))
    manifest=out/'plot_manifest.json'
    record=records[0] if len(records)==1 else build_bundle(dict(case=out.name,charts=charts),records)
    write_manifest(manifest,record,overwrite=True)
    violations=validate_manifest(manifest)
    if violations:raise ValueError('; '.join(violations))

if __name__=='__main__':main()
