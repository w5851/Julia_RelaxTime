# Ward/Møller Mott threshold response

This is the numerical source for the transport manuscript's Fig. 7. The
calculation uses Julia 1.12.5 and the core at
`096199b96d34045bf8e2d9611c2d786d6f92fb8e` (PR #352):
variational PNJL with a 96×24 background, Ward-matched one-line integrals
and couplings, medium-frame Møller flux, and direction-averaged linear
two-line loops. The local transport integrations tighten the collision-energy
sampling relative to the production grid and keep the other quadratures fixed.

[Numerical manifest](numerical_manifest.json) binds the effective profile,
remote runs and driver hashes, the directly solved threshold roots, and
the compact tables. Each window retains 33 fine and 11 reference nodes.
`window_replay.csv` records the terminal shear correction from the cached
\(s\bar u\) rates on W. The displayed scalar secants use computed endpoints
at \(\Delta\xi=10^{-3}\); the negative kaon wing also uses wider secants.
The remaining pseudoscalar secants keep the native fine intervals. The
maximum scalar display slope difference across the four intervals with
matching reference endpoints is \(1.38\times 10^{-6}\). The former narrow
scalar residual was below the numerical resolution and is omitted.

The isolated scalar integral is the squared mixed \(\sigma/\sigma'\)
contribution to elastic \(u\bar u\to u\bar u\), with the production
final-state factors and a continuous Møller thermal weight. The two
numerical tiers were run on separate remote machines. They compare
stencil, loop, and thermal quadrature refinement with an independently
evaluated threshold coefficient. The integration domain is
\(0<x<0.01\,\mathrm{fm}^{-2}\), where \(x=\hat{s}-4M_u^2\).
The peak area uses 64 half widths on either side of the resonance.
This integral does not enter the species-rate assembly, so its data are unchanged.
These fits concern the isolated integral; panels (a)–(c) show complete
shear slopes including the other collision channels.

## Reproduction

Use the [remote workflow](../../../../../../../.github/workflows/ward-moller-mott-local-response.yml)
or run the Julia drivers with the pinned core in an isolated remote checkout.
The workflow specifies the command, shard counts, and reference budget
recovery. The numerical manifest identifies the actual run artifacts
selected for this case. Large thermal caches and worker logs are not tracked.

The plot can be regenerated from these tables without running Julia:

```powershell
$casePath = 'relaxtime/transport/mott_threshold/ward_moller_manuscript_20261121_v1'
python scripts/analysis/relaxtime/plot_ward_moller_mott_response.py --data "data/outputs/results/$casePath" --out 'D:/Temp/mott-figure-rebuild'
python scripts/plotting/validate_plot_artifact.py "data/outputs/figures/$casePath/plot_manifest.json"
```

The PDF and 600 dpi color/gray PNGs support manuscript review. Their hashes
and the figure coordinates are in the
[plot manifest](../../../../../figures/relaxtime/transport/mott_threshold/ward_moller_manuscript_20261121_v1/plot_manifest.json).
