# meson_density phase-f stage bundle

output_dir: `data\outputs\results\relaxtime\meson_density\phase_f_stage_v1`

source files:
- `stable`: `data\outputs\results\relaxtime\scan\meson_density_scan_208_220_step2.csv` -> `data\outputs\results\relaxtime\meson_density\phase_f_stage_v1\meson_density_scan_208_220_step2.csv` (7 data rows)
- `current_bu`: `data\outputs\results\relaxtime\scan\phase_shift_meson_density_scan_208_220_step2.csv` -> `data\outputs\results\relaxtime\meson_density\phase_f_stage_v1\phase_shift_meson_density_scan_208_220_step2.csv` (7 data rows)
- `gbu_reference`: `data\outputs\results\relaxtime\scan\phase_shift_meson_density_scan_gbu_reference_208_220_step2.csv` -> `data\outputs\results\relaxtime\meson_density\phase_f_stage_v1\phase_shift_meson_density_scan_gbu_reference_208_220_step2.csv` (7 data rows)
- `strict_bw_stage2`: `data\outputs\results\relaxtime\scan\strict_bw_meson_density_scan_stage2_208_220_step2_converged.csv` -> `data\outputs\results\relaxtime\meson_density\phase_f_stage_v1\strict_bw_meson_density_scan_stage2_208_220_step2_converged.csv` (7 data rows)

note:

- This directory is the first formal meson-density result bundle promoted out of `scan/` for downstream baseline/regression governance.
- It remains a stage-level bundle, not the final paper-production directory.
