# T150 最小复扫记录

- 命令：`julia --project=. scripts/relaxtime/run_gap_transport_scan.jl --tmin 150 --tmax 150 --tstep 1 --mubmin 0 --mubmax 0 --mubstep 1 --xi-list -0.5 --mode finite_15 --tau-p-nodes 28 --tau-angle-nodes 8 --tau-phi-nodes 8 --tau-n-sigma 6 --sigma-grid-n 128 --output "D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05.csv" --failed-points-output "D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05_failed_points.csv" --overwrite`
- 结果：`xi=-0.5` 成功，输出行存在（`T_MeV=150.0, muB_MeV=0.0, xi=-0.5`）
- sidecar：`D:\Desktop\Temp\relaxtime_t150_repro\transport_vs_xi_T150_muB0_xim05_failed_points.csv`（仅表头，无失败行）
- 结论：T150 连续性起点可恢复，且失败机制具备可见 sidecar 通道
