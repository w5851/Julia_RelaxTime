# scripts/debug/

调试与排障脚本目录（一次性或临时诊断用途）。

## 当前内容

- `debug_tau_point.jl`：定位单点弛豫时间计算问题
- `triage_integration_unit_files.jl`：排查 integration/unit 测试文件组织或引用问题
- `triage_pnjl_unit_files.jl`：排查 PNJL 单元测试文件组织或引用问题

## 约定

- `scripts/relaxtime/`：面向稳定工作流的数据/结果生成脚本。
- `scripts/debug/`：定位问题的临时脚本，不纳入 CI 常规路径。

```bash
julia --project=. scripts/debug/debug_tau_point.jl
```
