# `pnjl_mag` equilibrium replay v1

本目录是 `pnjl_mag@e1fc81d3c3c9d220c49972e54307b66a604cb9db` 在当前机器上的轻量 equilibrium 重放诊断。

## 口径

- Python/JAX/Optimistix 使用外部仓库自己的 `.venv` 和 `uv.lock`；
- 固定 `muB=0`、`eB={0.2, 0.4, 0.8} GeV^2`；
- 沿作者脚本相同的单 seed 温度 continuation：`T=300..50 MeV`，步长 `1 MeV`；
- `p_num=128`、Landau levels=`80`、`zeta_num=256`；
- 只提取 `T={50, 150, 240} MeV` 的 9 行；每个场值的 251 点、三个场值合计
  753 个 continuation 点不完整写入主项目；
- 每个点记录五维状态、`Omega` 和最大 gap residual，并与作者已提交的
  `data/orders_muB0.csv` 对照。

## 结论边界

本机重放的 9 个代表点与作者 CSV 的五维状态逐字段一致，且 residual 保持有限并很小。
这证明了依赖环境和作者单 seed continuation 路线在当前机器上可以重现。
它仍是 `diagnostic_only`：没有多 seed 分支集合、全局极小值搜索、Julia/外部 solver
分支等价或 production target admission。

生成命令：

```powershell
python scripts/analysis/pnjl/extract_pnjl_mag_equilibrium_replay.py
```

表格、hash、运行耗时和边界见 [`manifest.json`](manifest.json) 与
[`provenance.json`](provenance.json)。
