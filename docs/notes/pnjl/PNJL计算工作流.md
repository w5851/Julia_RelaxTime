# PNJL 计算工作流（历史导航）

> 状态：deprecated 兼容说明
>
> 原始版本：2025-12-30
>
> 当前导航更新：2026-07-10

本页原先描述早期 PNJL 相结构脚本、内部模块路径和旧输出目录。相关实现入口与产物合同已经迁移，旧正文不再作为当前执行依据；原始内容可通过 git history 追溯。

## 当前执行入口

需要运行 PNJL/RPNJL 相结构计算时，按以下顺序阅读：

1. [PNJL 相结构计算 SOP](../../guides/sop/workflows/pnjl_phase_structure.md)
2. [相结构 API 总览](../../api/models/phase/README.md)
3. [稳定脚本入口清单](../../guides/scripts/README.md)

当前稳定 CLI：

```text
scripts/pnjl/calculate_phase_structure.jl
```

当前统一实现入口：

```text
Models.run_phase_pipeline
```

默认配置：

```text
config/models/pnjl/phase_pipeline_default.toml
```

## 当前产物口径

新计算结果写入 `data/outputs/results/` 下的命名 case 目录，并以 `run_manifest.json` 中的 effective config、config hash、git commit 和 artifact paths 为实际运行口径。

历史 `data/reference` 结果不能仅因路径或文件名被视为当前正式 reference；任何晋升都必须经过现行 SOP 中的收敛、回归和人工审阅门槛。

## 历史说明边界

- 旧模块拆分、旧内部路径和旧待办状态只具有迁移溯源价值。
- 新代码和新文档不得从本页复制历史入口。
- 如需解释早期相结构方法选择，可结合 git history 与 `docs/notes/pnjl/一阶相变方法对比.md` 阅读，并明确标注为历史证据。
