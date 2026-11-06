---
name: plotting-sop
description: 使用 Julia_RelaxTime 的共享绘图工具从冻结结果制作、审阅和交付科研图像，按 figure_production_v2 校验来源、图面质量与交付阶段。
---

# 项目绘图 SOP

执行[图像生产工作流程](../../../docs/guides/sop/workflows/figure_production.md)。
该页维护输入资格、两阶段交付、尺寸/字体/图例要求与验收契约；本技能只组织执行。

1. 确定图族、冻结输入、显示变换、输出尺寸和用户已授权的交付阶段。新图使用
   `candidate_aps_v2.toml`；具备 strict 输入资格时使用 `strict_aps_v2.toml`。
2. 复用已有图族生成器，以及 `scripts/plotting/` 的样式、系列身份、manifest、
   质量和交付 helper。接口见[公共样式 API](../../../docs/api/plotting/PlotStyle.md)。
   显示变换写入 manifest；需要数值重算时按实际计算任务处理。
3. 图像和唯一 `plot_manifest.json` 放在 `data/outputs/figures/`；数据、图注和
   必需的交付记录放在 `data/outputs/results/`。多图 case 复用 bundle 的共享来源。
   已接受的 case 使用新版本；临时失败或未接受的草稿可在确认无消费者后清理重试。
4. 按 PNG 审阅与 PDF 交付契约推进。已有明确接受及交付授权时继续执行；
   尚缺作者视觉判断时提交实际 PNG。作者接受记录绑定实际图件，后处理质量检查
   不能代替科学输入资格。
5. 运行 `validate_plot_artifact.py` 核对实际输出及来源。历史快照使用
   `--snapshot --code-ref <full-sha>`。仅检查本批输入、图件和受影响的合同；
   共享代码变化才运行对应公共测试，图族逻辑变化运行该图族测试。

汇报图件位置、已完成阶段、验证结果及确实剩余的交付项。更新现有入口即可，
不为每次重试另建状态说明、源码副本或清理回执。
