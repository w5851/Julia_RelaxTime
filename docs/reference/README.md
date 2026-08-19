# 参考资料总入口

`docs/reference/` 保存公式、领域知识、数值方法与实现能力之间的可追溯参考资料。需要先判断“当前仓库已经实现什么、如何计算、证据到哪一层”时，应先阅读能力总账，再按主题进入公式或 API 文档。

公式块采用 GitHub 兼容的 `$$ ... $$` 分隔符。若某个 Markdown 阅读器只显示原始 TeX 文本，说明该阅读器没有启用数学排版插件；公式源仍然保留在文档中，可在 GitHub 或支持 KaTeX/MathJax 的预览器中查看排版结果。

## 首读：实现能力与方法总账

- [已实现计算能力与方法追踪清单](implemented_capabilities.md)：按计算流程拆分各条路线，登记公式、信源、源码/API 入口、测试证据与验证边界。

## 公式索引

- [Models 公式索引](formula/models/README.md)：NJL/PNJL/rPNJL、磁场、旋转、GasLiquid/RMF 与统一平均场热力学。
- [relaxtime 公式索引](formula/relaxtime/README.md)：极化函数、传播子、Mott、介子密度/热力学、散射与输运。
- `formula/`：逐主题公式页及其文献映射。

## 其他参考资料

- `domain-knowledge/`：单位、运动学、数值策略与领域推导笔记。
- `mathematica/`：解析积分、符号推导与归一化核对记录。
- `usage/`：面向具体计算任务的使用说明。

## 建议阅读顺序

1. 先读能力总账，确认路线、成熟度和非目标边界。
2. 再读对应的 `docs/api/` 入口，确认调用合同和输出字段。
3. 需要复核公式时进入 `formula/models/` 或 `formula/relaxtime/`。
4. 需要运行或晋升数值产物时，再读取 `docs/guides/` 与相应 SOP。

## 维护规则

新增稳定入口、公式语义、数值策略或验证等级时，先更新权威专题文档，再同步 [能力总账](implemented_capabilities.md) 的路线、入口、证据和边界。不要只在总览表新增一个名称，也不要把内部回归或诊断候选写成外部验证或 formal production。
