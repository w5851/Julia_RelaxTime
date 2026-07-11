---
name: experiment-logbook-append
description: 以 append-only 方式记录 Julia_RelaxTime 参数扫描、补点、rerun、异常复测和分析输入输出，保留可复现实验事实且不重写历史结论。仅用于科研实验日志；开发任务执行台账使用 doc-exec-log-append。
---

# Experiment Logbook Append

## Hard rules

- 只追加，不重写、整理或升级历史结论。
- 每条记录自包含命令、配置、输入、输出、观察和下一步。
- 把推测标为推测，不写成已验证机制。
- 输出路径必须具体到文件或目录。

## Destinations

- 分析实验：`docs/analysis/<domain>/`
- 研究/开发混合任务：`docs/dev/active/` 下的专用实验日志
- 临时 rerun：与目标分析包同目录

## Entry contract

- 日期和实验目的
- 实际命令与环境差异
- 参数集、preset 或配置文件
- 结果、图表和诊断路径
- 关键观察与异常
- 下一步补点、回归、绘图或复核动作

## Workflow

1. 复用同一主题的明确日志目标。
2. 收集足以复现实验的事实。
3. 追加一个独立记录块，不读取或改写旧结论。
4. 核对命令、输入和输出路径的一一对应。
5. 以具体下一步结束记录。

使用 `analysis-artifact` 将成熟实验证据整理成正式分析包；使用 `formal-production-artifact` 晋升正式数值产物。
