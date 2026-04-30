---
name: experiment-logbook-append
description: 以 append-only 方式记录 Julia_RelaxTime 的实验、扫描、补点、复测与异常诊断，不重写历史结论。适用于参数扫描、阈值邻域补点、rerun 记录、输出路径留痕。关键词：experiment logbook, append-only, scan, rerun, analysis
---

# experiment-logbook-append

## Purpose

为研究实验而不是一般开发流程维护可追溯日志，重点记录命令、输入、输出路径、关键观察和下一步动作。

本 Skill 与 `doc-exec-log-append` 相邻，但场景更偏科研实验、扫描和分析记录。

## Apply When

- 记录参数扫描、补点、rerun、异常点复测。
- 记录图表生成前的实验输入与输出路径。
- 需要把一次实验结果以 append-only 方式写入日志。

## Prefer Another Skill When

- 目标是开发任务执行台账：优先 `doc-exec-log-append`。
- 目标是归档已完成任务文档：优先 `doc-archive`。

## Hard Rules

- 只追加，不回写历史结论。
- 每条记录必须自包含：命令、配置、输出、观察、下一步。
- 输出路径必须具体到文件或目录。
- 若记录依赖特定脚本或配置，必须写出其路径。

## Recommended Destinations

- 研究分析笔记：`docs/analysis/<domain>/`
- 研究/开发混合任务：`docs/dev/active/` 下专门日志
- 临时 rerun 说明：与相关分析文档同目录

## Standard Entry Template

- 日期：`YYYY-MM-DD`
- 任务：一句话说明本次实验目的
- 命令：实际执行命令
- 输入：
  - 参数集 / preset / 配置文件
- 输出：
  - 结果文件
  - 图表文件
- 观察：
  - 关键数值变化
  - 异常点或与预期不一致处
- 下一步：
  - 是否需要补点、回归、画图或复核

## Standard Workflow

1. Choose the correct log target
   - 优先复用同一主题的现有日志。
2. Gather reproducibility facts
   - 命令、配置、输出路径、时间范围、关键参数。
3. Append one self-contained block
   - 不改写历史正文，不整理旧批次。
4. End with the next action
   - 明确这条记录之后该做什么。

## Quality Checks

- 读者不回看上下文也能复现实验。
- 没有把“推测”写成“已证实结论”。
- 命令、输入和输出路径能一一对应。
