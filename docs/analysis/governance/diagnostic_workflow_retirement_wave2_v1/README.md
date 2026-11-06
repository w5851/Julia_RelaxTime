# 历史手动 workflow 定义

[PR #298](https://github.com/w5851/Julia_RelaxTime/pull/298) 在 2026-08-31
退役了本目录的 16 个一次性 workflow。原始定义保留在 `definitions/`，
[`retired_workflows.csv`](retired_workflows.csv) 索引其用途、来源、历史运行及 hash。
它们用于读取历史结果的执行条件。

当时保留的 Maxwell-local expansion 入口及其专用执行脚本，在桌面分支收尾中
一并退役。历史实现与原退役记录可从提交
`5946aaf73a4d26f0f392bc7d43b49419f0ee6f26` 恢复。

当前测试覆盖仍使用的收集器、输入校验和数值方法，不再逐条断言退役 YAML 的
文字、入口数量或历史清理布尔值。后续确有计算需求时，从所需算法和输入范围
建立任务即可。
