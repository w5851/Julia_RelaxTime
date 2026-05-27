## 背景与动机

- 解决了什么问题？为什么要改？

## 改动内容

- [ ] 功能/修复：
- [ ] 重构：
- [ ] 文档：
- [ ] 性能：

## 验证方式

- [ ] 分层 smoke / 单文件测试（请注明命令）
- [ ] 治理检查（如 `check_docs_consistency.jl`、`check_script_entrypoints.jl`、`check_models_entry_contract.jl`）
- [ ] 运行脚本或产物验证（请注明命令与输出位置）
- [ ] `Pkg.test()` / core/full profile（如适用）
- [ ] 其他：

命令/日志（可选）：

```text
(在这里粘贴关键命令与摘要输出)
```

## 影响范围与兼容性

- 是否改变稳定入口、脚本契约、`Models` 公共 API 或前后端数据契约？如果是，是否更新了 `docs/api/` / `docs/guides/` / `docs/dev/`？
- 是否涉及单位/公式变更？如果是，请给出参考或推导链接。
- 是否可能引入数值漂移？如果是，是否运行或更新了对应 regression/validation 覆盖？

## Checklist

- [ ] 代码风格与命名符合 `.github/copilot-instructions.md`
- [ ] 对外行为改变已更新文档
- [ ] 新增/修改逻辑有对应测试或说明
- [ ] 未改变 mixed-meson 治理与 non-fixedmu 统一 joint-solve 语义，或已在 PR 中明确说明并补充回归证据
