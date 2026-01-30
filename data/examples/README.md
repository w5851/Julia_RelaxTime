示例数据与输出存放说明

本目录 `data/examples/` 用于保存用于演示与测试的 CSV 文件及由脚本生成的示例图像（PDF/PNG）。为便于管理，约定将数据与生成的图像分别放在子目录中：

- 示例数据（CSV 等）：`data/examples/results/`
- 示例图像（PDF/PNG 等导出结果）：`data/examples/figures/`

约定：
- 示例 CSV 文件命名建议为 `*_example.csv` 或 `*_sample.csv`，放在 `data/examples/results/`。
- 输出子目录建议为 `data/examples/figures/<script_name>/`，用于存放由脚本生成的 PDF/PNG 文件。

当前步骤：
- 若要复制当前临时示例到此目录，请运行：

```bash
mkdir -p data/examples/results data/examples/figures
cp data/tmp_test_scan.csv data/examples/results/tmp_test_scan.csv
# 可选：将 docs/dev/active/sample_output 下的 PDF/PNG 复制到 data/examples/figures
```

建议在仓库中保留一套小型示例（CSV + 生成的 PDF），以便新开发者快速验证脚本行为与 CI 测试。