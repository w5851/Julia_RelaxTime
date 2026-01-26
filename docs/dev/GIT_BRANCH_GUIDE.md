# Git 分支使用指南

## 当前分支结构

### main 分支 (生产环境)
- **用途**: 稳定的生产代码和重要文档
- **内容**: 只包含归档的调查结论文档
- **提交**: 只提交经过验证的、需要长期保留的内容

### investigation/ssbar-uubar-convergence-2026-01-27 分支
- **用途**: 保存完整的调查过程
- **内容**: 
  - 所有测试脚本 (step1_compare_sigma_s.jl, step2_test_convergence.jl 等)
  - 所有分析文档 (SYSTEMATIC_INVESTIGATION_PLAN.md 等)
  - 所有中间结果和验证脚本
- **提交**: 包含完整的调查历史

---

## 常用操作

### 1. 查看所有分支

```bash
git branch
```

### 2. 切换到调查分支查看完整内容

```bash
git checkout investigation/ssbar-uubar-convergence-2026-01-27
```

### 3. 切换回 main 分支

```bash
git checkout main
```

### 4. 查看分支之间的差异

```bash
# 查看文件列表差异
git diff main investigation/ssbar-uubar-convergence-2026-01-27 --name-only

# 查看具体内容差异
git diff main investigation/ssbar-uubar-convergence-2026-01-27
```

### 5. 从调查分支复制特定文件到 main

```bash
# 在 main 分支上执行
git checkout investigation/ssbar-uubar-convergence-2026-01-27 -- path/to/file
```

### 6. 推送分支到远程仓库

```bash
# 推送 main 分支
git push origin main

# 推送调查分支
git push origin investigation/ssbar-uubar-convergence-2026-01-27
```

---

## 工作流程示例

### 场景 1: 需要查看完整的调查过程

```bash
# 1. 切换到调查分支
git checkout investigation/ssbar-uubar-convergence-2026-01-27

# 2. 查看所有文件
ls scripts/relaxtime/

# 3. 运行测试脚本
julia scripts/relaxtime/step2_test_convergence.jl

# 4. 完成后切换回 main
git checkout main
```

### 场景 2: 需要将调查分支的某个文件添加到 main

```bash
# 1. 确保在 main 分支
git checkout main

# 2. 从调查分支复制文件
git checkout investigation/ssbar-uubar-convergence-2026-01-27 -- scripts/relaxtime/INVESTIGATION_COMPLETE.md

# 3. 提交
git add scripts/relaxtime/INVESTIGATION_COMPLETE.md
git commit -m "添加调查完成总结"
```

### 场景 3: 继续在调查分支上工作

```bash
# 1. 切换到调查分支
git checkout investigation/ssbar-uubar-convergence-2026-01-27

# 2. 进行修改和测试
# ... 你的工作 ...

# 3. 提交改动
git add .
git commit -m "进一步的分析和测试"

# 4. 切换回 main
git checkout main
```

---

## 分支管理最佳实践

### 命名规范

- **功能分支**: `feature/功能名称`
- **调查分支**: `investigation/调查主题-日期`
- **修复分支**: `fix/问题描述`
- **实验分支**: `experiment/实验名称`

### 提交信息规范

- **main 分支**: 简洁明了,说明主要改动
- **调查分支**: 详细记录,包含背景和发现

### 何时合并分支

- **不合并**: 如果调查分支包含大量临时文件和实验代码
- **选择性合并**: 只将重要的结论和有用的脚本合并到 main
- **完全合并**: 如果调查分支的所有内容都有价值

---

## 当前项目的分支策略

### main 分支
- ✅ 归档的调查结论 (docs/dev/archived/)
- ✅ 稳定的源代码
- ✅ 文档和指南
- ❌ 临时测试脚本
- ❌ 中间分析文档

### investigation 分支
- ✅ 所有测试脚本
- ✅ 所有分析文档
- ✅ 中间结果
- ✅ 实验代码
- ✅ 完整的调查历史

---

## 查看本次调查的内容

### 在 main 分支上 (当前)

```bash
ls docs/dev/archived/2026_01_27_*
```

输出:
- `2026_01_27_SSBAR_UUBAR_INVESTIGATION_COMPLETE.md` - 完整结论
- `2026_01_27_STEP2_CONVERGENCE_TEST_RESULTS.md` - 收敛性测试
- `2026_01_27_WHY_ONLY_SSBAR_UUBAR_DIFFERS.md` - 差异原因

### 在调查分支上

```bash
git checkout investigation/ssbar-uubar-convergence-2026-01-27
ls scripts/relaxtime/
```

输出: 117 个文件,包括所有测试脚本和分析文档

---

## 推送到 GitHub

### 推送 main 分支 (推荐)

```bash
git push origin main
```

这会将归档的三个文档推送到远程仓库。

### 推送调查分支 (可选)

```bash
git push origin investigation/ssbar-uubar-convergence-2026-01-27
```

这会将完整的调查过程推送到远程仓库,供将来参考。

**建议**: 两个分支都推送,这样:
- main 分支保持整洁
- 调查分支保留完整历史
- 可以随时切换查看

---

## 删除分支 (谨慎操作)

### 删除本地分支

```bash
# 确保不在要删除的分支上
git checkout main

# 删除分支
git branch -d investigation/ssbar-uubar-convergence-2026-01-27
```

### 删除远程分支

```bash
git push origin --delete investigation/ssbar-uubar-convergence-2026-01-27
```

**注意**: 删除前请确保已经保存了所有需要的内容!

---

*创建时间: 2026-01-27*
*用途: 帮助理解和管理 Git 分支*
