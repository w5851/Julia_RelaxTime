# Git分支策略

## 概述

本文档定义QCD模型库项目的Git分支管理策略，包括分支命名规范、工作流程、合并策略和发布流程。

**策略类型**：Feature Branch Workflow（特性分支工作流）

---

## 1. 分支类型

### 1.1 主要分支

#### `main` - 主分支
- **用途**：生产就绪代码，始终可发布
- **保护**：受保护，禁止直接推送
- **合并来源**：仅接受来自`develop`的PR（发布时）
- **标签**：所有发布版本在此打标签（v1.0.0, v1.1.0等）

#### `develop` - 开发分支
- **用途**：集成分支，包含下一版本的最新开发代码
- **保护**：受保护，禁止直接推送
- **合并来源**：接受来自特性分支、修复分支的PR
- **CI要求**：所有测试必须通过才能合并

### 1.2 临时分支

#### Feature Branches - 特性分支
- **命名**：`feature/<feature-name>`
- **示例**：`feature/qcd-model-refactoring`, `feature/rpnjl-implementation`
- **来源**：从`develop`创建
- **合并到**：`develop`
- **生命周期**：特性完成后删除

#### Bugfix Branches - 错误修复分支
- **命名**：`fix/<bug-description>`
- **示例**：`fix/vandermonde-negative-value`, `fix/convergence-tolerance`
- **来源**：从`develop`创建
- **合并到**：`develop`
- **生命周期**：修复完成后删除

#### Hotfix Branches - 紧急修复分支
- **命名**：`hotfix/<version>-<description>`
- **示例**：`hotfix/v1.2.1-critical-nan-bug`
- **来源**：从`main`创建
- **合并到**：`main` 和 `develop`
- **生命周期**：修复完成后删除
- **用途**：修复生产环境的紧急问题

#### Release Branches - 发布分支
- **命名**：`release/v<version>`
- **示例**：`release/v1.2.0`
- **来源**：从`develop`创建
- **合并到**：`main` 和 `develop`
- **生命周期**：发布完成后删除
- **用途**：准备发布，只允许bug修复和文档更新

---

## 2. 分支命名规范

### 2.1 命名格式

```
<type>/<description>
```

**类型前缀**：
- `feature/` - 新功能
- `fix/` - Bug修复
- `hotfix/` - 紧急修复
- `release/` - 发布准备
- `docs/` - 纯文档更新
- `refactor/` - 代码重构
- `test/` - 测试相关
- `perf/` - 性能优化

**描述部分**：
- 使用小写字母
- 使用连字符分隔单词
- 简洁但描述性强
- 避免使用特殊字符

### 2.2 命名示例

✅ **好的命名**：
```
feature/qcd-model-architecture
feature/rpnjl-eight-quark-term
fix/polyakov-potential-nan
fix/integration-convergence
hotfix/v1.2.1-memory-leak
release/v2.0.0
docs/api-documentation-update
refactor/thermodynamics-module
test/property-based-tests
perf/optimize-integral-calculation
```

❌ **不好的命名**：
```
feature/new-stuff              # 太模糊
fix/bug                        # 不描述性
my-branch                      # 缺少类型前缀
feature/QCD_Model_Refactoring  # 使用大写和下划线
fix/issue-123                  # 应该描述问题而非issue号
```

---

## 3. 工作流程

### 3.1 开发新功能

```bash
# 1. 确保develop是最新的
git checkout develop
git pull origin develop

# 2. 创建特性分支
git checkout -b feature/my-new-feature

# 3. 开发和提交
git add .
git commit -m "feat: implement core functionality"
git commit -m "test: add unit tests for new feature"
git commit -m "docs: update API documentation"

# 4. 定期同步develop的更新
git fetch origin
git rebase origin/develop

# 5. 推送到远程
git push origin feature/my-new-feature

# 6. 创建Pull Request到develop
# （在GitHub/GitLab界面操作）

# 7. PR合并后，删除本地分支
git checkout develop
git pull origin develop
git branch -d feature/my-new-feature
```

### 3.2 修复Bug

```bash
# 1. 从develop创建修复分支
git checkout develop
git pull origin develop
git checkout -b fix/bug-description

# 2. 修复和测试
git add .
git commit -m "fix: resolve issue with calculation"
git commit -m "test: add regression test"

# 3. 推送并创建PR
git push origin fix/bug-description
# 创建PR到develop
```

### 3.3 紧急修复（Hotfix）

```bash
# 1. 从main创建hotfix分支
git checkout main
git pull origin main
git checkout -b hotfix/v1.2.1-critical-bug

# 2. 修复问题
git add .
git commit -m "fix: critical bug in production"

# 3. 推送并创建PR到main
git push origin hotfix/v1.2.1-critical-bug
# 创建PR到main

# 4. 合并到main后，也要合并到develop
git checkout develop
git pull origin develop
git merge hotfix/v1.2.1-critical-bug
git push origin develop

# 5. 在main上打标签
git checkout main
git pull origin main
git tag -a v1.2.1 -m "Hotfix release v1.2.1"
git push origin v1.2.1
```

### 3.4 发布流程

```bash
# 1. 从develop创建release分支
git checkout develop
git pull origin develop
git checkout -b release/v2.0.0

# 2. 更新版本号和文档
# 编辑 Project.toml, CHANGELOG.md等
git add .
git commit -m "chore: prepare release v2.0.0"

# 3. 最后的bug修复（如果需要）
git commit -m "fix: minor issue in release"

# 4. 合并到main
git checkout main
git pull origin main
git merge --no-ff release/v2.0.0
git tag -a v2.0.0 -m "Release version 2.0.0"
git push origin main
git push origin v2.0.0

# 5. 合并回develop
git checkout develop
git pull origin develop
git merge --no-ff release/v2.0.0
git push origin develop

# 6. 删除release分支
git branch -d release/v2.0.0
git push origin --delete release/v2.0.0
```

---

## 4. Commit消息规范

### 4.1 Conventional Commits格式

```
<type>(<scope>): <subject>

<body>

<footer>
```

**类型（type）**：
- `feat`: 新功能
- `fix`: Bug修复
- `docs`: 文档更新
- `style`: 代码格式（不影响功能）
- `refactor`: 重构
- `perf`: 性能优化
- `test`: 测试相关
- `chore`: 构建/工具配置

**范围（scope）**（可选）：
- `pnjl`: PNJL模型相关
- `rpnjl`: rPNJL模型相关
- `relaxtime`: 弛豫时间模块
- `integration`: 积分工具
- `api`: API接口
- `docs`: 文档

**主题（subject）**：
- 简洁描述（≤50字符）
- 使用祈使句（"add"而非"added"）
- 不以句号结尾

### 4.2 Commit消息示例

✅ **好的commit消息**：
```
feat(pnjl): implement PNJLModel with abstract interface

Add PNJLModel struct that inherits from AbstractIsotropicModel.
Implement all required methods: calculate_mass_vec, calculate_chiral,
polyakov_potential, and dispersion_relation.

Closes #123
```

```
fix(integration): handle convergence failure gracefully

Add fallback mechanism when integration does not converge.
Increase node count automatically up to max_nodes limit.
Log warning when using fallback.

Fixes #456
```

```
docs(api): update data contracts documentation

Add detailed parameter validation rules and conversion examples.
Include unit conversion formulas and edge case handling.
```

```
test(rpnjl): add property tests for eight-quark term

Verify that eight-quark term calculation satisfies physical constraints.
Test with random parameter generation using Supposition.jl.
```

❌ **不好的commit消息**：
```
update code                    # 太模糊
fixed bug                      # 没有说明什么bug
WIP                           # 不应该推送WIP commit
Added new feature.            # 使用过去时，有句号
```

### 4.3 关联Issue

在commit消息或PR描述中关联Issue：

```
Closes #123          # 关闭issue
Fixes #456           # 修复issue
Resolves #789        # 解决issue
Relates to #101      # 相关issue（不关闭）
```

---

## 5. Pull Request流程

### 5.1 创建PR

**PR标题格式**：
```
<type>: <description>
```

**PR描述模板**：
```markdown
## 描述
简要说明此PR的目的和改动内容。

## 改动类型
- [ ] 新功能
- [ ] Bug修复
- [ ] 重构
- [ ] 文档更新
- [ ] 性能优化
- [ ] 测试

## 相关Issue
Closes #123

## 改动清单
- 实现了XXX功能
- 修复了YYY问题
- 更新了ZZZ文档

## 测试
- [ ] 所有现有测试通过
- [ ] 添加了新的测试
- [ ] 手动测试通过

## 检查清单
- [ ] 代码遵循项目规范
- [ ] 添加/更新了文档
- [ ] 添加/更新了测试
- [ ] 更新了CHANGELOG.md
- [ ] CI检查全部通过

## 截图/示例（如适用）
```

### 5.2 PR审查要求

**必需审查者**：
- 至少1名核心开发者审查
- 涉及算法变更：需要物理学家审查
- 涉及性能关键代码：需要性能专家审查

**审查检查点**：
- [ ] 代码质量和可读性
- [ ] 测试覆盖率充分
- [ ] 文档完整且准确
- [ ] 性能影响可接受
- [ ] 无明显的安全问题
- [ ] 符合项目架构设计

### 5.3 合并策略

**推荐策略**：Squash and Merge（压缩合并）

**原因**：
- 保持主分支历史清晰
- 每个特性一个commit
- 便于回滚和cherry-pick

**何时使用Merge Commit**：
- 发布分支合并到main
- 长期特性分支的重要里程碑

**何时使用Rebase**：
- 个人特性分支同步develop更新
- 清理本地commit历史

---

## 6. 分支保护规则

### 6.1 `main`分支保护

- ✅ 要求PR审查（至少1人）
- ✅ 要求CI通过
- ✅ 要求分支是最新的
- ✅ 禁止强制推送
- ✅ 禁止删除
- ✅ 要求签名commit（可选）

### 6.2 `develop`分支保护

- ✅ 要求PR审查（至少1人）
- ✅ 要求CI通过
- ✅ 要求分支是最新的
- ✅ 禁止强制推送
- ❌ 允许管理员绕过（紧急情况）

### 6.3 特性分支

- ❌ 无保护规则
- ✅ 允许强制推送（rebase时）
- ✅ 完成后自动删除

---

## 7. 版本标签规范

### 7.1 语义化版本

格式：`v<major>.<minor>.<patch>`

**示例**：
- `v1.0.0` - 首次正式发布
- `v1.1.0` - 新增功能（向后兼容）
- `v1.1.1` - Bug修复
- `v2.0.0` - 破坏性变更

**版本号递增规则**：
- **Major**：不兼容的API变更
- **Minor**：向后兼容的新功能
- **Patch**：向后兼容的Bug修复

### 7.2 预发布版本

格式：`v<major>.<minor>.<patch>-<pre-release>`

**示例**：
- `v2.0.0-alpha.1` - Alpha版本
- `v2.0.0-beta.1` - Beta版本
- `v2.0.0-rc.1` - Release Candidate

### 7.3 标签创建

```bash
# 创建带注释的标签
git tag -a v1.2.0 -m "Release version 1.2.0

- Add rPNJL model implementation
- Improve numerical stability
- Update documentation
"

# 推送标签
git push origin v1.2.0

# 推送所有标签
git push origin --tags
```

---

## 8. 常见场景

### 8.1 同步远程更新

```bash
# 更新本地develop
git checkout develop
git pull origin develop

# 将develop的更新合并到特性分支
git checkout feature/my-feature
git rebase develop

# 如果有冲突，解决后继续
git add .
git rebase --continue

# 强制推送（因为rebase改变了历史）
git push origin feature/my-feature --force-with-lease
```

### 8.2 撤销错误的commit

```bash
# 撤销最后一次commit（保留改动）
git reset --soft HEAD~1

# 撤销最后一次commit（丢弃改动）
git reset --hard HEAD~1

# 修改最后一次commit消息
git commit --amend -m "new message"

# 交互式rebase清理历史
git rebase -i HEAD~3
```

### 8.3 Cherry-pick特定commit

```bash
# 从其他分支挑选commit
git checkout target-branch
git cherry-pick <commit-hash>

# 挑选多个commit
git cherry-pick <commit1> <commit2>

# 挑选一个范围
git cherry-pick <start-commit>^..<end-commit>
```

### 8.4 解决合并冲突

```bash
# 1. 尝试合并/rebase
git merge develop
# 或
git rebase develop

# 2. 查看冲突文件
git status

# 3. 手动编辑冲突文件，解决冲突标记
# <<<<<<< HEAD
# 当前分支的内容
# =======
# 要合并的内容
# >>>>>>> develop

# 4. 标记为已解决
git add <resolved-file>

# 5. 继续合并/rebase
git merge --continue
# 或
git rebase --continue
```

---

## 9. CI/CD集成

### 9.1 CI检查

每个PR必须通过以下CI检查：

1. **代码格式检查**
   - Julia代码风格检查
   - 文档格式检查

2. **测试套件**
   - 单元测试
   - 集成测试
   - 属性测试
   - 回归测试

3. **性能基准**
   - 关键路径性能测试
   - 与基线对比

4. **文档构建**
   - API文档生成
   - 示例代码验证

5. **依赖检查**
   - 依赖图生成
   - 循环依赖检测

### 9.2 自动化操作

- **PR创建时**：自动运行CI检查
- **PR更新时**：重新运行CI检查
- **合并到develop**：自动部署到测试环境
- **合并到main**：自动创建GitHub Release
- **标签推送**：自动构建和发布文档

---

## 10. 最佳实践

### 10.1 Commit频率

- ✅ 频繁commit（逻辑单元完成即commit）
- ✅ 每个commit应该是可编译的
- ✅ 每个commit应该通过测试
- ❌ 不要commit半成品代码到共享分支

### 10.2 分支生命周期

- ✅ 特性分支应该短命（<1周）
- ✅ 定期同步develop的更新
- ✅ 完成后立即删除
- ❌ 不要长期保留未合并的分支

### 10.3 代码审查

- ✅ 及时响应审查意见
- ✅ 小的PR更容易审查（<500行）
- ✅ 提供清晰的PR描述
- ❌ 不要在PR中混入无关改动

### 10.4 冲突处理

- ✅ 优先使用rebase保持历史清晰
- ✅ 及时解决冲突，不要拖延
- ✅ 冲突解决后重新测试
- ❌ 不要盲目接受一方的改动

---

## 11. 故障排查

### 11.1 分支落后太多

```bash
# 方案1：Rebase（推荐）
git checkout feature/my-feature
git fetch origin
git rebase origin/develop

# 方案2：Merge
git merge origin/develop
```

### 11.2 误删分支

```bash
# 查找分支的最后commit
git reflog

# 恢复分支
git checkout -b feature/my-feature <commit-hash>
```

### 11.3 推送被拒绝

```bash
# 原因：远程分支有更新
# 解决：先拉取再推送
git pull --rebase origin feature/my-feature
git push origin feature/my-feature
```

---

## 12. 检查清单

开始新特性前：
- [ ] 从最新的develop创建分支
- [ ] 分支命名符合规范
- [ ] 了解相关的Issue和需求

开发过程中：
- [ ] 频繁commit，消息清晰
- [ ] 定期同步develop更新
- [ ] 运行本地测试
- [ ] 更新相关文档

创建PR前：
- [ ] 所有测试通过
- [ ] 代码已自我审查
- [ ] 文档已更新
- [ ] CHANGELOG已更新
- [ ] PR描述完整

---

## 13. 参考资料

- Git Flow: https://nvie.com/posts/a-successful-git-branching-model/
- Conventional Commits: https://www.conventionalcommits.org/
- Semantic Versioning: https://semver.org/
- GitHub Flow: https://guides.github.com/introduction/flow/
