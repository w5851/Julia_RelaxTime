# Changelog

All notable changes to the QCD Model Library project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- 架构设计规范文档 (`docs/dev/active/2026_02_01_QCD模型库架构设计规范.md`)
- 开发实施计划文档 (`docs/dev/active/2026_02_01_QCD模型库开发实施计划.md`)
- 错误处理指南 (`docs/guides/error_handling.md`)
- 数据契约定义 (`docs/api/data_contracts.md`)
- Git分支策略文档 (`docs/dev/branching_strategy.md`)
- 依赖关系分析文档 (`docs/architecture/qcd_refactoring_dependencies.md`)
- 代码审查指南 (`docs/dev/code_review_guidelines.md`)
- 配置管理指南 (`docs/guides/configuration.md`)
- 架构决策记录（ADR）框架 (`docs/decisions/`)
- ADR-0001: 使用多重分派实现模型多态

### Changed
- 更新需求文档，增加需求13-15（错误处理、输入验证、日志监控）

### Planned
- 核心抽象接口实现 (`src/core/interface.jl`)
- 自定义异常类型定义 (`src/core/exceptions.jl`)
- 基础工具库 (`src/models/base/`)
- PNJL模型重构 (`src/models/isotropic/pnjl.jl`)
- rPNJL模型实现 (`src/models/anisotropic/rpnjl.jl`)
- 向后兼容层 (`src/compatibility/legacy_api.jl`)
- 模型工厂 (`src/models/factory.jl`)

---

## [0.1.0] - 2026-01-26

### Added
- 参数结构体系统 (`src/ParameterTypes.jl`)
  - `QuarkParams` 结构体
  - `ThermoParams` 结构体
  - 参数格式转换工具
- RelaxTime模块的结构体支持
  - `RelaxationTime.jl`
  - `AverageScatteringRate.jl`
  - `TotalCrossSection.jl`
  - `DifferentialCrossSection.jl`
  - `ScatteringAmplitude.jl`
  - `TotalPropagator.jl`
  - `ParticleSymbols.jl`
- 完整的测试套件 (`tests/unit/struct_migration/`)
  - 单元测试
  - 属性测试（Property-Based Tests）
  - 边缘情况测试
  - 集成测试
  - 向后兼容性测试
- 性能基准测试 (`tests/perf/relaxtime/benchmark_struct_vs_namedtuple.jl`)
- 文档
  - 参数类型API文档 (`docs/api/PARAMETER_TYPES_API.md`)
  - 参数结构体迁移指南 (`docs/guides/PARAMETER_STRUCT_MIGRATION.md`)
  - 使用示例 (`docs/guides/examples/`)

### Changed
- RelaxTime模块函数签名支持多种参数格式（Struct、NamedTuple、Dict）
- 内部实现使用NamedTuple作为统一格式

### Performance
- 结构体参数与NamedTuple性能相当（差异<5%）
- 保持了向后兼容性，无性能回归

---

## [0.0.1] - 2025-12-01

### Added
- 初始项目结构
- PNJL模型核心实现
  - 热力学计算 (`src/pnjl/core/Thermodynamics.jl`)
  - 数值积分 (`src/pnjl/core/Integrals.jl`)
  - Gap方程求解器 (`src/pnjl/solver/`)
- 扫描功能
  - T-μ扫描 (`src/pnjl/scans/TmuScan.jl`)
  - T-ρ扫描 (`src/pnjl/scans/TrhoScan.jl`)
  - 双分支扫描 (`src/pnjl/scans/DualBranchScan.jl`)
- 弛豫时间计算模块 (`src/relaxtime/`)
  - 散射振幅计算
  - 总截面计算
  - 平均散射率计算
  - 弛豫时间计算
- 基础测试套件 (`tests/unit/`)
- 基础文档 (`docs/`)

---

## 版本说明

### 版本号格式

本项目使用[语义化版本](https://semver.org/)：`MAJOR.MINOR.PATCH`

- **MAJOR**：不兼容的API变更
- **MINOR**：向后兼容的新功能
- **PATCH**：向后兼容的bug修复

### 变更类型

- **Added**：新增功能
- **Changed**：现有功能的变更
- **Deprecated**：即将移除的功能
- **Removed**：已移除的功能
- **Fixed**：bug修复
- **Security**：安全相关的修复
- **Performance**：性能改进

### Breaking Changes标记

破坏性变更使用 **⚠️ BREAKING** 标记：

```markdown
### Changed
- ⚠️ BREAKING: 函数签名变更，旧代码需要更新
```

---

## 贡献指南

更新CHANGELOG时请遵循以下规则：

1. **及时更新**：每个PR合并时更新CHANGELOG
2. **放在Unreleased下**：新变更放在`[Unreleased]`部分
3. **使用正确的类型**：选择合适的变更类型（Added、Changed等）
4. **清晰描述**：简洁但完整地描述变更内容
5. **链接PR/Issue**：如适用，链接到相关的PR或Issue
6. **标记Breaking Changes**：使用⚠️ BREAKING标记破坏性变更

### 示例

```markdown
## [Unreleased]

### Added
- 新功能描述 (#123)

### Changed
- ⚠️ BREAKING: API变更描述 (#124)
- 非破坏性变更描述 (#125)

### Fixed
- Bug修复描述 (#126)
```

---

## 发布流程

发布新版本时：

1. 将`[Unreleased]`下的内容移动到新版本号下
2. 添加发布日期
3. 创建新的空`[Unreleased]`部分
4. 更新版本链接
5. 提交并打标签

### 示例

```markdown
## [Unreleased]

## [1.0.0] - 2026-02-15

### Added
- 功能列表...

### Changed
- 变更列表...
```

---

## 参考资料

- [Keep a Changelog](https://keepachangelog.com/)
- [Semantic Versioning](https://semver.org/)
- 分支策略：`docs/dev/branching_strategy.md`
- 代码审查指南：`docs/dev/code_review_guidelines.md`

---

**维护者**：QCD开发团队  
**最后更新**：2026-02-01
