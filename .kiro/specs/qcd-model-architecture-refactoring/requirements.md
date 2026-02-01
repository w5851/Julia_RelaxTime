# 需求文档

## 引言

本文档规定了将QCD模型库架构从单体PNJL实现重构为多态、可扩展架构的需求，该架构支持多种QCD模型（PNJL、rPNJL、磁场PNJL等），同时保持向后兼容性和数值精度。

重构将在模型接口、基础工具和具体实现之间建立清晰的分离，通过Julia的多重分派系统实现新模型的轻松添加和代码复用的最大化。

## 术语表

- **QCD_System**: 正在重构的量子色动力学模型库系统
- **PNJL_Model**: Polyakov环扩展的Nambu-Jona-Lasinio模型（现有实现）
- **rPNJL_Model**: 带有八夸克相互作用项的排斥PNJL模型
- **Abstract_Interface**: 定义所有模型必需方法的Julia抽象类型
- **Base_Utilities**: 共享的数值工具（积分、Polyakov计算、安全数学运算）
- **Concrete_Model**: 具体的模型实现（例如PNJLModel、RPNJLModel）
- **Legacy_API**: 当前脚本使用的现有公共接口
- **Compatibility_Layer**: 将Legacy_API调用转换为新架构的代码
- **Vandermonde_Term**: Polyakov势中的数学项（式3.27-3.29），具有b2(T)温度依赖性
- **Eight_Quark_Term**: rPNJL质量方程中的排斥相互作用项（式3.31）
- **Property_Test**: 验证生成输入上通用属性的自动化测试
- **Numerical_Equivalence**: 重构代码产生与原始实现相同结果的要求

## 需求

### 需求1：核心架构

**用户故事：** 作为开发者，我希望有一个带有抽象接口的多态架构，以便我可以在不修改现有代码的情况下添加新的QCD模型。

#### 验收标准

1. THE QCD_System SHALL 定义AbstractQCDModel类型作为所有模型的根抽象类型
2. THE QCD_System SHALL 定义从AbstractQCDModel继承的领域抽象类型（AbstractIsotropicModel、AbstractAnisotropicModel、AbstractMagneticModel）
3. THE QCD_System SHALL 要求所有Concrete_Models实现以下方法：vacuum_contribution、thermal_contribution、polyakov_potential、dispersion_relation、calculate_mass_vec、calculate_chiral
4. WHEN Concrete_Model未实现必需方法时，THEN THE QCD_System SHALL 抛出带有描述性消息的MethodError
5. THE QCD_System SHALL 以分层结构组织代码：core/interface.jl → models/{domain}/types.jl → models/{domain}/{model}.jl

### 需求2：基础工具和代码复用

**用户故事：** 作为开发者，我希望有共享的数值工具，以便我可以避免跨模型的代码重复。

#### 验收标准

1. THE QCD_System SHALL 提供接受能量函数回调的数值积分Base_Utilities
2. THE QCD_System SHALL 提供所有模型共享的Polyakov对数计算Base_Utilities
3. THE QCD_System SHALL 提供安全数值运算的Base_Utilities（safe_log、safe_vandermonde）
4. THE QCD_System SHALL 提供夸克和反夸克分布函数的Base_Utilities
5. WHEN Concrete_Model需要数值积分时，THEN THE Concrete_Model SHALL 使用Base_Utilities而不是实现自定义积分

### 需求3：PNJL模型重构

**用户故事：** 作为物理学家，我希望将现有的PNJL实现重构到新架构中，以便我可以受益于改进的结构同时保持数值精度。

#### 验收标准

1. THE QCD_System SHALL 将现有PNJL实现封装为继承自AbstractIsotropicModel的PNJLModel类型
2. WHEN PNJLModel计算热力学量时，THEN THE 结果SHALL与原始实现数值等价（在浮点精度范围内）
3. THE PNJLModel SHALL 实现所有必需的Abstract_Interface方法
4. THE PNJLModel SHALL 使用Base_Utilities进行积分和Polyakov计算
5. THE QCD_System SHALL 提供Compatibility_Layer，允许Legacy_API调用无需修改即可工作

### 需求4：rPNJL模型实现

**用户故事：** 作为物理学家，我希望实现带有八夸克项和修改的Polyakov势的rPNJL模型，以便我可以研究QCD中的排斥相互作用。

#### 验收标准

1. THE QCD_System SHALL 定义继承自AbstractAnisotropicModel的RPNJLModel类型
2. THE RPNJLModel SHALL 支持额外参数：g1、g2、kappa、T0，并进行适当的单位转换
3. WHEN RPNJLModel计算夸克质量时，THEN THE 计算SHALL包含方程3.31中指定的Eight_Quark_Term
4. WHEN RPNJLModel计算Polyakov势时，THEN THE 计算SHALL包含方程3.27-3.29中指定的Vandermonde_Term
5. THE Vandermonde_Term计算SHALL包含b2(T)温度依赖性
6. WHEN Vandermonde_Term产生负值时，THEN THE QCD_System SHALL 记录警告而不自动修正值

### 需求5：安全数值计算

**用户故事：** 作为开发者，我希望有安全的数值运算，以便我可以优雅地处理边缘情况而不会崩溃。

#### 验收标准

1. THE QCD_System SHALL 提供处理零和负输入的safe_log函数
2. THE QCD_System SHALL 提供处理负判别式的safe_vandermonde函数
3. WHEN safe_vandermonde遇到负值时，THEN THE 函数SHALL记录警告并返回计算值
4. WHEN safe_log遇到零或负输入时，THEN THE 函数SHALL返回安全的回退值或抛出描述性错误
5. THE Base_Utilities SHALL 在所有关键计算中使用安全数值函数

### 需求6：向后兼容性

**用户故事：** 作为用户，我希望现有脚本继续工作，以便我不需要重写分析代码。

#### 验收标准

1. THE QCD_System SHALL 提供保留Legacy_API的Compatibility_Layer
2. WHEN 脚本调用Legacy_API函数时，THEN THE Compatibility_Layer SHALL 将调用转换为新架构
3. WHEN 脚本使用旧参数格式时，THEN THE QCD_System SHALL 接受新旧格式
4. THE QCD_System SHALL 维护所有现有的公共函数导出
5. WHEN 使用Compatibility_Layer时，THEN THE 结果SHALL与原始实现数值等价

### 需求7：性能要求

**用户故事：** 作为用户，我希望重构的代码保持性能，以便我的模拟不会变慢。

#### 验收标准

1. WHEN 执行关键计算路径时，THEN THE 性能SHALL相比原始实现不降低超过5%
2. THE QCD_System SHALL 对性能关键函数使用@inline提示
3. THE QCD_System SHALL 为重复计算返回闭包以启用JIT优化
4. THE QCD_System SHALL 利用Julia的多重分派进行类型特定优化
5. WHEN 运行性能测试时，THEN THE 结果SHALL被记录并与基线测量进行比较

### 需求8：测试和验证

**用户故事：** 作为开发者，我希望有全面的测试，以便我可以验证正确性并捕获回归。

#### 验收标准

1. THE QCD_System SHALL 包含所有Base_Utilities的单元测试
2. THE QCD_System SHALL 包含通用正确性属性的Property_Tests
3. THE QCD_System SHALL 包含比较重构PNJL与原始实现的回归测试
4. THE QCD_System SHALL 包含完整热力学计算的集成测试
5. WHEN 实现rPNJL时，THEN THE QCD_System SHALL 包含验证Eight_Quark_Term和Vandermonde_Term计算的测试
6. THE QCD_System SHALL 包含数值边缘情况的测试（零温度、极端化学势等）

### 需求9：模型工厂模式

**用户故事：** 作为用户，我希望有一种简单的方法来创建模型实例，以便我可以轻松地在模型之间切换。

#### 验收标准

1. THE QCD_System SHALL 提供接受模型名称符号的create_model工厂函数
2. WHEN 使用:PNJL调用create_model时，THEN THE 函数SHALL返回PNJLModel实例
3. WHEN 使用:rPNJL调用create_model时，THEN THE 函数SHALL返回RPNJLModel实例
4. THE create_model函数SHALL接受模型参数的关键字参数
5. WHEN 使用未知模型名称调用create_model时，THEN THE 函数SHALL抛出描述性错误

### 需求10：文档和示例

**用户故事：** 作为新开发者，我希望有清晰的文档和示例，以便我可以理解如何使用和扩展架构。

#### 验收标准

1. THE QCD_System SHALL 为所有Abstract_Interface方法提供API文档
2. THE QCD_System SHALL 提供展示如何实现新Concrete_Model的示例
3. THE QCD_System SHALL 提供展示如何使用工厂模式的示例
4. THE QCD_System SHALL 提供从Legacy_API过渡到新架构的迁移指南
5. THE QCD_System SHALL 记录Property_Tests中使用的所有正确性属性

### 需求11：可扩展性

**用户故事：** 作为研究人员，我希望轻松添加新模型，以便我可以探索不同的QCD理论。

#### 验收标准

1. WHEN 开发者创建新Concrete_Model时，THEN THE 开发者SHALL只需定义类型并实现必需方法
2. THE QCD_System SHALL 在添加新模型时不需要修改核心基础设施
3. WHEN 添加新Concrete_Model时，THEN THE 模型SHALL自动受益于Base_Utilities
4. THE QCD_System SHALL 支持领域特定优化的可选方法覆盖
5. WHEN 定义新模型类型时，THEN THE Julia编译器SHALL为该特定类型生成优化代码

### 需求12：参数管理

**用户故事：** 作为用户，我希望有灵活的参数处理，以便我可以使用基于结构体和基于字典的参数。

#### 验收标准

1. THE QCD_System SHALL 支持QuarkParams和ThermoParams结构体类型
2. THE QCD_System SHALL 支持NamedTuple参数格式
3. THE QCD_System SHALL 支持Dict{Symbol, Any}参数格式
4. WHEN Concrete_Model接收参数时，THEN THE 模型SHALL接受任何支持的参数格式
5. THE QCD_System SHALL 提供参数格式之间的转换工具

### 需求13：错误处理与降级策略

**用户故事：** 作为开发者，我希望有统一的错误处理机制，以便我可以优雅地处理异常情况并提供有用的错误信息。

#### 验收标准

1. THE QCD_System SHALL 定义自定义异常类型层次结构（ModelParameterError、NumericalConvergenceError、PhysicalConstraintError、ModelConfigurationError、NumericalInstabilityError）
2. THE QCD_System SHALL 为所有自定义异常提供描述性错误消息，包含问题描述、上下文信息和建议操作
3. WHEN 数值积分不收敛时，THEN THE QCD_System SHALL 实施降级策略（增加节点数或切换积分方法）
4. WHEN 求解器不收敛时，THEN THE QCD_System SHALL 尝试调整初值或放宽容差，最多尝试3次
5. WHEN Vandermonde项产生负值时，THEN THE QCD_System SHALL 记录WARN级别日志但继续计算
6. THE QCD_System SHALL 在各层明确异常处理责任：基础工具层捕获数值库异常，模型接口层转换为模型异常，高层API提供用户友好消息
7. THE QCD_System SHALL 提供异常情况的测试用例，覆盖所有定义的异常类型

### 需求14：输入验证与安全性

**用户故事：** 作为用户，我希望系统能够验证我的输入参数，以便我可以及早发现错误并避免无效计算。

#### 验收标准

1. THE QCD_System SHALL 定义所有参数的约束条件（范围、单位、物理意义）
2. THE QCD_System SHALL 在模型构造时验证参数（推荐）或在计算前验证（备选）
3. WHEN 参数超出有效范围时，THEN THE QCD_System SHALL 抛出ModelParameterError并说明约束条件
4. THE QCD_System SHALL 验证物理约束（如Polyakov loop在[0,1]范围内，夸克质量非负）
5. THE QCD_System SHALL 支持自定义验证函数用于复杂约束
6. THE QCD_System SHALL 确保模型实例的线程安全性（推荐使用immutable设计）
7. THE QCD_System SHALL 在日志中避免记录完整的大型数组，使用统计信息代替
8. THE QCD_System SHALL 提供参数验证的单元测试，覆盖边界值和无效值

### 需求15：日志与监控

**用户故事：** 作为开发者和用户，我希望有清晰的日志记录，以便我可以理解系统行为、调试问题和监控性能。

#### 验收标准

1. THE QCD_System SHALL 使用Julia标准Logging库，支持DEBUG、INFO、WARN、ERROR四个级别
2. THE QCD_System SHALL 在DEBUG级别记录详细的内部状态和中间计算结果
3. THE QCD_System SHALL 在INFO级别记录正常操作的关键步骤（如扫描开始、求解器收敛）
4. THE QCD_System SHALL 在WARN级别记录可恢复的异常情况和性能问题
5. THE QCD_System SHALL 在ERROR级别记录严重错误但程序可继续的情况
6. THE QCD_System SHALL 使用结构化日志格式：`@info "Operation" key1=value1 key2=value2`
7. THE QCD_System SHALL 为关键操作提供性能监控点（计算时间、迭代次数、收敛残差）
8. THE QCD_System SHALL 允许用户配置日志级别和输出目标（控制台、文件）
9. THE QCD_System SHALL 在文档中说明各日志级别的使用场景和典型消息

## 实施约束

### 技术约束

1. **Julia版本**：THE QCD_System SHALL 兼容Julia 1.9+（当前项目环境）
2. **依赖包**：THE QCD_System SHALL NOT 引入新的重型依赖包（如特殊积分库）
3. **内存使用**：WHEN 重构完成后，THEN THE 内存使用SHALL NOT 显著增加（增长<20%）
4. **编译时间**：WHEN 重构完成后，THEN THE 首次编译时间增加SHALL <30%
5. **数值精度**：THE QCD_System SHALL 使用Float64作为默认精度，支持可选的更高精度

### 业务约束

1. **时间线**：所有P0（优先级0）需求必须在3周内完成
2. **团队技能**：假设团队成员熟悉Julia多重分派但可能需要架构指导
3. **现有代码**：THE QCD_System SHALL 保持现有科学结果的完全可复现性
4. **测试覆盖率**：核心功能的测试覆盖率必须≥80%
5. **文档同步**：代码修改必须同步更新相关文档

### 性能约束

1. **关键路径性能**：热力学计算的关键路径性能降低必须<5%
2. **内存分配**：重复计算中的内存分配应最小化（使用闭包和预分配）
3. **类型稳定性**：所有性能关键函数必须保持类型稳定
4. **编译器优化**：必须使用@inline、@inbounds等提示以启用编译器优化

## 交付物

### 代码交付物

1. **核心接口** (`src/core/interface.jl`) - 完整的抽象类型和方法签名定义
2. **异常定义** (`src/core/exceptions.jl`) - 自定义异常类型层次结构
3. **领域抽象类型** (`src/models/{isotropic,anisotropic,magnetic}/types.jl`) - 领域特定的抽象类型和默认实现
4. **PNJL模型** (`src/models/isotropic/pnjl.jl`) - 重构的PNJL模型实现
5. **rPNJL模型** (`src/models/anisotropic/rpnjl.jl`) - 新的rPNJL模型实现
6. **基础工具** (`src/models/base/`) - 所有共享的数值工具模块
   - `integrals.jl` - 通用积分器
   - `polyakov.jl` - Polyakov势计算
   - `safe_math.jl` - 安全数值运算
   - `distributions.jl` - 分布函数
   - `validation.jl` - 参数验证工具
   - `mass_chiral.jl` - 质量和手征凝聚计算
7. **兼容层** (`src/compatibility/legacy_api.jl`) - 向后兼容的API包装器
8. **工厂模式** (`src/models/factory.jl`) - 模型创建工厂函数

### 文档交付物

1. **API文档** (`docs/api/`) - 从DocString生成的完整API文档
   - `data_contracts.md` - 参数格式和转换规则
2. **架构设计文档** (`docs/dev/architecture/`) - 详细的架构设计说明
   - `qcd_refactoring_dependencies.md` - 依赖关系分析
3. **开发指南** (`docs/dev/`) - 开发流程和规范
   - `branching_strategy.md` - Git分支策略
   - `code_review_guidelines.md` - 代码审查指南
4. **使用指南** (`docs/guides/`) - 用户和开发者指南
   - `error_handling.md` - 错误处理指南
   - `configuration.md` - 配置管理指南
   - `migration_guide.md` - 从旧API迁移到新架构的指南
5. **性能基准报告** (`docs/performance/benchmark_report.md`) - 重构前后的性能对比
6. **示例代码** (`docs/examples/`) - 展示如何使用和扩展架构的示例
7. **架构决策记录** (`docs/decisions/`) - ADR文档记录关键决策
8. **变更日志** (`CHANGELOG.md`) - 版本变更历史

### 测试交付物

1. **单元测试套件** (`tests/unit/models/`) - 所有模块的单元测试，包含覆盖率报告
2. **属性测试** (`tests/property/`) - 基于属性的测试验证通用正确性
3. **回归测试** (`tests/regression/pnjl_refactoring/`) - 验证数值等价性的回归测试
4. **集成测试** (`tests/integration/`) - 完整工作流的端到端测试
5. **性能基准** (`tests/perf/models/`) - 性能基准测试和结果数据
6. **边缘情况测试** (`tests/edge_cases/`) - 极端参数和边界条件测试
7. **异常处理测试** (`tests/unit/models/test_exceptions.jl`) - 验证所有异常类型和降级策略
8. **参数验证测试** (`tests/unit/models/test_validation.jl`) - 验证参数约束和验证逻辑

### 验证交付物

1. **测试覆盖率报告** - 显示≥80%核心功能覆盖率
2. **回归测试结果** - 证明数值等价性（误差<1e-12）
3. **性能基准数据** - 证明性能降低<5%
4. **代码审查记录** - 所有关键模块的审查记录
