# Phase-guided transport 绘图案例合同

本页保存 v13 的当前审查布局，以及 v11/v12 的历史布局。公共 APS v2 规则见
[绘图 SOP](../../../guides/sop/workflows/figure_production.md)，逐图 panel specs、caption 与
冻结输入仍以各版本 manifest 为准。历史产物字节和审核资格不随新布局改变。

## v13 图内共享／分置图例布局

v13 从 v12 同一冻结点表生成独立 PNG 审查 case。布局入口为
`scripts/analysis/relaxtime/build_phase_guided_publication_clean_v13.py --png-review`。

- 两张主图的参数 key 放在 (a) 左上空白区，标题 `α_T`，条目仅
  `1.0 / 1.1 / 1.2`；端点 key 放在 (c) 左上空白区，标题 `First-order`，
  圆圈为 `restored`、方块为 `broken`。前者全图通用，后者仅解释
  `mu_B=900 MeV, alpha_T=1.0` 的实际一阶分支端点。
- 参数数学标题 13 pt；三个数字条目及一阶 key 的标题／条目 11 pt。
  列标题和行标签 13 pt，复合图刻度和 panel 编号 11 pt。参数标题的
  数学下标仍须实测达到 2 mm；不使用 1.5 mm 紧凑字形例外。
- 保留 171.45 mm 宽度和 Figure 1/2 的 7.1/5.9 in 高度，回收顶部公共
  图例占用的空间。Figure 1 全 log-y、Figure 2 全 linear-y；独立 y 范围、
  行内精度、三个 x 数字标签、全部源曲线顶点及 gap 均保留。
- 四面板局部图沿用 v12 的 (a)/(g)/(h)/(i) 和全部线性 y 窗口，参数 key
  在 (a)；端点在窗口外，不添加端点 key，也不把 viewport 当成曲线终点。
- 72 张单图按原生双栏尺寸绘制。确定性地尝试声明的四角组合，按真实
  曲线／端点几何选择无挡位置；同一 axes 可保留两份独立 key，必须逐份
  检查。全部组合不合适时允许简短顶部备选，理由和实际位置写入 manifest。
  mode-A 使用 α_T key；固定 T 的 mode-B 使用 `μ_B (MeV)` 标题与
  `0 / 450 / 900` 条目。颜色／线型仍绑定各自原 series，不混淆两种参数。
  端点 key 的适用曲线在 `legend_placements[].applies_to` 逐项登记。
- 新 policy 为 `declared_geometry_checked`。每份 key 的数量、host、标题、
  条目和位置声明与实测逐项一致；图例和 axes 相交本身允许，曲线／端点
  遮挡、图例互相覆盖、图内溢出、文字裁切／重叠均不允许。
- 75 对彩色／灰度 PNG 先交付审查；同尺寸灰度转换、逐图缩放区间、
  最终 600 dpi 几何、线型追踪和端点检查保留。单图也不能默认缩成单栏。
  PNG 已由作者接受，见 [接受记录](publication_clean_v13_png_acceptance_v1.json)；
  对应 75 份矢量 PDF 已使用冻结 renderer 导出，见 [PDF 图包](phase_guided_transport_publication_clean_v13_pdf/README.md)。
  论文资格及 current-layer 晋升分别处理。

v12 更新前的 27 份源码／合同字节保存在
`publication_clean_v12_code_snapshot_v1.zip`，archive hash 登记于
`config/plotting/historical_snapshots.toml`。其中 base commit 仅作上下文，
不把未提交字节伪装成提交内容；原 v12 manifest 和图像保持不变。
文献样本的 [调查报告](literature_figure_style_audit_20261005/README.md) 也作为
冻结输入保留；其比例用于选择试排顺序，不用于禁止顶部或其他图外布局。

## v12 最终尺寸审查布局（历史）

v12 是新 sibling 审查 case，v11 的下面各节继续作为历史布局合同保存。
当前 v12 使用 `build_phase_guided_publication_clean_v12.py --png-review`：

- 原生推荐宽度 171.45 mm；Figure 1/2 保留 v11 的 7.1/5.9 in 高度。
  图例和列标题 13 pt，主刻度与 panel 编号 11 pt，行轴标签 13 pt。
  这些字号需逐图实测，不使用 1.5 mm 紧凑字形例外。
- 两组公共图例置于图外顶部：参数 key 在第一行，第二行用
  `First-order (restored)`／`First-order (broken)`。端点语义只适用于
  `mu_B=900 MeV, alpha_T=1.0`，不能说全部 panel 都具有一阶端点。
- Figure 1 仍为全 log-y、Figure 2 为全线性轴；保留独立 y 范围、
  三个 x 数字标签、行内小数位、原曲线顶点和分支 gap。
- 另给四面板低值视图，对应 Figure 1 的 (a)/(g)/(h)/(i)，使用线性
  y 轴和明确范围：前两幅 0.58--0.66 fm，(h) 0.375--0.43 fm，
  (i) 0.27--0.41 fm。离开窗口的曲线仍保留在完整主图；一阶端点
  在这些窗口外。局部图不能用于确认新增物理结构或数值收敛。
- 所有 72 张单图仍按 6.75 x 4.6 in 绘制，改用图外 key，避免自动
  `best` 图例在个别 mode-B 曲线／端点上遮挡。逐图报告缩放区间及单栏
  复用资格；不允许将原生双栏字高的通过结论用于半宽插图。
- 每张图给出同尺寸彩色及灰度 PNG，共 75 对；一个总 manifest 保留
  逐图记录。灰度转换不改变数据。作者视觉接受、PDF 与论文资格待后续阶段。

图注见 v12 case 的 `caption_handoff.md`，其中明确图外 key、局部图与
原图 panel 的对应关系，以及继承的 display-only 数值边界。机器检查
和人工审查范围分别留证，不能用总 manifest 或测试通过替代人工审图。

## v11 图例与端点语义（历史）

对于论文/组内博士论文风格的密集输运复合图，图例可以采用一次性的
case-level in-panel key，但必须遵守以下收缩规则：参数曲线只保留
`$\alpha_T=1.0$`、`$\alpha_T=1.1$`、`$\alpha_T=1.2$`；端点项使用
`First-order (restored)` 与 `First-order (broken)`；若最终 panel 的曲线几何
不允许容纳完整短语，可以退为 marker 旁的 `restored` 与 `broken`，但必须
在 caption 中明确它们是 first-order transition 的 chirally restored/broken
endpoints。分开的端点 key 使用共同标题 `First-order`，其下分别写圆圈
`restored` 和方块 `broken`，避免重复相同文字。圆圈和方块的完整物理含义放入 caption，而不是在 panel 内重复
长句。标题与条目左对齐、同字号、常规字重；使用公开的 legend
`alignment="left"`，不操作 `_legend_box` 等私有属性。可将参数图例与端点
图例分放到两个有空白的 panel，各只出现一次；
不能为固定的“一处五项图例”遮盖数据。子图编号统一放在左上方；若该位置
仍有数据，可以全图统一移到框外左上方，并保留行间距和标题间隔。不可只把
承载图例的编号移到低值曲线区域来规避图例。把复合图轴标签、
刻度、标题及 legend 字号作为 case override 记录（通常约 10--11 pt）；不能
因为压缩图例而牺牲最终尺寸下的可读性。图例 bbox 与所有实际曲线路径必须
通过几何检查，不能仅检查文字与坐标框是否重叠；`legend_curve_overlap_count`
必须为 0，除非 manifest 明确列出并由作者接受该例外。圆圈/方块等端点
还必须通过 `legend_landmark_overlap_count=0` 检查。caption 分别限定适用范围：
参数颜色/线型 key 全图通用；端点 key 仅解释实际具有一阶端点 marker 的曲线。
本输运图族的端点对应 `mu_B=900 MeV, alpha_T=1.0`，不能笼统写两个 key
均适用于所有 panel。

## v11 坐标轴与数字显示（后续版本沿用）

输运 v7--v10 保留历史 `1st-order` 措辞；新 v11 审查层使用 `First-order`。
完整的 chirally restored/broken、endpoint、温度及分支断线解释在 caption handoff。
复合图默认只显示三个 x 轴数字刻度（左右端点和零点）；底层 tick 策略和
case override 必须写入 manifest。Figure 1 指四行弛豫时间图，Figure 2 指三行
输运系数图，不按外部 review 的图号推断目标。v11 的 Figure 1 全部 12 个
panel 统一 log-y，匹配的 mode-A 弛豫时间单图也使用 log-y；Figure 2 保持
线性轴和已审查的布局。独立 y 范围保留，caption 必须说明对数轴及各 panel
范围不同；不能把 log 变换后的视觉曲率解释为绝对增长率或饱和证据。这是
本正值输运图族的 case-level 规则，不是对其他图族强制使用 log 的要求。
标题单位使用显式可见空格，轴标签统一使用 `$\tau_u\;(\mathrm{fm})$`
一类的数学变量+正体单位格式。复合图记录 figsize、hspace/wspace、legend
host、legend 字号和 `legend_curve_overlap_count`。

### 多面板数字刻度一致性

- 线性轴同一物理量的一行统一小数位，按该行最细的主刻度间距确定；这只是
  显示格式，不表示数据不确定度或有效精度。本图族的 `eta/s`、`zeta/s`
  为两位小数，`sigma/T` 为三位；不同物理量之间不强求相同小数位。
- log-y 优先在 `1,2,5` 乘以十的整数次幂处标出普通数字，如
  `0.5,1,2,5,10,20`；统一去除冗余尾零，不强制写成 `1.00`，不使用含混的
  offset。窄范围可补充其他整洁的实值刻度，争取每个 panel 约 3--5 个
  标签；宽范围允许稀疏采样，具体刻度与格式规则写入 manifest。
- 普通数字标签不改变对数位置。log 小刻度保持 `2--9` locator，去除与主
  刻度重复的位置，不标小刻度数字。改变轴类型、formatter 或 legend 宽度后，
  必须重新检查文字拥挤、裁切以及曲线/端点遮挡。

线性复合图 y 轴默认约 4--5 个主刻度，不能以降低数字字号掩盖刻度线过密。
`sigma_over_T` 采用三位小数时，主刻度必须落在 `0.001` 的整数倍上，避免
把 `0.0175` 等实际刻度标成 `0.018`。保留前置零；如使用公共倍率，
必须是真实的乘法尺度（如 `10^{-3}`）且记录 transform，不能把重复的
`0.` 提取成含混的公共前缀。

## 入口与历史核验

- 当前审查入口为上述 v13 builder；v11/v12 入口用于历史追溯，已存在的 case 均拒绝覆盖。
- 已接受阶段用 `scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --check` 按登记的历史 commit 核验。
- `publication_clean_v11_code_snapshot_v1.zip` 保存该提交的 21 份代码/合同文件与原阶段验收记录；registry 记录 archive hash。浅克隆或 squash 后不依赖旧 Git object 可用性；现存数据/图像仍直接验证，不从源码 archive 替代。
- 阶段验收、正式矢量交付、数值生产与 current-layer 晋升分别依赖各自的作者授权。
