# Phase-guided transport 绘图案例合同

本页保存 v11 图族的已审查布局规则，移自 figure_production_v2 SOP。公共 APS v2 规则见
[绘图 SOP](../../../guides/sop/workflows/figure_production.md)，逐图 panel specs、caption 与
冻结输入仍以各版本 manifest 为准。本次迁移不改变 renderer、产物字节或审核资格。

## 图例与端点语义

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

## 坐标轴与数字显示

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

- 新审查运行使用 `scripts/analysis/relaxtime/build_phase_guided_publication_clean_v11.py --png-review`，输出新 case，保留原产物。
- 已接受阶段用 `scripts/analysis/relaxtime/formalize_phase_guided_publication_clean_v11_stage.py --check` 按登记的历史 commit 核验。
- `publication_clean_v11_code_snapshot_v1.zip` 保存该提交的 21 份代码/合同文件与原阶段验收记录；registry 记录 archive hash。浅克隆或 squash 后不依赖旧 Git object 可用性；现存数据/图像仍直接验证，不从源码 archive 替代。
- 阶段验收、正式矢量交付、数值生产与 current-layer 晋升分别依赖各自的作者授权。
