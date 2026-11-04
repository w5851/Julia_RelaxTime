# Full-RS publication_clean_v14

本包直接读取全网格验收通过的 1819 个原始点和 12 个同名列。
无平滑、插值、删点或局部替换，关闭绘制路径简化；470 个有限 τ 比值提示保留。
仅断开背景审计确认的一阶相变相邻原始点，空心圆/方块表示恢复/破缺侧最近采样点。
它们不是外推的共存端点。所有数值检查在 quote-remote，纯 Python 绘图在本机。

figure1：列为 μB=0/450/900 MeV，行为 τu/τs/τū/τs̄，纵轴统一 log。
figure2：同样三列，行为 η/s、ζ/s、σ/T，纵轴 linear。
两图实线/虚线/点划线依次为 αT=1.0/1.1/1.2，T=αT Tphase(μB,ξ=0)。
mode-A 一阶断线只适用于 μB=900 MeV、αT=1.0 的 ξ=−0.003 到 +0.003。
温度原始精确映射、每条曲线输入行和每幅图轴范围均见 plot_manifest.json。

figure4：列为 T=120/160/200 MeV，行为 η/s、ζ/s、σ/T；前两行 log，第三行 linear。
实线/虚线/点划线依次为 μB=0/450/900 MeV。图例保留用户指定的面板内布局。
一阶断线只适用于 T=120 MeV、μB=900 MeV 的 ξ=−0.13 到 −0.12。
建议 FIG4 图注：Transport ratios at fixed temperatures, shown as functions of ξ.
Columns correspond to T=120, 160 and 200 MeV. Line styles identify μB=0, 450 and
900 MeV. Open circles and squares mark the nearest sampled restored/broken points
on either side of the reviewed first-order discontinuity. All raw samples are retained.

低值细节图保留线性窗口：0.57–0.66、0.57–0.66、0.375–0.43、0.27–0.41 fm。
前两窗比历史 v13 下延 0.01 fm，以容纳本次原始数据；离开窗口不表示物理端点。
其余 72 张单图沿用 v13 的面板、轴类型、图例角落、线型与蓝/红/绿色。
坐标范围随原始值适配。FIG4 尺寸 6.75×7.4 in，其他图保留 v13 尺寸。

τ 单位 fm；η/ζ 单位 fm⁻³；σ 单位 fm⁻¹；η/s、ζ/s、σ/T 与 ξ 无量纲。
μB=3μq；hadron/quark 标签指当前 PNJL 模型的手征破缺/恢复分支。
图的最小插入宽度由各自 placement_limits 指定，不能默认缩成单栏。
PNG 接受、数值晋升及 PDF 验收记录保存在同名 results case；论文文件并未自动替换。
