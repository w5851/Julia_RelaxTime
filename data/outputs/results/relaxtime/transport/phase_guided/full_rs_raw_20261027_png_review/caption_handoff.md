# Full-RS 原始 PNG 图注与参数

本包从已验收的 1819 个原始点直接绘制：mode-A 910 点，mode-B 909 点。
12 个量均读取 CSV 同名列；不重新计算比值，不平滑、插值、拟合、删点或替换值。
470 个有限 τ 比值提示点全部保留。相邻原始点以直线相连；关闭路径顶点简化。
图像为本机纯 Python 绘制，没有调用求解器或补做本机数值验证。

图包与 v13 对应：两张 mode-A 主图、四面板低值细节图、72 张 A/B 单图；
另有逐张同尺寸灰度 PNG。只借用历史结构与排版，未读取旧 clean_value 点表、
平滑表、替换表或旧 gap 表。

Figure 1 建议图注：

> Raw full-RS relaxation times versus the anisotropy parameter ξ. Columns show
> μB=0, 450 and 900 MeV; rows show τu, τs, τū and τs̄. All vertical axes are
> logarithmic with independent ranges. Solid, dashed and dash-dotted curves
> denote αT=1.0, 1.1 and 1.2, with T=αT Tphase(μB, ξ=0).
> For μB=900 MeV and αT=1.0, the line is disconnected between the original
> samples at ξ=−0.003 and +0.003. Open circles and squares identify the
> nearest sampled points on the chirally restored and broken sides of the
> reviewed first-order discontinuity; they are not extrapolated coexistence endpoints.
> All original samples are retained.

Figure 2 建议图注：

> Raw full-RS transport ratios η/s, ζ/s and σ/T versus ξ, with the same
> columns and parameter encoding as Figure 1. Vertical axes are linear and
> use independent ranges. The first-order markers have the same limited
> scope as in Figure 1. Ratios are read directly from the accepted source CSV.

低值细节图建议图注：

> Linear low-value views of Figure 1 panels (a), (g), (h) and (i), using
> exactly the same raw vertices. Their vertical windows are 0.57–0.66,
> 0.57–0.66, 0.375–0.43 and 0.27–0.41 fm. The first two windows extend
> 0.01 fm below the v13 windows to retain the new αT=1.2 low-value branch.
> Curves outside these windows continue in the full-range main figure;
> leaving a viewport is not a physical endpoint.

单图：mode-A 的参数图例为 αT；mode-B 固定 T=120/160/200 MeV，
参数图例为 μB=0/450/900 MeV。mode-B 的 T=120 MeV、μB=900 MeV
只在 ξ=−0.13 与 −0.12 两个原始点之间断线，并保留两侧点及相态 marker。
其余相态字段变化（包括 crossover）不触发断线。单图保持全部原始 y 范围；
mode-A 弛豫时间统一 log-y，其余单图仅在正值动态范围至少 20 时选 log-y。
每图具体轴范围、输入行和分段记录在 manifest。

自然单位：τ 为 fm；η、ζ 为 fm⁻³；σ 为 fm⁻¹；η/s、ζ/s、σ/T 和 ξ 无量纲。
μB=3μq。源表 hadron/quark 对应此 PNJL 计算的手征破缺/恢复分支，
不意味着另算了一套强子气体输运。分段依据来自本次 quote-remote 背景审计。

温度映射（原始字符串，显示位数不表示不确定度）：

| μB (MeV) | αT=1.0 | αT=1.1 | αT=1.2 |
| ---: | ---: | ---: | ---: |
| 0.0 | 200.0898076157829 | 220.0987883773612 | 240.10776913893946 |
| 450.0 | 182.72820247485038 | 201.00102272233542 | 219.27384296982044 |
| 900.0 | 125.7360002645551 | 138.3096002910106 | 150.8832003174661 |

版面：原生宽 171.45 mm（6.75 in），600 dpi，普通及数学字体来源逐图记录。
推荐原生宽度查看；逐图 placement_limits 给出最小插入宽度，不默认允许缩成单栏。
彩色与灰度各 75 张。PNG 作者接受、矢量 PDF、论文资格与当前图层晋升分别处理。
本次交付阶段为 png_review，manuscript_eligible=false，current_publication_layer=false。
图像可用于检查局部起伏；排版通过或相邻变化筛查通过不证明所有非物理尖峰已消除。
