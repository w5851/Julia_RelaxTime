# publication_clean_v15：η 阈值局部加密后的原始网格

本包直接绘制已正式验收的 1867 个数值点（mode A 942、mode B 925），读取源表的
六种 τ、η、ζ、σ、η/s、ζ/s、σ/T 共 12 个原始列，不在绘图时重算物理量。
1810 个 v14 原值原样保留；9 个原坐标采用本轮重算值，另加入 48 个新 ξ 点。
三条加密曲线是 mode A 的 μB=0/450 MeV、αT=1.0，以及 mode B 的 T=200 MeV、μB=0。
这57点含39个 s 表 rtol=10⁻⁵、18个 rtol=10⁻⁶ 的已验证结果；其余生产参数与 v14 相同。
精确参数、残差要求的实际覆盖、逐点选源和接缝审查见数值包的 manifest 与 review。
图中没有平滑、插值、拟合、删点或值混合；关闭路径简化，相邻实算点以直线相连。

新增点不伪造历史求根迭代数或相态追踪信息。48个新点的 phase 标签为 unknown；
9个原坐标保留旧相态元数据。两处已审计一阶断线及470个有限τ比值提示点保持不变。
断线只在 mode A μB=900、αT=1.0 的 ξ=−0.003/+0.003，以及 mode B T=120、μB=900
的 ξ=−0.13/−0.12。空心圆／方块为两侧最近的手征恢复／破缺采样点，不是外推共存端点。

全部76张彩色图与76张同尺寸灰度图复用 v14 renderer、面板、轴类型、单位、图例和线型：
72张单图、Figure 1/2/4 三张主组合图，以及一张四面板低值细节图。
Figure 1：列 μB=0/450/900 MeV，行 τu/τs/τū/τs̄，纵轴 log。
Figure 2：同样三列，行 η/s、ζ/s、σ/T，纵轴 linear。
这两图实线／虚线／点划线对应 αT=1.0/1.1/1.2，T=αT Tphase(μB,ξ=0)。
Figure 4：列 T=120/160/200 MeV，行 η/s、ζ/s、σ/T；前两行 log、第三行 linear。
其参数线型对应 μB=0/450/900 MeV，保留20261021版面板内图例。
低值细节窗口仍为0.57–0.66、0.57–0.66、0.375–0.43、0.27–0.41 fm。
τ单位 fm；η/ζ 为 fm⁻³；σ为 fm⁻¹；ξ、η/s、ζ/s、σ/T 无量纲；μB=3μq。
Figure 4 为6.75×7.4 in；其余尺寸沿用v14。最小插入宽度以各图 placement_limits 为准。

建议补充图注：The accepted full-RS grid includes additional computed samples around three
reverse eta-Mott thresholds. Nine original-coordinate values are replaced by the verified
recomputations at the refined integration precision; all other original values are retained.
The plotted curves connect the accepted samples without smoothing or fitted interpolation.
The open markers retain their limited meaning as the nearest sampled points across the two
reviewed first-order discontinuities.

全部数值回放与派生审计在 quote-remote 完成；本机仅进行文件核验及 Python 绘图。
本包为 PNG 审阅阶段：manuscript_eligible=false、current_publication_layer=false，
作者对新PNG的接受与后续PDF交付另行记录；数值产物的正式资格不代替视觉接受。
