# v13 图注与插入尺寸交接

两张主图均在 171.45 mm（6.75 in）宽度直接绘制。不要默认缩为单栏；
逐图 `rendering.placement_limits` 给出由字形、线宽、端点和 PNG 有效 DPI
共同约束的插入区间。推荐保持原生宽度，最终整篇稿件另行测量。

mode-A 单图与主图的参数 key 为 α_T；mode-B 单图固定 T，参数 key 为
μ_B (MeV)，条目 0 / 450 / 900。各单图的一阶 key 仅解释 manifest
`legend_placements[].applies_to` 列出的实际分支端点，不将 mode-B 改标为 α_T。

Figure 1（弛豫时间）建议图注：

> Relaxation times as functions of the anisotropy parameter $\xi$. Columns
> correspond to $\mu_B=0$, 450, and 900 MeV; rows show $\tau_u$, $\tau_s$,
> $\tau_{\bar u}$, and $\tau_{\bar s}$ in fm. All vertical axes are logarithmic,
> with independent ranges. The common color and line-style key in panel (a)
> applies throughout: solid, dashed and dash-dotted lines correspond
> to $\alpha_T=1.0$, 1.1 and 1.2. For $\mu_B=900$ MeV and $\alpha_T=1.0$,
> the First-order key in panel (c) identifies open circles and squares as
> the endpoints of the chirally restored and
> chirally broken branches at the first-order transition, respectively.
> Disconnected segments retain the branch gap. Low-value details of panels
> (a), (g), (h) and (i) are provided separately with linear vertical axes.

Figure 2（输运系数）建议图注：

> Transport ratios as functions of $\xi$. Columns correspond to $\mu_B=0$,
> 450 and 900 MeV; rows show $\eta/s$, $\zeta/s$ and $\sigma/T$. Vertical axes
> are linear with independent ranges. The common line-style key in panel (a)
> and the First-order key in panel (c) use the conventions of Figure 1. The endpoint
> key applies only to $\mu_B=900$ MeV, $\alpha_T=1.0$ curves.

局部视图建议图注：

> Low-value views of Figure 1 panels (a), (g), (h) and (i), using the same
> display values and line styles. The vertical axes are linear and restricted
> to the indicated ranges; curves leaving the view continue in Figure 1.
> First-order endpoints are outside these windows. Panel labels refer to the
> corresponding main-figure panels. The parameter key in panel (a) applies
> to all four views; these views introduce no new data.

注意：这里的 chirally broken quark branch 不表示另行计算了强子输运。
显示曲线继承 v5 已披露的 display adjustments；没有新增平滑或数值收敛证据。
论文讨论不得把不同纵轴范围下的视觉斜率当成可直接比较的增长率。

精确冻结温度映射（MeV；表格精度不代表不确定度）：

| panel | alpha_T=1.0 | alpha_T=1.1 | alpha_T=1.2 |
| --- | ---: | ---: | ---: |
| muB0.0 | 200.088901 | 220.097792 | 240.106682 |
| muB450.0 | 182.763066 | 201.039372 | 219.315679 |
| muB900.0 | 125.737258 | 138.310984 | 150.884710 |

PNG 审查层：作者接受待定，PDF 待后续交付，manuscript_eligible=false，current=v5。
