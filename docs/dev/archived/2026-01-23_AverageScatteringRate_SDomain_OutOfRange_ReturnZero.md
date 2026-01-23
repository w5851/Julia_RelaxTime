---
title: AverageScatteringRate s-domain out-of-range returns zero
archived: true
original: docs/dev/active/任务6.md
archived_date: 2026-01-23
---

计算平均散射率时，当质心能量超出(s_bo,s_up)区间时，不应该截断到边界，而应该直接归零，因为超出区间意味着质心能量不符合物理意义，归零更符合实际情况,对应作出修改并同步到文档的提示中
