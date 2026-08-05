# 审计记录

- source run：`30980094983`
- calculation SHA：`3217bed3635574f00c04cbee75e843b4c49451db`
- workflow head SHA：`df1fcdece7bc0888dd57c465c0015828743596c5`
- curve rows：`1293`；metric rows：`13`；frontier rows：`5`
- targeted additions：`12`
- source job manifest SHA256：`ee2c74821d11d4b79cde932ccffc247e4602b7246697ebd6fe8ab4f937ffced5`
- source aggregate manifest SHA256：`e91e3cd3fb783213f7d6e5ec179e8c25f90763c20a292aa88326c2dbdc661af1`
- solver called：`false`；reference write：`false`
- contract errors：`[]`

作者裁决只接受 `rho_h` 的有界零密度端点极限，不把该区间伪装成严格正密度根，也不
替代后续 production integration 与 targeted/full shadow gate。
