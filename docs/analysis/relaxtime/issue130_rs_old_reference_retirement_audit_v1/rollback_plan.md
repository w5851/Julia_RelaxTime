# RS old-reference retirement rollback plan

本文件只定义后续 implementation PR 的回退顺序，不执行删除。

1. 保持 current `prod_v2` registry entry、分析默认值和 figure path 不变。
2. 若要退出旧 `prod_v1` canonical path，先把 raw/figure tree 做 byte-preserving versioned snapshot，并写入 source path、tree SHA-256、manifest SHA-256。
3. 更新 resolver/registry/documentation，使显式 legacy fallback 指向该 snapshot；默认解析仍必须指向 `prod_v2`。
4. 在 solver-free smoke 中验证 default、explicit legacy 和 rollback 三条路径。
5. 若任一检查失败，恢复 registry/path pointers；不改 raw bytes。
6. 物理删除必须另开 PR，携带 deletion allowlist、历史依赖审计和作者的单独授权。

回退原则：先恢复 pointer，再恢复 active path，最后才考虑文件操作；任何阶段都不重算数值。
