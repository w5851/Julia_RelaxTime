# Restore commands for the deleted PNJL legacy snapshot

The pre-delete tree is preserved by the path-retirement merge
`9aa4c313901ca0c91e851f58514e3df9aa124df4`. Restoring it is an explicit
historical-recovery operation, not a runtime rollback or fallback.

```text
git switch --detach 9aa4c313901ca0c91e851f58514e3df9aa124df4
git switch -c restore/issue130-pnjl-legacy-phase-reference-v1
git restore --source=9aa4c313901ca0c91e851f58514e3df9aa124df4 -- data/reference/pnjl/legacy_phase_reference_v1
```

After a restore, re-run the physical-deletion validator against the disposable
branch and verify every row in `deletion_allowlist.csv`. Do not restore these
files in the main worktree or re-enable a runtime legacy path without a new
reviewed contract.
