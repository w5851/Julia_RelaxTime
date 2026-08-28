# Restore commands for the deleted RS `prod_v1` snapshot

The deletion proposal does not rewrite Git history. The pre-delete snapshot is
available at path-retirement merge commit
`74b53b47ebcca2b292cee72f70a70a84b0d2eea5` as long as that commit remains
reachable.

To restore the complete snapshot into a disposable branch/worktree:

```text
git switch --detach 74b53b47ebcca2b292cee72f70a70a84b0d2eea5
git switch -c restore/issue130-rs-legacy-prod-v1
git restore --source=74b53b47ebcca2b292cee72f70a70a84b0d2eea5 -- data/outputs/results/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1 data/outputs/figures/relaxtime/transport/phase_guided/legacy_prod_v1_snapshot_v1
```

The allowlist and tree hashes in this package must be rechecked after any
restore. Restoring the files is not a runtime switch and must not be done in
the main worktree without a new review.
