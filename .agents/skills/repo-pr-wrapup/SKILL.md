---
name: repo-pr-wrapup
description: Complete Julia_RelaxTime branch, selective staging, validation, commit, push, PR creation or PR update after implementation is ready. Use for repository wrap-up and merge preparation; unresolved GitHub review comments must first be handled by gh-address-comments or an equivalent review workflow.
---

# Repository PR Wrap-up

## Core rules

1. Determine whether the target is an existing PR or a new PR.
2. If unresolved review comments require code changes, hand off to `gh-address-comments` before wrap-up.
3. Use `codex/{topic}` only when a task branch is needed and the user did not specify another name.
4. Inspect `git status --short` and stage only files belonging to the requested scope.
5. Use a concise commit subject matching the staged scope; consult recent history if the convention is unclear.
6. Run tests and governance checks proportional to changed behavior.
7. Treat the PR title as the intended squash-merge subject.
8. Check active-document completion and archive obligations when `docs/dev/active/` is in scope.

## Workflow

### Inspect state

- Check the branch, worktree, staged diff, open PR, and unresolved reviews.
- Preserve unrelated user changes and generated artifacts outside scope.
- Route review-comment implementation before continuing wrap-up.

### Select and validate

- Stage explicit files or tightly scoped path groups; do not default to `git add -A`.
- Review `git diff --cached --stat` and the cached diff.
- Choose focused unit, integration, regression, validation, and governance checks.
- Record skipped checks and risk; never loosen numerical tolerances to pass.

### Commit subject

- Use a familiar prefix such as `docs:`, `fix:`, `refactor:`, `ci:`, or `feat:`.
- Keep the subject concise, intent-driven, and single-line; match the actual staged change.

### Prepare the PR

Use a merge-ready title and keep these body sections:

1. `变更范围`
2. `用户影响`
3. `实现方式`
4. `验证项`
5. `已知非目标`

Update an existing PR instead of creating a duplicate.

### Active documents

- Archive a completed task through `doc-archive` when appropriate.
- Keep unfinished documents active and align status with evidence.
- Never mark unverified work complete.

## Final report

Report the branch, change scope, validation results, PR link and any remaining work. If publication was not executed, state the exact step reached.

Do not commit, push, or create a PR when the user requested review or preparation only.
