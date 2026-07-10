---
name: doc-archive
description: Archive completed docs/dev/active task documents into docs/dev/archived with required metadata by using scripts/dev/archive_docs.jl and validating the result. Do not use for unfinished tasks, ordinary content edits, or custom nonstandard moves.
---

# Document Archive

## Contract

- Archive only completed active documents.
- Use `scripts/dev/archive_docs.jl`; do not reproduce its move/frontmatter logic manually.
- Preserve the original content and record `title`, `archived`, `original`, and `archived_date`.
- Follow `docs/dev/README.md` naming and metadata rules.

## Workflow

1. Confirm that the task is complete and its status matches implementation and validation evidence.
2. Run a dry run for the explicit active filename:

```powershell
julia --project=. scripts/dev/archive_docs.jl --dry-run <filename.md>
```

3. Execute the archive, optionally passing `--date YYYY-MM-DD`.
4. Confirm the source left `docs/dev/active/`, the destination exists under `docs/dev/archived/`, and all required frontmatter fields are present.
5. Run the archive format check:

```powershell
julia --project=. scripts/dev/archive_docs.jl --check
```

6. Report the destination path and validation result.

## Guardrails

- Do not overwrite an existing archive without explicit authorization.
- Do not archive a partially complete document to make active governance pass.
- Do not search archived history for current facts unless the task requires provenance.
