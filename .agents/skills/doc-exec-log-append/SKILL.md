---
name: doc-exec-log-append
description: Append one self-contained execution record to an existing or explicitly requested development log without rewriting history. Use only when the user asks for a log entry or the selected active task explicitly requires an execution ledger; do not trigger from ordinary implementation progress alone.
---

# Execution Log Append

## Invocation contract

- Require either an explicit user request or an explicit logging obligation in the selected task document.
- Do not infer a logging obligation merely because commands, artifacts, or test results exist.
- Pass `--log-file` or `--task-file` explicitly; do not rely on automatic latest-task discovery.
- Create a missing log only when the request or task contract authorizes creation.

## Hard rules

- Append only; never rewrite, reorder, or summarize historical entries.
- Do not read the historical body unless checking a user-requested duplicate/conflict or repairing a damaged format.
- Keep every entry self-contained: goal, changes, validation commands, artifacts, result, and mainline mapping.
- Keep task definition and DoD in the task document, not in the execution log.
- Record facts only; do not mark unverified work complete.

## Required input

- target `log_file` or `task_file`;
- `batch_id`;
- `goal`;
- code or document change summary;
- validation commands;
- artifact paths;
- result;
- mainline mapping.

## Workflow

1. Confirm the logging authorization and explicit target.
2. Collect reproducibility facts without reading the old log body.
3. Run `scripts/dev/append_exec_log.jl` with `--dry-run` and the explicit target.
4. Inspect the rendered block for missing fields or unintended paths.
5. Run the same command without `--dry-run` to append exactly one block.
6. Report the target file, batch id, command count, artifact count, and append result.

Use the script's `--help` output as the CLI contract instead of copying its full option list into this skill.

## Handoff

- Use `doc-implementation` when the request is to advance the task mainline.
- Use `experiment-logbook-append` for parameter scans, reruns, or scientific experiment records.
- Use `doc-archive` when a completed active document should be archived.
