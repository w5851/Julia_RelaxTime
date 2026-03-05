---
name: vscode-symbol-tools
description: Prefer semantic symbol tools for precise navigation and refactoring. Use vscode_listCodeUsages and vscode_renameSymbol instead of grep-style text search when changing code symbols.
---

# Skill: vscode-symbol-tools

## When to apply
- User asks to rename a function, variable, class, method, or type.
- User asks where a symbol is used, referenced, or implemented.
- Refactoring needs high precision across files.

## Tool policy
- First choice for usages: `vscode_listCodeUsages`.
- First choice for rename: `vscode_renameSymbol`.
- Do not use `grep`/text search for symbol refactors unless the symbol tool is unavailable.

## Required workflow
1. Locate one real occurrence of the symbol in a file.
2. Use `vscode_listCodeUsages` to inspect impact before edits.
3. Use `vscode_renameSymbol` to apply semantic rename.
4. Run tests/lint relevant to changed scope when available.

## Example requests
- `Use vscode_renameSymbol and rename fib to fibonacci.`
- `Use vscode_listCodeUsages to find all references of thermalConductivity.`

## Notes
- Symbol tools depend on language server support.
- If rename cannot be resolved semantically, report why and fall back carefully.