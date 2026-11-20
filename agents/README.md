# Portable agent pack

Rules and skills for any IDE. Cursor also mirrors rules to `.cursor/rules/`.

Do not depend on a machine-local `~/.cursor` after clone.

- `rules/` — always-on project conventions
- `skills/` — `cfa-viz`, `code-review`, `debug`. Only `cfa-viz` is mirrored to `.cursor/skills/` (a symlink). `code-review` and `debug` were copied at scaffold time; the templates they link (`review-template.md`, `debug-report-template.md`, `_shared/artifact-registry-template.md`) and the `artifact-registry` rule are not in this repository.
