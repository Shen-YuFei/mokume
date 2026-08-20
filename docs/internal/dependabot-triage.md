# Dependabot PR triage

Internal triage of the 12 open Dependabot PRs on `bigbio/mokume`. This file
lives under `docs/internal/` and is excluded from the built site by
`mkdocs.yml` (`exclude_docs: internal/*`), so it is a repo-internal note, not a
published page.

State was refetched with `gh pr list --repo bigbio/mokume --state open` and
`gh pr view <n> --repo bigbio/mokume --json number,title,mergeable,mergeStateStatus`.

## Summary

All 12 PRs are Dependabot bumps, all target `base=main`, and all report
`MERGEABLE`. On the current refetch the CI picture has narrowed: only the two
coupled Rust bumps `#32` (arrow 58 -> 59) and `#35` (parquet 58 -> 59) are still
`UNSTABLE` with real `FAILURE` checks; the other ten (the five GitHub Actions
PRs `#27`-`#31`, the three Python PRs `#33`/`#34`/`#38`, and the two low-risk
Rust PRs `#36`/`#37`) now show `CLEAN` / all-`SUCCESS`. This matches the
hypothesis that the widespread earlier instability came from a shared,
pre-existing failure on `main` rather than from the bumps themselves: once
`main` is green, most of these PRs go green on rebase/re-run. So the standing
recommendation is unchanged: keep `main` CI green first, then re-trigger any PR
that predates the fix so its status reflects the current tree, and merge in the
order below. Note also that when the packaging/Rust layer is re-architected (the
`pyproject.toml` / `hatchling` build backend and the Rust workspace under
`rust/`), some of these bumps will need a rebase or a re-trigger to pick up the
new tree, and a couple of the Python PRs can conflict or become redundant --
concretely, PR `#34` only edits the `numpy<2.1.0` pin in `python/pyproject.toml`
but the same pin also lives in `python/requirements.txt`, so that PR leaves the
two files out of sync and may need a follow-up edit or become moot if the
requirements files are consolidated.

## Triage table

| PR  | area           | bump                                | risk | verdict                          | how-to-resolve |
|-----|----------------|-------------------------------------|------|----------------------------------|----------------|
| #27 | GitHub Actions | upload-artifact 4 -> 7              | low  | merge-now                        | CI already CLEAN. Bumps the `actions/upload-artifact@v4` uses in `wheels.yml` and `rust-publish.yml`; v7 keeps v4 naming semantics, no code change. |
| #28 | GitHub Actions | checkout 4 -> 7                    | low  | merge-now                        | Ubiquitous `actions/checkout` bump, no behavior change relevant here. Merge once main is green. |
| #29 | GitHub Actions | setup-python 5 -> 6               | low  | merge-now                        | `actions/setup-python` bump; default Python resolution unchanged. Merge. |
| #30 | GitHub Actions | download-artifact 4 -> 8         | med  | verify-cargo-or-pytest-then-merge | download-artifact v4 changed artifact naming/merge behavior; verify the `actions/download-artifact@v4` uses in `wheels.yml` (line ~120) and `rust-publish.yml` (lines ~128/133) still resolve the same artifact names produced by the matching upload steps. Run the `wheels` and `rust-publish` jobs (or dry-run) before merge. |
| #31 | GitHub Actions | actions/labeler 5.0.0 -> 6.1.0   | low  | merge-now                        | Labeler v6 tolerates the existing labeler config; only affects PR labeling, not build. Merge. |
| #32 | Rust (cargo)   | arrow 58.3.0 -> 59.0.0           | high | needs-code-changes               | COUPLED with #35 (see notes). arrow 59 is a major bump with likely breaking API. Bump arrow AND parquet together in `mokume-io` and `mokume-golden-tests`, then `cargo build`/`cargo test` those crates and fix compile errors; must also pass `cargo clippy -D warnings` and `cargo doc -D warnings` per `rust.yml`. |
| #33 | Python (pip)   | pymdown-extensions >=10.0 -> >=11.0.1 | low | verify-cargo-or-pytest-then-merge | Docs-only dep in `python/requirements-docs.txt`. Verify with `mkdocs build --strict`; no runtime/pytest impact. Merge if strict build passes. |
| #34 | Python (pip)   | numpy <2.1.0 -> <2.3.0          | med  | verify-cargo-or-pytest-then-merge | Loosens the `numpy<2.1.0` cap in `python/pyproject.toml` (line 35). Confirm why it was pinned (coexists with `pyopenms`, a numpy-ABI-sensitive dep), then run pytest against numpy 2.2. Also update the duplicate pin in `python/requirements.txt` (line 4) or the two files drift. |
| #35 | Rust (cargo)   | parquet 58.3.0 -> 59.0.0        | high | needs-code-changes               | COUPLED with #32 (see notes). Same 58 -> 59 major bump in the same two crates (`mokume-io`, `mokume-golden-tests`). Do not merge alone -- must move in lockstep with arrow 59. Failing checks: `test`, `pre-commit`, and all five `wheel *` jobs. |
| #36 | Rust (cargo)   | regex 1.12.3 -> 1.12.4          | low  | merge-now                        | Patch-level bump in `mokume-pipeline`, no API change. CI already SUCCESS. Merge. |
| #37 | Rust (cargo)   | hdf5-metno 0.12.6 -> 0.13.0     | low  | verify-cargo-or-pytest-then-merge | Minor bump of the `hdf5-metno` (features `static`) dep in `mokume-command`. Low risk but a 0.x minor can shift API; `cargo test --workspace` in `mokume-command` to confirm, then merge. CI currently SUCCESS. |
| #38 | Python (pip)   | mkdocs-material >=9.0 -> >=9.7.6 | low | verify-cargo-or-pytest-then-merge | Docs-only dep in `python/requirements-docs.txt`. Verify with `mkdocs build --strict` only; no code impact. Merge if strict build passes. |

## Key notes

- **arrow #32 + parquet #35 are one coupled change.** Both bump `58.3.0 ->
  59.0.0` and both crates (`rust/crates/mokume-io/Cargo.toml`,
  `rust/crates/mokume-golden-tests/Cargo.toml`) declare arrow and parquet at the
  same version, so they must move together to keep the arrow/parquet ecosystem
  version-aligned. arrow 59 is a major release with a likely breaking Rust API.
  Resolution: bump both deps to `59.0.0` in both crates in a single change, run
  `cargo build`/`cargo test` on `mokume-io` and `mokume-golden-tests`, fix the
  compile errors, and confirm the whole workflow in `.github/workflows/rust.yml`
  passes (`cargo test --workspace --all-targets`, `cargo clippy --workspace
  --all-targets -- -D warnings`, and `cargo doc --workspace -D warnings`, in both
  debug and release). Merging either PR alone will leave arrow and parquet on
  mismatched majors and keep CI red -- these are the only two PRs currently
  failing.

- **download-artifact v4 -> v8 (#30) has breaking artifact-naming behavior.**
  The v4 line of `actions/download-artifact` changed how artifacts are named and
  merged versus v3, and jumping to v8 must be checked against the upload steps it
  pairs with. Inspect the `download-artifact` uses in
  `.github/workflows/wheels.yml` and `.github/workflows/rust-publish.yml` and
  make sure the artifact names/patterns they download still match what the
  `upload-artifact` steps in the same workflows produce. Verify by running (or
  dry-running) the `wheels` and `rust-publish` jobs before merge.

- **numpy <2.1 -> <2.3 (#34) loosens the pin.** The `numpy<2.1.0` cap is set in
  `python/pyproject.toml` and duplicated in `python/requirements.txt`; it likely
  exists because of a numpy-ABI-sensitive dependency (`pyopenms` is in the tree).
  Before merging, confirm the original reason for the cap and run pytest with
  numpy 2.2 installed. If it passes, also raise the cap in
  `python/requirements.txt` so the two files stay consistent -- PR #34 only
  edits pyproject.

- **Docs deps (#33 pymdown-extensions, #38 mkdocs-material)** live in
  `python/requirements-docs.txt` and only affect the documentation build. Verify
  each with `mkdocs build --strict` and nothing else; no pytest or cargo needed.

- **regex #36 and hdf5-metno #37 are low-risk.** #36 is a patch bump
  (`1.12.3 -> 1.12.4`) in `mokume-pipeline` with no API surface change. #37 is a
  minor 0.x bump (`0.12.6 -> 0.13.0`) of `hdf5-metno` in `mokume-command`; a `cargo
  test --workspace` covering `mokume-command` is enough to confirm before merge. Both
  are already CI-green.

## Suggested merge order

1. Get `main` CI green (resolves the shared pre-existing failure that caused the
   earlier broad instability), then re-trigger any PR opened before the fix.
2. Merge the CI-green low-risk PRs: #27, #28, #29, #31, #36, and (after a
   `mkdocs build --strict`) #33 and #38.
3. Verify then merge the medium-risk PRs: #30 (artifact naming), #34 (numpy),
   #37 (hdf5-metno).
4. Handle the coupled high-risk pair #32 + #35 last, as a single arrow+parquet
   59 upgrade with the required Rust code changes.
