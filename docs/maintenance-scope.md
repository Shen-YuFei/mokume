# Maintenance scope: Rust-first

mokume ships its computation in two codebases — a **Rust compute kernel**
(`rust/crates/`) and a **pure-Python package** (`python/mokume/`, `pip install
mokume`) — and both currently implement the same four computation commands:
`features2proteins`, `features2peptides`, `peptides2protein`, and
`correct-batches`. This page states which implementation leads, what each is
maintained for, and where new work goes, so the overlap does not drift into two
subtly different behaviors.

!!! abstract "The policy in one line"
    **Rust is the leading implementation. The pure-Python computation package is
    added value — a readable reference and parity check — not the place new
    computation lands first.**

## What this covers

This page is about the **computation implementations only**: the quantification,
normalization, imputation, batch-correction, and differential-expression logic
behind the four commands above.

It does **not** govern the Python periphery — the pipeline API, shared
post-processing, plotting, reporting, [TissueMap](periphery/tissuemap.md), and
the `agentic` optimizer. Those are Python-only by design and are out of scope
here; see [Architecture](architecture.md) and [CLI vs Wheel](cli-vs-wheel.md).

## The rule

1. **The Rust kernel is authoritative.** For every computed quantity, the Rust
   kernel defines the correct behavior, supported options, and validation. The
   CLI binary and the `mokume-rs` wheel are two front doors onto it (see
   [Architecture](architecture.md)).
2. **New computation is written in Rust first.** When a feature touches the
   computation commands — for example [native SDRF + MSstats input
   support](https://github.com/bigbio/mokume/pull/74) — it is implemented and
   tested in the Rust I/O, CLI, and pipeline crates. It does **not** need a
   matching pure-Python change to ship.
3. **The pure-Python computation package is added value.** It is kept public and
   usable so that people can plug individual functions into Python-based
   pipelines, and so it can serve as a readable reference and a parity check for
   the Rust kernel. A Python equivalent of new computation is built **only** when
   someone needs it for their pipeline or when a test requires it as a parity
   reference — and that can be done later, or agentically.
4. **Parity is checked, not assumed.** Where both implementations exist, results
   agree within floating-point tolerance (the kernel computes in `f32`); this is
   asserted for `features2proteins` in
   `rust/tests/test_rust_python_equivalence.py`.

## Support table

Both codebases implement the four computation commands today. The table records
which implementation leads and how each is maintained.

| Computation command | Rust kernel (`rust/`) | Pure-Python package (`python/mokume/`) |
| --- | --- | --- |
| `features2proteins` | ✅ Leading — authoritative | ✅ Added value · parity-checked |
| `features2peptides`  | ✅ Leading — authoritative | ✅ Added value · best-effort |
| `peptides2protein`   | ✅ Leading — authoritative | ✅ Added value · best-effort |
| `correct-batches`    | ✅ Leading — authoritative (native ComBat) | ✅ Added value · best-effort |

Legend:

- **Leading — authoritative:** new computation lands here first; defines correct
  behavior, options, and validation.
- **Added value · parity-checked:** maintained as a usable reference and covered
  by an automated equivalence test against the Rust kernel.
- **Added value · best-effort:** kept public and usable; updated when a user
  needs the function in a Python pipeline or a test needs a parity reference,
  not as a release gate.

### Backend selection

Within the pure-Python pipeline API, `RuntimeConfig.backend` selects, **for
`features2proteins` only**, whether the features-to-proteins step computes in
pure Python or delegates to the compiled Rust kernel; both paths then share the
Python post-processing layer. The other three commands have implementations in
both codebases but are not routed through `RuntimeConfig.backend`. See the
hybrid stage boundary in [Architecture](architecture.md#the-pure-python-package-and-its-compute-backends).

## Practical guidance

- **Contributing computation changes?** Implement them in the Rust crates and
  their tests. A pure-Python counterpart is optional and can follow later.
- **Need a computation function in a Python pipeline?** Prefer the `mokume-rs`
  wheel (`mokume.features2proteins(...)`, `mokume.peptides2protein(...)`, …),
  which runs the Rust kernel in-process. The pure-Python package remains
  available when you specifically want the readable Python implementation or to
  compose individual functions.
- **Filing an issue about behavior differences?** Treat the Rust kernel's result
  as the reference and the pure-Python result as the one to reconcile.
