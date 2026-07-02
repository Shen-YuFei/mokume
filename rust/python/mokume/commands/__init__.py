"""Pure-Python periphery commands for the mokume wheel.

These modules draw the plots / reports / tissue atlas that are not part of the
Rust compute kernel. The numbers are produced by the kernel (the ``mokume``
extension or the ``mokume`` CLI binary); these commands read those tables and
render the figures, so the compute stays single-sourced in Rust.

Each module exposes ``main(argv)`` and is runnable as
``python -m mokume.commands.<name>``; the ergonomic wrappers in
:mod:`mokume` (``mokume.tsne_visualization`` etc.) build the argv from keyword
arguments.
"""
