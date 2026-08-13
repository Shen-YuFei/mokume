"""Pure-Python periphery commands for the mokume wheel.

These modules provide plots, reports, and tissue-atlas analysis outside the Rust
compute kernel. Plotting and reporting consume tables written by the
``mokume._mokume`` extension or the standalone ``mokume`` CLI; TissueMap derives
downstream analysis from QPX data, and explicit fallbacks cover operations the
kernel does not provide.

Each module exposes ``main(argv)`` and is runnable as
``python -m mokume.commands.<name>``; the ergonomic wrappers in
:mod:`mokume` (``mokume.tsne_visualization`` etc.) build the argv from keyword
arguments.
"""
