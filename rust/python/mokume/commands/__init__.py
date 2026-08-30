"""Pure-Python periphery commands for the mokume wheel.

These modules provide plots, reports, and tissue-atlas analysis outside the Rust
compute kernel. Plotting and reporting consume tables written by the
``mokume._mokume`` extension; TissueMap derives downstream analysis from QPX
data, and explicit fallbacks cover operations the kernel does not provide.

Each public module exposes ``main(argv)`` through the wheel's unified
``mokume`` console command; the ergonomic wrappers in :mod:`mokume`
(``mokume.tsne_visualization`` etc.) build the argv from keyword arguments.
"""
