# Key Concepts

Understanding the core concepts behind mokume helps you choose the right methods for your proteomics experiment.

## Quantification

Protein quantification transforms peptide-level mass spectrometry measurements into protein-level abundance estimates. mokume supports multiple approaches, each with different trade-offs.

[:octicons-arrow-right-24: Quantification Methods](quantification.md)

## Normalization

Normalization corrects systematic biases between runs and samples so that intensity differences reflect true biological variation rather than technical artifacts.

[:octicons-arrow-right-24: Normalization](normalization.md)

## Batch Correction

When samples are processed in multiple batches (different days, instruments, or labs), batch effects can dominate biological signal. ComBat-based correction removes these while preserving biology.

[:octicons-arrow-right-24: Batch Correction](batch-correction.md)

## IRS Normalization

For multi-plex TMT experiments, Internal Reference Scaling uses shared reference channels across plexes to make protein intensities comparable.

[:octicons-arrow-right-24: IRS Normalization](irs.md)

## Preprocessing Filters

Quality control filters remove low-quality features before quantification. mokume provides a comprehensive, configurable filter system.

[:octicons-arrow-right-24: Preprocessing Filters](preprocessing.md)
