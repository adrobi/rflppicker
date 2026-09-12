# Reference genome data

This directory is intentionally kept free of large reference sequences in the
GitHub source repository. The local development copy contains
`ARS-UI_Ramb_v2.0.fa` (about 2.66 GB), which is ignored by Git.

To run the application after cloning:

1. Put the required FASTA file in this directory, or keep it elsewhere and
   enter its full path in the FASTA field of the application.
2. Make sure the chromosome mapping file in `names/` matches the FASTA contig
   identifiers.
3. If the FASTA has no `.fai` index, `pyfaidx` will create one next to it.

The FASTA itself must be obtained from the source/licence that applies to your
project. Do not commit restricted or very large reference datasets to GitHub.
