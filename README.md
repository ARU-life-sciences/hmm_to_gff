# Various HMM output to GFF

Generating GFF3 output from HMM (like) programs. My use case if fairly specific, and I may make more general later on.

## Supported formats

- `nhmmscan --tblout` with the `fam` file originating from oatkDB.
- `cmscan --tblout` with the `cm` file originating from Rfam.
- `hmmscan --tblout` with the `hmm` file originating from Pfam, and the input fasta from ORFfinder.

If you do not have these specific outputs, do not expect that this will work for you!
