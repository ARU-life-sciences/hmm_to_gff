# Format the output of hmmscan to a GFF file

Specifically we are interested in the output from ORFfinder, where the pipeline is:

```
fasta -> ORFfinder -> hmmscan -> GFF
```
We use the Pfam HMM database to achieve this.

Or when we are running cmscan (Infernal):

```
fasta -> cmscan -> GFF
```

Example of metadata which indicates this script will work on the file. No guarantees, no tests (yet).

```
#
# Program:         hmmscan
# Version:         3.4 (Aug 2023)
# Pipeline mode:   SCAN
# Query file:      ../data/Acaena_ovalifolia_orfs.fasta
# Target file:     ../pfam_hmm/Pfam-A.hmm
# Options:         /software/team301/hmmer-3.4/src/hmmscan --tblout ../data/pfam_hmmscans/Acaena_ovalifolia.tbl --cpu 16 ../pfam_hmm/Pfam-A.hmm ../data/Acaena_ovalifolia_orfs.fasta
# Current dir:     /lustre/scratch123/tol/teams/blaxter/users/mb39/ARU/mito_structural_variation/annotation/src
# Date:            Fri Feb  7 19 13 12 2025
# [ok]
```

```
# Program:         cmscan
# Version:         1.1.4 (Dec 2020)
# Pipeline mode:   SCAN
# Query file:      ../data/Acaena_ovalifolia.fasta
# Target file:     ../rfam_introns/group_II_intron.cm
# Option settings: /software/team301/infernal-1.1.4-linux-intel-gcc/binaries/cmscan --tblout test.tbl ../rfam_introns/group_II_intron.cm ../data/Acaena_ovalifolia.fasta 
# Current dir:     /lustre/scratch123/tol/teams/blaxter/users/mb39/ARU/mito_structural_variation/annotation/src
# Date:            Tue Feb 11 14:09:12 2025
# [ok]
```
