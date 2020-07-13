# SkipGuide
Prediction of CRISPR-Cas9 mediated Exon Skipping

Paper: TBA

## Installation
`pip install cython`

`pip install git+https://github.com/gifford-lab/skipguide.git`

## Example Usage
```
from skipguide import SkipGuide

sg = SkipGuide()

intron = 'GTAAGTTATCACCTTCGTGGCTACAGAGTTTCCTTATTTGTCTCTGTTGCCGGCTTATATGGACAAGCATATCACAGCCATTTATCGGAGCGCCTCCGTACACGCTATTATCGGACGCCTCGCGAGATCAATACGATTACCAGCTGCCCTCGTCGACCCAGGTAGCCTGGCGTGACCCCCTCCCGCTGCCCCAG'
exon = 'TTCTTCTCAGATGTGCGGGAGGCCTGATTACACATATAGACACGCGAGCAGCCATCTTTTATAGAATGGGTAGAACCCGTCCTAAGGACTCAGATTGAGCATCGTTTGCTTCTCGAGTACTACCTGGTACAGATGTCTCTTCAAACAG'

seq = intron + exon
splice_acceptor_site = len(intron)
cutsite = len(intron)
gRNA_orientation = '-'

PSI = skipguide.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)
```
