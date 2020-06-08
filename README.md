# SkipGuide
Info and Doc TBA

## Installation
`pip install cyvcf2 cython`

`pip install mmsplice==1.0.3`

`pip install scikit-learn==0.20.0`

`pip install git+https://github.com/gifford-lab/skipguide.git`

## Example Usage
```
intron = 'GTAAGTTATCACCTTCGTGGCTACAGAGTTTCCTTATTTGTCTCTGTTGCCGGCTTATATGGACAAGCATATCACAGCCATTTATCGGAGCGCCTCCGTACACGCTATTATCGGACGCCTCGCGAGATCAATACGATTACCAGCTGCCCTCGTCGACCCAGGTAGCCTGGCGTGACCCCCTCCCGCTGCCCCAG'

exon = 'TTCTTCTCAGATGTGCGGGAGGCCTGATTACACATATAGACACGCGAGCAGCCATCTTTTATAGAATGGGTAGAACCCGTCCTAAGGACTCAGATTGAGCATCGTTTGCTTCTCGAGTACTACCTGGTACAGATGTCTCTTCAAACAG'

seq = intron + exon
splice_acceptor_site = len(intron)
cutsite = len(intron)
gRNA_orientation = '-'

PSI = skipguide.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)
```
