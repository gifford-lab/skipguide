# SkipGuide
Prediction of CRISPR-Cas9 mediated Exon Skipping

Paper: TBA

## Installation (Python Package)
```shell
pip install cython
```

```shell
pip install git+https://github.com/gifford-lab/skipguide.git
```

## Example Usage
Please refer to [`skipguide/skipguide.py`](skipguide/skipguide.py) for documentation.

```python
from skipguide import SkipGuide

sg = SkipGuide()

intron = 'GTAAGTTATCACCTTCGTGGCTACAGAGTTTCCTTATTTGTCTCTGTTGCCGGCTTATATGGACAAGCATATCACAGCCATTTATCGGAGCGCCTCCGTACACGCTATTATCGGACGCCTCGCGAGATCAATACGATTACCAGCTGCCCTCGTCGACCCAGGTAGCCTGGCGTGACCCCCTCCCGCTGCCCCAG'
exon = 'TTCTTCTCAGATGTGCGGGAGGCCTGATTACACATATAGACACGCGAGCAGCCATCTTTTATAGAATGGGTAGAACCCGTCCTAAGGACTCAGATTGAGCATCGTTTGCTTCTCGAGTACTACCTGGTACAGATGTCTCTTCAAACAG'

seq = intron + exon
splice_acceptor_site = len(intron)
cutsite = len(intron)
gRNA_orientation = '-'

# The predicted percent spliced in of the exon, which measures the fraction of transcripts containing the exon.
# One minus this value gives the predicted exon skipping frequency.
PSI = sg.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)
```
