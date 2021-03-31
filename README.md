# SkipGuide
Prediction of CRISPR-Cas9 mediated Exon Skipping

Paper: [Machine learning based CRISPR gRNA design for therapeutic exon skipping](https://doi.org/10.1371/journal.pcbi.1008605)

## Installation (Python Package)
```shell
pip install cython
```

```shell
pip install git+https://github.com/gifford-lab/skipguide.git
```

inDelphi requires version scikit-learn version 0.20.0. Although the MMSplice package requires scikit-learn version 0.19.2, it'll still work with version 0.20.0. Make sure version 0.20.0 is installed:
```shell
pip install scikit-learn==0.20.0 --no-deps
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
