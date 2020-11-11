import unittest

from skipguide import SkipGuide

class TestSkipGuide(unittest.TestCase):
    def test_skipguide(self):
        sg = SkipGuide()
        intron = 'GTAAGTTATCACCTTCGTGGCTACAGAGTTTCCTTATTTGTCTCTGTTGCCGGCTTATATGGACAAGCATATCACAGCCATTTATCGGAGCGCCTCCGTACACGCTATTATCGGACGCCTCGCGAGATCAATACGATTACCAGCTGCCCTCGTCGACCCAGGTAGCCTGGCGTGACCCCCTCCCGCTGCCCCAG'
        exon = 'TTCTTCTCAGATGTGCGGGAGGCCTGATTACACATATAGACACGCGAGCAGCCATCTTTTATAGAATGGGTAGAACCCGTCCTAAGGACTCAGATTGAGCATCGTTTGCTTCTCGAGTACTACCTGGTACAGATGTCTCTTCAAACAG'

        seq = intron + exon
        splice_acceptor_site = len(intron)
        cutsite = len(intron)
        gRNA_orientation = '-'

        PSI = sg.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)
        print(PSI)
        PSI2 = sg.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)
        print(PSI2)

        self.assertTrue(PSI < 0.15)

if __name__ == '__main__':
    unittest.main()