import warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=DeprecationWarning)
warnings.filterwarnings('ignore', category=UserWarning)

from mmsplice import MMSplice, predict_all_table
from mmsplice.vcf_dataloader import SplicingVCFDataloader

from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

import os
import sys
import tempfile
import pickle
import numpy as np

import importlib
inDelphi = importlib.import_module(
    ".models.inDelphi.inDelphi", package='skipguide')


class SkipGuide():
    """Prediction of CRISPR-Cas9 mediated Exon Skipping

    SkipGuide is a machine learning method that predicts the 
    skipping level of an exon given a Streptococcus pyogenes Cas9 
    guide RNA that targets its splice acceptor site.

    Examples
    --------
    >>> from skipguide import SkipGuide
    >>> sg = SkipGuide()
    >>> intron = 'GTAAGTTATCACCTTCGTGGCTACAGAGTTTCCTTATTTGTCTCTGTTGCCGGCTTATATGGACAAGCATATCACAGCCATTTATCGGAGCGCCTCCGTACACGCTATTATCGGACGCCTCGCGAGATCAATACGATTACCAGCTGCCCTCGTCGACCCAGGTAGCCTGGCGTGACCCCCTCCCGCTGCCCCAG'
    >>> exon = 'TTCTTCTCAGATGTGCGGGAGGCCTGATTACACATATAGACACGCGAGCAGCCATCTTTTATAGAATGGGTAGAACCCGTCCTAAGGACTCAGATTGAGCATCGTTTGCTTCTCGAGTACTACCTGGTACAGATGTCTCTTCAAACAG'
    >>> seq = intron + exon
    >>> splice_acceptor_site = len(intron)
    >>> cutsite = len(intron)
    >>> gRNA_orientation = '-'
    >>> PSI = sg.predict(seq, cutsite, splice_acceptor_site, gRNA_orientation)  # predicted percent spliced in of the exon, i.e. fraction of transcripts containing the exon
    >>> PSI
    0.12883051446158653

    References
    -------
    TBA
    """

    def __init__(self, cell_type='mESC', verbose=1):
        """
        Parameters
        ----------
        cell_type : {'mESC', 'U2OS', 'HEK293', 'HCT116', 'K562'}, default='mESC'
            The cell type. Although inDelphi accepts many cell types, the MetaSplice component
            of SkipGuide is tuned and evaluated on mESC cells, so predictions on other cell 
            types may not be reliable.
        verbose : int
            The verbosity. 0 = silent, 1 = progress print out.
        """
        self.verbose = verbose

        # Initialize inDelphi
        self._reset_inDelphi()
        inDelphi.init_model(celltype=cell_type)
        self._inDelphi_input_length = 140
        self._inDelphi_cutsite = int(self._inDelphi_input_length / 2)

        if self.verbose:
            print("Initializing MetaSplice...")

        # Initialize MetaSplice
        MetaSplice_model_path = os.path.join(os.path.dirname(
            __file__), 'models/MetaSplice/MetaSplice_model.p')
        with open(MetaSplice_model_path, 'rb') as handle:
            self._MetaSplice = pickle.load(handle)

    def predict(self, seq, cutsite, splice_acceptor_site, gRNA_orientation='+'):
        """Predicts the percent spliced in of the exon.

        Precent spliced in of an exon measures the fraction of transcripts containing the exon.

        Parameters
        ----------
        seq : str
            The sequence of the intron-exon splice acceptor region.
        cutsite : int
            The 0-index position of the cutsite, such that seq[:cutsite] and seq[cutsite:] describe the cut products.
        splice_acceptor_site : int
            The 0-index position of the intron-exon boundary, such that 
            seq[:splice_acceptor_site] gives the intron sequence and 
            seq[splice_acceptor_site:] gives the exon sequence.
        gRNA_orientation : {'-', '+'}, default='+'
            The guide RNA orientation. '|' means cutsite in the following examples.
            '+' means the NGG PAM is on the seq strand, e.g.:
                gRNA:    CATCGCGGCTCCCCTCC|AGA
                seq:  ...CATCGCGGCTCCCCTCC|AGA(NGG)...
            '-' means the NGG PAM is on the strand complement to the provided seq, e.g.:
                gRNA:         GCA|GGTGGACGTGCTGCGGG
                seq:  ...(CCN)GCA|GGTGGACGTGCTGCGGG...

        Returns
        -------
        float
            The predicted percent spliced in of the exon. One minus this value gives the predicted exon skipping frequency.
        """

        seq = seq.upper()

        # ========== inDelphi predict (P(p_i), for product p_i) ==========
        inDelphi_offset = cutsite - self._inDelphi_cutsite
        inDelphi_cutsite = self._inDelphi_cutsite
        inDelphi_seq = seq[inDelphi_offset:(
            inDelphi_offset + self._inDelphi_input_length)]

        if gRNA_orientation == '-':
            inDelphi_seq = self._reverse_complement(inDelphi_seq)
            inDelphi_cutsite = len(inDelphi_seq) - inDelphi_cutsite

        pred_df, stats = inDelphi.predict(inDelphi_seq, inDelphi_cutsite)
        pred_df = inDelphi.add_mhless_genotypes(pred_df, stats)

        # ========== MetaSplice predict (P(R|p_i), R = exon retention) ==========
        # MMSplice Intron, Acceptor, Exon scores
        if self.verbose:
            print("MMSplice Scoring...")
        fasta = tempfile.NamedTemporaryFile(mode='w')
        gtf = tempfile.NamedTemporaryFile(mode='w')
        vcf = tempfile.NamedTemporaryFile(mode='w')

        self._write_fasta(fasta, seq)
        self._write_gtf(gtf, seq, splice_acceptor_site)
        inDelphi_freqs, inDelphi_genotypes = self._write_vcf(
            vcf, seq, gRNA_orientation, inDelphi_offset, inDelphi_seq, inDelphi_cutsite, pred_df)

        dl = SplicingVCFDataloader(gtf.name, fasta.name, vcf.name)
        mmsplice = MMSplice()
        predictions = predict_all_table(
            mmsplice, dl, progress=False, pathogenicity=False, splicing_efficiency=False)

        fasta.close()
        gtf.close()
        vcf.close()

        MMSplice_scores = []
        for _, row in predictions.iterrows():
            MMSplice_scores.append(row[[
                'alt_acceptorIntron',
                'alt_acceptor',
                'alt_exon'
            ]].values)

        # SpliceAI predictions
        # MMSplice's keras model loading does not play well with SpliceAI. Need to reload SpliceAI
        # every predict call...
        # TODO: introduce a predict method that takes in multiple sequences, and have SpliceAI and MMSplice
        # process them in batches and recombine afterwards to speed this up.
        if self.verbose:
            print("SpliceAI Scoring...")
        self._load_spliceai()
        context = 10000
        max_seq_len = context + len(seq) + 1
        x = np.zeros((len(inDelphi_genotypes), max_seq_len, 4))
        for i, repair_genotype in enumerate(inDelphi_genotypes):
            input_sequence = 'N'*(context // 2) + repair_genotype
            input_sequence += 'N'*(max_seq_len - len(input_sequence))
            x[i, :, :] = one_hot_encode(input_sequence)[None, :]
        splice_ai_preds = np.mean([self._SpliceAI_models[m].predict(x, verbose=self.verbose) for m in range(5)], axis=0)
        SpliceAI_logit_scores = self._logit(np.max(splice_ai_preds[:, :, 1], axis=1).reshape(-1, 1))
        
        # MetaSplice combines MMSplice and SpliceAI scores into a single P(R|p_i) prediction
        if self.verbose:
            print("MetaSplice Predicting...")
        X = np.hstack((np.array(MMSplice_scores), SpliceAI_logit_scores))
        y_preds = self._expit(self._MetaSplice.predict(X))

        # SkipGuide final weighted average (sum(P(p_i)P(R|p_i) over all i) = P(R) = PSI)
        return np.average(y_preds, weights=inDelphi_freqs)

    def _load_spliceai(self):
        paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
        self._SpliceAI_models = [load_model(resource_filename('spliceai', x)) for x in paths]

    def _reset_inDelphi(self):
        inDelphi.init_flag = False
        inDelphi.nn_params = None
        inDelphi.nn2_params = None
        inDelphi.normalizer = None
        inDelphi.rate_model = None
        inDelphi.bp_model = None
        inDelphi.CELLTYPE = None

    def _reverse_complement(self, seq):
        rev = {'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G'}
        seq = seq.strip().upper()
        return ''.join([rev[seq[i]] if seq[i] in rev else seq[i] for i in range(len(seq) - 1, -1, -1)])

    def _expit(self, x):
        return 1. / (1. + np.exp(-np.array(x)))

    def _logit(self, x):
        x = np.clip(x, 0.00001, 0.99999)
        return np.log(x) - np.log(1 - x)

    def _write_fasta(self, f, seq):
        f.write('>0\n')
        f.write(seq + '\n')
        f.seek(0)

    def _write_gtf(self, f, seq, splice_acceptor_site):
        # Gene
        f.write('\t'.join([
            '0', 'artificial', 'gene', '1', str(len(seq)), '.', '+', '.',
            'gene_id "g0"; transcript_id ""; gene_name "g0";\n'
        ]))

        # Transcript
        f.write('\t'.join([
            '0', 'artificial', 'transcript', '1', str(len(seq)), '.', '+', '.',
            'gene_id "g0"; transcript_id "t0"; gene_name "g0";\n'
        ]))

        # Exon
        f.write('\t'.join([
            '0', 'artificial', 'exon', str(
                splice_acceptor_site + 1), str(len(seq)), '.', '+', '.',
            'gene_id "g0"; transcript_id "t0"; gene_name "g0"; exon_id "e0";\n'
        ]))

        f.seek(0)

    def _write_vcf(self, f, seq, g_orientation, inDelphi_offset, inDelphi_seq, inDelphi_cutsite, inDelphi_preds):
        inDelphi_freqs = []
        inDelphi_genotypes = []

        f.write('##fileformat=VCFv4.0\n')
        f.write('##contig=<ID=0,length={0}>\n'.format(len(seq)))
        f.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF',
                           'ALT', 'QUAL', 'FILTER', 'INFO\n']))

        # Insertions
        insertion_pred_df = inDelphi_preds.loc[inDelphi_preds['Category'] == 'ins']
        for _, row in insertion_pred_df.iterrows():
            inserted_base = row['Inserted Bases'] if g_orientation == '+' else self._reverse_complement(
                row['Inserted Bases'])
            frequency = row['Predicted frequency']
            inDelphi_freqs.append(frequency)

            pos = inDelphi_offset + \
                (inDelphi_cutsite if g_orientation ==
                 '+' else (len(inDelphi_seq) - inDelphi_cutsite))
            ref = seq[pos - 1]
            alt = seq[pos - 1] + inserted_base

            f.write('\t'.join([
                '0', str(pos), '.', ref, alt, '.', '.', '.', '\n'
            ]))

            genotype = seq[0:pos] + inserted_base + seq[0:pos]
            inDelphi_genotypes.append(genotype)

        # Deletions
        deletion_pred_df = inDelphi_preds.loc[inDelphi_preds['Category'] == 'del']
        for _, row in deletion_pred_df.iterrows():
            deletion_size = row['Length']
            frequency = row['Predicted frequency']
            inDelphi_freqs.append(frequency)

            genotype_pos = row['Genotype position']
            deletion_genotype = inDelphi_seq[(
                inDelphi_cutsite + genotype_pos - deletion_size):(inDelphi_cutsite + genotype_pos + 1)]
            if g_orientation == '+':
                pos = inDelphi_offset + inDelphi_cutsite + genotype_pos - deletion_size
            else:
                deletion_genotype = self._reverse_complement(deletion_genotype)
                pos = inDelphi_offset + \
                    (len(inDelphi_seq) - (inDelphi_cutsite + genotype_pos))

            ref = seq[pos - 1] + deletion_genotype
            alt = seq[pos - 1]

            f.write('\t'.join([
                '0', str(pos), '.', ref, alt, '.', '.', '.', '\n'
            ]))

            genotype = seq[0:pos] + seq[(pos + deletion_size):]
            inDelphi_genotypes.append(genotype)

        f.seek(0)

        return inDelphi_freqs, inDelphi_genotypes
