"""External-sequence helpers for the SNAF surface/B-antigen viewer.

Both functions here are used ONLY by the interactive Dash B-antigen viewer, not by
the batch CLI pipeline. They have been de-REST-ed:

  * run_emboss  -> pure-Python global alignment via Biopython PairwiseAligner
                   (replaces the EMBOSS-needle EBI REST web service + emboss.py).
  * run_ensembl -> still web (Ensembl REST) because there is no bundled offline
                   source for arbitrary Ensembl IDs, but it now degrades gracefully
                   (returns '' + a logged note) instead of raising.
"""
import os
import logging

logger = logging.getLogger(__name__)


def run_ensembl(ens):
    """Fetch a reference sequence (ENST cDNA or ENSP peptide) from Ensembl REST.
    VIEWER-ONLY. Degrades gracefully to '' on any failure so the viewer still renders."""
    try:
        import requests
        r = requests.get("https://rest.ensembl.org/sequence/id/{}?".format(ens),
                         headers={"Content-Type": "text/plain"}, timeout=20)
        return r.text
    except Exception as e:
        logger.warning("run_ensembl: could not fetch %s from Ensembl (%s); returning empty.", ens, e)
        return ''


def _read_fasta_seq(path_or_seq):
    """Accept a fasta file path or a raw sequence string; return the bare sequence."""
    try:
        if os.path.exists(str(path_or_seq)):
            seq = []
            with open(path_or_seq) as f:
                for line in f:
                    if not line.startswith('>'):
                        seq.append(line.strip())
            return ''.join(seq)
    except Exception:
        pass
    return str(path_or_seq)


def run_emboss(asequence, bsequence, python_executable=None):
    """Global (Needleman-Wunsch) protein alignment of two sequences -- PURE PYTHON via
    Biopython PairwiseAligner, replacing the EMBOSS-needle EBI REST web service.
    `asequence`/`bsequence` may be fasta file paths or raw sequence strings.
    `python_executable` is accepted for signature compatibility and ignored.
    Returns a formatted alignment string, or '' (with a note) if Biopython is missing."""
    a = _read_fasta_seq(asequence)
    b = _read_fasta_seq(bsequence)
    try:
        from Bio.Align import PairwiseAligner, substitution_matrices
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        try:
            aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
            aligner.open_gap_score = -10
            aligner.extend_gap_score = -0.5
        except Exception:
            pass  # fall back to default match/mismatch scoring
        return str(aligner.align(a, b)[0])
    except Exception as e:
        logger.warning("run_emboss: pairwise alignment unavailable (%s); returning empty.", e)
        return ''
