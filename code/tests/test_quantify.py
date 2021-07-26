import random
from collections import Counter
from Bio import SeqIO
from tnseq import sequence
from tnseq import tnseq
from tnseq import extract_barcodes
from tnseq import quantify_barcodes
import json
import subprocess
import shlex


def test_quantify_load_fq_barcodes():
    in_fq_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/test100.fasta"
    cnter = quantify_barcodes.quantify_load_fq_barcodes(in_fq_file)
    expectedCounter = Counter({'ATTGTTCTATGCCTGCC': 2, 'ATCCATCAAGCAAACGC': 2, 'CGAAAAGATATACCAAT': 2, 'TCGCGATGTAATATATA': 2,
                               'GTCACGCGCCCGCGCCA': 1, 'AGACCGATGCGTCAAGG': 1, 'GAAAGCCCGAGATCGAT': 1, 'TCGCCGAGCGATTTTTA': 1,
                               'AGGAAACCAAATATAAA': 1, 'TGCAAAGTGAGGGCTAC': 1, 'AGAGCACACGACCGGTG': 1, 'TACATCGCTAGTGAACT': 1,
                               'GTGCGCTATGCAACGAC': 1, 'TGGCACGAAAGCAAGAT': 1, 'CTGCATATCTAGCGCGG': 1, 'CCAAGAGCACGAGACCA': 1,
                               'GCGGTTCGGAGAGACAG': 1, 'ACGGATCCGACCCGGGG': 1, 'TGGAAACCCTATCGCTA': 1, 'AGATGCTCGCGTGCGGT': 1,
                               'CACATTGGAGTCTACTG': 1, 'CGCAATGGGGCGTGGGG': 1, 'CTCCCAGCTTATCAGAG': 1, 'GACCAATACGCGGGGAG': 1,
                               'CCGTCTGGCTGAGACTG': 1, 'TCATTAAAGGCTTGGTG': 1, 'ATAAGCTTTAATAGTTA': 1, 'AGGAGTAGCCGAGTGCT': 1,
                               'AGTATGCAGACATGTGT': 1, 'TAAGGACGCTACCGTAC': 1, 'GCAAAGTAATCACAGAG': 1, 'CAATATCATCCTCGACT': 1,
                               'CCCCAAGGCTGCATCCA': 1, 'CTTTATACGCGAACTCG': 1, 'TCGCGCCTAGAGCTATG': 1, 'GCGCACGACGGTACGCC': 1,
                               'CAAAAATGATCGCGTAG': 1, 'ACGAAGACCAGCTATGG': 1, 'GATTTCAGACTACCTGC': 1, 'GAGGGTGATTCACGACT': 1,
                               'CATAAGCAGGCCAACAG': 1, 'CGACTTGTAGATCTCTG': 1, 'GATGGGCCTGAGGCAGA': 1, 'TAATGGATGAATTTGTG': 1,
                               'AGCCCAGCATACAGAGA': 1, 'CCCCATGACAGCGTTAT': 1, 'TTCCCAGCAAAATACGG': 1, 'GCCTCCACACGAGTACA': 1,
                               'CGTTATGTACATGTTCC': 1, 'CCATCGCCTGTAACGAT': 1, 'GTAACGCTGAACAAAAA': 1, 'AAAACCTCCCTGCCCAT': 1,
                               'TTGAAGACGGTACCTGT': 1, 'CGACCATTGGTAGACAC': 1, 'GCCGTAATCCTTAGAAG': 1, 'TCAGAGCGTATTCATCC': 1,
                               'TATCAGCGCGCGACTTA': 1, 'GAATTTCAGCTGGCAAG': 1, 'CAAACCTGATTAAATCA': 1, 'CATAGCTCTTGTTACAC': 1,
                               'GCGAACTTTATGGCAGA': 1, 'AGTACGGGTACAATGCG': 1, 'CGATCACCCTCAGTATG': 1, 'CACGTATCCTCGTGGAC': 1,
                               'AGGGTACACGCAGCGCG': 1, 'GAGGCCATGGGCAGTCC': 1, 'GCGGTGACAGGATCGGA': 1, 'CCTGAAAACCTTCACTC': 1,
                               'ACGCGCCAGACTTACGC': 1, 'ATGGCGCCTTCCGACGT': 1, 'TCCTTTAGGGGCGAATG': 1, 'AGACACTCCATTTAAAG': 1,
                               'CTCCCAAACACGAGAAT': 1, 'TTACGAGCTCATGAGCA': 1, 'GGTGGCACCTATAACAA': 1, 'CTCAAGGCCCGACGGGC': 1,
                               'TATTAAGGTAATGTAGA': 1, 'TACCGTAAATGCAAATA': 1, 'GGTGACTAGCCGGTAGT': 1, 'GCAAACTGTTACATAAG': 1,
                               'CAGCCTGTCTGCGACAT': 1, 'TAAACTACCTTTGACCA': 1, 'TACATACTGATGCCCCT': 1, 'AAGTAACCAGTCGAAGA': 1,
                               'AAGTATATGGTGGTATA': 1, 'AAGAGCCGGCACGCAAC': 1, 'CACAAGCAGTCAAACAT': 1, 'CGATACCCGTAACGCGT': 1,
                               'ACCCGCCTGAGCTAACA': 1, 'TCCAGACCTTTGCGCGA': 1, 'CGTATGTATCGGCCAAG': 1, 'TAGCATGGGGGGCTGAA': 1,
                               'AACCTGAGAACCGCTCT': 1, 'AGAATCCTCAAACTATA': 1, 'GATAGCTTGATGACGCA': 1, 'ATGGACTACTGCACGCG': 1})
    print(cnter==expectedCounter)

def test_quantify_read_barcode_map_files():
    outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"

    barcode_2_pos, barcode_2_abundance = quantify_barcodes.quantify_read_barcode_map_files(outMap)
    with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2P.json') as jf:
        expectedB2P = json.load(jf)
        expectedB2P = {key: tuple(val) for key, val in expectedB2P.items()}
    with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2A.json') as jf2:
        expectedB2A = json.load(jf2)
    assert barcode_2_pos == expectedB2P
    assert barcode_2_abundance == expectedB2A



def capture(command_str):
    command = shlex.split(command_str)
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,)
    out, err = proc.communicate()
    return out, err, proc

def to_str(bytes_or_str):
    if isinstance(bytes_or_str, bytes):
        value = bytes_or_str.decode('utf-8')
    else:
        value = bytes_or_str
    return value

def test_quantify_extract_annotated_correct():
    with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2P.json') as jf:
        expectedB2P = json.load(jf)
        barcode_2_position = {key: tuple(val) for key, val in expectedB2P.items()}
    with open('/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/B2A.json') as jf2:
        barcode_2_abundance = Counter(json.load(jf2))
    print(barcode_2_abundance.most_common()[0:10])
    in_fq_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/test100.fasta"
    cnter = quantify_barcodes.quantify_load_fq_barcodes(in_fq_file)
    max_edit_distance = 2
    expectedOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/quantifyExtract.sorted.out"
    testOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/quantifyExtractTest.out"
    quantify_barcodes.quantify_extract_annotated_correct(barcode_2_position, barcode_2_abundance, cnter, testOutFile, max_edit_distance)
    out, err, proc = capture(f'sort -k1 {testOutFile}')
    with open(expectedOutFile, 'r') as o:
        expectedLines = o.read()
    assert to_str(out) == expectedLines


def test_quantify():
    in_fq_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/test100.fasta"
    map_file = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
    testOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/quantifyQuantOut.tsv"
    expectedOutFile = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/seqData/quantifyExtract.sorted.out"
    tp2 = 'GTGTATAAGAGACAG'
    bc2tp2 = 13
    bcLen = 17
    before = True
    max_edit_distance = 2
    quantify_barcodes.quantify(in_fq_file, map_file, testOutFile, tp2, bc2tp2, bcLen, before, max_edit_distance)
    out, err, proc = capture(f'sort -k1 {testOutFile}')
    with open(expectedOutFile, 'r') as o:
        expectedLines = o.read()
    assert to_str(out) == expectedLines


