import subprocess
import shlex
from pathlib import Path

# def get_test_data_map():
#     r1 = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_map/RB-TnSeq_mapping_library_1.1.fq.gz"
#     r2 = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline/input_map/RB-TnSeq_mapping_library_1.2.fq.gz"
#     Sdb = "/nfs/nas22/fs2202/biol_micro_sunagawa/Projects/DEV/TNSEQ_DEV/data/processed/tnseq/tnseq_pipeline//blastdb/Salmonella_genome_FQ312003.1_SL1344.fasta"
#     outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
#     expectedMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/outMap.tsv"
#     return r1, r2, Sdb, outMap, expectedMap



def get_test_data_quant():
    test_files = ["dnaid2023_13.fasta.gz", "dnaid2023_131.fasta.gz", "dnaid2023_33.fasta.gz"]
    file_root = "/nfs/home/ansintsova/test_data/processed/demux/dnaid2023"
    library_rool = ""
    faFiles = [Path(file_root)/f for f in test_files]
    libraries = []
    outMap = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testMap.tsv"
    expectedQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/outQuant.tsv"
    outQuant = "/nfs/home/ansintsova/TNSEQ_DEV/code/package/tests/data/testQuant.tsv"
    return r, outMap, expectedQuant, outQuant

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


# def test_cli_map():  # todo Pytest this
#     r1, r2, Sdb, outMap, expectedMap = get_test_data_map()
#     cmd_str = f'tnseq map -r1 {r1} -r2 {r2} -o {outMap} -db {Sdb} -f '
#     _, err, proc = capture(cmd_str)
#     if proc.returncode !=0:
#         print(to_str(err))
#     assert proc.returncode == 0
#     cmd_str2 = f'cmp --silent  {outMap} {expectedMap}'
#     out, err, proc = capture(cmd_str2)
#     assert proc.returncode == 0



def test_cli_quantify():

    r, outMap, expectedQuant, outQuant = get_test_data_quant()
    cmd_str = f'tnseq quantify -r {r}  -o {outQuant} -b {outMap} -f'
    _, err, proc = capture(cmd_str)
    assert proc.returncode == 0
    out, err, proc = capture(f'sort -k1 {outQuant}')
    with open(expectedQuant, 'r') as o:
       expectedLines = o.read()
    assert to_str(out) == expectedLines


if __name__ == "__main__":
    test_cli_map()