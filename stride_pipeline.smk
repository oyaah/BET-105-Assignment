import os

configfile: "config.yaml"

PDB_PATH = config.get("pdb_dir", "pdbs")

# Read PDB list from a batch file (one ID per line)
BATCH_FILE = config.get("batch_file", "batch.txt")

with open(BATCH_FILE) as fh:
    PDB_IDS = [line.strip() for line in fh if line.strip()]

rule all:
    input:
        expand("stride_out/{pdb}.ss.out", pdb=PDB_IDS)

rule decompress_pdb:
    input:
        gz = PDB_PATH + "/{pdb}.pdb.gz"
    output:
        pdb = temp("unzipped_pdbs/{pdb}.pdb")
    shell:
        "gzip -dc {input.gz} > {output.pdb}"

rule run_stride:
    input:
        pdb = "unzipped_pdbs/{pdb}.pdb"
    output:
        ss = "stride_out/{pdb}.ss.out"
    shell:
        "stride {input.pdb} > {output.ss}"
