import glob
import os

configfile: "config.yaml"

AA = config["target_aa"]
PDB_DIR = config.get("pdb_dir", "pdbs")

# discover PDB files
PDB_IDS = [os.path.basename(f).replace(".pdb.gz", "")
           for f in glob.glob(os.path.join(PDB_DIR, "*.pdb.gz"))]

rule all:
    input:
        "final/angle_plot.png"

# --- step 1: decompress and run stride ---

rule decompress_pdb:
    input:
        gz = PDB_DIR + "/{pdb}.pdb.gz"
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

# --- step 2: extract tripeptide contexts ---

rule extract_context:
    input:
        ss = "stride_out/{pdb}.ss.out"
    output:
        tsv = "contexts/context_for_{pdb}.tsv"
    params:
        aa = AA
    shell:
        "python3 scripts/parse_stride.py {input.ss} {output.tsv} {params.aa}"

# --- step 3: compute angles ---

rule calculate_angles:
    input:
        contexts = expand("contexts/context_for_{pdb}.tsv", pdb=PDB_IDS)
    output:
        tsv = "final/angles.tsv"
    params:
        aa = AA
    shell:
        "python3 scripts/calc_angles.py contexts {output.tsv} {params.aa}"

# --- step 4: plot ---

rule plot:
    input:
        tsv = "final/angles.tsv"
    output:
        png = "final/angle_plot.png"
    shell:
        "python3 scripts/plot_angles.py"
