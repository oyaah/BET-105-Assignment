import glob

configfile: "config.yaml"

AA = config["target_aa"]

PDBS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("stride_out/*.ss.out")]

rule all:
    input:
        "final/angles.tsv"

rule calculate_angles:
    input:
        expand("contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDBS, aa=[AA])
    output:
        "final/angles.tsv"
    params:
        aa = AA
    shell:
        "python3 scripts/calc_angles.py contexts {output} {params.aa}"
