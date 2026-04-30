import glob

configfile: "config.yaml"

AA = config["target_aa"]

PDBS = [f.split("/")[-1].replace(".ss.out", "") for f in glob.glob("stride_out/*.ss.out")]

rule all:
    input:
        expand("contexts/context_for_{aa}_in_{pdb}.tsv", pdb=PDBS, aa=[AA])

rule extract_context:
    input:
        ss = "stride_out/{pdb}.ss.out"
    output:
        tsv = "contexts/context_for_{aa}_in_{pdb}.tsv"
    params:
        aa = "{aa}"
    shell:
        "python3 scripts/parse_stride.py {input.ss} {output.tsv} {params.aa}"
