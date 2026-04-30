# parse stride output and extract tripeptide contexts for a target amino acid
import sys
import re
import csv
import os

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

def read_stride(filepath):
    residues = []
    with open(filepath) as f:
        for line in f:
            if not line.startswith("ASG"):
                continue
            cols = line.split()
            if len(cols) < 10:
                continue
            m = re.match(r"\d+", cols[3])
            if m is None:
                continue
            try:
                residues.append({
                    "name3": cols[1],
                    "ch": cols[2],
                    "rnum": int(m.group()),
                    "seq_idx": int(cols[4]),
                    "ss": cols[5],
                    "ss_full": cols[6],
                    "phi": float(cols[7]),
                    "psi": float(cols[8]),
                    "asa": float(cols[9]),
                })
            except (ValueError, IndexError):
                continue
    return residues


def build_triplets(residues, target_3letter, pdb_name):
    rows = []
    for k in range(1, len(residues) - 1):
        cur = residues[k]
        if cur["name3"] != target_3letter:
            continue

        prev = residues[k - 1]
        nxt = residues[k + 1]

        # triplet sequence like "MRR" and secondary structure like "HHH"
        seq_triple = "".join(
            THREE_TO_ONE.get(r["name3"], "X") for r in [prev, cur, nxt]
        )
        ss_triple = prev["ss"] + cur["ss"] + nxt["ss"]

        pos_label = (
            f"{cur['ch']}:{prev['rnum']},"
            f"{cur['ch']}:{cur['rnum']},"
            f"{cur['ch']}:{nxt['rnum']}"
        )

        # output all 3 residues of the triplet
        for r in [prev, cur, nxt]:
            rows.append([
                r["name3"], r["ch"], r["rnum"], r["seq_idx"],
                r["ss"], r["ss_full"], r["phi"], r["psi"], r["asa"],
                THREE_TO_ONE.get(r["name3"], "X"),
                seq_triple, ss_triple, pdb_name, pos_label,
            ])
    return rows


if __name__ == "__main__":
    stride_file = sys.argv[1]
    out_file = sys.argv[2]
    aa_target = sys.argv[3]

    pdb_name = os.path.basename(stride_file).split(".")[0]

    residue_list = read_stride(stride_file)
    triplet_rows = build_triplets(residue_list, aa_target, pdb_name)

    with open(out_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        for row in triplet_rows:
            writer.writerow(row)
