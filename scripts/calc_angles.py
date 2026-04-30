# calculate angles between CA->centroid vectors for tripeptide contexts
# only considers HHH (all-helix) triplets
# classifies by size of the LEFT neighbor residue

import os
import sys
import gzip
import csv
import numpy as np
from Bio.PDB import PDBParser
from tqdm import tqdm

# size classification from prof's instructions
RESIDUE_SIZE = {
    "G": "Tiny",   "A": "Tiny",
    "V": "Small",  "P": "Small",  "S": "Small",  "T": "Small",  "C": "Small",
    "N": "Intermediate", "L": "Intermediate", "I": "Intermediate", "D": "Intermediate",
    "K": "Large",  "E": "Large",  "M": "Large",  "Q": "Large",  "H": "Large",
    "R": "Bulky",  "F": "Bulky",  "Y": "Bulky",  "W": "Bulky",
}

PDB_DIR = os.environ.get("PDB_DIR", "pdbs")


# signed angle between two vectors given a reference axis
# code provided by prof
def signed_angle_3d(v1, v2, axis):
    v1 = np.array(v1, dtype=float)
    v2 = np.array(v2, dtype=float)
    axis = np.array(axis, dtype=float)
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    axis /= np.linalg.norm(axis)
    cross = np.cross(v1, v2)
    dot = np.dot(v1, v2)
    angle = np.arctan2(np.dot(cross, axis), dot)
    return np.degrees(angle)


def load_pdb(pdb_id):
    path = os.path.join(PDB_DIR, f"{pdb_id}.pdb.gz")
    if not os.path.isfile(path):
        return None
    parser = PDBParser(QUIET=True)
    with gzip.open(path, "rt") as f:
        return parser.get_structure(pdb_id, f)


def get_residue(struct, chain_id, resnum):
    for model in struct:
        try:
            chain = model[chain_id]
        except KeyError:
            continue
        for res in chain:
            if res.id[1] == resnum:
                return res
    return None


def get_ca(res):
    if res is not None and "CA" in res:
        return res["CA"].get_coord()
    return None


def get_centroid(res):
    # centroid of sidechain atoms only (exclude backbone N, CA, C, O)
    backbone = {"N", "CA", "C", "O"}
    coords = [a.get_coord() for a in res if a.get_name() not in backbone]
    if not coords:
        return None
    return np.mean(coords, axis=0)


def process_all(ctx_dir, target_aa, out_path):
    files = [f for f in os.listdir(ctx_dir) if f.endswith(".tsv")]
    results = []

    for fname in tqdm(files, desc="processing"):
        fpath = os.path.join(ctx_dir, fname)
        if os.path.getsize(fpath) == 0:
            continue

        try:
            with open(fpath) as f:
                rows = list(csv.reader(f, delimiter="\t"))
        except:
            continue

        if not rows:
            continue

        pdb_id = rows[0][12]
        struct = load_pdb(pdb_id)
        if struct is None:
            continue

        # each triplet = 3 consecutive rows
        i = 0
        while i + 2 < len(rows):
            rp, rc, rn = rows[i], rows[i+1], rows[i+2]
            i += 3

            # only want target AA in center
            if rc[0] != target_aa:
                continue
            # only helix-helix-helix
            if rc[11] != "HHH":
                continue

            ch = rc[1]
            try:
                nrp, nrc, nrn = int(rp[2]), int(rc[2]), int(rn[2])
            except (ValueError, IndexError):
                continue

            resp = get_residue(struct, ch, nrp)
            resc = get_residue(struct, ch, nrc)
            resn = get_residue(struct, ch, nrn)
            if not all([resp, resc, resn]):
                continue

            ca_p = get_ca(resp)
            ca_c = get_ca(resc)
            ca_n = get_ca(resn)
            sc_p = get_centroid(resp)
            sc_c = get_centroid(resc)

            if any(x is None for x in [ca_p, ca_c, ca_n, sc_p, sc_c]):
                continue

            # vectors: CA -> sidechain centroid
            v1 = sc_p - ca_p
            v2 = sc_c - ca_c
            axis = ca_c - ca_p  # along the backbone

            if np.linalg.norm(axis) < 1e-9:
                continue

            try:
                angle = signed_angle_3d(v1, v2, axis)
            except:
                continue

            # classify by left (previous) residue size
            left_aa = rp[9]
            size = RESIDUE_SIZE.get(left_aa)
            if size is None:
                continue

            results.append([pdb_id, left_aa, size, f"{angle:.4f}"])

    # save angles
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["pdb", "left_aa", "size_class", "angle"])
        for r in results:
            w.writerow(r)

    # also save list of PDBs that contributed data
    pdb_list_path = os.path.join(os.path.dirname(out_path), "pdb_list.txt")
    used = sorted(set(r[0] for r in results))
    with open(pdb_list_path, "w") as f:
        for p in used:
            f.write(p + "\n")

    print(f"saved {len(results)} angles to {out_path}")
    print(f"saved {len(used)} pdb names to {pdb_list_path}")


if __name__ == "__main__":
    ctx_dir = sys.argv[1]
    out = sys.argv[2]
    aa = sys.argv[3]
    process_all(ctx_dir, aa, out)
