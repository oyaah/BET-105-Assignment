#!/bin/bash
# runs the whole pipeline end to end
# usage: bash scripts/run_all.sh /path/to/pdbs ARG

PDB_DIR="${1:-pdbs}"
TARGET_AA="${2:-ARG}"
BATCH_SIZE=1000

export PDB_DIR

mkdir -p unzipped_pdbs stride_out contexts final

echo "generating pdb list..."
find "$PDB_DIR" -name "*.pdb.gz" -exec basename {} .pdb.gz \; > all_pdbs.txt
total=$(wc -l < all_pdbs.txt | tr -d ' ')
echo "found $total pdb files"

# step 1: run stride in batches (otherwise snakemake DAG takes forever)
echo "running stride..."
batch_num=0
while IFS= read -r line; do
    echo "$line" >> "batch_${batch_num}.txt"
    count=$(wc -l < "batch_${batch_num}.txt" | tr -d ' ')
    if [ "$count" -ge "$BATCH_SIZE" ]; then
        echo "  batch $batch_num..."
        snakemake -s stride_pipeline.smk --cores 4 --keep-going \
            --config batch_file="batch_${batch_num}.txt" pdb_dir="$PDB_DIR" 2>/dev/null
        rm "batch_${batch_num}.txt"
        batch_num=$((batch_num + 1))
    fi
done < all_pdbs.txt

if [ -f "batch_${batch_num}.txt" ]; then
    echo "  final batch..."
    snakemake -s stride_pipeline.smk --cores 4 --keep-going \
        --config batch_file="batch_${batch_num}.txt" pdb_dir="$PDB_DIR" 2>/dev/null
    rm "batch_${batch_num}.txt"
fi

echo "stride done: $(find stride_out/ -name '*.ss.out' | wc -l | tr -d ' ') / $total"

# step 2: extract tripeptide contexts
echo "extracting contexts..."
count=0
for ss in stride_out/*.ss.out; do
    pdb_id=$(basename "$ss" .ss.out)
    outf="contexts/context_for_${TARGET_AA}_in_${pdb_id}.tsv"
    [ -f "$outf" ] && continue
    python3 scripts/parse_stride.py "$ss" "$outf" "$TARGET_AA" 2>/dev/null
    count=$((count + 1))
    [ $((count % 5000)) -eq 0 ] && echo "  $count done..."
done
echo "contexts: $(find contexts/ -name '*.tsv' | wc -l | tr -d ' ')"

# step 3: calculate angles
echo "calculating angles..."
python3 scripts/calc_angles.py contexts final/angles.tsv "$TARGET_AA"

# step 4: plot
echo "plotting..."
python3 scripts/plot_angles.py

rm -f all_pdbs.txt batch_*.txt
echo "done!"
