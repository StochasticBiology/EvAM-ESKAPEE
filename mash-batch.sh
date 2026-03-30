#!/bin/bash

# exit if critical commands are missing
command -v ncbi-genome-download >/dev/null 2>&1 || { echo "ncbi-genome-download not found"; exit 1; }
command -v mash >/dev/null 2>&1 || { echo "mash not found"; exit 1; }

# input list of assemblies
IDFILE="gca_ids.txt"

# genomes per batch
BATCHSIZE=400

# directories
TEMPDIR="temp"
SKETCHDIR="sketches"

mkdir -p "$TEMPDIR"
mkdir -p "$SKETCHDIR"

# split accession list
split -l $BATCHSIZE "$IDFILE" batch_

for chunk in batch_*; do

    echo "Processing batch $chunk"

    # convert newline list to comma separated
    IDS=$(paste -sd, "$chunk")

    # download genomes from GenBank
    ncbi-genome-download bacteria \
        -s genbank \
        --assembly-accessions "$IDS" \
        -F fasta \
        --retries 3 \
        -o "$TEMPDIR"

    # sketch each genome
    find "$TEMPDIR" -name "*.fna.gz" | while read fasta; do

        base=$(basename "$fasta" .fna.gz)
        sketch="$SKETCHDIR/$base.msh"

        # skip if sketch already exists
        if [ -f "$sketch" ]; then
            echo "Sketch exists, skipping $base"
            continue
        fi

        echo "Sketching $base"
        mash sketch -s 1000 -o "$SKETCHDIR/$base" "$fasta"

    done

    # remove downloaded genomes
    rm -rf "$TEMPDIR"/*

done

echo "All batches finished"

# assumes all .msh files are in 'sketches/' folder
mash triangle sketches/*.msh > mash_distances.tab
