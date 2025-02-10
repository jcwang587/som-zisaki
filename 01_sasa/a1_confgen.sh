#!/bin/bash

# Start a local job server
/apps/prod/COMPCHEM/schrodinger/schrodinger/jsc local-server-start

# Directory containing .mae files
mae_dir="/home/wanjiach/github/som-zisaki/data/mae"
conformer_dir="/home/wanjiach/github/som-zisaki/data/conformer"

# Create the maegz_c directory if it doesn't exist
mkdir -p "$conformer_dir"
rm -f "$conformer_dir"/*

# Copy .mae files to maegz_dir
for mae_file in "$mae_dir"/*.mae; do
    base_name=$(basename "$mae_file" .mae)
    cp "$mae_file" "$conformer_dir/${base_name}.mae"
done

# Find all .mae files and run confgen on each
for mae_file in "$conformer_dir"/*.mae; do
    echo "Processing $mae_file..."
    /apps/prod/COMPCHEM/schrodinger/schrodinger/confgen "$mae_file" -n 1000 -m 1000 -optimize -WAIT

    generated_maegz=$(basename "$mae_file" .mae)-out.maegz
    if [ -f "$generated_maegz" ]; then
        mv "$generated_maegz" "$conformer_dir/"
    else
        echo "Warning: $generated_maegz not found"
    fi
    echo "Finished processing $mae_file"
done

