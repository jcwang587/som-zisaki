# Script to combine BDE and SASA data

import os
import shutil
from schrodinger.structure import StructureReader
import pandas as pd


def combine_bde_mae_sasa(bde_mae_file, sasa_file, output_file):
    structure_bde = StructureReader.read(bde_mae_file)
    structure_sasa = StructureReader.read(sasa_file)

    bde_dict = {}
    # Get the SASA property r_user_sasa
    for atom in structure_sasa.atom:
        if "r_user_sasa" in atom.property:
            bde_dict[atom.index] = atom.property["r_user_sasa"]

    # Add the SASA property r_user_sasa to the BDE structure
    for atom in structure_bde.atom:
        if atom.index in bde_dict:
            atom.property["r_user_sasa"] = bde_dict[atom.index]

    # Export the combined structure
    structure_bde.write(output_file)


def combine_bde_csv_sasa(bde_csv_file, sasa_csv_file, output_file):
    df_bde = pd.read_csv(bde_csv_file)
    structure_sasa = StructureReader.read(sasa_csv_file)

    bde_dict = {}
    # Add the BDE property to the BDE structure
    for index, row in df_bde.iterrows():
        # get the atom index from the column "Atom" remove the "C"
        atom_index = int(row["Atom"].replace("C", ""))
        bde_dict[atom_index] = row["BDE (kcal/mol)"]

    # Add the BDE property to the SASA structure
    for atom in structure_sasa.atom:
        if atom.index in bde_dict:
            atom.property["r_user_CH-BDE"] = bde_dict[atom.index]

    # Export the combined structure
    structure_sasa.write(output_file)




output_dir = "./output"
shutil.rmtree(output_dir, ignore_errors=True)
os.makedirs(output_dir, exist_ok=True)

bde_csv_dir = "./bde_csv"
bde_csv_files = os.listdir(bde_csv_dir)
bde_csv_files.sort()

sasa_csv_dir = "./sasa_mae"
sasa_csv_files = os.listdir(sasa_csv_dir)
sasa_csv_files.sort()

for bde_csv_file, sasa_csv_file in zip(bde_csv_files, sasa_csv_files):
    bde_csv_path = f"{bde_csv_dir}/{bde_csv_file}"
    sasa_csv_path = f"{sasa_csv_dir}/{sasa_csv_file}"

    bde_file = os.path.basename(bde_csv_path)
    sasa_file = os.path.basename(sasa_csv_path)

    structure = StructureReader.read(sasa_csv_path)
    output_file = f"{output_dir}/{structure.title}_bde_sasa.mae"
    combine_bde_csv_sasa(bde_csv_path, sasa_csv_path, output_file)
    print(f"Combined {bde_file} and {sasa_file} into {output_file}")    