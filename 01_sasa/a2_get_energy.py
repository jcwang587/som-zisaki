import os
import numpy as np

from schrodinger.application.macromodel.apps import get_energy
from schrodinger.structure import StructureReader
from schrodinger.structure import StructureWriter


def calculate_sasa(conformer_maegz_path, energy_sdf_path):
    # Read the structure
    energy_list = []

    # Read all structures into a list
    structures = []
    with StructureReader(conformer_maegz_path) as reader:
        for idx, st in enumerate(reader):
            structures.append(st)
            energy = get_energy(st)["Total energy kJ/mol"]
            energy_list.append(energy)

    energy_list = np.array(energy_list, dtype=float)

    # Assign the energy as a property to the structure
    for idx, st in enumerate(structures):
        st.property["r_user_macromodel_energy"] = energy_list[idx]

    with StructureWriter(energy_sdf_path) as writer:
        for st in structures:
            writer.append(st)

    print(f"Exported {len(structures)} structures to {energy_sdf_path}")



if __name__ == "__main__":
    # Set the directory for maegz files
    maegz_dir = "../data/conformer_maegz"
    energy_dir = "../data/energy_sdf"

    # Get the list of maegz files
    maegz_files = [f for f in os.listdir(maegz_dir) if f.endswith(".maegz")]

    # Calculate the SASA for each maegz file
    for maegz_file in maegz_files:
        calculate_sasa(
            conformer_maegz_path=os.path.join(maegz_dir, maegz_file),
            energy_sdf_path=os.path.join(energy_dir, maegz_file.replace(".maegz", ".sdf")),
        )
