import os
import numpy as np

from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import calculate_sasa_by_atom


def calculate_sasa(conformer_path):
    # Read the structure
    energy_list = []
    sasa_dict = {}

    # Read all structures into a list
    structures = []
    with StructureReader(conformer_path) as reader:
        for idx, st in enumerate(reader):
            structures.append(st)
            # get the energy from the strucural property
            energy = st.property["r_user_macromodel_energy"]
            energy_list.append(energy)

    energy_list = np.array(energy_list, dtype=float)
    boltzmann_factor = np.exp(-energy_list / (8.314 * 298.15))

    # Get the hydrogen atoms index in schrodinger
    structure = structures[0]
    total_atoms = structure.atom_total

    # Calculate the SASA for each atom across all conformers, loop with schrodinger index
    for atom in range(1, total_atoms + 1):
        atom_sasa_list = []
        for st in structures:
            sasa = calculate_sasa_by_atom(
                st=st, atoms=[atom], probe_radius=2, resolution=0.1
            )
            atom_sasa_list.append(sasa[atom - 1])

        # Calculate the Boltzmann average SASA for the current atom
        atom_sasa_list = np.array(atom_sasa_list, dtype=float)
        atom_sasa = np.sum(atom_sasa_list * boltzmann_factor) / np.sum(boltzmann_factor)
        sasa_dict[atom] = atom_sasa

    # Print the Boltzmann average SASA for each atom
    for atom_idx, sasa in sasa_dict.items():
        print(f"Atom {atom_idx}: Boltzmann average SASA = {sasa:.2f}")

    # Calculate the SASA for each atom for the lowest-energy conformer
    # Find the lowest-energy conformer
    lowest_energy_index = np.argmin(energy_list)
    lowest_energy_conformer = structures[lowest_energy_index]

    # Add the SASA property to the structure
    for atom in lowest_energy_conformer.atom:
        atom.property["r_user_sasa"] = sasa_dict[atom.index]

    # Write the structure with the SASA property
    lowest_energy_conformer.write(conformer_path.split("-out.sdf")[0] + "_sasa.mae")


def main():
    sdf_conformer_energy_dir = "../data/energy_sdf"

    for file in os.listdir(sdf_conformer_energy_dir):
        if file.endswith(".sdf"):
            sdf_path = os.path.join(sdf_conformer_energy_dir, file)
            calculate_sasa(sdf_path)


if __name__ == "__main__":
    main()
