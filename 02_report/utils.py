from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdDepictor
from rdkit.Chem import rdCoordGen
import os
import cairosvg

from schrodinger.structure import StructureReader
from schrodinger.rdkit.rdkit_adapter import to_rdkit
from PIL import Image


def draw_molecule(structure: StructureReader, name="labeled"):

    mol = to_rdkit(structure)

    # Remove hydrogens
    mol = Chem.RemoveHs(mol)

    # Compute 2D coordinates
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)

    # Remove chiral information
    Chem.RemoveStereochemistry(mol)

    # Calculate the size based on the number of heavy atoms
    num_atoms = mol.GetNumHeavyAtoms()
    size = min(1000, num_atoms * 30)

    # Draw the molecule with highlighted atoms without bonds
    mol_draw = Draw.rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(size, size)

    opts = drawer.drawOptions()

    for i in range(mol_draw.GetNumAtoms()):
        opts.atomLabels[i] = f"{mol_draw.GetAtomWithIdx(i).GetSymbol()}<sub>{i+1}</sub>"

    drawer.drawOptions().prepareMolsBeforeDrawing = False
    drawer.drawOptions().maxFontSize = 30

    drawer.DrawMolecule(mol_draw)

    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()

    # Save the SVG image
    with open(f"{name}.svg", "w") as f:
        f.write(svg)

    cairosvg.svg2png(
        url=f"{name}.svg",
        write_to=f"{name}.png",
    )

    # Open the PNG image
    img = Image.open(f"{name}.png").convert("RGBA")

    # Make white areas transparent
    datas = img.getdata()
    new_data = []
    for item in datas:
        if item[0] > 200 and item[1] > 200 and item[2] > 200:
            new_data.append((255, 255, 255, 0))
        else:
            new_data.append(item)

    img.putdata(new_data)

    img_cropped = img.crop(img.getbbox())

    # Save the resized image with transparency
    img_cropped.save(f"{name}.png")

    # Remove the SVG file
    os.remove(f"{name}.svg")
