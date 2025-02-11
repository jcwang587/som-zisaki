from pathlib import Path
import os, warnings

warnings.filterwarnings("ignore")


from schrodinger.structure import StructureReader
from schrodinger.structutils import build

from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image
from reportlab.platypus import Spacer, Table, TableStyle
from reportlab.lib import colors

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from PIL import Image as ImagePIL
from PIL import ImageChops

from rdkit.Chem import rdDepictor
rdDepictor.SetPreferCoordGen(True)

from utils import draw_molecule


def create_image(st, out_dir=Path().resolve(), name="labeled"):
    """
    Create image function compatible with the Schrodinger python
    """

    draw_molecule(st, name)
    out_png = out_dir / f"{name}.png"

    return out_png


def check_missing_sites(structure: StructureReader):
    """Check if the structure has missing sites"""
    for atom in structure.atom:
        if atom.property.get("r_user_CH-BDE") is None:
            print(f"Missing site: {atom.element} {atom.index}")


def autocrop(image: str, bgcolor: str = "white"):
    im = ImagePIL.open(image)
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = ImagePIL.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        image_crop = im.crop(bbox)
        image_crop.save(image)
        width, length = image_crop.size
        return width, length


def create_risk_scale_png(
    medium: float, high: float, filename: Path = Path().resolve() / "risk_scale.png"
):
    """Create risk figure based on medium and high cut-offs"""

    fig = plt.figure()
    ax = fig.add_subplot(111)
    rect1 = Rectangle((50, 50), 50, 40, color="yellow")

    if high > medium:
        rect = Rectangle((0, 50), 50, 40, color="green")
        rect2 = Rectangle((100, 50), 50, 40, color="red")
        ax.add_patch(rect)
        ax.add_patch(rect1)
        ax.add_patch(rect2)
        ax.annotate(str(medium), xy=(44, 94), xytext=(44, 94), fontsize=14)
        ax.annotate(str(high), xy=(94, 94), xytext=(94, 94), fontsize=14)
        ax.annotate("Low", xy=(18, 38), xytext=(18, 38), fontsize=14)
        ax.annotate("Medium", xy=(58, 38), xytext=(58, 38), fontsize=14)
        ax.annotate("High", xy=(115, 38), xytext=(115, 38), fontsize=14)

    else:
        rect = Rectangle((100, 50), 50, 40, color="green")
        rect2 = Rectangle((0, 50), 50, 40, color="red")
        ax.add_patch(rect2)
        ax.add_patch(rect1)
        ax.add_patch(rect)
        ax.annotate(str(high), xy=(44, 94), xytext=(44, 94), fontsize=14)
        ax.annotate(str(medium), xy=(94, 94), xytext=(94, 94), fontsize=14)
        ax.annotate("High", xy=(18, 38), xytext=(18, 38), fontsize=14)
        ax.annotate("Medium", xy=(58, 38), xytext=(58, 38), fontsize=14)
        ax.annotate("Low", xy=(115, 38), xytext=(115, 38), fontsize=14)

    ax.annotate("kcal/mol", xy=(-40, 94), xytext=(-40, 94), fontsize=14)
    plt.xlim(-40, 175)
    plt.ylim(-10, 175)
    plt.axis("off")

    plt.savefig(str(filename), dpi=300)

    autocrop(str(filename))


class CoxidesAnalysis:

    def __init__(
        self,
        dir_report: Path,
        path_structure: Path,
        title: str,
        medium_coff: float,
        high_coff: float,
    ):
        self.dir_report = dir_report
        self.path_structure = path_structure
        self.title = title
        self.medium_coff = medium_coff
        self.high_coff = high_coff

        structure = StructureReader.read(self.path_structure)
        build.add_hydrogens(structure)
        self.structure = structure
        
        self.dft = "B3LYP"
        self.basis = "6-31G(d,p)"
        self.risk_scale = self.dir_report / "C-oxidation_risk_scale.png"
        self.img = create_image(self.structure, name=f"{self.dir_report}/{self.title}")
        self.data = self.GenerateDataList()
        self.missing_sites = self.MissingSites()


    def MissingSites(self):
        """Check if the structure has missing sites"""

        # Get all the carbon atoms
        carbons = [atom for atom in self.structure.atom if atom.element == "C"]

        # Get the carbon atoms with no hydrogen atoms
        no_hydrogen_carbons = [
            atom for atom in carbons if any(neighbor.element == "H" for neighbor in atom.bonded_atoms)
        ]

        missing_sites = []
        for atom in no_hydrogen_carbons:
            if atom.property.get("r_user_CH-BDE") is None:
                missing_sites.append(f"{atom.element}{atom.index}, ")
        
        # Change the last comma to a period
        if len(missing_sites) > 0:
            missing_sites[-1] = missing_sites[-1].replace(",", ".")

        return missing_sites



    def GenerateDataList(self):
        """
        Take information from self.structure to create a Data list for report
        generation
        """

        print(f"Creating Reports for {self.title}")
        data = [["Atom", "BDE (kcal/mol)", "Propensity"]]

        for atom in self.structure.atom:
            try:
                bde = float(atom.property["r_user_CH-BDE"])
                if bde < self.high_coff:
                    propensity = "High"
                if bde <= self.medium_coff and bde >= self.high_coff:
                    propensity = "Moderate"
                if bde > self.medium_coff:
                    propensity = "Low"
                row = [
                    str(atom.element) + str(atom.index),
                    str(atom.property["r_user_CH-BDE"]),
                    propensity,
                ]
                data.append(row)
            except KeyError:
                pass

        # Check if SASA property is available
        if "r_user_sasa" in self.structure.atom[1].property:
            data = [["Atom", "BDE (kcal/mol)", "Propensity", "SASA (Å²)"]]
            for atom in self.structure.atom:
                try:
                    bde = float(atom.property["r_user_CH-BDE"])
                    if bde < self.high_coff:
                        propensity = "High"
                    if bde <= self.medium_coff and bde >= self.high_coff:
                        propensity = "Moderate"
                    if bde > self.medium_coff:
                        propensity = "Low"

                    hydrogen_atoms = []
                    # Find the hydrogen atoms attached to the carbon atoms
                    for bond in atom.bond:
                        if bond.atom1.element == "H" or bond.atom2.element == "H":
                            hydrogen_atoms.append(
                                bond.atom1 if bond.atom1.element == "H" else bond.atom2
                            )

                    hydrogen_sasa = 0
                    # Get the SASA of the hydrogen atoms
                    for hydrogen_atom in hydrogen_atoms:
                        if "r_user_sasa" in hydrogen_atom.property:
                            hydrogen_sasa += hydrogen_atom.property["r_user_sasa"]

                    hydrogen_sasa = (
                        hydrogen_sasa / len(hydrogen_atoms)
                        if len(hydrogen_atoms) > 0
                        else 0
                    )

                    row = [
                        str(atom.element) + str(atom.index),
                        str(atom.property["r_user_CH-BDE"]),
                        propensity,
                        f"{hydrogen_sasa:.2f}",
                    ]
                    data.append(row)
                except KeyError:
                    pass

        return data

    def CreatePDFReport(self):
        """Create pdf report from a data dictionary and structural image"""
        # Define style
        style = TableStyle(
            [
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("FONTNAME", (0, 0), (-1, 0), "Courier-Bold"),
                ("FONTSIZE", (0, 0), (-1, 0), 10),
                ("BACKGROUND", (0, 1), (-1, -1), colors.beige),
                ("BOX", (0, 0), (-1, -1), 2, colors.black),
                ("GRID", (0, 1), (-1, -1), 2, colors.black),
            ]
        )

        # Get other styles
        styles = getSampleStyleSheet()
        styleNormal = styles["Normal"]
        styleHeading = styles["Heading1"]
        styleHeading2 = styles["Heading2"]
        styleHeading3 = styles["Heading3"]
        styleHeading.alignment = 1

        # Process data
        pil_image = ImagePIL.open(self.img)
        image_size = pil_image.size

        pil_st = Image(
            str(self.img),
            width=image_size[0] / 200 * inch,
            height=image_size[1] / 200 * inch,
        )
        pil_risk = Image(str(self.risk_scale), width=3.0 * inch, height=0.8 * inch)

        # Create PDF
        story = []

        # Append HEAD of document
        story.append(
            Paragraph("C-oxidation BDE Energy Report for: " + self.title, styleHeading)
        )
        story.append(Spacer(inch, 0.25 * inch))
        story.append(
            Paragraph(
                "This report covers the results for bond dissociation enthalpies (BDE) and "
                f"solvent accessible surface area (SASA) calculations performed for {self.title}. "
                "Oxidation propensity is established using C-H BDE. The lower the C-H BDE values the "
                "higher the propensity for C-oxidation. Details for the density functional theory (DFT) "
                "calculations and overall workflow are explained at the end of this document.",
                styleNormal,
            )
        )

        # Append molecule picture
        story.append(Paragraph("BDE and SASA", styleHeading2))
        story.append(Spacer(inch, 0.15 * inch))
        story.append(pil_st)
        story.append(Spacer(inch, 0.25 * inch))

        # Append BDE table
        main_table = Table(self.data, colWidths=85, rowHeights=18)
        main_table.setStyle(style)
        story.append(main_table)

        # Append missing sites
        story.append(Paragraph("Missing Sites:", styleHeading3))
        

        if isinstance(self.missing_sites, list):
            if len(self.missing_sites) > 0:
                self.missing_sites = "\n".join(self.missing_sites)
                story.append(Paragraph(self.missing_sites, styleNormal))
            else:
                story.append(Paragraph("None", styleNormal))

        # Risk scale
        story.append(Paragraph("Risk Scale:", styleHeading3))
        story.append(pil_risk)

        # Append calculation details
        story.append(Paragraph("Calculation Details", styleHeading2))
        story.append(
            Paragraph(
                "Conformational search calculations were performed only "
                "for the base ground state molecule. The lowest energy conformer was selected "
                "to generate radicals and run optimization DFT calculations. DFT calculations "
                f"were performed using Gaussian with {self.dft} level of theory and {self.basis} "
                "basis set. The BDE protocol was adapted from: <i>Lienard, P., Gavartin, J., "
                "Boccardi, G., & Meunier, M. (2015). Predicting drug substances autoxidation. "
                "Pharmaceutical research, 32, 300-310.</i>",
                styleNormal,
            )
        )

        # Save PDF
        pdf_file = self.dir_report / f"{self.title}_CH-BDE_report.pdf"
        doc = SimpleDocTemplate(
            str(pdf_file),
            pagesize=letter,
            title=f"{self.title} BDE Report",
            author="BDE 2.0",
        )
        doc.build(story)

        # Remove the image file
        os.remove(self.img)

        return pdf_file


if __name__ == "__main__":

    C_HIGH_COFF = 88
    C_MEDIUM_COFF = 94

    mae_dir = Path("./output")
    mae_files = mae_dir.glob("*.mae")

    # Generate the risk scale image
    create_risk_scale_png(
        medium=C_MEDIUM_COFF,
        high=C_HIGH_COFF,
        filename=mae_dir / "C-oxidation_risk_scale.png",
    )

    # Generate the pdf report using the CoxidesAnalysis
    for mae_file in mae_files:
        mae_structure = StructureReader.read(mae_file)

        qmox_analysis = CoxidesAnalysis(
            dir_report=mae_dir,
            path_structure=Path(mae_file),
            title=mae_structure.title,
            medium_coff=C_MEDIUM_COFF,
            high_coff=C_HIGH_COFF,
        )
        qmox_analysis.CreatePDFReport()

    # Remove the risk scale image
    os.remove(mae_dir / "C-oxidation_risk_scale.png")
