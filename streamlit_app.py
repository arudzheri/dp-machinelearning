import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from PIL import Image
import io

# Title and description
st.title("Chemical Molecule Explorer")
st.write("Enter a SMILES string to visualize a molecule and see some basic properties.")

# Sidebar for input
st.sidebar.header("Molecule Input")
smiles = st.sidebar.text_input("Enter SMILES notation:", "CCO")  # Example SMILES for ethanol

# Function to draw molecule
def draw_molecule(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=(300, 300))
        return img
    else:
        return None

# Function to get molecular properties
def get_molecule_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        properties = {
            "Molecular Weight": Descriptors.MolWt(mol),
            "LogP (octanol/water partition coefficient)": Descriptors.MolLogP(mol),
            "Number of Rotatable Bonds": Descriptors.NumRotatableBonds(mol),
            "Number of Hydrogen Bond Donors": Descriptors.NumHDonors(mol),
            "Number of Hydrogen Bond Acceptors": Descriptors.NumHAcceptors(mol),
        }
        return properties
    else:
        return None

# Draw the molecule and show properties
if smiles:
    # Display molecule structure
    st.subheader("Molecule Structure")
    molecule_image = draw_molecule(smiles)
    if molecule_image:
        st.image(molecule_image, caption="Molecule Structure", use_column_width=True)
    else:
        st.error("Invalid SMILES string. Please enter a valid notation.")

    # Display molecular properties
    st.subheader("Molecular Properties")
    properties = get_molecule_properties(smiles)
    if properties:
        for prop, value in properties.items():
            st.write(f"**{prop}:** {value}")
    else:
        st.error("Could not compute properties. Check the SMILES input.")
