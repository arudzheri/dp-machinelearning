# Install the required libraries (uncomment these if needed)
# !pip install streamlit rdkit-pypi

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

st.title('🤖 Machine Learning- Chemical Molecule App')

st.info('This is app that builds a machine learning model for chemical molecules!')

def main():
    st.title("Molecule Viewer")

    # Input the SMILES string here (example placeholder)
    smiles = "COc1ccccc1C(=O)Nc2ccc(cc2)C(=O)c3ccccc3O"  # Replace with actual SMILES of your molecule
    
    # Generate the molecule from the SMILES string
    molecule = Chem.MolFromSmiles(smiles)
    
    # Draw the molecule
    if molecule:
        st.write("### Molecular Structure")
        image = Draw.MolToImage(molecule)
        st.image(image)
    else:
        st.write("Invalid SMILES string. Please enter a valid one.")

if __name__ == "__main__":
    main()
