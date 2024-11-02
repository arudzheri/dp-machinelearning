# Install the required libraries (uncomment if needed)
# !pip install streamlit rdkit-pypi

import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

st.title('🤖 Machine Learning - Chemical Molecule App')

st.info('This is an app that builds a machine learning model for chemical molecules!')

def main():
    st.title("Molecule Viewer")

    # Input the SMILES string here (example placeholder)
    smiles = st.text_input("Enter a SMILES string:", "COc1ccccc1C(=O)Nc2ccc(cc2)C(=O)c3ccccc3O")
    
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
