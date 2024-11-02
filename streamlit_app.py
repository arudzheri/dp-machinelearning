import streamlit as st
from chembl_webresource_client.new_client import new_client as ch
from streamlit_ketcher import st_ketcher

# Utility functions
def name_to_molecule(name: str):
    """Fetches the molfile and ChEMBL ID for a given molecule name."""
    columns = ['molecule_chembl_id', 'molecule_structures']
    result = ch.molecule.filter(molecule_synonyms__molecule_synonym__iexact=name).only(columns)
    if result:
        best_match = result[0]
        return best_match["molecule_structures"]["molfile"], best_match["molecule_chembl_id"]
    return None, None

def find_similar_molecules(smiles: str, threshold: int):
    """Fetches molecules similar to the given SMILES string."""
    columns = ['molecule_chembl_id', 'similarity', 'pref_name', 'molecule_structures']
    try:
        return ch.similarity.filter(smiles=smiles, similarity=threshold).only(columns)
    except Exception as e:
        st.error(f"Error in similarity search: {e}")
        return None

# Streamlit app layout
st.title("Chemical Molecule Similarity Search")

# Input options
with st.sidebar:
    st.subheader("Molecule Input")
    option = st.radio("Choose input type", ("Draw Molecule", "Search by Name"))
    
    if option == "Search by Name":
        molecule_name = st.text_input("Enter molecule name:")
        if st.button("Fetch Molecule"):
            molfile, chembl_id = name_to_molecule(molecule_name)
            if molfile:
                smiles = st_ketcher(molfile=molfile)
                st.write(f"ChEMBL ID: {chembl_id}")
            else:
                st.warning("Molecule not found.")
    else:
        smiles = st_ketcher()

# Similarity search settings
similarity = st.slider("Similarity threshold:", min_value=60, max_value=100, value=80)

# Similarity search results
if smiles:
    with st.expander("Raw SMILES Data"):
        st.markdown(f"```{smiles}```")
    
    similar_molecules = find_similar_molecules(smiles, similarity)
    
    if similar_molecules is None:
        st.warning("No similar molecules found.")
    else:
        st.subheader("Similar Molecules")
        for mol in similar_molecules:
            chembl_id = mol.get('molecule_chembl_id')
            similarity_score = mol.get('similarity')
            pref_name = mol.get('pref_name') or "None"
            molfile = mol['molecule_structures'].get('molfile')
            
            st.write(f"**ChEMBL ID:** [{chembl_id}](https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id})")
            st.write(f"**Preferred Name:** {pref_name}")
            st.write(f"**Similarity:** {similarity_score}%")
            
            if molfile:
                st.image(st_ketcher(molfile=molfile))
            st.markdown("---")

