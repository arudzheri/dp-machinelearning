import streamlit as st
from chembl_webresource_client.new_client import new_client as ch
from streamlit_ketcher import st_ketcher

# Display app info in an expander at the top
with st.expander("ℹ️ About this App"):
    st.write("""
        This application allows you to search for chemical molecules by name or draw them manually. It finds molecules similar to your input based on a similarity threshold you set.
        - **Choose input type**: You can either search by name or draw a molecule structure.
        - **Similarity Threshold**: Adjust the threshold to control how close the results are to your input.
    """)

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

def display_similar_molecule(molecule):
    """Displays the details of a similar molecule."""
    chembl_id = molecule.get('molecule_chembl_id')
    similarity_score = molecule.get('similarity')
    pref_name = molecule.get('pref_name') or "None"
    molfile = molecule['molecule_structures'].get('molfile')

    st.markdown(f"### Molecule: {pref_name} (Similarity: {similarity_score}%)")
    st.write(f"**ChEMBL ID:** [{chembl_id}](https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id})")
    st.write(f"**Preferred Name:** {pref_name}")
    st.write(f"**Similarity:** {similarity_score}%")
    st.write(f"**Molecular Weight:** {molecule.get('molecular_weight', 'None')}")

    if molfile:
        st.image(st_ketcher(molfile=molfile))
    st.markdown("---")

# Streamlit app layout
st.title("🧪 Chemical Molecule Similarity Search")

# Sidebar and Input Options
with st.sidebar:
    st.subheader("Molecule Input")
    option = st.radio("Choose input type", ("Draw Molecule", "Search by Name"))

    smiles = None  # Ensure smiles is defined
    
    if option == "Search by Name":
        molecule_name = st.text_input("Enter molecule name:")
        if st.button("Fetch Molecule"):
            if molecule_name.strip():
                with st.spinner("Searching for molecule..."):
                    molfile, chembl_id = name_to_molecule(molecule_name)
                    if molfile:
                        smiles = st_ketcher(molfile=molfile)
                        st.write(f"**ChEMBL ID:** {chembl_id}")
                    else:
                        st.warning("Molecule not found. Please check the name and try again.")
            else:
                st.warning("Please enter a valid molecule name.")
    else:
        smiles = st_ketcher()

# Similarity Search Settings
st.subheader("🔍 Similarity Search Settings")
st.write("Adjust the similarity threshold below to control the strictness of matching. A higher threshold will find molecules more similar to your input.")
similarity = st.slider("Similarity threshold:", min_value=60, max_value=100, value=80)

# Displaying Raw SMILES Data
if smiles:
    st.subheader("🧬 Raw SMILES Data")
    with st.expander("Click to view SMILES data for the selected molecule"):
        st.markdown(f"```{smiles}```")

    # Similarity Search Results
    st.subheader("🔗 Similar Molecules")
    st.write("Below are the molecules that are similar to your input based on the chosen similarity threshold. Each result includes the molecule’s ChEMBL ID, name, similarity score, and structure.")
    
    similar_molecules = find_similar_molecules(smiles, similarity)
    
    if similar_molecules:
        for mol in similar_molecules:
            col1, col2 = st.columns(2)
            with col1:
                display_similar_molecule(mol)
            with col2:
                if mol['molecule_structures'].get('molfile'):
                    st.image(st_ketcher(molfile=mol['molecule_structures']['molfile']))
    else:
        st.warning("No similar molecules found.")
else:
    st.info("Please provide a molecule name or draw a molecule to proceed.")

