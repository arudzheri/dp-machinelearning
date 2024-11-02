import streamlit as st
from streamlit_ketcher import st_ketcher

smiles = st_ketcher()

st.title('🤖 Machine Learning - Chemical Molecule App')

st.info('This is an app that builds a machine learning model for chemical molecules!')

molecule = st.text_input("Molecule", DEFAULT_MOL)
smile_code = st_ketcher(molecule)
st.markdown(f"Smile code: ``{smile_code}``")

with editor_column:
    smiles = st_ketcher(st.session_state.molfile)
    similarity = st.slider("Similarity threshold:", min_value=60, max_value=100)
    with st.expander("Raw data"):
        st.markdown(f"```{smiles}```")
    with results_column:
        similar_molecules = utils.find_similar_molecules(smiles, similarity)
        if not similar_molecules:
            st.warning("No results found")
        else:
            table = utils.render_similarity_table(similar_molecules)
            similar_smiles = utils.get_similar_smiles(similar_molecules)
            st.markdown(f'<div style="overflow:scroll; height:600px; padding-left: 80px;">{table}</div>',
                        unsafe_allow_html=True)

def find_similar_molecules(smiles: str, threshold: int):
    columns = ['molecule_chembl_id', 'similarity', 'pref_name', 'molecule_structures']
    try:
        return ch.similarity.filter(smiles=smiles, similarity=threshold).only(columns)
    except Exception as _:
        return None
