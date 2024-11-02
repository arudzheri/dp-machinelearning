import streamlit as st
import pandas as pd

st.title('🤖 Machine Learning- Chemical Molecule App')

st.info('This is app that builds a machine learning model for chemical molecules!')

molecule = st.text_input("Molecule", DEFAULT_MOL)
smile_code = st_ketcher(molecule)
st.markdown(f"Smile code: ``{smile_code}``")
