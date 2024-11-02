import streamlit as st
from PIL import Image

st.title('🤖 Machine Learning- Chemical Molecule App')

st.info('This is app that builds a machine learning model for chemical molecules!')

# Load and display the image
image_path = "/mnt/data/Screenshot 2024-11-02 142900.png"
image = Image.open(image_path)

st.title("Chemical Structure Image")
st.image(image, caption="Chemical Structure", use_column_width=True)
