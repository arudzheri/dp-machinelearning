import streamlit as st

# App Title
st.title("User Profile App")

# Collecting User Input
name = st.text_input("Enter your name:", "")
age = st.number_input("Enter your age:", min_value=0, max_value=120, step=1)
color = st.color_picker("Pick your favorite color:", "#00f900")
gender = st.radio("Select your gender:", ("Male", "Female", "Other"))

# Button to Submit
if st.button("Submit"):
    # Check if name is provided
    if name:
        # Display the customized message
        st.write(f"Hello, {name}!")
        st.write(f"You're {age} years old.")
        st.write(f"Your favorite color is {color}.")
        st.write(f"You identify as {gender}.")
    else:
        st.error("Please enter your name.")
