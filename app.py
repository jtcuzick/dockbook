import os
import json
import streamlit as st

st.set_page_config(page_title="DockBook", layout="centered")

NOTEBOOK_FILE = "notebooks.json"

# Load and Save Functions
def load_notebooks():
    if os.path.exists(NOTEBOOK_FILE):
        with open(NOTEBOOK_FILE, "r") as f:
            return json.load(f)
    return []

def save_notebooks(notebooks):
    with open(NOTEBOOK_FILE, "w") as f:
        json.dump(notebooks, f)

# Session State Initialization
if "has_visited" not in st.session_state:
    st.session_state.has_visited = False
if "notebooks" not in st.session_state:
    st.session_state.notebooks = load_notebooks()
if "current_notebook" not in st.session_state:
    st.session_state.current_notebook = None
if "creating_notebook" not in st.session_state:
    st.session_state.creating_notebook = False

# Global Button Style
st.markdown("""
    <style>
    div.stButton > button {
        background-color: #4a90e2;
        color: white;
        padding: 0.5em 1em;
        font-size: 1em;
        font-weight: 500;
        border: none;
        border-radius: 6px;
        box-shadow: 0 2px 6px rgba(0, 0, 0, 0.1);
        transition: background-color 0.2s ease;
    }
    div.stButton > button:hover {
        background-color: #3a78c2;
    }
    .dockbook-hr {
        border: none;
        height: 1px;
        background: #000;
        margin: 0.4em 0;
    }
    </style>
""", unsafe_allow_html=True)

# Welcome Page
if not st.session_state.has_visited:
    st.markdown(
        """
        <div style='text-align: center; margin-top: 10em; font-family: "Segoe UI", sans-serif;'>
            <h1 style='font-size: 3em;'>Welcome to <span style="color:#4a90e2;">DockBook</span></h1>
            <p style='font-size: 1.3em; color: #000000;'>Your interactive notebook for molecular docking.</p>
        </div>
        """,
        unsafe_allow_html=True
    )

    col1, col2, col3 = st.columns([1, 0.7, 1])
    with col2:
        if st.button("🚀 GET TO WORK", key="welcome_start"):
            st.session_state.has_visited = True
            st.rerun()

    st.stop()

# Main Page
st.title("DockBook Files")

# Create new notebook
if st.button("➕  New DockBook", key="main_new_btn"):
    st.session_state.creating_notebook = True

if st.session_state.creating_notebook:
    notebook_input = st.text_input("Enter name", key="main_input")
    if st.button("Create", key="main_create_btn"):
        name = notebook_input.strip()
        if name and name not in st.session_state.notebooks:
            st.session_state.notebooks.append(name)
            save_notebooks(st.session_state.notebooks)
        st.session_state.current_notebook = name
        st.session_state.creating_notebook = False
        st.switch_page("pages/DockingNotebook.py")

# Notebook List Display
for i, nb in enumerate(st.session_state.notebooks):
    st.markdown("<hr class='dockbook-hr' />", unsafe_allow_html=True)
    col_open, col_name, col_del = st.columns([0.16, 0.60, 0.11])
    
    with col_open:
        if st.button("Open", key=f"main_open_{nb}"):
            st.session_state.current_notebook = nb
            st.switch_page("pages/DockingNotebook.py")
    
    with col_name:
        st.markdown(f"""
            <div style="
                font-family: 'Source Sans Pro', sans-serif;
                font-size: 1.5em;
                font-weight: 600;
                color: #262730;
                padding-left: 0.25em;
            ">{nb}</div>
        """, unsafe_allow_html=True)
    
    with col_del:
        if st.button("Delete", key=f"main_delete_{nb}"):
            st.session_state.notebooks.remove(nb)
            save_notebooks(st.session_state.notebooks)
            st.rerun()
    
    st.markdown("<hr class='dockbook-hr' />", unsafe_allow_html=True)







