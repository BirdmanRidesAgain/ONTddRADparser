import streamlit as st
import time
import subprocess

with st.sidebar:
    with st.echo():
        st.write("This code will be printed to the sidebar.")

    with st.spinner("Loading..."):
        time.sleep(5)
    st.success("Done!")

st.markdown("# :red[ON]:orange[Td]:yellow[dR]:green[AD]:blue[pa]:violet[rs]:gray[er] 🌈")
st.subheader("A GUI frontend for [ONTddRADparser](https://github.com/BirdmanRidesAgain/ONTddRADparser)")


st.set_page_config(
    page_title="ONTddRADparser",
    page_icon=":rainbow:",
)

st.markdown('## Set parameters')

# Initialize touched fields tracking
if "touched_fields" not in st.session_state:
    st.session_state.touched_fields = set()

# Define input parameters as a list of dicts for easy expansion
input_params = [
    {
        "label": "Path to fasta",
        "key": "fasta",
        "default": ""
    },
    {
        "label": "Path to demux file",
        "key": "demux",
        "default": ""
    },
    {
        "label": "Minimum alignment percent for a 'fuzzy' match",
        "key": "fuzzy_aln_percent",
        "default": ".9"
    },
    {
        "label": "Minimum alignment percent for an 'exact' match",
        "key": "exact_aln_percent",
        "default": "1"
    },
    {
        "label": "Num\. bp around boundaries of long element where short element is valid",
        "key": "buffer",
        "default": "9"
    }
]

# Create callback functions to mark fields as touched
def make_touch_callback(field_key):
    def callback():
        st.session_state.touched_fields.add(field_key)
    return callback

for param in input_params:
    st.text_input(
        param["label"], 
        value=param["default"], 
        key=param["key"],
        on_change=make_touch_callback(param["key"])
    )


# Validation functions
def check_all_params_valid():
    """Checks if all parameters are valid (for button state)"""
    fasta = st.session_state.get("fasta", "")
    demux = st.session_state.get("demux", "")
    
    # Check fasta and demux are non-empty strings
    if not fasta or not isinstance(fasta, str) or not fasta.strip():
        return False
    if not demux or not isinstance(demux, str) or not demux.strip():
        return False
    
    # Check fuzzy_aln_percent is a float between 0 and 1
    try:
        fuzzy_aln_percent = float(st.session_state.get("fuzzy_aln_percent", ""))
        if not (0 <= fuzzy_aln_percent <= 1):
            return False
    except (ValueError, TypeError):
        return False
    
    # Check exact_aln_percent is a float between 0 and 1
    try:
        exact_aln_percent = float(st.session_state.get("exact_aln_percent", ""))
        if not (0 <= exact_aln_percent <= 1):
            return False
    except (ValueError, TypeError):
        return False
    
    # Check buffer is an integer
    try:
        int(st.session_state.get("buffer", ""))
    except (ValueError, TypeError):
        return False
    
    return True

def get_error_message_for_touched_fields():
    """Returns error message only for fields that have been touched"""
    touched = st.session_state.touched_fields
    
    # Check fasta
    if "fasta" in touched:
        fasta = st.session_state.get("fasta", "")
        if not fasta or not isinstance(fasta, str) or not fasta.strip():
            return "Please provide a valid path to fasta file"
    
    # Check demux
    if "demux" in touched:
        demux = st.session_state.get("demux", "")
        if not demux or not isinstance(demux, str) or not demux.strip():
            return "Please provide a valid path to demux file"
    
    # Check fuzzy_aln_percent
    if "fuzzy_aln_percent" in touched:
        try:
            fuzzy_aln_percent = float(st.session_state.get("fuzzy_aln_percent", ""))
            if not (0 <= fuzzy_aln_percent <= 1):
                return "fuzzy_aln_percent must be between 0 and 1"
        except (ValueError, TypeError):
            return "fuzzy_aln_percent must be a valid number between 0 and 1"
    
    # Check exact_aln_percent
    if "exact_aln_percent" in touched:
        try:
            exact_aln_percent = float(st.session_state.get("exact_aln_percent", ""))
            if not (0 <= exact_aln_percent <= 1):
                return "exact_aln_percent must be between 0 and 1"
        except (ValueError, TypeError):
            return "exact_aln_percent must be a valid number between 0 and 1"
    
    # Check buffer
    if "buffer" in touched:
        try:
            int(st.session_state.get("buffer", ""))
        except (ValueError, TypeError):
            return "buffer must be a valid integer"
    
    return ""

# Validate parameters
is_valid = check_all_params_valid()
error_message = get_error_message_for_touched_fields()

def run_ONTddRADparser():
    '''Builds and runs the command to run ONTddRADparser locally from user inputs.'''
    fasta = st.session_state.get('fasta', '')
    demux = st.session_state.get('demux', '')
    buffer = st.session_state.get('buffer', '')
    fuzzy_aln_percent = st.session_state.get('fuzzy_aln_percent', '')
    exact_aln_percent = st.session_state.get('exact_aln_percent', '')
    
    command = [
        'python3', './ONTddRADparser.py',
        '-f', fasta,
        '-d', demux,
        '-b', buffer,
        '-fa', fuzzy_aln_percent,
        '-ea', exact_aln_percent
    ]
    
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        st.success("ONTddRADparser completed successfully!", icon="✅")
        if result.stdout:
            st.text(result.stdout)
    except subprocess.CalledProcessError as e:
        st.error(f"Error running ONTddRADparser: {e}")
        if e.stderr:
            st.text(e.stderr)
    except Exception as e:
        st.error(f"Unexpected error: {e}")

# Add button at the bottom
st.markdown("---")  # Add a separator
if is_valid:
    if st.button("Run Analysis", type="primary", use_container_width=True):
        run_ONTddRADparser()
else:
    st.button("Run Analysis", disabled=True, use_container_width=True)
    if error_message:
        st.error(error_message)
