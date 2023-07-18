import streamlit as st
import pandas as pd
from io import StringIO
from HD import parse_seq, filter_df, read_file, validate
from Bio import SeqIO
import os
import tempfile


# about app
st.title('HairpinDetector')
st.subheader('Poorly writen app to detect hairpins and analyse them.')
# upload file 
files = st.file_uploader('Upload your fasta file here', accept_multiple_files=True, type = ['fasta','fa'])
# parameters of hairpin
loop_length = st.number_input('Length of a loop',0,10)
stem_length = st.number_input('Length of a stem',0,10)
threshold = st.number_input('GC count in a stem',0,10)

for file in files:

    # creating temporary file, need to be overwriten.
    with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
        tmp_filepath = tmp_file.name
        tmp_file.write(file.read())
    # Read the temporary file using SeqIO
    sequences = SeqIO.parse(tmp_filepath, "fasta")

    # Iterate over the sequences and display their information
    for seq in sequences:
        st.write(f"ID: {seq.id}")
    if stem_length is not None and loop_length is not None and threshold is not None:
        hairpins_df = parse_seq(seq.seq, stem_length,loop_length)
        filtered_df = filter_df(hairpins_df,threshold)
        st.write(filtered_df)


    # Remove the temporary file
    os.remove(tmp_filepath)
    # bytes_data = file.read()
    # data = str(bytes_data).replace("\\n", "").replace("b'",'\>')
    # st.write("filename:", file.name)
    # st.write(data)
    # st.write(read_file(data))





st.write()