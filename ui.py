import streamlit as st
import pandas as pd
from io import StringIO
from src.HD import parse_seq, filter_df, validate
from Bio import SeqIO
import os
import tempfile

def main():
    st.set_page_config(
    page_title="HairpinDetector",
    page_icon="🧊",
    initial_sidebar_state="expanded")
# about app
    st.title('HairpinDetector')
    st.text('Poorly writen app to detect hairpins and analyse them.')
    # upload file 
    files = st.file_uploader('Upload your fasta file here', accept_multiple_files=True, type = ['fasta','fa'])
    # parameters of hairpin
    columns = st.columns([1, 1, 1])
    loop_length = columns[0].number_input('Length of a loop',0,10)
    stem_length = columns[1].number_input('Length of a stem',0,10)
    threshold = columns[2].number_input('GC count in a stem',0,10)



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
        if st.button('Search') and stem_length is not None and loop_length is not None and threshold is not None:
            # applying function to find hairpins and filter them 
            hairpins_df = parse_seq(seq.seq, stem_length,loop_length)
            filtered_df = filter_df(hairpins_df,threshold)
            # display df
            st.write(filtered_df)


            st.download_button(
                label="Download data as CSV",
                data= filtered_df.to_csv().encode('utf-8'),
                file_name=f'{seq.id}_l{loop_length}_s{stem_length}_t{threshold}.csv',
                mime='text/csv',
            )


        # Remove the temporary file
        os.remove(tmp_filepath)
    


if __name__ == '__main__':
    main()