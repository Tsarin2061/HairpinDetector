# TsarinLev@gmail.com
import streamlit as st
import pandas as pd
from src.HD import parse_seq, filter_df, validate
from Bio import SeqIO
import tempfile

def main():
    # page set up
    
    st.set_page_config(
    page_title = "HairpinDetector",
    initial_sidebar_state = "expanded",
    )
    st.title('HairpinDetector')
    st.text('Detects hairpins formed by perfect palindromes')
    

    # uploading data to analyse/ 0 - data itself, 1 - way we analyse it(text, multi/fasta).
    input_info = input_options()     
    data = input_info[0]
    way = input_info[1]

    # parameters of the hairpin
    columns = st.columns([1, 1, 1,1])
    loop_length = columns[0].number_input('Length of a loop',0,10)
    stem_length = columns[1].number_input('Length of a stem',0,10)
    threshold = columns[2].number_input('GC count in a stem',0,10)
    

    analyse(way,data,loop_length,stem_length,threshold)
    st.runtime.legacy_caching.clear_cache()



def input_options():

    input_format = st.radio(
    "Choose input format:",
    ('Text', 'Fasta files'),
    horizontal=True,
    )

    if input_format == 'Text':
        sequence = st.text_input(
            'Enter the sequence', max_chars = 400
            ).upper()
        return sequence, 'text'
    elif input_format == 'Fasta files':
        sequence = st.file_uploader(
            'Upload your fasta file', accept_multiple_files = True, type = ['fasta','fa']
            )
        return sequence, 'fasta'
    

def analyse(way,data,loop_length, stem_length, threshold):
    if way == 'fasta':
        for file in data:
        # creating temporary file, need to be overwriten.
            with tempfile.NamedTemporaryFile(delete=False) as tmp_file:
                tmp_filepath = tmp_file.name
                tmp_file.write(file.read())
            # Read the temporary file using SeqIO
            sequences = SeqIO.parse(tmp_filepath, "fasta")

            # Iterate over the sequences and display their information
            for seq in sequences:
                st.write(f"ID: {seq.id}")
            col1,col2 = st.columns(2)
            with col1:
                st.text('Search hairpins \nwith exact characteristics')
                search_button = st.button('Search')
            with col2:
                st.text("Search all hairpins \nby using values as tresholds")
                searchall_button = st.button("Search all")
            if search_button and stem_length is not None and loop_length is not None and threshold is not None:
                # applying function to find hairpins and filter them 
                with st.spinner("Searching..."):
                    hairpins_df = parse_seq(seq.seq, stem_length,loop_length)
                    filtered_df = filter_df(hairpins_df,threshold)
                # display df
                st.write(filtered_df)
            if searchall_button:
                # applying function to find hairpins and filter them 
                with st.spinner("Searching..."):
                    hairpins_df = full_search(seq.seq)
                # display df
                st.write(hairpins_df)

                # st.download_button(
                #     label="Download data as CSV",
                #     data = filtered_df.to_csv().encode('utf-8'),
                #     file_name=f'{seq.id}_l{loop_length}_s{stem_length}_t{threshold}.csv',
                #     mime='text/csv',
                # )
            

    elif way == 'text':
        seq = data
        if validate(seq) == False:
            return st.write("This string must consist exclusively of nucleotides(A, T, G, C)!")
        col1,col2 = st.columns(2)
        with col1:
            st.text('Search hairpins \nwith precise structure')
            search_button = st.button('Search')
        with col2:
            st.text("Search all hairpins \nwhich could occur")
            searchall_button = st.button("Search all")
        if search_button and stem_length is not None and loop_length is not None and threshold is not None:
            # applying function to find hairpins and filter them
            with st.spinner("Searching..."):
                hairpins_df = parse_seq(seq, stem_length,loop_length)
                filtered_df = filter_df(hairpins_df,threshold)
            # display df
            st.write(filtered_df)
        if searchall_button and stem_length is not None and loop_length is not None and threshold is not None:
            # applying function to find hairpins and filter them 
            with st.spinner("Searching..."):
                hairpins_df = full_search(seq)
            # display df
            st.write(hairpins_df)



def full_search(data,loop_len = 15, stem_len = 15):
    df = pd.DataFrame()
    for i in range(loop_len):
        for j in range(4,stem_len):
            try:
                iter_df = parse_seq(seq = data,loop_length = i,ir_length = j)
                df = pd.concat([df,iter_df],axis = 0)
            except IndexError:
                pass
    df = df.drop_duplicates(subset="Hairpin_region")
    df.reset_index(inplace=True)
    df.drop(columns="index",inplace=True)
    return df


if __name__ == '__main__':
    main()