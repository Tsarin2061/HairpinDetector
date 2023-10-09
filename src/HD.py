#!/usr/bin/env python3
# Author: Lev Tsarin

import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os
import pandas as pd
from concurrent.futures import ThreadPoolExecutor
import time
from multiprocessing import Pool



def parse_arg():
    """fucntion that parse aguments

    Returns:
        list : list of arguments that will be rewriten in variables in main()
    """
    parser = argparse.ArgumentParser(
        prog="HairpinFinder",
        usage="%(prog)s[options]",
        description="Hi there! This is a program for detection of hairpins in nucleotide sequence.",
        epilog="Contact me here: TsarinLev@gmail.com",
    )

    parser.add_argument("--path", "-p", metavar="", type=str, help="path to fasta file")

    parser.add_argument(
        "--stem_length", "-s", metavar="", type=int, help="lenght of inverted repeats"
    )

    parser.add_argument(
        "--loop_length",
        "-l",
        metavar="",
        type=int,
        help="length of loop seq, amount of nucleotides between inverted repeats",
    )

    parser.add_argument(
        "--threshold_GC",
        "-t",
        metavar="",
        type=int,
        help="setting up MINIMUM GCs in inverted repeats (pcs)",
    )
    parser.add_argument(
        "--output_name",
        "-o",
        metavar="",
        type=str,
        help="path&name of output file (without extension)",
    )
    parser.add_argument(
        "--search_all",
        "-all",
        metavar="",
        type=bool,
        help="option to parse for all hairpins that are within thresholds",
    )

    args = parser.parse_args()
    return (
        args.path,
        args.stem_length,
        args.loop_length,
        args.threshold_GC,
        args.output_name,
        args.search_all
    )


def main():
    # Parse command-line arguments
    fasta_file, stem_length, loop_length, threshold_GC, output_names, search_all = parse_arg()
    
    if search_all:
        # Read sequences from the input file
        print("reading file")
        seq = read_file(fasta_file)
        print(f"finished seq{type(seq)}")
        
        # Process each sequence separately and save results to individual CSV files
        start = time.time()
        for id, subseq in seq.items():
            print(f"Working on {id}...")
            subseq = str(subseq).replace('N','')

            final_df = full_search(Seq(subseq), stem_length, loop_length).reset_index().drop(columns='index')
            # Save the results to a CSV file with a unique name based on the sequence ID
            final_df.to_csv(f'{output_names}_{id}.csv')
            print(f"{id} done!\nNEXT!")     
        end = time.time()
        print(f"{(end - start)/60} min")
        # Merge all individual CSV files into one
        # directory_path = os.path.dirname(output_names)
        # contents = os.listdir(directory_path)
        
        # dfs = []
        # for file in contents:
        #     if os.path.isfile(f"{directory_path}/{file}"):
        #         df_to_merge = pd.read_csv(f"{directory_path}/{file}")
        #         dfs.append(df_to_merge)

        # merged_df = pd.concat(dfs, ignore_index=True)
        # # Save the merged results to a single CSV file
        # merged_df.to_csv(f'{directory_path}/all_merged.csv')
        print(f"\nResults are stored in {output_names}.csv\nPlease use different output name to avoid overwriting data!")
    else:
        # Read sequences from the input file
        seq = read_file(fasta_file)
        
        # Process sequences and filter the results based on threshold
        if isinstance(seq, dict):
            dfs = []
            
            for id, subseq in seq.items():
                final_df = filter_df(parse_seq(subseq, stem_length, loop_length), threshold_GC).reset_index().drop(columns='index')
                final_df = final_df.assign(ID=id)
                final_df = final_df[["ID"] + [col for col in final_df.columns if col != 'ID']]
                dfs.append(final_df)
            merged_df = pd.concat(dfs, ignore_index=True)
        else:
            final_df = filter_df(parse_seq(seq, stem_length, loop_length), threshold_GC).reset_index().drop(columns='index')

        if len(final_df) == 0:
            sys.exit("Hairpins were not found!")
        else:
            # Save the results to a CSV file
            final_df.to_csv(f'{output_names}.csv')
            print(f"\nResults are stored in {output_names}.csv\nPlease use different output name to avoid overwriting data!")




def read_file(path):
    dict = {}
    records = list(SeqIO.parse(path,"fasta"))
    if len(records) == 1:
        return records[len(records)-1].seq
    else:
        for i in records:
            dict[i.id] = i.seq
        return dict


def validate(s):
    """Validation of input

    Args:
        s (string): nucleotides sequence
    Returns:
        Boolean :   False - if there are some artifacts,
                    True - if sequence consist  of only nucleotides
    """
    nucleotides = ["A", "T", "C", "G"]
    for i in s:
        if i.upper() not in nucleotides:
            return False
    else:
        return True
    
def parse_seq(seq, ir_length, loop_length):
    """Inverted repeat detection

    Args:
        seq (string): nucleotide sequence
       ir_length(integer): lenght of inverted repeat we are looking for

    Returns:
        seqs_df : a banch of information about hairpins.
    """
    seqs_df = pd.DataFrame(columns =['Coordinates','IR1', 'IR2','loop_seq' ,'Hairpin_region', 'Adjacent_region(30nt)','AR_coordinates','loop_len','stem_len'])
    print('seq object heeereee')
    start = 0
    end = ir_length
    ir2_length = ir_length
    # list with coordinates and sequences
    while end != len(seq):
        act_loop_len = (ir2_length + 1) - end
        if seq[start:end] != Seq(seq[ir2_length+1 :ir2_length + ir_length + 1]).reverse_complement():
            ir2_length += 1
            if act_loop_len > loop_length:
                # if inversted repead was not found - we looking for another one
                start += 1
                end += 1
                ir2_length = ir_length + start
            else:
                pass
        else:
            # checks if loop meet requirments and correctness of inverted repeats.
            if (
                act_loop_len >= loop_length 
                and seq[start - 1] != Seq(seq[ir2_length + ir_length + 1]).reverse_complement()
            ):
                if seq[start:end+1] != Seq(seq[ir2_length: ir2_length + ir_length + 1]).reverse_complement():
                    # if we found inverted repeat - we add it to list and look for next one
                    seqs_df.loc[len(seqs_df)] =[
                        # coordinates
                            f"{start}-{ir2_length + ir_length +1}",
                        # IRs
                            str(seq[start:end]),
                            str(seq[ir2_length+1: ir2_length + ir_length + 1]),
                        # Loop
                            str(seq[end:ir2_length]),
                        # entire hairpin
                            str(seq[start : ir2_length +ir_length+ 1]),
                        # extended hairpin region
                            str(seq[start-15:   ir2_length +ir_length+ 1+15]),
                        # coordinates #2
                            f"{start-15}-{ir2_length+ir_length+1+15}",
                            f"{act_loop_len}",
                            f"{ir_length+1}"
                            
                        ]
                elif seq[start:end+1] == Seq(seq[ir2_length : ir2_length + ir_length + 1]).reverse_complement():
                    seqs_df.loc[len(seqs_df)] =[
                        # coordinates
                            f"{start}-{ir2_length + ir_length +1}",
                        # IRs
                            str(seq[start:end+1]),
                            str(seq[ir2_length: ir2_length + ir_length + 1]),
                            str(seq[end:ir2_length]),
                        # entire hairpin
                            str(seq[start : ir2_length +ir_length+ 1]),
                        # extended hairpin region
                            str(seq[start-15:   ir2_length +ir_length+ 1+15]),
                        # coordinates #2
                            f"{start-15}-{ir2_length+ir_length+1+15}",
                            f"{act_loop_len-1}",
                            f"{ir_length+1}"
                        ]

                start += 1
                end += 1
                ir2_length = ir_length + start

            else:
              ir2_length += 1
              pass
        if end == len(seq):
            return seqs_df
        


def process_chunk(params):
    data, i_range, j_range, k_range, threshold = params
    df_chunk = pd.DataFrame()
    for i in i_range:
        for j in j_range:
            for k in k_range:
                try:
                    print('PARSING')
                    iter_df = parse_seq(seq=data, ir_length=j, loop_length=i)
                    print("filtering...")
                    filtered = filter_df(iter_df, threshold_GC=k)
                    print(f"recorded:\n {filtered}\n NEXT!")
                    df_chunk = pd.concat([df_chunk, filtered], axis=0)
                except IndexError:
                    pass
    return df_chunk

def full_search(data, loop_len=15, stem_len=15, threshold=0, num_processes=5):
    print('started full search')
    df = pd.DataFrame()
    params_list = []
    
    if loop_len <= 0 or num_processes <= 0:
        print("Error: loop_len and num_processes must be greater than zero.")
        return df  # Return an empty DataFrame to avoid further errors.

    # Split the work into smaller chunks for parallel processing
    chunk_size = loop_len // num_processes
    if chunk_size == 0:
        chunk_size = 1  # Set a minimum chunk size of 1 if loop_len is small compared to num_processes
    
    for i in range(0, loop_len, chunk_size):
        i_range = range(i, min(i + chunk_size, loop_len))
        params_list.append((data, i_range, range(4, stem_len), range(round(stem_len / 2), stem_len), threshold))

    with Pool(num_processes) as pool:
        print('processing chunk!')
        results = pool.map(process_chunk, params_list)

    for result in results:
        df = pd.concat([df, result], axis=0)


    df = df.drop_duplicates(subset="Hairpin_region")
    df.reset_index(inplace=True)
    df.drop(columns="index", inplace=True)
    return df






def compare_dna_sequences(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Input sequences must have the same length")
    
    mismatches = 0
    for base1, base2 in zip(seq1, seq2):
        if base1 != base2:
            mismatches += 1
            
    return mismatches


def filter_df(df_to_filter,threshold_GC):
    """Filtering list by GC%

    Args:
     df_to_filter(DataFrame): contain cordinates and sequences of  all iverted repeats

    Returns:
         : filtered by threshold
    """
    try:
        mask = (df_to_filter['IR1'].str.count('C') + df_to_filter['IR1'].str.count('G')) >= threshold_GC

    # Apply the mask to filter the dataframe
        filtered_df = df_to_filter[mask]
        return filtered_df
    except KeyError:
        pass


if __name__ == "__main__":
    main()
