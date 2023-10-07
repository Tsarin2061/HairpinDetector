#!/usr/bin/env python3
# Author: Lev Tsarin

import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import os
import pandas as pd


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
    # defining global variables with CL arguments as a values
    # global fasta_file, stem_length, loop_length, threshold_GC, output_names
    fasta_file, stem_length, loop_length, threshold_GC, output_names, search_all = parse_arg()
    if search_all == True:
        seq = read_file(fasta_file)
        for id,subseq in seq.items():
            final_df = full_search(subseq, stem_length,loop_length).reset_index().drop(columns = 'index')
            final_df = final_df.assign(ID = id)
            final_df = final_df[["ID"] + [col for col in final_df.columns if col != 'ID']]
            final_df.to_csv(f'{output_names}_{id}.csv')

        directory_path = os.path.dirname(output_names)
        contents = os.listdir(directory_path)
        print(directory_path)
        dfs = []
        for file in contents:
            if os.path.isfile(f"{directory_path}/{file}"):
                df_to_merge = pd.read_csv(f"{directory_path}/{file}")
                dfs.append(df_to_merge)

        merged_df = pd.concat(dfs, ignore_index=True)
        merged_df.to_csv(f'{directory_path}/all_merged.csv')
        print(f"\n  Results are stored in {output_names}.csv\n  Please use different output name to avoid overwriting data!")
    else:
        # creating df instead of lists
        seq = read_file(fasta_file)
        if type(seq) != dict:
            final_df = filter_df(parse_seq(seq, stem_length,loop_length),threshold_GC).reset_index().drop(columns = 'index')
        else:
            dfs = []
            for id,subseq in seq.items():
                final_df = filter_df(parse_seq(subseq, stem_length,loop_length),threshold_GC).reset_index().drop(columns = 'index')
                final_df = final_df.assign(ID = id)
                final_df = final_df[["ID"] + [col for col in final_df.columns if col != 'ID']]
                dfs.append(final_df)
            merged_df = pd.concat(dfs, ignore_index=True)
            print(merged_df)
        if len(final_df) == 0:
            sys.exit("Hairpins were not found!")
        else:
            print(final_df)
            final_df.to_csv(f'{output_names}.csv')
            f"\n  Results are stored in {output_names}.csv\n  Please use different output name to avoid overwriting data!"



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
    seq = Seq(seq)
    start = 0
    end = ir_length
    ir2_length = ir_length
    # list with coordinates and sequences
    list = []
    while start != len(seq) - 2 * ir_length + loop_length:
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
                            f"{loop_length}",
                            f"{ir_length}"
                            
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
                            f"{loop_length}",
                            f"{ir_length}"     
                        ]

                start += 1
                end += 1
                ir2_length = ir_length + start

            else:
              ir2_length += 1
              pass
        if start == len(seq) - 2 * ir_length + loop_length:
            return seqs_df
        

def full_search(data,loop_len = 15, stem_len = 15, threshold = 0):
    df = pd.DataFrame()
    for i in range(loop_len):
        for j in range(4,stem_len):
            for k in range(round(stem_len/2),stem_len):
                try:
                    iter_df = parse_seq(seq = data,loop_length = i,ir_length = j)
                    filtered = filter_df(df,threshold_GC=k)
                    df = pd.concat([df,iter_df],axis = 0)
                except IndexError:
                    pass
    df = df.drop_duplicates(subset="Hairpin_region")
    df.reset_index(inplace=True)
    df.drop(columns="index",inplace=True)
    return df

def parse_imperfect_palindormes(sequence, palindrome_len, threshold, max_distance):
    seq = Seq(sequence)
    palindrome_list = []
    
    for i in range(len(seq) - palindrome_len + 1):
        position1_start = i
        position1_end = i + palindrome_len
        position2_start = i + palindrome_len
        position2_end = i + palindrome_len * 2
        gap = 0
        
        while position2_end <= len(seq) and gap <= max_distance:
            first = seq[position1_start:position1_end]
            second = str(Seq(seq[position2_start:position2_end]).reverse_complement())
            mismatches = compare_dna_sequences(str(first), str(second))
            
            if mismatches <= threshold:
                palindrome_list.append(seq[position1_start:position2_end])
            
            position2_start += 1
            position2_end += 1
            gap += 1
    
    return str(palindrome_list[0])



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
