#!/usr/bin/env python3
# Author: Lev Tsarin

import argparse
import sys
from Bio.Seq import Seq
from Bio import SeqIO
import pandas as pd
from src.viz import seq_visualisation




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

    args = parser.parse_args()
    return (
        args.path,
        args.stem_length,
        args.loop_length,
        args.threshold_GC,
        args.output_name,
    )


def main():
    # defining global variables with CL arguments as a values
    # global fasta_file, stem_length, loop_length, threshold_GC, output_names
    fasta_file, stem_length, loop_length, threshold_GC, output_names = parse_arg()
    # creating df instead of lists
    seq = read_file(fasta_file)    
    final_df = filter_df(parse_seq(seq, stem_length,loop_length),threshold_GC).reset_index().drop(columns = 'index')

    if validate(seq):
        if len(final_df) == 0:
            sys.exit("Hairpins were not found!")
        else:
            print(final_df)
            final_df.to_excel(f'{output_names}.xlsx')
            f"\n  Results are stored in {output_names}.xlsx\n  Please use different output name to avoid overwriting data!"
    else:
        print("\n  Your fasta file contains an errors. Check the sequence please!")

    #vizualization part

    # for index,seq in final_df.iterrows():
    #     print(seq['Adjacent_region(30nt)'])
    #     seq_visualisation(seq['Adjacent_region(30nt)'], name = seq['AR_coordinates'])

def read_file(path):
    # reads fasta file
    record = SeqIO.read(path, "fasta")
    return record.seq


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
    seqs_df = pd.DataFrame(columns =['Coordinates','IR1', 'IR2', 'Hairpin_region', 'Adjacent_region(30nt)','AR_coordinates'])
    seq = Seq(seq)
    start = 0
    end = ir_length
    ir2_length = ir_length
    # list with coordinates and sequences
    list = []
    for i in range(len(seq) ** 2):
        act_loop_len = (ir2_length + 1) - end
        if seq[start:end] != Seq(seq[ir2_length + 1 :ir2_length +ir_length+ 1]).reverse_complement():
            ir2_length += 1
            if act_loop_len >= loop_length:
                # if inversted repead was not found - we looking for another one
                start += 1
                end += 1
                ir2_length = ir_length+ start
            else:
                pass
        else:
            # checks if loop meet requirments and correctness of inverted repeats.
            if (
                act_loop_len == loop_length
                and seq[start - 1] != Seq(seq[ir2_length + ir_length + 1]).reverse_complement()
                and seq[end] != Seq(seq[ir2_length]).reverse_complement()
            ):
                # if we found inverted repeat - we add it to list and look for next one
                seqs_df.loc[len(seqs_df)] =[
                    # coordinates
                        f"{start}-{ir2_length + ir_length +1}",
                    # IRs
                        str(seq[start:end]),
                        str(seq[ir2_length + 1 : ir2_length + ir_length + 1]),
                    # entire hairpin
                        str(seq[start : ir2_length +ir_length+ 1]),
                    # extended hairpin region
                        str(seq[start-15:   ir2_length +ir_length+ 1+15]),
                    # coordinates #2
                        f"{start-15}-{ir2_length+ir_length+1+15}",
                        

                    ]            
                start += 1
                end += 1
                ir2_length = ir_length + start

            else:
              ir2_length += 1
              pass
        if start == len(seq) - 2 * ir_length + loop_length:
            return seqs_df


def filter_df(df_to_filter,threshold_GC):
    """Filtering list by GC%

    Args:
     df_to_filter(DataFrame): contain cordinates and sequences of  all iverted repeats

    Returns:
         : filtered by threshold
    """
    mask = (df_to_filter['IR1'].str.count('C') + df_to_filter['IR1'].str.count('G')) >= threshold_GC

    # Apply the mask to filter the dataframe
    filtered_df = df_to_filter[mask]
    return filtered_df
if __name__ == "__main__":
    main()
