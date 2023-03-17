#!/usr/bin/env python3
# Author: Lev Tsarin

import argparse
import xlsxwriter
import sys
from prettytable import PrettyTable
from Bio.Seq import Seq
from Bio import SeqIO



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
    global fasta_file, stem_length, loop_length, threshold_GC, output_names
    fasta_file, stem_length, loop_length, threshold_GC, output_names = parse_arg()
    seq = read_fasta(fasta_file)
    if validate(seq):
        if len(filter_list(parse_seq(seq, stem_length))) == 0:
            sys.exit("Hairpins were not found!")
        else:
            table = PrettyTable()
            table.field_names = ["Coordinates","Inverted repeat1","Hairpin region" ,'Inverted repeat2']
            for i in filter_list(parse_seq(seq, stem_length)):
                table.add_row(i)
            write_table(filter_list(parse_seq(seq, stem_length)))
            print(table)
        print(
            f"\n  Results are stored in {output_names}.xlsx\n  Please use different output name to avoid overwriting data!"
        )
    else:
        print("\n  Your fasta file contains an errors. Check the sequense please!")


def read_fasta(path):
    # reads fasta file
    record = SeqIO.read(path, "fasta")
    return record.seq


def write_table(data):
    wb = xlsxwriter.Workbook(f"{output_names}.xlsx")
    ws = wb.add_worksheet()
    cell = wb.add_format()
    cell.set_text_wrap()
    ws.add_table(
        f"A1:D{len(data)+1}",
        {
            "data": data,
            "columns": [
                {"header": "Coordinates"},
                {"header": "Inverted repeat 1"},
                {"header": "Hairpin region"},
                {"header": "Inverted repeat 2"},
            ],
        },
    )
    wb.close()


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


def parse_seq(seq, l):
    """Inverted repeat detection

    Args:
        seq (string): nucleotide sequence
        l (integer): lenght of inverted repeat we are looking for

    Returns:
        list : [coordinates of hairpin, sequence, inverted sequece]
    """
    seq = Seq(seq)
    start = 0
    end = l
    g = l
    # list with coordinates and sequences
    list = []
    for _ in range(len(seq) ** 2):
        act_loop_len = (g + 1) - end
        if seq[start:end] != Seq(seq[g + 1 : g + l + 1]).reverse_complement():
            g += 1
            if act_loop_len >= loop_length:
                # if inversted repead was not found - we looking for another one
                start += 1
                end += 1
                g = l + start
            else:
                pass
        else:
            # checks if loop meet requirments and correctness of inverted repeats.
            if (
                act_loop_len == loop_length
                and seq[start - 1] != Seq(seq[g + l + 1]).reverse_complement()
                and seq[end] != Seq(seq[g]).reverse_complement()
            ):
                # if we found inverted repeat - we add it to list and look for next one
                list.append(
                    [
                        f"{start}-{g+l+1}",
                        str(seq[start:end]),
                        str(seq[start : g + l + 1]),
                        str(seq[g + 1 : g + l + 1]),
                    ]
                )
                start += 1
                end += 1
                g = l + start

            else:
                g += 1
                pass
        if start == len(seq) - 2 * l + loop_length:
            return list


def filter_list(list_to_filter):
    """Filtering list by GC%

    Args:
        list_to_filter (list): contain cordinates and sequenses of  all iverted repeats

    Returns:
        list : filtering by threshold
    """
    filtered_list = []
    for sample in list_to_filter:
        if sample[1].count("G") + sample[1].count("C") >= threshold_GC:
            filtered_list.append(sample)

    return filtered_list


if __name__ == "__main__":
    main()
