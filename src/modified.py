import pandas as pd

def parse_seq(seq, ir_length):
    """Inverted repeat detection

    Args:
        seq (string): nucleotide sequence
       ir_length(integer): lenght of inverted repeat we are looking for

    Returns:
        list : [coordinates of hairpin, sequence, inverted sequece]
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
         if seq[start:end] != Seq(seq[g + 1 :ir2_length +ir_length+ 1]).reverse_complement():
            ir2_length += 1
            if act_loop_len >= loop_length:
                # if inversted repead was not found - we looking for another one
                start += 1
                end + = 1
                ir2_length = ir_length+ start
            else:
                pass
        else:
            # checks if loop meet requirments and correctness of inverted repeats.
            if (
                act_loop_len == loop_length
                and seq[start - 1] != Seq(seq[g + ir_length + 1]).reverse_complement()
                and seq[end] != Seq(seq[g]).reverse_complement()
            ):
                # if we found inverted repeat - we add it to list and look for next one
                seqs_df.loc[len(seqs_df)] =[
                    # coordinates
                        f"{start}-{g+l+1}",
                    # IRs
                        str(seq[start:end]),
                        str(seq[g + 1 : ir2_l ength +ir_length+ 1]),
                    # entire hairpin
                        str(seq[start : ir2_lengt h +ir_length+ 1]),
                    # extended hairpin region
                        str(seq[start-15:   ir2_length +ir_length+ 1+15]),
                    # coordinates #2
                        f"{start-15}-{g+l+1+15}",
                        

                    ]            
                start += 1
                end + = 1
                ir2_length = ir_length + start

            else:
              ir2_length += 1
                pass
        if start == len(seq) - 2 * ir_length + loop_length:
            return seqs_df