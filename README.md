# HairpinDetector
 Poorly written tiny script for detection of hairpins in DNA sequence using CLI.  

## Dependencies
*Operating Systems:*  Linux, MacOS

*Required:* 
* Python 3.10.7+
* Biopython 1.81
* XlsxWriter 3.0.9
* Prettytable 3.6.0 


## Setup
Find a location where you want to locate the script and clone the script using git:

    git clone git@github.com:Tsarin2061/HairpinDetector.git
Installation of requirements:

*You should install python manually, other requirements will be installed with the following command:

    pip install -r requirements.txt


## Input/output
The script can accept **fasta** files (so far only with one record) as input.
The output of the script creates CLI table and also an excel file with coordinates, sequences of inverted repeats, and the entire region of the hairpin. 

Highly recommended to create a separate folder for results.


## Example of usage
The following code has to be executed via the terminal in the relative directory:

    ./HD.py -p <path to file> -s <length of stem>  -l <length of loop> -t <amount of GCs in inverted repeat> -o <path/to/file>

To see more details about arguments:
    
    ./HD.py -h


## Data to test&train
There is a test file called **NM_006231.4.fasta**. You can use it to test the script.

The following code:

    ./HD.py -p NM_006231.4.fasta -l 4 -s 5 -t 0 -o /home/bio/folder1/folder2/results

Should prompt the following output:

    +-------------+------------------+----------------+------------------+
    | Coordinates | Inverted repeat1 | Hairpin region | Inverted repeat2 |
    +-------------+------------------+----------------+------------------+
    |   774-788   |      GGAAA       | GGAAATGCTTTTCC |      TTTCC       |
    |  3035-3049  |      TGGCT       | TGGCTCTGTAGCCA |      AGCCA       |
    |  4527-4541  |      GGGAT       | GGGATCTTCATCCC |      ATCCC       |
    |  5990-6004  |      GGCAG       | GGCAGCCTCCTGCC |      CTGCC       |
    |  7200-7214  |      CACCA       | CACCACACCTGGTG |      TGGTG       |
    |  7229-7243  |      CTAGG       | CTAGGTGCCCCTAG |      CCTAG       |
    +-------------+------------------+----------------+------------------+

        Results are stored in /home/bio/folder1/folder2/results.xlsx
        Please use a different output name to avoid overwriting data!

        
You can also set your parameters!

## TO-DO list
* visualization of hairpin - **WIP**
* reading fasta files with few records  - **WIP**
* ... - **WIP**


