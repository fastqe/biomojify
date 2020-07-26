'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) FASTQE-CONVERT_AUTHOR, 22 Jul 2020 
License     : MIT 
Maintainer  : FASTQE-CONVERT_EMAIL 
Portability : POSIX

The program reads one or more input FASTA files and convets them to emoji.
'''

from argparse import ArgumentParser
from math import floor
import sys
import logging
import pkg_resources
from Bio import SeqIO
from fastqe import fastqe_map as emaps
from pyemojify import emojify 
from Bio.SeqIO import QualityIO
import binascii
import gzip


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_FASTA_FILE_ERROR = 3
DEFAULT_MIN_LEN = 0
DEFAULT_VERBOSE = False
PROGRAM_NAME = "biomojify"

# #PyCharm testing command line processing
# sys.argv = [
#     __file__,
# #    '--bin',
# #    '--long','3000',
# # #   '--output', 'testouput.txt',
#     'fastq',
#     '../functional_tests/test_data/test.fastq',
# #    'test/test.fastq',
# #    'test/test_wiki.fq',
# ]

local_seq_emoji_map = {
    'A': '🥑',  # avocado not in pyemojify, trying a failthrough which works for the noemoji mode
    'C': ':corn:',
    'T': ':tomato:',
    'G': ':grapes:',
    'N': ':question:'
}



try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args(error=False):
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Read one or more FASTA or FASTQ files, and convert them to emoji.😀'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    subparsers = parser.add_subparsers(help='sub-command help')
    
    # FASTA processing
    parser_fasta = subparsers.add_parser('fasta', help='fasta help')
    parser_fasta.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to include in stats (default {})'.format(
            DEFAULT_MIN_LEN))
    parser_fasta.add_argument('fasta_files',
                        nargs='*',
                        metavar='FASTA_FILE',
                        type=str,
                        help='Input FASTA files')
    parser_fasta.set_defaults(func=convert_fasta)

    #TODO add FASTQ parser and convert both sequence and quality      
    # FASTQ processing
    parser_fastq = subparsers.add_parser('fastq', help='fastq help')
    parser_fastq.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to convert (default {})'.format(
            DEFAULT_MIN_LEN))
    parser_fastq.add_argument('fastq_files',
                              nargs='*',
                              metavar='FASTQ_FILE',
                              type=str,
                              help='Input FASTQ files')
    parser_fastq.set_defaults(func=convert_fastq)

    if(error):
        parser.print_help()
        return
    else:
        return parser.parse_args()


class FastaStats(object):
    '''Compute various statistics for a FASTA file:

    num_seqs: the number of sequences in the file satisfying the minimum
       length requirement (minlen_threshold).
    num_bases: the total length of all the counted sequences.
    min_len: the minimum length of the counted sequences.
    max_len: the maximum length of the counted sequences.
    average: the average length of the counted sequences rounded down
       to an integer.
    '''
    #pylint: disable=too-many-arguments
    def __init__(self,
                 num_seqs=None,
                 num_bases=None,
                 min_len=None,
                 max_len=None,
                 average=None):
        "Build an empty FastaStats object"
        self.num_seqs = num_seqs
        self.num_bases = num_bases
        self.min_len = min_len
        self.max_len = max_len
        self.average = average

    def __eq__(self, other):
        "Two FastaStats objects are equal iff their attributes are equal"
        if type(other) is type(self):
            return self.__dict__ == other.__dict__
        return False

    def __repr__(self):
        "Generate a printable representation of a FastaStats object"
        return "FastaStats(num_seqs={}, num_bases={}, min_len={}, max_len={}, " \
            "average={})".format(
                self.num_seqs, self.num_bases, self.min_len, self.max_len,
                self.average)

    def from_file(self, fasta_file, minlen_threshold=DEFAULT_MIN_LEN):
        '''Compute a FastaStats object from an input FASTA file.

        Arguments:
           fasta_file: an open file object for the FASTA file
           minlen_threshold: the minimum length sequence to consider in
              computing the statistics. Sequences in the input FASTA file
              which have a length less than this value are ignored and not
              considered in the resulting statistics.
        Result:
           A FastaStats object
        '''
        num_seqs = num_bases = 0
        min_len = max_len = None
        for seq in SeqIO.parse(fasta_file, "fasta"):
            this_len = len(seq)
            if this_len >= minlen_threshold:
                if num_seqs == 0:
                    min_len = max_len = this_len
                else:
                    min_len = min(this_len, min_len)
                    max_len = max(this_len, max_len)
                num_seqs += 1
                num_bases += this_len
        if num_seqs > 0:
            self.average = int(floor(float(num_bases) / num_seqs))
        else:
            self.average = None
        self.num_seqs = num_seqs
        self.num_bases = num_bases
        self.min_len = min_len
        self.max_len = max_len
        return self

    def pretty(self, filename):
        '''Generate a pretty printable representation of a FastaStats object
        suitable for output of the program. The output is a tab-delimited
        string containing the filename of the input FASTA file followed by
        the attributes of the object. If 0 sequences were read from the FASTA
        file then num_seqs and num_bases are output as 0, and min_len, average
        and max_len are output as a dash "-".

        Arguments:
           filename: the name of the input FASTA file
        Result:
           A string suitable for pretty printed output
        '''
        if self.num_seqs > 0:
            num_seqs = str(self.num_seqs)
            num_bases = str(self.num_bases)
            min_len = str(self.min_len)
            average = str(self.average)
            max_len = str(self.max_len)
        else:
            num_seqs = num_bases = "0"
            min_len = average = max_len = "-"
        return "\t".join([filename, num_seqs, num_bases, min_len, average,
                          max_len])


def convert_fasta(options):
    '''Convert FASTA file to emoji. If no FASTA files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    if options.fasta_files:
        for fasta_filename in options.fasta_files:
            logging.info("Processing FASTA file from %s", fasta_filename)
            try:
                if fasta_filename.endswith(".gz"):
                    fasta_file = gzip.open(fasta_filename, 'rt')
                else:
                    fasta_file = open(fasta_filename)

            except IOError as exception:
                exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
            else:
                with fasta_file:
                    #stats = FastaStats().from_file(fasta_file, options.minlen)
                    for seq in SeqIO.parse(fasta_file, "fasta"):
                        print(emojify(":arrow_forward:") + " " + seq.id)
                        #print(">"+seq.id)
                        original = seq.seq
                        bioemojify = "".join([emojify(local_seq_emoji_map.get(s,":heart_eyes:")) for s in original])
                        print(bioemojify)
    else:
        logging.info("Processing FASTA file from stdin")
        #stats = FastaStats().from_file(sys.stdin, options.minlen)
        if (binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b'):
            # print("zipped")
            stdin_file = gzip.open(sys.stdin.buffer, 'rt')
        else:
            stdin_file = sys.stdin
        for seq in SeqIO.parse(stdin_file, "fasta"):
                         print(emojify(":arrow_forward:") + " " + seq.id)
                         #print(">"+seq.id)
                         original = seq.seq
                         bioemojify = "".join([emojify(local_seq_emoji_map.get(s,":heart_eyes:")) for s in original])
                         print(bioemojify)

def convert_fastq(options):
    '''Convert FASTQ file to emoji. If no FASTQ files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    if options.fastq_files:
        for fastq_filename in options.fastq_files:
            logging.info("Processing FASTA file from %s", fastq_filename)
            try:
                if fastq_filename.endswith(".gz"):
                    fastq_file = gzip.open(fastq_filename, 'rt')
                else:
                    fastq_file = open(fastq_filename)

            except IOError as exception:
                exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
            else:
                with fastq_file:
                    for seq in SeqIO.parse(fastq_file, "fastq"):
                        print(emojify(":arrow_forward:")+"  "+seq.id)
                        #print(">"+seq.id)
                        original = seq.seq
                        bioemojify = "".join([emojify(local_seq_emoji_map.get(s,":heart_eyes:")) for s in original])
                        original_qual = QualityIO._get_sanger_quality_str(seq)
                        bioemojify_qual = "".join([emojify(emaps.fastq_emoji_map.get(s,":heart_eyes:")) for s in original_qual])
                        print(bioemojify+"\n"+bioemojify_qual)
#                        print(*zip([a for a in bioemojify if a != " "],[b for b in bioemojify_qual if b != " "]))
    else:
        logging.info("Processing FASTA file from stdin")
        #stats = FastaStats().from_file(sys.stdin, options.minlen)
        if (binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b'):
            # print("zipped")
            stdin_file = gzip.open(sys.stdin.buffer, 'rt')
        else:
            stdin_file = sys.stdin

        for seq in SeqIO.parse(stdin_file, "fasta"):
                         print(">"+seq.id)
                         original = seq.seq
                         bioemojify = " ".join([emojify(local_seq_emoji_map.get(s,":heart_eyes:")) for s in original])
                         print(bioemojify)





def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt="%Y-%m-%dT%H:%M:%S%z")
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    try:
        options.func(options)
    except:
        parse_args(error=TRUE)
# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
