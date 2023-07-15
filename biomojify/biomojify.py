'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Andrew Lonsdale, 22 Jul 2020 
License     : BSD 
Maintainer  : Andrew Lonsdale
Portability : POSIX

The program reads one or more input FASTA files and convets them to emoji.
'''

from argparse import ArgumentParser, Namespace
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
# from . import biomojify_map
import ast
import vcf
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq



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

# only suitable for english speakers
prot_seq_emoji_map = {
    'A': ':green_apple:',
    'B': ':banana:',
    'C': ':cherries:',
    'D': ':doughnut:',
    'E': ':elephant:',
    'F': ':fries:',
    'G': ':grapes:',
    'H': ':hamburger:',
    'I': ':icecream:',
    'J': ':jeans:',
    'K': ':key:',
    'L': ':lemon:',
    'M': ':mushroom:',
    'N': ':nose:',
    'O': ':octopus:',
    'P': ':pineapple:',
    'Q': ':princess:',
    'R': ':rabbit:',
    'S': ':strawberry:',
    'T': ':tomato:',
    'U': ':umbrella:',
    'V': ':volcano:',
    'W': ':watermelon:',
    'X': ':x:',
    'Y': ':sailboat:',
    '*': ':hand:',
    '-': ':wavy_dash:',
}

prot_scale = "ABCDEFGHIJKLMNOPQRSTUVWXY*-"




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
    subparsers = parser.add_subparsers(help='sub-command help',dest='subparser_name')

    # FASTA processing
    parser_fasta = subparsers.add_parser('fasta', help='fasta --help')
    parser_fasta.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to include in stats (default {})'.format(
            DEFAULT_MIN_LEN))
    parser_fasta.add_argument('--custom',
                              metavar='CUSTOM_DICT',
                              type=str,
                              help='use a mapping of custom emoji to nucleotides in CUSTOM_DICT (' + emojify(":yellow_heart:") + emojify(
                                  ":blue_heart:") + ')')
    parser_fasta.add_argument('fasta_files',
                              nargs='*',
                              metavar='FASTA_FILE',
                              type=str,
                              help='Input FASTA files')
    parser_fasta.set_defaults(func=convert_fasta)

    # FASTA protein processing
    parser_fasta_protein = subparsers.add_parser('fasta_protein', help='fasta_protein --help')
    parser_fasta_protein.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to include in stats (default {})'.format(
            DEFAULT_MIN_LEN))
    parser_fasta_protein.add_argument('--custom',
                              metavar='CUSTOM_DICT',
                              type=str,
                              help='use a mapping of custom emoji to proteins in CUSTOM_DICT (' + emojify(":yellow_heart:") + emojify(
                                  ":blue_heart:") + ')')
    parser_fasta_protein.add_argument('fasta_files',
                              nargs='*',
                              metavar='FASTA_FILE',
                              type=str,
                              help='Input FASTA files')
    parser_fasta_protein.set_defaults(func=convert_fasta_protein)

    #TODO add FASTQ parser and convert both sequence and quality      
    # FASTQ processing
    parser_fastq = subparsers.add_parser('fastq', help='fastq --help')
    parser_fastq.add_argument(
        '--minlen',
        metavar='N',
        type=int,
        default=DEFAULT_MIN_LEN,
        help='Minimum length sequence to convert (default {})'.format(
            DEFAULT_MIN_LEN))
    parser_fastq.add_argument('--bin',
                        action='store_true',
                        help='use binned scores (' + emojify(":no_entry_sign:") + emojify(":skull:")
                             + emojify(":poop:") + emojify(":warning:") + " " + emojify(":smile:") + emojify(
                            ":laughing:") + emojify(":sunglasses:") + emojify(":heart_eyes:") + ")")
    parser_fastq.add_argument('--custom',
                              metavar='CUSTOM_DICT',
                              type=str,
                              help='use a mapping of custom emoji to nucleotides in CUSTOM_DICT (' + emojify(":yellow_heart:") + emojify(
                                  ":blue_heart:") + ')')
    parser_fastq.add_argument('--custom_qual',
                              metavar='CUSTOM_DICT',
                              type=str,
                              help='use a mapping of custom emoji to quality scores in CUSTOM_DICT (' + emojify(":moneybag:") + emojify(
                                  ":snake:") + ')')
    parser_fastq.add_argument('fastq_files',
                              nargs='*',
                              metavar='FASTQ_FILE',
                              type=str,
                              help='Input FASTQ files')
    parser_fastq.set_defaults(func=convert_fastq)
    
    
   

    # file  processing template
    parser_vcf = subparsers.add_parser('vcf', help='vcf --help')
    parser_vcf.add_argument('vcf_files',
                              nargs='*',
                              metavar='VCF_FILE',
                              type=str,
                              help='(experimental) Input VCF files')
    parser_vcf.set_defaults(func=convert_vcf)



    
    # 
    # # file  processing template
    # parser_filetype = subparsers.add_parser('filetype', help='filetype help')
    # parser_filetype.add_argument(
    #     '--minlen',
    #     metavar='N',
    #     type=int,
    #     default=DEFAULT_MIN_LEN,
    #     help='Minimum length sequence to include in stats (default {})'.format(
    #         DEFAULT_MIN_LEN))
    # parser_filetype.add_argument('--custom',
    #                           metavar='CUSTOM_DICT',
    #                           type=str,
    #                           help='use a mapping of custom emoji to proteins in CUSTOM_DICT (' + emojify(":yellow_heart:") + emojify(
    #                               ":blue_heart:") + ')')
    # parser_filetype.add_argument('fasta_files',
    #                           nargs='*',
    #                           metavar='FASTA_FILE',
    #                           type=str,
    #                           help='Input FASTA files')
    # parser_filetype.set_defaults(func=convert_filetype)




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



def convert_vcf(options):
    '''Convert VCF file to emoji ''' 
    print("\t".join(["CHROM","POS","ID","REF","ALT","QUAL","FILTER"]))
    if options.vcf_files:
        for vcf_filename in options.vcf_files:
            logging.info("Processing VCF file from %s", vcf_filename)
            try:
                if vcf_filename.endswith(".gz"):
                    vcf_file = gzip.open(vcf_filename, 'rt')
                else:
                    vcf_file = open(vcf_filename)

            except IOError as exception:
                exit_with_error(str(exception), EXIT_FILE_IO_ERROR)
            else:
                with vcf_file:
                    for record in vcf.Reader(vcf_file):
                        print("\t".join([str(a) for a in [record.CHROM,
                                                          record.POS,
                                                          record.ID,
                                                          "".join([a for a in map(get_vcf_emoji, record.REF)]),
                                                          ",".join([get_vcf_emoji(str(rec)) for rec in record.ALT]),
                                                          get_vcf_qual(record.QUAL),
                                                          get_vcf_filter(record.FILTER),
                                                         # record.INFO,
                                                         # record.FORMAT,
                                                         # record.samples
                                                          ]]))
    else:
        logging.info("Processing vcf file from stdin")
        if (binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b'):
            # print("zipped")
            stdin_file = gzip.open(sys.stdin.buffer, 'rt')
        else:
            stdin_file = sys.stdin
        with stdin_file as vcf_file:
            for record in vcf.Reader(vcf_file):
                print("\t".join([str(a) for a in [record.CHROM,
                                                  record.POS,
                                                  record.ID,
                                                  "".join([a for a in map(get_vcf_emoji, record.REF)]),
                                                  ",".join([get_vcf_emoji(str(rec)) for rec in record.ALT]),
                                                  get_vcf_qual(record.QUAL),
                                                  get_vcf_filter(record.FILTER),
                                                  # record.INFO,
                                                  # record.FORMAT,
                                                  # record.samples
                                                  ]]))


def get_vcf_emoji(orig_c, map_dict=local_seq_emoji_map, default=":heart_eyes:"):
     if (orig_c == "None"):
        return(emojify((":x:")))
     #print("orig:",orig_c,"\n")
     return "".join([emojify(map_dict.get(e, ":heart_eyes:")) for e in orig_c])

def get_vcf_qual(quality):
    '''Map a quality value to an emoji'''

    # Hack to do this quickly - use same trick as FASTQE and convert from value to a PHRED encoding then map
    #TODO make this better
    #
    if quality == None:
        bioemojify_qual = emojify(":question:")
    else:
        fake_seq = 'N'
        record_qual = SeqRecord(Seq(fake_seq), id="test", name="lookup",
                                description="example",
                                letter_annotations={'phred_quality': [int(quality)]})
        mapping_dict_qual_use = emaps.fastq_emoji_map_binned
        original_qual = QualityIO._get_sanger_quality_str(record_qual)
        #print(original_qual)
        bioemojify_qual = "".join([emojify(mapping_dict_qual_use.get(s, ":heart_eyes:")) for s in original_qual])

    return(bioemojify_qual)

def get_vcf_filter(filter_val):

   filt_emoji = ""
   if filter_val == None:
       filt_emoji = emojify(":question:")
   elif filter_val == []:
       filt_emoji = emojify(":thumbsup:")
   else:
       filt_emoji = emojify(":thumbsdown:")+":"+ str(",".join(filter_val))

   return(filt_emoji)

# 
# def convert_filetype(options):
#     return


def convert_fasta_protein(options):
    convert_fasta(options, mapping_dict=prot_seq_emoji_map)

    return


def convert_fasta(options, mapping_dict=local_seq_emoji_map):
    '''Convert FASTA file to emoji. If no FASTA files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''

    if options.custom:
        with open(options.custom) as f:
            mapping_dict_use =ast.literal_eval(f.read())
    else:
        mapping_dict_use=mapping_dict

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
                        bioemojify = "".join([emojify(mapping_dict_use.get(s,":heart_eyes:")) for s in original])
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
                         bioemojify = "".join([emojify(mapping_dict_use.get(s,":heart_eyes:")) for s in original])
                         print(bioemojify)

def convert_fastq(options):
    '''Convert FASTQ file to emoji. If no FASTQ files are specified on the command line then
    read from the standard input (stdin).

    Arguments:
       options: the command line options of the program
    Result:
       None
    '''

    if options.custom:
        with open(options.custom) as f:
            mapping_dict_use = ast.literal_eval(f.read())
    else:
        mapping_dict_use = local_seq_emoji_map

    if options.custom_qual:
        with open(options.custom_qual) as f:
            mapping_dict_qual_use = ast.literal_eval(f.read())
    elif options.bin:
        mapping_dict_qual_use = emaps.fastq_emoji_map_binned
    else:
        mapping_dict_qual_use = emaps.fastq_emoji_map

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
                        bioemojify = "".join([emojify(mapping_dict_use.get(s,":heart_eyes:")) for s in original])
                        original_qual = QualityIO._get_sanger_quality_str(seq)
                        bioemojify_qual = "".join([emojify(mapping_dict_qual_use.get(s,":heart_eyes:")) for s in original_qual])
                        print(bioemojify+"\n"+bioemojify_qual)
#                        print(*zip([a for a in bioemojify if a != " "],[b for b in bioemojify_qual if b != " "]))
    else:
        logging.info("Processing FASTQ file from stdin")
        #stats = FastaStats().from_file(sys.stdin, options.minlen)
        if (binascii.hexlify(sys.stdin.buffer.peek(1)[:2]) == b'1f8b'):
            # print("zipped")
            stdin_file = gzip.open(sys.stdin.buffer, 'rt')
        else:
            stdin_file = sys.stdin

        for seq in SeqIO.parse(stdin_file, "fastq"):
                        print(emojify(":arrow_forward:")+"  "+seq.id)
                        #print(">"+seq.id)
                        original = seq.seq
                        bioemojify = "".join([emojify(mapping_dict_use.get(s,":heart_eyes:")) for s in original])
                        original_qual = QualityIO._get_sanger_quality_str(seq)
                        bioemojify_qual = "".join([emojify(mapping_dict_qual_use.get(s,":heart_eyes:")) for s in original_qual])
                        print(bioemojify+"\n"+bioemojify_qual)





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

def run_biomojify(files,subparser_name="fastq",minlen=0,version=False,log=None,bin=False,custom=None, custom_qual=None):
    options = Namespace(files=files, bin=bin, version=version,custom=custom,custom_qual=custom_qual,log=log, minlen=minlen)
    if options.version:
        print(PROGRAM_NAME,PROGRAM_VERSION)
        return
    options.func(options)
                   

def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    try:
        options.func(options)
    except Exception as e:
        s = str(e)
        print("\n\n")
        print(s)
        print("\n\n")
        parse_args(error=True)
# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
