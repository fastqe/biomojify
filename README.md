# Overview 

`biomojify` reads one or more input FASTA or FASTQ files, and converts them to emoji.

## FASTA files
For a DNA sequence: 
```
$ cat test.fasta 
>SEQUENCE_1
ATAGTCAGTACGTAGTCGATGCTAGCTAGTAGGGGGGCTGATGATGTAGCTCGATCGTACGTACGTACGCTGAGTCAGTG
CACGTACGCTGCATGCCNTAGNCTAAAAGNCTANGCTAGCNTANGCTGACTNAGNTGACTGCNTCGTCGNATCATGTACG
TAGCGAGCTTTTTTTAGTGTACGTAGTACTACCCCCCCCCGTACGTACGTACGTCACGTGCTGACTNNNACGATCGTAGT
AGCTGACTGATGCT
```

`biomojify` will convert the ATCG to 🥑🍅🌽🍇:

```
$ biomojify fasta test.fasta 
▶️ SEQUENCE_1
🥑🍅🥑🍇🍅🌽🥑🍇🍅🥑🌽🍇🍅🥑🍇🍅🌽🍇🥑🍅🍇🌽🍅🥑🍇🌽🍅🥑🍇🍅🥑🍇🍇🍇🍇🍇🍇🌽🍅🍇🥑🍅🍇🥑🍅🍇🍅🥑🍇
🌽🍅🌽🍇🥑🍅🌽🍇🍅🥑🌽🍇🍅🥑🌽🍇🍅🥑🌽🍇🌽🍅🍇🥑🍇🍅🌽🥑🍇🍅🍇🌽🥑🌽🍇🍅🥑🌽🍇🌽🍅🍇🌽🥑🍅🍇🌽🌽❓
🍅🥑🍇❓🌽🍅🥑🥑🥑🥑🍇❓🌽🍅🥑❓🍇🌽🍅🥑🍇🌽❓🍅🥑❓🍇🌽🍅🍇🥑🌽🍅❓🥑🍇❓🍅🍇🥑🌽🍅🍇🌽❓🍅🌽🍇🍅
🌽🍇❓🥑🍅🌽🥑🍅🍇🍅🥑🌽🍇🍅🥑🍇🌽🍇🥑🍇🌽🍅🍅🍅🍅🍅🍅🍅🥑🍇🍅🍇🍅🥑🌽🍇🍅🥑🍇🍅🥑🌽🍅🥑🌽🌽🌽🌽🌽
🌽🌽🌽🌽🍇🍅🥑🌽🍇🍅🥑🌽🍇🍅🥑🌽🍇🍅🌽🥑🌽🍇🍅🍇🌽🍅🍇🥑🌽🍅❓❓❓🥑🌽🍇🥑🍅🌽🍇🍅🥑🍇🍅🥑🍇🌽🍅🍇
🥑🌽🍅🍇🥑🍅🍇🌽🍅

```

Protein sequences and custom emoji dictionaries are also supported.

```
▶️ Protein
🍄🍅🐘🍦🍅🍏🍏🍄🌋🔑🐘🍋🐰🐘🍓🍅🍇🍏🍇🍄🍄🍩🍒🔑👃🍏🍋🍓🐘🍅👃🍇🍩🍟🍩🔑🍏🌋👸🍋🍋🐰🐘🔑🍇🍋🍇🔑🍏
🍏🔑🔑🍏🍩🐰🍋🍏🍏🐘🍇🍋🌋🍓🌋🔑🌋🍓🍩🍩🍟🍅🍦🍏🍏🍄🐰🍍🍓⛵🍋🍓⛵🐘🍩🍋🍩🍄🍅🍟🌋🐘👃🐘⛵🔑🍏🍋🌋
🍏🐘🍋🐘🔑🐘👃🐘🐘🐰🐰🐰🍋🔑🍩🍍👃🔑🍍🐘🍔🔑🍦🍍👸🍟🍏🍓🐰🔑👸🍋🍓🍩🍏🍦🍋🔑🐘🍏🐘🐘🔑🍦🔑🐘🐘🍋🔑
🍏👸🍇🔑🍍🐘🔑🍦🍉🍩👃🍦🍦🍍🍇🔑🍄👃🍓🍟🍦🍏🍩👃🍓👸🍋🍩🍓🔑🍋🍅🍋🍄🍇👸🍟⛵🌋🍄🍩🍩🔑🔑🍅🌋🐘👸🌋
🍦🍏🐘🔑🐘🔑🐘🍟🍇🍇🔑🍦🔑🍦🌋🐘🍟🍦🍒🍟🐘🌋🍇🐘🍇🍋🐘🔑🔑🍅🐘🍩🍟🍏🍏🐘🌋🍏🍏👸🍋
```


## FASTQ files

For a FASTQ file, both sequence and quality information are converted to emoji:

```
$ cat test.fq
@ Sequence
GTGCCAGCCGCCGCGGTAGTCCGACGTGGC
+
GGGGGGGGGGGGGGGGGGGGGG!@#$%&%(
```
```
$ biomojify fastq test.fq
▶️  Sequence
🍇🍅🍇🌽🌽🥑🍇🌽🌽🍇🌽🌽🍇🌽🍇🍇🍅🥑🍇🍅🌽🌽🍇🥑🌽🍇🍅🍇🍇🌽
😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁😁🚫😄👺💔🙅👾🙅💀
```

Binned quality scores, and custom emoji dictionaries for either sequence or quality are also supported ([fastqe format](https://github.com/fastqe/fastqe#custom)).

```
▶️  Sequence 
💜💛💜💛💚💙💜💚💚💜💚💚💜💚💜💜💛💙💙💛💙💚💜💜💙💜💜💜💛💜💚💙💙💜💚💜💛💛💙💛💚💚💜💜💙💛💛💛💙💛
😎😎😎😎😎😆😆😆😎⚠️😆😎😎😆😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😎😆😆😎😎😎😎😎😎😎😆😄😆😎😎😎💩
```



## VCF files

Converting VCF files is in early development:

```
$biomojify vcf test.vcf
CHROM  POS      ID         REF       ALT          QUAL  FILTER
20     14370    rs6054257  🍇        🥑           😄    👍
20     17330    None       🍅        🥑           💀    👎:q10
20     1110696  rs6040355  🥑        🍇,🍅        😍    👍
20     1230237  None       🍅        ❌             😍    👍
20     1234567  microsat1  🍇🍅🌽    🍇,🍇🍅🌽🍅  😍    👍
Z      99999    deletion   🍅🍅🍅🍅  ❌             ❓      ❓
```

## Help message

`biomojify` can display usage information on the command line via the `-h` or `--help` argument:


```
$ biomojify --help
usage: biomojify [-h] [--version] [--log LOG_FILE] {fasta,fastq} ...

Read one or more FASTA files, and convert them to emoji.😀

positional arguments:
  {fasta,fasta_protein,fastq,vcf}
                        sub-command help
    fasta               fasta --help
    fasta_protein       fasta_protein --help
    fastq               fastq --help
    vcf                 vcf --help

optional arguments:
  -h, --help      show this help message and exit
  --version       show program's version number and exit
  --log LOG_FILE  record program progress in LOG_FILE
```


# Development

## New file formats

To add a new file format to `biomojify`, create a fork on GitHub with the name of the format you wish to to add, and
submt a pull request with a subparser: 
```
parser_filetype = subparsers.add_parser('filetype', help='filetype help')
parser_filetype.add_argument('filetype_files',
                           nargs='*',
                           metavar='FILETYPE_FILE',
                           type=str,
                           help='Input FILETYPE files')
parser_filetype.set_defaults(func=convert_filetype)
```

and with a function name that refers to a new function you create that will process that file type:

```
def convert_filetype(options):

    for file in options.filetype_files:
        # print emoji format 
```

That's it! 

## Contributors

`biomojify` extends the ideas and work in `fastqe` and shares many of the same contributors: https://github.com/fastqe/fastqe#contributors

# Licence

This program is released as open source software under the terms of [BSD License](https://raw.githubusercontent.com/fastqe/biomojify/master/LICENSE).

# Installing

You can install biomojify via `pip`, directly from the source code or build and run it from within Docker container (experimental).


## Pip

`pip install biomojify`

## Installing directly from source code

Clone this repository: 
```
$ git clone https://github.com/fastqe/biomojify
```

Move into the repository directory:
```
$ cd biomojify
```

Python 3 is required for this software.

Fastqe-convert can be installed using `pip` in a variety of ways (`$` indicates the command line prompt):

1. Inside a virtual environment:
```
$ python3 -m venv biomojify_dev
$ source biomojify_dev/bin/activate
$ pip install -U /path/to/biomojify
```
2. Into the global package database for all users:
```
$ pip install -U /path/to/biomojify
```
3. Into the user package database (for the current user only):
```
$ pip install -U --user /path/to/biomojify
```



## Logging

If the ``--log FILE`` command line argument is specified, biomojify will output a log file containing information about program progress. The log file includes the command line used to execute the program, and a note indicating which files have been processes so far. Events in the log file are annotated with their date and time of occurrence. 

```
$ biomojify --log bt.log file1.fasta file2.fasta 
```
```
$ cat bt.log
2016-12-04T19:14:47 program started
2016-12-04T19:14:47 command line: /usr/local/bin/biomojify --log bt.log file1.fasta file2.fasta
2016-12-04T19:14:47 Processing FASTA file from file1.fasta
2016-12-04T19:14:47 Processing FASTA file from file2.fasta
```


## Exit status values

`biomojify` returns the following exit status values:

* 0: The program completed successfully.
* 1: File I/O error. This can occur if at least one of the input FASTA files cannot be opened for reading. This can occur because the file does not exist at the specified path, or biomojify does not have permission to read from the file. 
* 2: A command line error occurred. This can happen if the user specifies an incorrect command line argument. In this circumstance biomojify will also print a usage message to the standard error device (stderr).

# Running within the Docker container

The following section describes how to run biomojify within the Docker container. It assumes you have Docker installed on your computer and have built the container as described above. 
The container behaves in the same way as the normal version of biomojify, however there are some Docker-specific details that you must be aware of.

The general syntax for running biomojify within Docker is as follows:
```
$ docker run -i biomojify CMD
```
where CMD should be replaced by the specific command line invocation of biomojify. Specific examples are below.

Display the help message:
```
$ docker run -i biomojify biomojify -h
```
Note: it may seem strange that `biomojify` is mentioned twice in the command. The first instance is the name of the Docker container and the second instance is the name of the biomojify executable that you want to run inside the container.

Display the version number:
```
$ docker run -i biomojify biomojify --version
```

Read from a single input FASTA file redirected from standard input:
```
$ docker run -i biomojify biomojify fasta < file.FASTA 
```

Read from multuple input FASTA files named on the command line, where all the files are in the same directory. You must replace `DATA` with the absolute file path of the directory containing the FASTA files:  
```
$ docker run -i -v DATA:/in biomojify biomojify fasta /in/file1.fasta /in/file2.fasta /in/file3.fasta
```
The argument `DATA:/in` maps the directory called DATA on your local machine into the `/in` directory within the Docker container.

Logging progress to a file in the directory OUT: 
```
$ docker run -i -v DATA:/in -v OUT:/out biomojify biomojify fasta --log /out/logfile.txt /in/file1.fasta /in/file2.fasta /in/file3.fasta
```
Replace `OUT` with the absolute path of the directory to write the log file. For example, if you want the log file written to the current working directory, replace `OUT` with `$PWD`.
As above, you will also need to replace `DATA` with the absolite path to the directory containing your input FASTA files.

# Testing

## Unit tests

## Test suite

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[biomojify issue tracker](https://github.com/fastqe/biomojify/issues)


We are also on Gitter along with `fastqe`: https://gitter.im/fastqe/community#

