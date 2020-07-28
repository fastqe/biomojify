#!/usr/bin/env python

from setuptools import setup

LONG_DESCRIPTION = """The program reads one or more input FASTA or FASTQ files and converts them to emoji."""


setup(
    name="biomojify",
    version="0.2.0.0",
    author="Andrew Lonsdale",
    author_email="andrew.lonsdale@lonsbio.com.au",
    packages=["biomojify"],
    package_dir={"biomojify": "biomojify"},
    entry_points={
        "console_scripts": ["biomojify = biomojify.biomojify:main"]
    },
    url="https://github.com/fastqe/biomojify",
    license="LICENSE",
    description=("Convert FASTQ and FASTA files to emoji."),
    long_description=(LONG_DESCRIPTION),
    install_requires=["biopython","fastqe","pyvcf"],
)
