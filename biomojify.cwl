#!/usr/bin/env cwl-runner

baseCommand: [biomojify]
class: CommandLineTool
cwlVersion: v1.0
id: biomojify
inputs:
  - doc: Input FASTA files
    id: fasta_file
    inputBinding: {position: 0}
    type: File
  - doc: Minimum length sequence to include in stats (default 0)
    id: min_len
    inputBinding: {prefix: --minlen}
    type: long?
  - doc: record program progress in LOG_FILE
    id: log
    inputBinding: {prefix: --log}
    type: string?
outputs:
  - doc: Stats file
    id: stats
    type: stdout
requirements:
  - class: DockerRequirement
    dockerPull: biomojify
