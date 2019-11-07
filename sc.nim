#
# Author: Daniel E. Cook
#
import sugar
import argparse
import strformat
import colorize
import tables
import hts
import tables
import strutils
import sequtils
import terminal
import asyncfile
import zip/gzipfiles

import src/fq_meta
import src/fq_count
import src/fq_dedup

import src/insert_size

import src/vcf2fasta
import src/vcf2tsv
import src/vcf2json
import src/vcf_window

import src/utils/helpers

from posix import signal, SIG_PIPE, SIG_IGN
signal(SIG_PIPE, SIG_IGN)

const VERSION = "0.0.1"

proc get_vcf(vcf: string): string =
    if vcf == "STDIN":
        return "-"
    return vcf


import tables

var p = newParser("sc"):
    flag("--debug", help="Debug")
    help(fmt"Sequence data utilities (Version {VERSION})")
    command("fq-meta", group="FASTQ"):
        help("Output metadata for FASTQ")
        arg("fastq", nargs = -1, help="List of FASTQ files")
        option("-n", "--lines", help="Number of sequences to sample (n_lines) for qual and index/barcode determination", default = "100")
        flag("--header", help="Output the header")
        flag("-s", "--symlinks", help="Follow symlinks")
        run:
            if opts.fastq.len == 0:
                quit_error("No FASTQ specified", 3)
            if opts.fastq.len > 0:
                for fastq in opts.fastq:
                    fq_meta.fq_meta(fastq, parseInt(opts.lines), opts.symlinks, opts.header)
    command("fq-count", group="FASTQ"):
        help("Counts lines in a FASTQ")
        flag("--header", help="Output just header")
        arg("fastq", nargs = -1, help = "Input FASTQ")
        run:
            if opts.header:
                echo fq_count_header
            if opts.fastq.len == 0:
                quit_error("No FASTQ specified", 3)
            if opts.fastq.len > 0:
                for fastq in opts.fastq:
                    fq_count.fq_count(fastq)
    command("fq-dedup", group="FASTQ"):
        help("Removes exact duplicates from FASTQ Files")
        arg("fastq", nargs = 1, help = "Input FASTQ")
        run:
            fq_dedup.fq_dedup(opts.fastq)
    
    command("json", group="VCF"):
        help("Convert a VCF to JSON")
        arg("vcf", nargs = 1, help="VCF to convert to JSON")
        arg("region", nargs = -1, help="List of regions")
        option("-i", "--info", help="comma-delimited INFO fields; Use 'ALL' for everything", default="")
        option("-f", "--format", help="comma-delimited FORMAT fields; Use 'ALL' for everything", default="")
        option("-s", "--samples", help="Set Samples", default="ALL")
        flag("-p", "--pretty", help="Prettify result")
        flag("-a", "--array", help="Output as a JSON array instead of individual JSON lines")
        flag("-z", "--zip", help="Zip sample names with FORMAT fields (e.g. {'sample1': 25, 'sample2': 34})")
        flag("-n", "--annotation", help="Parse ANN Fields")
        flag("--pass", help="Only output variants where FILTER=PASS")
        flag("--debug", help="Debug")
        run:
            to_json(get_vcf(opts.vcf), opts.region, opts.samples, opts.info, opts.format, opts.zip, opts.annotation, opts.pretty, opts.array, opts.pass)
    
    command("insert-size", group="BAM"):
        help("Calculate insert-size metrics")
        flag("--header", help="Output the header")
        option("-d", "--dist", default="", help = "Output raw distribution(s)")
        arg("bam", nargs = -1, help = "Input BAM")
        run:
            if opts.header:
                echo insert_size.header
            if opts.bam.len == 0:
                quit_error("No BAM specified", 3)
            if opts.bam.len > 0:
                for bam in opts.bam:
                    insert_size.cmd_insert_size(bam, opts.dist)

    command("vcf2tsv", group="VCF"):
        help("Converts a VCF to TSV or CSV")
        flag("--header", help="Output the header")
        flag("--long", help="Output in long format instead of wide")
        arg("vcf", nargs = 1, help = "Input FASTQ")
        arg("regions", nargs = -1, help = "Regions to subset on")
        run:
            if opts.vcf.len == 0:
                quit_error("No VCF specified", 3)
            elif opts.vcf.len > 0:
                vcf2tsv(opts.vcf, opts.long, opts.regions)

    command("window", group="VCF"):
        help("Generate windows from a VCF for parallel execution")
        arg("vcf", nargs = 1, help = "Input VCF")
        arg("width", nargs = 1)
        run:
            vcf_window(opts.vcf)
    
    command("fasta", group="VCF"):
        help("Convert a VCF to a FASTA file")
        arg("vcf", nargs = 1, help="VCF to convert to JSON")
        arg("region", nargs = -1, help="List of regions or bed files")
        option("-s", "--samples", help="Set Samples", default="ALL")
        option("-r", "--reference", help="Output full reference sequence")
        flag("-f", "--force", help="Force output even if genotypes are not phased")
        flag("-c", "--concat", help="Combine chromosomes into a single sequence")
        flag("-m", "--merge", help="Merge samples into a single file and send to stdout")
        run:
            to_fasta(get_vcf(opts.vcf), opts.region, opts.samples, opts.force)
    # command("index-swap", group="BAM"):
    #     arg("BAM", nargs= -1, help="List of BAMs or CRAMs to examine")
    #     option("-s", "--sites", help="List of sites to check (required)")
    #     option("-f", "--fasta", help="Reference for use with CRAM files")
    #     option("-t", "--threads", default="1", help="Threads")
    #     run:
    #         index_swaps(opts.BAM, opts.sites, opts.fasta, parseInt(opts.threads))

# Check if input is from pipe
var input_params = commandLineParams()
if getFileInfo(stdin).id.device==0:
    if input_params.find("-") > -1:
       input_params[input_params.find("-")] = "STDIN"
    else:
        input_params.add("STDIN")

if input_params.len <= 1:
    input_params.add("-h")
    p.run(input_params)
else:
    try:
        p.run(input_params)
    except UsageError as E:
        input_params.add("-h")
        stderr.write_line "Error".bgWhite.fgRed & fmt": {E.msg}".fgRed
        if input_params.find("--debug") > -1:
            p.run(input_params)
    except Exception as E:
        if commandLineParams().find("--debug") > -1:
            stderr.write_line "Error".bgWhite.fgRed & fmt": {E.msg}".fgRed
            raise
        else:
            if E.msg != "errno: 32 `Broken pipe`":
                quit_error(E.msg)

proc ctrlc() {.noconv.} =
  echo "Ctrl+C fired!"
  quit_error("G")

setControlCHook(ctrlc)