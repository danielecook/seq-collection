#
# Author: Daniel E. Cook
#
import sugar
import argparse
import strformat
import colorize
import tables
import streams
import hts
import json
import tables
import strutils
import sequtils
import terminal
import asyncfile
import zip/gzipfiles
import src/fqmeta
import src/utils
from constants import ANN_header

from posix import signal, SIG_PIPE, SIG_IGN
signal(SIG_PIPE, SIG_IGN)

proc get_vcf_it(vcf: VCF): iterator(): Variant =
  return iterator(): Variant =
    for i in vcf:
      yield i

iterator variants*(vcf:VCF, regions: seq[string]): Variant =
    ## iterator over region or just the variants.
    if regions.len == 0:
        for v in vcf: yield v
    for region in regions:
        if fileExists(region):
            ## must be in bed format.
            for l in region.lines:
                if l[0] == '#' or l.strip().len == 0: continue
                var toks = l.strip().split(seps={'\t'})
                for v in vcf.query(&"{toks[0]}:{parseInt(toks[1]) + 1}-{toks[2]}"):
                    yield v
        else:
            for v in vcf.query(region): yield v

var fmt_zip = newJObject()
var fmt_arr = newJArray()
let null_json = newJNull()

proc out_fmt[T](record: T, fmt_field: FormatField, zip: bool, samples: seq[string]): JsonNode =
    ## For fascilitating formatting FORMAT fields
    var 
        idx_start: int
        idx_end: int
        rec_out: JsonNode
    for idx in 0..<samples.len:
        if fmt_field.n_per_sample == 1:
            if record[idx] < 0:
                rec_out = %* null_json
            else:
                rec_out = %* record[idx]
        else:
            idx_start = idx * fmt_field.n_per_sample
            idx_end = idx * fmt_field.n_per_sample + fmt_field.n_per_sample - 1
            var rec_arr = newJArray()
            for i in idx_start..idx_end:
                if record[i] < 0:
                    rec_arr.add(%* nil)
                else:
                    rec_arr.add(%* record[i])
            rec_out = %* rec_arr
        if zip:
            fmt_zip.add(samples[idx], rec_out)
        else:
            fmt_arr.add(rec_out)
    if zip:
        return fmt_zip
    else:
        return fmt_arr


proc to_json(vcf: string, region_list: seq[string], sample_set: string, info: string, format: string, zip: bool, annotation: bool, pretty: bool, array: bool) =
    var v:VCF

    ## Format Fields
    let info_keep = filterIt(info.split({',', ' '}), it.len > 0)
    let format_keep = filterIt(format.split({',',' '}), it.len > 0)
    
    ## Output fields
    var field_float = newSeq[float32](4)
    var field_int = newSeq[int32](4)
    var field_string = new_string_of_cap(4)

    doAssert open(v, vcf)
    if sample_set != "ALL":
        let samples_keep = filterIt(sample_set.split({',', ' '}), it.len > 0)
        set_samples(v, samples_keep)
    var samples = v.samples

    if array:
        echo "["
    for rec in variants(v, region_list):
        var info = rec.info
        var format = rec.format
        # Fetch INFO Fields
        var j_info = newJObject()
        let output_all_info = ("ALL" in info_keep or annotation)
        if output_all_info or info_keep.len >= 1:
            for info_field in rec.info.fields:
                if annotation and info_field.name == "ANN":
                    if info_field.name == "ANN":
                        discard info.get(info_field.name, field_string)
                        var ann_record_set = newJArray()
                        for ANN in field_string.split(","):
                            var ann_record = newJObject()
                            var ann_split = ANN.split("|")
                            for ann_col in 0..<ANN_header.len:
                                ann_record.add(ANN_header[ann_col], newJString(ann_split[ann_col]))
                            ann_record_set.add(ann_record)
                        j_info.add("ANN", ann_record_set)
                    elif info_field.name == "BCSQ":
                        discard
                elif output_all_info or info_field.name in info_keep:
                    if info_field.n == 1:
                        if info_field.vtype == BCF_TYPE.FLOAT:
                            ## Single-field Float
                            discard info.get(info_field.name, field_float)
                            j_info.add(info_field.name, %* field_float[0])
                        elif info_field.vtype in [BCF_TYPE.INT8, BCF_TYPE.INT16, BCF_TYPE.INT32]:
                            ## Single-field Int
                            discard info.get(info_field.name, field_int)
                            j_info.add(info_field.name, %* field_int[0])
                    elif info_field.vtype == BCF_TYPE.FLOAT:
                        ## Multi-field Float
                        discard info.get(info_field.name, field_float)
                        j_info.add(info_field.name, %*field_float)
                    elif info_field.vtype in [BCF_TYPE.INT8, BCF_TYPE.INT16, BCF_TYPE.INT32]:
                        ## Multi-field Int
                        discard info.get(info_field.name, field_int)
                        j_info.add(info_field.name, %* field_int)
                    elif info_field.vtype == BCF_TYPE.CHAR:
                        discard info.get(info_field.name, field_string)
                        j_info.add(info_field.name, %* field_string)
                    elif info_field.vtype == BCF_TYPE.NULL:
                        j_info.add(info_field.name, %* true)
        
        ## Fetch FORMAT fields
        var j_format = newJObject()
        let output_all_format = ("ALL" in format_keep)
        if output_all_format or format_keep.len >= 1:
            for format_field in rec.format.fields:
                if output_all_format or format_field.name in format_keep:
                    if format_field.vtype == BCF_TYPE.FLOAT:
                        discard format.get(format_field.name, field_float)
                        j_format.add(format_field.name, out_fmt(field_float, format_field, zip, samples))
                    elif format_field.vtype in [BCF_TYPE.INT8, BCF_TYPE.INT16, BCF_TYPE.INT32]:
                        discard format.get(format_field.name, field_int)
                        j_format.add(format_field.name, out_fmt(field_int, format_field, zip, samples))
            
            
        var jnode = newJObject()
        var variant = newJObject()
        var json_out = %* { "CHROM": $rec.CHROM,
                    "POS": rec.POS,
                    "ID": $rec.ID,
                    "REF": rec.REF,
                    "ALT": rec.ALT,
                    "QUAL": rec.QUAL,
                    "FILTER": rec.FILTER.split(";")}
        if info_keep.len > 0:
            json_out.add("INFO", j_info)
        if format_keep.len > 0:
            json_out.add("FORMAT", j_format)
    
        if pretty:
            stdout.write $json_out.pretty()
        else:
            stdout.write $json_out
        if array:
            stdout.write ",\n"
        else:
            stdout.write "\n"
    if array:
        echo "]"

proc `$`*(f:FileStream): string {.inline.} =
    echo "<filestream>"

proc to_fasta(vcf: string, region_list: seq[string], sample_set: string, force: bool) =
    # vcf - VCF input
    # region_list
    # sample_set
    # force - Force even if is not phased
    var v:VCF
    doAssert open(v, vcf)
    if sample_set != "ALL":
        let samples_keep = filterIt(sample_set.split({',', ' '}), it.len > 0)
        set_samples(v, samples_keep)

    var gts = new_seq[int32](10)

    var allele_set: seq[string]
    var seq_add: cstring
    var sample: string
    # sample → chromosome_n → genotype

    var n_samples = v.samples.len
    # Initialize streams
    echo n_samples
    var sequences = new_seq[seq[FileStream]](n_samples)
    var chrom_set: seq[FileStream]
    for n_sample in 0..<n_samples:
        sample = v.samples[n_sample]
        chrom_set = @[]
        for n_chrom in 0..<2:
            echo fmt"{n_sample} {n_chrom} - {v.samples[n_sample]}"
            chrom_set.add(newFileStream(fmt"{sample}_{n_chrom}.fa", fmWrite))
        sequences[n_sample] = chrom_set

    var n_variant = 0
    var n_chrom = 0
    var n_sample = 0
    var allele_out: string
    var gt_set: seq[string]
    for rec in variants(v, region_list):
        allele_set = sequtils.concat(@[rec.REF], rec.ALT)
        # Record ploidy
        n_sample = 0
        for g in rec.format.genotypes(gts):
            n_chrom = 0
            for a in g:
                if a.phased == false and force == false:
                    quit_error("Genotypes are not phased", 99)
                else:
                    if a.value() != -1:
                        allele_out = allele_set[a.value()]
                    else:
                        allele_out = "N"
                    echo fmt"{allele_out} → {n_sample} → {n_chrom}"
                    sequences[n_sample][n_chrom].write(allele_out)
                n_chrom.inc()
            n_sample.inc()

proc check_file(fname: string): bool =
    if not fileExists(fname):
        print_error(fmt"{fname} does not exist or is not readable")
        return false
    return true

proc get_vcf(vcf: string): string =
    if vcf == "STDIN":
        return "-"
    return vcf

# Convert a stream into an iterator of lines
proc linesIterator(stream: Stream): iterator(): string =
  result = iterator(): string =
    while not stream.atEnd:
      yield stream.readLine()

import tables

proc fq_meta(fasta: string, sample_n = 20) =
    let stream: Stream =
        if fasta[^3 .. ^1] == ".gz":
            newGZFileStream(fasta)
        else:
            newFileStream(fasta, fmRead)
    if stream == nil:
        quit_error("Unable to open file: " & fasta, 2)
    
    var
        qual_min, qual_max: int
        line: string
        barcodes = newSeq[string](sample_n)
        i = 0
    
    while not stream.atEnd() and i < sample_n * 4:
        line = stream.readLine()
        if i %% 4 == 0:
            try:
                echo line
                barcodes[i.div(4)] = line.split(" ")[1].split(":")[^1]
            except IndexError:
                discard
        
        # Quality scores
        if i %% 4 == 3:
            (qual_min, qual_max) = qual_min_max(line, qual_min, qual_max)

        i.inc()
    
    stream.close()

    var fastq_scores = fastq_types.filterIt(qual_min >= it.minimum and qual_max <= it.maximum)
    var most_comm_barcode = barcodes.newCountTable().largest()[0]
    echo "MOST COMMON: ", most_comm_barcode
    echo fastq_scores, "    ", qual_max
    var fastq_scores_name = fastq_scores.mapIt(it.name).join(";")
    var fastq_scores_phred = fastq_scores.mapIt(it.phred).join(";")
    echo fastq_scores_name
    echo fastq_scores_phred


var p = newParser("sc"):
    help("Sequence data utilities")
    command("fq-meta", group="FASTQ"):
        help("Output metadata for FASTQ")
        arg("fastq", nargs = 1, help="fastq file")
        option("-n", "--lines", help="Number of sequences to sample for qual and index/barcode analysis", default = "20")
        flag("header", help="Output just header")
        run:
            fq_meta(opts.fastq, parseInt(opts.lines))
    command("json", group="VCF"):
        help("Convert a VCF to JSON")
        arg("vcf", nargs = 1, help="VCF to convert to JSON")
        arg("region", nargs = -1, help="List of regions")
        option("-i", "--info", help="comma-delimited INFO fields; Use 'ALL' for everything", default="")
        option("-f", "--format", help="comma-delimited FORMAT fields; Use 'ALL' for everything", default="")
        option("-s", "--samples", help="Set Samples", default="ALL")
        flag("-p", "--pretty", help="Prettify result")
        flag("-a", "--array", help="Output as a JSON array instead of ind. JSON lines")
        flag("-z", "--zip", help="Zip sample names with FORMAT fields (e.g. {'sample1': 25, 'sample2': 34})")
        flag("-n", "--annotation", help="Parse ANN Fields")
        flag("--debug", help="Debug")
        run:
            to_json(get_vcf(opts.vcf), opts.region, opts.samples, opts.info, opts.format, opts.zip, opts.annotation, opts.pretty, opts.array)
            quit()
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
    command("filter", group="BAM"):
        run:
            echo "G"

# Check if input is from pipe
var input_params = commandLineParams()
if terminal.isatty(stdin) == false and input_params[input_params.len-1] == "-":
    echo "ADD STIDN"
    input_params[input_params.len-1] = "STDIN"
elif terminal.isatty(stdin) == false:
    if input_params.find("-") > -1:
       input_params[input_params.find("-")] = "STDIN"
    else:
        input_params.add("STDIN")

if commandLineParams().len == 0:
    stderr.write p.help()
    quit()
else:
    try:
        p.run(input_params)
    except UsageError as E:
        input_params.add("-h")
        p.run(input_params)
    #except Exception as E:
    #    quit_error(E.msg)

