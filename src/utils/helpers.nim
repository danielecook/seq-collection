import os
import hts
import math
import streams
import zip/gzipfiles
import strutils
import strformat
import sequtils
import colorize
import sequtils
import tables
import nre

proc ending*(s: string, endings: seq[string]): bool = 
    for i in endings:
        if s.endswith(i):
            return true
    return false

proc is_fasta*(s: string): bool = 
    return s.toLower().ending(@[".fa.gz", ".fa", ".fasta", ".fasta.gz"])

proc is_bam*(s: string): bool = 
    return s.toLower().ending(@[".sam", ".bam", ".cram"])

proc is_vcf*(s: string): bool = 
    return s.toLower().ending(@[".vcf", ".vcf.gz", ".bcf"])

proc error_msg*(msg: string, error_code = 1) =
    stderr.write_line fmt"Error {error_code}: {msg}".fgRed

proc quit_error*(msg: string, error_code = 1) =
    error_msg(msg, error_code)
    quit(error_code)

proc warning_msg*(msg: string) = 
    stderr.write_line fmt"Warning: {msg}".fgYellow

proc assert_file*(fname: string): string =
    # Checks that a file exists and throws an error if not
    if not fileExists(fname) and fname != "STDIN":
        quit_error fmt"{fname} does not exist or is not readable"
    return fname

proc check_file*(fname: string): bool =
    if not fileExists(fname):
        warning_msg fmt"{fname} does not exist or is not readable"
        return false
    return true

proc check_file_list*(files: seq[string]): bool = 
    var missing_file = false
    for f in files:
        var check = check_file(f)
        if check == false:
            missing_file = true
    return missing_file

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

# Position iterator

type Position* = ref object
    chrom*: string
    pos*: int

proc pos0*(p: Position): int =
    # Returns the 0-based coordinate
    # of a position
    return p.pos - 1

proc `$`*(p: Position): string =
    return fmt"<{p.chrom}:{p.pos}>"

iterator iter_pos*(pos_in: string): Position =
    # Parses VCF, BED, or single string and outputs
    # as a position list.
    # ***Always*** returns 1-based coordinates
    var v:VCF
    var chrom: string
    var pos: string
    var bed_offset: int

    if ":" in pos_in and ("/" in pos_in) == false:
        # Read in string argument
        (chrom, pos) = pos_in.split(":")
        yield Position(chrom: chrom,
                       pos: pos.parseInt())
    else:
        var (_, name, ext) = splitFile(pos_in.to_lower())
        # VCF
        if ext in ["bcf", "vcf", "vcf.gz"]:
            doAssert open(v, pos_in)
            for line in v:
                yield Position(chrom: $line.CHROM, 
                               pos: line.POS.int)
        
        # TODO: Bed handling needs to be reworked.
        # Bedfiles
        elif name.ends_with(".bed") or name.endswith(".bed.gz"):
            # If bed file, use 0-based coordinates
            bed_offset =
                if name.ends_with(".bed") or name.endswith(".bed.gz"):
                    0
                else:
                    1
    
        let stream: Stream =
            if pos_in[^3 .. ^1] == ".gz":
                newGZFileStream(pos_in)
            else:
                newFileStream(pos_in, fmRead)

        var n = 0
        for line in stream.lines:
            # TODO: This needs to be cleaned up
            n += 1
            var curr_line = line.strip(chars = {'\t', ':', ' '})
            curr_line = curr_line.strip(chars = {'\t', ':', ' '})
            try:
                (chrom, pos) = nre.split(curr_line, re"[\t: ]+")[0..1]
            except IndexError:
                # If it is the first line, assume it is a header
                if n == 1:
                    continue
                else:
                    warning_msg fmt"""Invalid line: {n} in "{pos_in}" > {line}"""
                    continue
            try:
                yield Position(chrom: chrom,
                               pos: pos.parseInt() - bed_offset)
            except ValueError:
                # If it is the first line, assume it is a header
                if n == 1:
                    continue
                else:
                    warning_msg fmt"""Invalid line: {n} in "{pos_in}" > {line}"""
                    continue


proc is_numeric(s: string): bool =
    return all(s, isDigit)

proc fix_chr(s: string): string = 
    if s.to_lower().startswith("chr") and s.len > 3:
        return s[3..s.len-1].to_lower()
    return s.to_lower()

const CHROM_VALS = {"x": 1, "y": 2, "m": 3}.toTable

proc genome_cmp*(x, y: Position): int =
    # Sort genomic coordinates
    # First numericaly and then by X,Y,M
    let x_chr = x.chrom.fix_chr()
    let y_chr = y.chrom.fix_chr()
    if x_chr.is_numeric() and y_chr.is_numeric():
        if x_chr == y_chr:
            if x.pos < y.pos: return -1
            elif x.pos == y.pos: return 0
            else: return 1
        elif x_chr.parseInt() < y_chr.parseInt(): return -1
        else: return 1
    elif x_chr.is_numeric() and y_chr.is_numeric() == false:
        return -1
    elif x_chr.is_numeric() == false and y_chr.is_numeric():
        return 1
    else:
        # Both are non-numeric
        if x_chr in CHROM_VALS and y_chr in CHROM_VALS:
            if CHROM_VALS[x_chr] < CHROM_VALS[y_chr]: return -1
            elif CHROM_VALS[x_chr] < CHROM_VALS[y_chr]:
                if x.pos < y.pos: return -1
                elif x.pos == y.pos: return 0
                else: return 1
            else:
                return 1
        else:
            if x_chr < y_chr: return -1
            elif x_chr == y_chr: return 0
            else: return 1


#====================#
# Headers and Output #
#====================#

proc output_header*(header: string, basename: bool, absolute: bool): string =
    var
        basename_str = ""
        absolute_str = ""
    if basename == true:
        basename_str = "basename"
    if absolute == true:
        absolute_str = "absolute"
    return [header, basename_str, absolute_str].filterIt(it.len > 0).join("\t")

proc get_absolute*(path: string): string =
    if symlinkExists(path) == true:
        return absolutePath(expandSymlink(path))
    return absolutePath(path)

proc output_w_fnames*(header: string, path: string, basename: bool, absolute: bool): string =
    # Attaches the basename or absolute if selected in a consistant manner
    var
        basename_str = ""
        absolute_str = ""
    if basename == true:
        basename_str = lastPathPart(path)
    if absolute == true:
        absolute_str = get_absolute(path)
    return [header, basename_str, absolute_str].filterIt(it.len > 0).join("\t")

#=========#
# Parsing #
#=========#

proc sci_parse_int*(s: string): int =
    # Parses ints with comma delimiters and scientific notation.
    if 'e' in s:
        var scientific_notation = s.split("e", maxsplit = 1)
        let coeff = scientific_notation[0].parseFloat()
        let exponent = scientific_notation[1].parseInt()
        return math.pow(coeff * 10.0, exponent.float).int
    return s.replace(",", "").parseInt()