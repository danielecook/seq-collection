import os
import hts
import math
import streams
import zip/gzipfiles
import strutils
import strformat
import colorize
import sequtils


proc error_msg*(msg: string, error_code = 1) =
    stderr.write_line fmt"Error {error_code}: {msg}".fgRed

proc quit_error*(msg: string, error_code = 1) =
    error_msg(msg, error_code)
    quit(error_code)

proc warning_msg*(msg: string) = 
    stderr.write_line fmt"Warning {msg}".fgYellow

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

proc `$`*(p: Position): string =
    return fmt"<{p.chrom}:{p.pos}>"

iterator iter_pos*(pos_in: string): Position =
    # Parses a VCF, Bedfile, or single string and outputs
    # as a position list.
    var v:VCF
    var chrom: string
    var pos: string
    var bed_offset: int

    if ":" in pos_in and ("/" in pos_in) == false:
        (chrom, pos) = pos_in.split(":")
        yield Position(chrom: chrom,
                       pos: pos.parseInt() - 1)
    else:
        var (_, name, ext) = splitFile(pos_in.to_lower())
        # VCF
        if ext in ["bcf", "vcf", "vcf.gz"]:
            doAssert open(v, pos_in)
            for line in v:
                yield Position(chrom: $line.CHROM, 
                            pos: line.POS.int - 1)
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
            n += 1
            #try:
            (chrom, pos) = line.strip(chars = {'\t', ':', ' '}).split({'\t', ':', ' '})
            yield Position(chrom: chrom,
                            pos: pos.parseInt() - bed_offset)
            #except IndexError:
            #    warning fmt"""Invalid line: {n} in "{pos_in}"; skipping"""
            
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
        let coeff = scientific_notation[0].parseInt()
        let exponent = scientific_notation[1].parseInt()
        return math.pow(coeff.float * 10.0, exponent.float).int
    return s.replace(",", "").parseInt()