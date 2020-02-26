import os
import hts
import strutils
import strformat
import colorize
import sequtils

proc quit_error*(msg: string, error_code = 1) =
    stderr.write_line fmt"Error {error_code}".bgWhite.fgRed & fmt": {msg}".fgRed
    quit(error_code)

proc print_error*(msg: string) =
    stderr.write_line "Error".bgWhite.fgRed & fmt": {msg}".fgRed

proc check_file*(fname: string): bool =
    if not fileExists(fname):
        print_error(fmt"{fname} does not exist or is not readable")
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

#================#
# Bedfile reader #
#================#
type region_t* = ref object
    chrom: string
    start: uint32
    stop: uint32
    name: string

proc bed_line_to_region*(line: string): region_t =
    var
        cse = line.strip().split('\t', 5)
    
    if len(cse) < 3:
        stderr.write_line("[seq-collection] skipping bad bed line:", line.strip())
        return nil
    var
        s = strutils.parse_int(cse[1])
        e = strutils.parse_int(cse[2])
        reg = region_t(chrom: cse[0], start: uint32(s), stop: uint32(e))
    doAssert s <= e, "[seq-collection] ERROR: start > end in bed line:" & line
    if len(cse) > 3:
        reg.name = cse[3]
    return reg

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