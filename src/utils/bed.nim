import streams
import zip/gzipfiles
import strutils
import strformat

type
  Region* = ref object
    chrom*: string
    start*: int
    stop*: int
    name*: string

proc bed_line_to_region(line: string): Region =
    # Adapted from https://github.com/brentp/hileup
    var
        cse = line.strip().split('\t', 5)
    if len(cse) < 3:
        stderr.write_line("[seq-collection] skipping bad bed line:", line.strip())
        return
    result = Region()
    result.start = parse_int(cse[1])
    result.stop = parse_int(cse[2])
    result.chrom = cse[0]
    if len(cse) > 3:
        result.name = cse[3]

proc len*(r: Region): int =
    return r.stop - r.start

proc `$`*(r: Region): string = 
    return fmt"{r.chrom}:{r.start}-{r.stop}"

iterator iter_bed*(bedfile: string): Region =
    let stream: Stream =
        if bedfile[^3 .. ^1] == ".gz":
            newGZFileStream(bedfile)
        else:
            newFileStream(bedfile, fmRead)
    
    for line in bedfile.lines:
        yield line.bed_line_to_region()