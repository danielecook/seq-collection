import hts
import json
import sequtils
import strutils
import strformat
import utils/helpers
import utils/vcf_header
import streams
import math
from constants import ANN_header

proc out_fmt[T](record: T, fmt_field: FormatField, zip: bool, samples: seq[string]): JsonNode =
    # For fascilitating formatting FORMAT fields
    var 
        idx_start: int
        idx_end: int
        rec_out: JsonNode
        fmt_zip = newJObject()
        fmt_arr = newJArray()
    for idx in 0..<samples.len:
        if fmt_field.n_per_sample == 1:
            rec_out = %* record[idx]
        else:
            idx_start = idx * fmt_field.n_per_sample
            idx_end = idx * fmt_field.n_per_sample + fmt_field.n_per_sample - 1
            var rec_arr = newJArray()
            for i in idx_start..idx_end:
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


proc `$`*(f:FileStream): string {.inline.} =
    echo "<filestream>"

proc calc_tajima*(vcf: string, region_list: seq[string]) =
    var v:VCF

    var
        n = v.samples.len
        a1, a2, b1, b2, c1, c2, e1, e2: float

    for i in 1..n:
        a1 += 1.0 / i.float
    for i in 1..n:
        a2 += 1.0 / pow(i.float, 2.float)
    b1 = (n.float + 1.0) / (3.0 * (n.float - 1.0))
    b2 = (2.0 * (pow(n.float, 2.0) + n.float + 3.0)) /
         (9.0 * n.float * (n.float - 1.0))
    c1 = b1 - (1.0 / a1)
    c2 = b2 - ((n.float + 2.0) / (a1 * n.float)) + (a2 / (pow(a1, 2.0)))
    e1 = c1 / a1
    e2 = c2 / ((a1*a1) + a2)


