import hts
import json
import sequtils
import strutils
import utils/helpers
import streams
import math
from constants import ANN_header, BCSQ_header


proc `%`(s: string): JsonNode =
  # Overload JsonNode to string
  # Important for converting
  # '.' to null values
  if s == ".":
    result = newJNull()
  else:
    result = newJString(s)
    
proc `%`(s: int): JsonNode =
  # Overload JsonNode from int
  # Important for converting
  # '.' to null values
  if s == int32.low:
    result = newJNull()
  else:
    result = newJInt(s)

  proc `%`(s: float32): JsonNode =
    # Overload JsonNode from float32
    # Important for converting
    # '.' to null values
    if s.classify == fcNaN:
      result = newJNull()
    else:
      result = newJFloat(s)      
    
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

proc to_json*(vcf: string, region_list: seq[string], sample_set: string, info: string, format: string, zip: bool, annotation: bool, pretty: bool, array: bool, pass: bool) =
    var v:VCF

    ## Format Fields
    let info_keep = filterIt(info.split({',', ' '}), it.len > 0)
    let format_keep = filterIt(format.split({',',' '}), it.len > 0)
    let output_all_format = ("ALL" in format_keep)

    # Custom Format Fields
    var gt: FormatField
    gt.name = "GT"
    gt.n_per_sample = 1
    var sgt: FormatField
    sgt.name = "SGT"
    sgt.n_per_sample = 1
    var tgt: FormatField
    tgt.name = "TGT"
    tgt.n_per_sample = 1

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
        if pass and rec.FILTER != "PASS":
            continue
        var info = rec.info
        var format = rec.format
        # Fetch INFO Fields
        var j_info = newJObject()
        let output_all_info = ("ALL" in info_keep or annotation)
        if output_all_info or info_keep.len >= 1:
            for info_field in rec.info.fields:
                if annotation and info_field.name == "ANN":
                    discard info.get(info_field.name, field_string)
                    var ann_record_set = newJArray()
                    for ANN in field_string.split(","):
                        var ann_record = newJObject()
                        var ann_split = ANN.split("|")
                        for ann_col in 0..<ANN_header.len:
                            ann_record.add(ANN_header[ann_col], newJString(ann_split[ann_col]))
                        ann_record_set.add(ann_record)
                    j_info.add("ANN", ann_record_set)
                elif annotation and info_field.name == "BCSQ":
                    discard info.get(info_field.name, field_string)
                    var ann_record_set = newJArray()
                    for ANN in field_string.split(","):
                        var ann_record = newJObject()
                        var ann_split = ANN.split("|")
                        for ann_col in 0..<ann_split.len:
                            ann_record.add(BCSQ_header[ann_col], newJString(ann_split[ann_col]))
                        ann_record_set.add(ann_record)
                    j_info.add("BCSQ", ann_record_set)
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
        if output_all_format or format_keep.len >= 1:
            for format_field in format.fields:
                if output_all_format or format_field.name in format_keep and format_field.name != "GT":
                    if format_field.vtype == BCF_TYPE.FLOAT:
                        discard format.get(format_field.name, field_float)
                        j_format.add(format_field.name, out_fmt(field_float, format_field, zip, samples))
                    elif format_field.vtype in [BCF_TYPE.INT8, BCF_TYPE.INT16, BCF_TYPE.INT32]:
                        discard format.get(format_field.name, field_int)
                        j_format.add(format_field.name, out_fmt(field_int, format_field, zip, samples))
            if "GT" in format_keep:
                var i_gt: seq[int]
                var gt_set: seq[seq[int]]
                for g in format.genotypes(field_int):
                    for a in g:
                        i_gt.add (if a.value() >= 0: a.value() else: int.low)
                    gt_set.add(i_gt)
                    i_gt.setLen 0
                j_format.add(gt.name, out_fmt(gt_set, gt, zip, samples))
            if "SGT" in format_keep:
                j_format.add(sgt.name,
                             out_fmt(format.genotypes(field_int).mapIt($it),
                                     sgt,
                                     zip,
                                     samples))
            if "TGT" in format_keep:
                var alleles = @[rec.REF].concat(rec.ALT)
                var tgt_set: seq[string]
                for g in format.genotypes(field_int):
                    var gt: string
                    for a in g:
                        gt = gt & (if a.value() >= 0: alleles[a.value()] else: ".") & (if a.phased: '|' else: '/')
                    gt.set_len(gt.len - 1)
                    tgt_set.add(gt)
                j_format.add(tgt.name, out_fmt(tgt_set, tgt, zip, samples))
                
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