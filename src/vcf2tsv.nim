import hts
import json
import sequtils
import strutils
import strformat
import utils/helpers
import utils/vcf_header
import streams
from constants import ANN_header

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
  if s == int.low:
    result = newJNull()
  else:
    result = newJInt(s)
    
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

proc vcf2tsv*(vcf: string, region_list: seq[string], sample_set: string, info: string, format: string, long: bool, annotation: bool, pass: bool) =
    var v:VCF

    var info_order: seq[string]
    var format_order: seq[string]

    ## Set Fields
    var info_keep = filterIt(info.split({',', ' '}), it.len > 0)
    var format_keep = filterIt(format.split({',',' '}), it.len > 0)

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

    # Parse header
    var header = parse_vcf_header($v.header)

    if sample_set != "ALL":
        let samples_keep = filterIt(sample_set.split({',', ' '}), it.len > 0)
        set_samples(v, samples_keep)
    var samples = v.samples

    # Set items to get
    if info_keep[0] == "ALL":
        info_keep.setLen(0)
        for i in header.info:
            info_keep.add(i.ID.string)

    if format_keep[0] == "ALL":
        format_keep.setLen(0)
        for i in header.format:
            format_keep.add(i.ID.string)

    var header_out = info_keep & format_keep

    for rec in variants(v, region_list):
        if pass and rec.FILTER != "PASS":
            continue
    #     var info = rec.info
    #     var format = rec.format

    #     # Output fields
    #     var info_out: seq[string]

    #     echo header.info

    #     if info_keep.len == 0:
                #info_keep.add i.ID
        # echo info_keep
        # if output_all_info or info_keep.len >= 1:
        #     for info_field in rec.info.fields:
        #         if (info_field.name in info_order) == false:
        #             info_order.add(info_field.name)
                
                # if annotation and info_field.name == "ANN":
                #     #if info_field.name == "ANN":
                #         # discard info.get(info_field.name, field_string)
                #         # var ann_record_set: seq[string]
                #         # for ANN in field_string.split(","):
                #         #     var ann_record = newJObject()
                #         #     var ann_split = ANN.split("|")
                #         #     for ann_col in 0..<ANN_header.len:
                #         #         ann_record.add(ANN_header[ann_col], newJString(ann_split[ann_col]))
                #         #     ann_record_set.add(ann_record)
                #         # j_info.add("ANN", ann_record_set)
                #     elif info_field.name == "BCSQ":
                #         discard
        
        ## Fetch FORMAT fields
        # var j_format = newJObject()
        # if output_all_format or format_keep.len >= 1:
        #     for format_field in format.fields:
        #         if output_all_format or format_field.name in format_keep and format_field.name != "GT":
        #             if format_field.vtype == BCF_TYPE.FLOAT:
        #                 discard format.get(format_field.name, field_float)
        #                 j_format.add(format_field.name, out_fmt(field_float, format_field, zip, samples))
        #             elif format_field.vtype in [BCF_TYPE.INT8, BCF_TYPE.INT16, BCF_TYPE.INT32]:
        #                 discard format.get(format_field.name, field_int)
        #                 j_format.add(format_field.name, out_fmt(field_int, format_field, zip, samples))
        #     if "GT" in format_keep:
        #         var i_gt: seq[int]
        #         var gt_set: seq[seq[int]]
        #         for g in format.genotypes(field_int):
        #             for a in g:
        #                 i_gt.add (if a.value() >= 0: a.value() else: int.low)
        #             gt_set.add(i_gt)
        #             i_gt.setLen 0
        #         j_format.add(gt.name, out_fmt(gt_set, gt, zip, samples))
        #     if "SGT" in format_keep:
        #         j_format.add(sgt.name,
        #                      out_fmt(format.genotypes(field_int).mapIt($it),
        #                              sgt,
        #                              zip,
        #                              samples))
        #     if "TGT" in format_keep:
        #         var alleles = @[rec.REF].concat(rec.ALT)
        #         var tgt_set: seq[string]
        #         for g in format.genotypes(field_int):
        #             var gt: string
        #             for a in g:
        #                 gt = gt & (if a.value() >= 0: alleles[a.value()] else: ".") & (if a.phased: '|' else: '/')
        #             gt.set_len(gt.len - 1)
        #             tgt_set.add(gt)
        #         j_format.add(tgt.name, out_fmt(tgt_set, tgt, zip, samples))
                
        var base_rec = @[$rec.CHROM,
                         $rec.POS,
                         $rec.ID,
                         $rec.REF,
                         $rec.ALT,
                         $rec.QUAL,
                         $rec.FILTER]
        echo base_rec
        # if info_keep.len > 0:
        #     json_out.add("INFO", j_info)
        # if format_keep.len > 0:
        #     json_out.add("FORMAT", j_format)
    
        # if pretty:
        #     stdout.write $json_out.pretty()
        # else:
        #     stdout.write $json_out
        # if array:
        #     stdout.write ",\n"
        # else:
        #     stdout.write "\n"