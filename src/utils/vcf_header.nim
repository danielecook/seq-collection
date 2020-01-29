import re
import strutils
import strformat
const contig_fields = @["ID", "length"]
const info_fields = @["ID", "Number", "Type", "Description"]
const format_fields = @["ID", "Number", "Type", "Description"]
const filter_fields = @["ID", "Description"]

type
    contig* = object
        ID*: string
        length*: int
    
    info* = object
        ID*: string
        Number*: string
        Type*: string
        Description*: string

    format* = object
        ID*: string
        Number*: string
        Type*: string
        Description*: string

    filter* = object
        ID*: string
        Description*: string

    header* = object
        contig*: seq[contig]
        info*: seq[info]
        format*: seq[format]
        filter*: seq[filter]


proc extract_contig(line: string): contig =
    var matches: array[2, string]
    if match(line, re".*ID=([^,]+),.*length=([0-9]+)", matches, 2):
        return contig(ID: matches[0], length: matches[1].parseInt())

proc extract_info(line: string): info =
    var m: array[1, string]
    var r: info
    for field in format_fields:
        if match(line, re($(fmt(".*{field}=([^,>]+).*"))), m, 1):
            #r[field] = m[0]
            case field
            of "ID": r.ID = m[0]
            of "Number": r.Number = m[0]
            of "Type": r.Type = m[0]
            of "Description": r.Description = m[0] 
    return r

proc extract_format(line: string): format =
    var m: array[1, string]
    var r: format
    for field in format_fields:
        if match(line, re($(fmt(".*{field}=([^,>]+).*"))), m, 1):
            case field
            of "ID": r.ID = m[0]
            of "Number": r.Number = m[0]
            of "Type": r.Type = m[0]
            of "Description": r.Description = m[0]
    return r

proc extract_filter(line: string): filter =
    var m: array[1, string]
    var r: filter
    for field in filter_fields:
        if match(line, re($(fmt(".*{field}=([^,>]+).*"))), m, 1):
            case field
            of "ID": r.ID = m[0]
            of "Description": r.Description = m[0] 
    return r

    

proc parse_vcf_header*(h_in: string): header =
    var h: header
    for line in h_in.splitLines():
        if line.startsWith("##contig"):
            h.contig.add extract_contig(line)
        elif line.startsWith("##INFO"):
            h.info.add extract_info(line)
        elif line.startsWith("##FORMAT"):
            h.format.add extract_format(line)
        elif line.startsWith("##FILTER"):
            h.filter.add extract_filter(line)
    return h