import regex
import strutils
import strformat
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
    var m: RegexMatch
    if regex.match(line, regex.re".*ID=(?P<ID>[^,]+),.*length=(?P<length>[0-9]+)", m, 2):
        return contig(ID: m.group("ID", line)[0], length: m.group("length", line)[0].parseInt())


proc extract_info(line: string): info =
    var m: RegexMatch
    var r: info
    for field in format_fields:
        if match(line, regex.re($(fmt(".*{field}=(?P<val>[^,>]+).*"))), m, 1):
            case field
            of "ID": r.ID = m.group("val", line)[0]
            of "Number": r.Number = m.group("val", line)[0]
            of "Type": r.Type = m.group("val", line)[0]
            of "Description": r.Description = m.group("val", line)[0]
    return r

proc extract_format(line: string): format =
    var m: RegexMatch
    var r: format
    for field in format_fields:
        if match(line, regex.re($(fmt(".*{field}=(?P<val>[^,>]+).*"))), m, 1):
            case field
            of "ID": r.ID = m.group("val", line)[0]
            of "Number": r.Number = m.group("val", line)[0]
            of "Type": r.Type = m.group("val", line)[0]
            of "Description": r.Description = m.group("val", line)[0]
    return r

proc extract_filter(line: string): filter =
    var m: RegexMatch
    var r: filter
    for field in filter_fields:
        if match(line, re($(fmt(".*{field}=(?P<val>[^,>]+).*"))), m, 1):
            case field
            of "ID": r.ID = m.group("val", line)[0]
            of "Description": r.Description = m.group("val", line)[0]
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