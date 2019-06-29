{.experimental.}

import hts
import json
import stats
import strutils
import algorithm
import strformat
import progress
from ../fq_meta import extract_read_info

type Site* = object
  ref_allele*: char
  alt_allele*: char
  chrom*: string
  position*: int

type Stat4* = object
  dp*: RunningStat
  gtdp*: RunningStat # depth of genotyped sites
  un*: RunningStat
  ab*: RunningStat


type count* = object
  nref: uint32
  nalt: uint32
  nother: uint32

type pair = tuple[a:string, b:string, rel:float64]

proc `%`*(p:pair): JsonNode =
  return %*{"a":p.a, "b":p.b, "rel":p.rel}

proc ab*(c:count): float {.inline.} =
  return c.nalt.float / (c.nalt + c.nref).float

template proportion_other(c:count): float =
  if c.nother == 0: 0'f else: c.nother.float / (c.nother + c.nref + c.nalt).float

proc alts*(c:count, min_depth:int): int8 {.inline.} =
  ## These numbers modified as we are only interested
  ## in contaminated "ref-like" sites and homozygous alts.
  ## give an estimate of number of alts from counts of ref and alt
  ## AB < 0.30 is called as hom-ref
  ## AB > 0.98 is hom-alt
  ## 0.30 <= AB <= 0.99 is het
  ##
  ## 0 = HOM REF
  ## 1 = HET
  ## 2 = HOM ALT
  ## 3 = CONTAMINATED REF
  ##
  if c.proportion_other > 0.04: return -1
  if int(c.nref + c.nalt) < min_depth:
    return -1
  if c.nalt == 0:
    return 0

  # ab befow is the allele freq of the alt allele
  var ab = c.ab
  # "Contaminated" sites or seq. errors on reference
  if ab > 0.0 and ab < 0.30:
    return 3
  elif ab > 0.98:
    return 2 # HOM ALT

  # Anything else is a 'het' and is ignored
  return 1

proc count_alleles(b:Bam, site:Site): count {.inline.} =
    # Incr. mapping quality for site from 10 to 30
  for aln in b.query(site.chrom, site.position, site.position + 1):
    if aln.mapping_quality < 30: continue
    var off = aln.start
    var qoff = 0
    var roff_only = 0
    for event in aln.cigar:
      var cons = event.consumes
      if cons.query:
        qoff += event.len
      if cons.reference:
        off += event.len
        if not cons.query:
          roff_only += event.len
      if off <= site.position: continue

      if not (cons.query and cons.reference): continue

      var over = off - site.position - roff_only
      if over > qoff: break
      if over < 0: continue
      doAssert qoff - over >= 0
      var base = aln.base_at(qoff - over)
      if base == site.ref_allele:
        result.nref += 1
      elif base == site.alt_allele:
        result.nalt += 1
      else:
        #stderr.write_line $event, " -> over:", over, " -> ", site, " -> ", aln.tostring
        result.nother += 1


proc get_alts(bam:Bam, sites:seq[Site], nalts: ptr seq[int8], alt_depth: ptr seq[int], alt_alleles: ptr seq[int], depth: ptr seq[int], stat: ptr Stat4, min_depth:int=6, progress: ptr ProgressBar): bool =
  ## count alternate alleles in a single bam at each site.
  for i, site in sites:
    if i mod 1000 == 0:
        # Update progress bar
        if progress.isComplete() == false:
            progress.tick(1000)
    var c = bam.count_alleles(site)
    alt_depth[i] = c.nalt.int
    depth[i] = int(c.nref + c.nalt + c.nother)
    stat.dp.push(int(c.nref + c.nalt))
    if c.nref > 0'u32 or c.nalt > 0'u32 or c.nother > 0'u32:
      stat.un.push(c.nother.float64 / float64(c.nref + c.nalt + c.nother))
    if c.nref.float > min_depth / 2 and c.nalt.float > min_depth / 2:
      stat.ab.push(c.ab)

    alt_alleles[i] = c.nalt.int
    nalts[][i] = c.alts(min_depth)
    if nalts[][i] != -1:
      stat.gtdp.push(int(c.nref + c.nalt))
    


proc get_bam_alts*(path:string, fai:string, sites:seq[Site], nalts: ptr seq[int8], alt_depth: ptr seq[int], alt_alleles: ptr seq[int], depth: ptr seq[int], stat: ptr Stat4, min_depth:int=6, progress: ptr ProgressBar): bool =
  var bam: Bam
  if not open(bam, path, fai=fai):
    quit "couldn't open :" & $path
  bam.load_index(path & ".bai")
  discard bam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 8191 - SAM_QUAL.int - SAM_QNAME.int - SAM_RNAME.int)
  result = bam.get_alts(sites, nalts, alt_depth, alt_alleles, depth, stat, min_depth, progress)
  bam.close()

proc siteOrder(a:Site, b:Site): int =
  if a.chrom == b.chrom:
    return cmp(a.position, b.position)
  return cmp(a.chrom, b.chrom)

proc checkSiteRef(s:Site, fai:var Fai) =
  var fa_allele = fai.get(s.chrom, s.position, s.position)[0].toUpperAscii
  if s.ref_allele != fa_allele:
    quit "reference base from sites file:" & s.ref_allele & " does not match that from reference: " & fa_allele

{.push checks: off, optimization:speed.}
proc toSite(toks: seq[string]): Site =
  result = Site()
  result.chrom = toks[0]
  result.position = parseInt(toks[1]) - 1
  result.ref_allele = toks[3][0]
  result.alt_allele = toks[4][0]


proc readSites*(path: string, fai: var Fai): seq[Site] =
  result = newSeqOfCap[Site](8192)
  var kstr = kstring_t(l:0, m:0, s:nil)
  var hf = hts_open(path.cstring, "r")

  while hts_getline(hf, cint(10), kstr.addr) > 0:
    var line  = $kstr.s
    if line[0] == '#': continue
    var sep = '\t'
    # handle ":" or tab. with ":", there is no id field.
    if line.count(sep) == 0:
      sep = ':'
    var toks = line.strip().split(sep)
    if sep == ':':
      toks.insert(".", 2)

    result.add(toSite(toks))
  if len(result) > 65535:
    stderr.write_line "warning:cant use more than 65535 sites"
  sort(result, siteOrder)
  
  # check reference after sorting so we get much faster access.
  if fai != nil:
      for i, r in result:
          if i mod 10000 == 0 and i > 0:
              stderr.write_line "[seq-collection] checked reference for " & $i & " sites"

          checkSiteRef(r, fai)
  fai = nil
  return result

proc bam_like(path:string): bool {.inline.} =
    return path.endsWith(".bam") or path.endsWith(".cram")

proc get_sample_names*(path: string): seq[string] =
  if path.bam_like:
    var bam: Bam
    open(bam, path)
    var txt = newString(bam.hdr.hdr.l_text)
    copyMem(txt[0].addr, bam.hdr.hdr.text, txt.len)
    for line in txt.split("\n"):
      if line.startsWith("@RG") and "\tSM:" in line:
        var t = line.split("\tSM:")[1].split("\t")[0].strip()
        return @[t]

proc get_sample_flowcells*(path: string): string =
    if path.bam_like:
        var bam: Bam
        open(bam, path)
        for line in bam:
            return extract_read_info(line.qname)[4]