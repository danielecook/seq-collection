import zip/zlib
import streams

# Simple wrapper for exposing GzFile as a Stream.
# Based on https://stackoverflow.com/a/33104286

type
  GZipStream* = object of StreamObj
    f: GzFile

  GzipStreamRef* = ref GZipStream

proc fsClose(s: Stream) =
  discard gzclose(GZipStreamRef(s).f)

proc fsReadData(s: Stream, buffer: pointer, bufLen: int): int =
  return gzread(GZipStreamRef(s).f, buffer, bufLen)

proc fsWriteData(s: Stream, buffer: pointer, bufLen: int) {.gcsafe.} =
  discard gzwrite(GZipStreamRef(s).f, buffer, bufLen)

proc fsAtEnd(s: Stream): bool =
  return gzeof(GZipStreamRef(s).f) != 0

proc newGZipStream*(f: GzFile): GZipStreamRef =
  new result
  result.f = f
  result.closeImpl = fsClose
  result.readDataImpl = fsReadData
  result.writeDataImpl = fsWriteData
  result.atEndImpl = fsAtEnd

proc newGZipStream*(filename: cstring, mode: cstring): GZipStreamRef =
  var gz = gzopen(filename, mode)
  if gz != nil: return newGZipStream(gz)

proc openGzRead*(fname: string): GZipStreamRef =
  let s = newGZipStream(c_string(fname), "r")
  if s == nil:
    stderr.writeLine("Error opening ", fname, " for reading; exiting!")
    quit(2)
  return s

proc openGzWrite*(fname: string): GZipStreamRef =
  let s = newGZipStream(c_string(fname), "w")
  if s == nil:
    stderr.writeLine("Error opening ", fname, " for writing; exiting!")
    quit(2)
  return s