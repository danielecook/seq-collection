# Package

version       = "0.0.1"
author        = "Daniel E. Cook"
description   = "seq-collection: Sequence data utilities"
license       = "MIT"

# Dependencies

requires "argparse >= 0.7.1", "hts >= 0.2.8", "colorize", "zip >= 0.2.1"
requires "https://github.com/danielecook/BitVector#b8cc21271c90cca96ed31f5d5383711dc96a8d3f"

bin = @["sc"]
skipDirs = @["test"]

task test, "run tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"