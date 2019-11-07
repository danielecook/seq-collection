# Package

version       = "0.0.1"
author        = "Daniel E. Cook"
description   = "seq-collection: Sequence data utilities"
license       = "MIT"

# Dependencies

requires "argparse >= 0.9.0", "hts >= 0.2.8", "colorize", "zip >= 0.2.1"
requires "https://github.com/danielecook/BitVector#b8cc21271c90cca96ed31f5d5383711dc96a8d3f"

bin = @["sc"]
skipDirs = @["test"]

task test, "run tests":
  exec "bash ./scripts/functional-tests.sh"
  #exec "nim c --threads:on -d:release --lineDir:on --debuginfo -r tests/all"