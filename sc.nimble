# Package

version       = "0.0.2"
author        = "Daniel E. Cook"
description   = "seq-collection: Sequence data utilities"
license       = "MIT"

# Dependencies

requires "colorize", "zip >= 0.2.1"
requires "https://github.com/danielecook/BitVector"
requires "hts >= 0.3.8"
requires "argparse >= 0.10.0"
requires "alea >= 0.1.4"
requires "regex"

bin = @["sc"]
skipDirs = @["test"]

task test, "run tests":
  exec "bash ./scripts/functional-tests.sh"
  #exec "nim c --threads:on -d:release --lineDir:on --debuginfo -r tests/all"
