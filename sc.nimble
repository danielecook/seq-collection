# Package

version       = "0.0.1"
author        = "Daniel E. Cook"
description   = "seq-collection: Sequence data utilities"
license       = "MIT"

# Dependencies

requires "argparse >= 0.7.1", "hts >= 0.2.8", "colorize"

bin = @["sc"]
skipDirs = @["test"]

task test, "run tests":
  exec "nim c --lineDir:on --debuginfo -r tests/all"