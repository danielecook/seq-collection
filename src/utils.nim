import os
import strformat
import colorize

proc quit_error*(msg: string, error_code = 1) =
    stderr.write_line fmt"Error {error_code}".bgWhite.fgRed & fmt": {msg}".fgRed
    quit(error_code)

proc print_error*(msg: string) =
    stderr.write_line "\nError".bgWhite.fgRed & fmt": {msg}\n".fgRed

proc check_file*(fname: string): bool =
    if not fileExists(fname):
        print_error(fmt"{fname} does not exist or is not readable")
        return false
    return true