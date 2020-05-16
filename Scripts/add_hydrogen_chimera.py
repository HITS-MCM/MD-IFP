from chimera import runCommand as rc
from chimera import replyobj
import os
import sys


def main():
    directory = sys.argv[1]
    os.chdir(directory)
    file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]
    if "merged.pdb" in file_names:
        print("merged.pdb already exists. You probably already converted everything. " +
              "Doing so again would introduce mistakes. Stopping.")
        rc("stop now")

    for fn in file_names:
        replyobj.status("Processing " + fn)  # show what file we're working on
        rc("open " + fn)
        rc("addh")
        rc("write 0 " + fn)
        rc("close all")
    rc("stop now")
    return


if __name__ == "__main__":
    main()
