from chimera import runCommand as rc
from chimera import replyobj
import os
import sys


def main():
    directory = sys.argv[1]
    ligand = sys.argv[2]
    os.chdir(directory)
    replyobj.status("Processing " + ligand)  # show what file we're working on
    rc("open " + ligand + ".pdb")
    rc("addcharge all")
    rc("write format mol2 0 " + ligand + ".mol2")
    rc("close all")
    rc("stop now")
    return


if __name__ == "__main__":
    main()
