#!/usr/bin/env python

import sys
from pwem.viewers import Chimera


def main():
    fileName = sys.argv[1]

    # Clean filename in case it comes annotated like 1@xxxmrc:mrcs
    fileName = fileName.split(":")[0]

    # if contains @
    if "@" in fileName:
        fileName = fileName.split("@")[1]

    Chimera().runProgram(args=fileName)


if __name__ == '__main__':
    main()
