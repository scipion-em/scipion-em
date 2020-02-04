#!/usr/bin/env python

import sys
from pwem.viewers import Chimera

def main():

    Chimera().runProgram(args=sys.argv[1])

if __name__ == '__main__':
    main()
