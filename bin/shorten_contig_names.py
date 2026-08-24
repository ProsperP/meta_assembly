#!/usr/bin/env python3


import sys


shorten = False

with open(sys.argv[1], 'rt') as infh:
    for line in infh:
        if not line.startswith('>'):
            print(line.rstrip())
        else:
            if shorten:
                print("_".join(line.rstrip().split("_")[:4]))
            elif len(line) > 20 and len(line.split("_")) > 5:
                print("_".join(line.rstrip().split("_")[:4]))
                shorten = True
            else:
                print(line.rstrip())
