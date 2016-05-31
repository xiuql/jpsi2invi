#!/usr/bin/env python
"""
Get samples list. 
"""

__author__ = "SHI Xin <shixin@ihep.ac.cn>"
__copyright__ = "Copyright (c) SHI Xin"
__created__ = "[2016-05-30 Mon 09:19]"

import sys
import os


def usage():
    sys.stdout.write('''
NAME
    get_samples.py 

SYNOPSIS

    ./get_samples.py  input_dir output_file

    ./get_samples.py  /bes3fs/offline/data/664p03/psip/dst ../samples/data_664p03_psip.txt

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    May 2016 
\n''')

    
def main():
    args = sys.argv[1:]
    if len(args) != 2:
        return usage()
    
    src = args[0]
    dst = args[1]
    print src, dst

    file_list = []
    for root, dirs, files in os.walk(src):
        for f in files:
            file_list.append(os.path.join(root, f))

    nfiles = len(file_list)

    fo = open(dst, 'w')
    fo.write('EventCnvSvc.digiRootInputFile = {\n')

    n = 0
    for f in file_list:
        n = n+1
        if n<nfiles:
            fo.write('"%s",\n' % f)
        else:
            fo.write('"%s"\n};\n' % f)

    fo.close()
    

    
if __name__ == '__main__':
    main()
