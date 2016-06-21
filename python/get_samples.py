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

    ./get_samples.py  input_dir output_file [size]

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    June 2016 
\n''')

    
def main():
    args = sys.argv[1:]
    if len(args) < 2:
        return usage()
    
    src = args[0]
    dst = args[1]
    sys.stdout.write('Scanning %s...\n' %src)

    file_list = []
    for root, dirs, files in os.walk(src):
        for f in files:
            file_list.append(os.path.join(root, f))

    if len(args) == 2:
        save_list_into_file(file_list, dst)
        return
    
    if len(args) == 3:
        size = args[2]
        sys.stdout.write('Size limit per list: %s\n' %size)
        groups = group_files_by_size(file_list, size)
        path, name = os.path.split(dst)
        dst_com = name.split('.')[0] 
        dst_ext = name.split('.')[1]
        n = 0 
        for file_list in groups:
            n = n+1 
            new_dst = os.path.join(path, '%s_%s-%s.%s' %
                                   (dst_com, size, n, dst_ext))
            save_list_into_file(file_list, new_dst)
        

def save_list_into_file(file_list, dst):
    nfiles = len(file_list)
    
    path, tail = os.path.split(dst)
    if path != '' and not os.access(path, os.F_OK) :
        sys.stdout.write('Creating dir %s ...\n'  % path)
        os.makedirs(path)
                
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
    sys.stdout.write('Saved as: %s\n' %dst)

    
def group_files_by_size(name_list, size_max='20G'):
    size_max = convert_size_from_str(size_max)
    groups = []
    group = []    
    size_sum = 0

    for name in name_list:
        size = os.path.getsize(name)
        if size_sum < size_max:
            group.append(name)
            size_sum += float(size)
        else:
            groups.append(group)
            group = []
            size_sum = 0
            group.append(name)            
            size_sum += float(size)

        if name == name_list[-1]:
            groups.append(group)

    return groups


def convert_size_from_str(size_str):
    c_1GB = 1024*1024*1024
    factor = eval(size_str.split('G')[0])
    return factor*c_1GB


if __name__ == '__main__':
    main()
