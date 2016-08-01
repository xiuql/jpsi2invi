#!/usr/bin/env python
"""
Check PBS jobs 
"""

__author__ = "SHI Xin <shixin@ihep.ac.cn>"
__copyright__ = "Copyright (c) SHI Xin"
__created__ = "[2016-06-02 Thu 09:42]" 

import sys
import os
from hurry.filesize import size 
from tools import BossLogFile, EventsLogFile 


def usage():
    sys.stdout.write('''
NAME
    chk_pbsjobs.py 

SYNOPSIS

    ./chk_pbsjobs.py  input_dir num_of_files

AUTHOR 
    SHI Xin <shixin@ihep.ac.cn> 

DATE
    July 2016 
\n''')

    
def main():
    args = sys.argv[1:]
    if len(args) < 2:
        return usage()
    
    src = args[0]
    num = int(args[1])
    jobs_created = set(range(1, num+1))

    log = src 
    logdir = src.split('/')[-1]
    if logdir == 'data':
        logfiletype = 'BossLogFile'
    elif logdir == 'events':
        logfiletype = 'EventsLogFile'
    else:
        raise NameError(logdir)
    
    log = log.replace(logdir, 'log/%s' %logdir) 
    
    sys.stdout.write('Scanning %s...\n' %src)

    file_list = []
    total_size = 0 
    for root, dirs, files in os.walk(src):
        for f in files:
            file_list.append(int(f.split('-')[-1].split('.')[0]))
            total_size = total_size + os.path.getsize(os.path.join(root,f))

    sys.stdout.write('Found %s files, with total size %s.\n' %(
        len(file_list), size(total_size)))
    
    if len(file_list) < num:
        jobs_missing = jobs_created.difference(file_list)
        jobs_missing = [str(li) for li in jobs_missing]
        sys.stdout.write('Missing jobs are: %s\n' % ','.join(jobs_missing))
        
    sys.stdout.write('Checking log files...\n')
    jobs_not_terminated = []
    for root, dirs, files in os.walk(log):
        for f in files:
            if logfiletype == 'BossLogFile': 
                l = BossLogFile( os.path.join(root, f) )
            elif logfiletype == 'EventsLogFile':
                l = EventsLogFile( os.path.join(root, f) )
            else:
                raise NameError(logfiletype) 

            if not l.terminated:
                sys.stdout.write('%s ... Not OK.\n' %f)
                job = f.split('-')[-1]
                jobs_not_terminated.append(job)
            else:
                sys.stdout.write('%s ... OK.\n' %f)

    sys.stdout.write('Non-terminated jobs are (%s): %s\n' % (
        len(jobs_not_terminated), ','.join(jobs_not_terminated)))
    
                
if __name__ == '__main__':
    main()
