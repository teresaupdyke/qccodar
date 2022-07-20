#!/usr/bin/python

"""
run_RenameReprocess
"""

import sys
import os
import re
import glob
import subprocess

def RenameReprocess(datadir, fns, new_str):
    """
    """
    # [outdir, fn] = os.path.split(ofn)
    # if not os.path.exists(outdir):
    #     os.makedirs(outdir)

    cwd = os.getcwd()
    # Rename the first file (faster than copying)
    os.chdir(datadir)
    if len(fns)>0:
        print "... renaming %s to %s " % (fns[0], new_str)
        os.rename(fns[0], new_str)
    else:
        print "... No files to rename"

    # merge copy with "ditto" from other reprocessed time-stamped subdirs
    # into the desired name (usually one month)
    if len(fns)>1:
        for fn in fns[1:]:
            # ditto [ <options> ] src [ ... src ] dst
            cmdstr = "/usr/bin/ditto -v "+fn+" "+new_str
            print "Merging with: "+cmdstr
            subprocess.call(cmdstr, shell=True)
            rmstr = "rm -rf "+fn
            print "Removing old with: "+rmstr
            subprocess.call(rmstr, shell=True)
            

    os.chdir(cwd)
    return

if __name__ == '__main__':
    # 
    datadir = sys.argv[1]
    pattern_str = sys.argv[2]
    new_str = sys.argv[3]
    
    # datadir = './'
    # datadir = '/Users/codar/Documents/reprocess_DUCK/'
    # pattern_str = 'Reprocess_151215*' 
    # pattern_str = 'Reprocess_*'
    # new_str = '2014_09'

    cwd = os.getcwd()
    os.chdir(datadir)
    fns = glob.glob(pattern_str)
    fns.sort()
    print datadir
    print fns
    print new_str
    os.chdir(cwd)

    print "Run RenameReprocess ... "
    try:
        RenameReprocess(datadir, fns, new_str)
    except:
        pass

