#!/usr/bin/env python
## Last modified: Time-stamp: <2022-07-22 15:27:15 codar>
""" Runs the merge step (RadialShorts_qcd to Radials_qcd) for multiple files in <datadir>
    This tool is helpful if want to change some parameters for the merge

Usage: 
    run_do_merge.py <datadir> <patterntype> [configfile]

Options:
    <datadir>     Data directory to process: e.g. /Codar/SeaSonde/Data, or 
                                                  /path/to/files/reprocess_SITE/YYYY_MM
    <patterntype> Pattern tpype: IdealPattern or MeasPattern
    [configfile]  Path to qccodar configuration file:  e.g. /path/to/config/qccodar.plist or
                                                  /path/to/config/qccodar_SITE_YYYY_MM.plist

If no configfile is specified, default settings will be used (specific for long-range systems).

"""

import sys
import os

from qccodar.app import load_configs
from qccodar.codarutils import *
from qccodar.qcutils import *


def manual(datadir, pattern, configfile):
    """ Manual mode runs do_merge for all files in datadir """
    # load config
    qccodar_values = load_configs(configfile)
    print(qccodar_values)

    # get file list of RadialShorts
    # depending on system and desired time span for merge, change the target time for file search in configfile
    fns = recursive_glob(os.path.join(datadir, 'RadialShorts_qcd', pattern),
                         'RDL'+ qccodar_values['merge']['shorts_minute_filter'] + '.ruv')

    print('qccodar (manual) -- merge RadialShorts_qcd to Radials_qcd: ...')

    # run LLUVMerger for each
    for fullfn in fns:
        print('###################')
        print('... input: %s' % fullfn)
        fn = os.path.basename(fullfn)
        ofn = do_merge(datadir, fn, pattern, qccodar_values)
        print('... output: %s' % ofn)


if __name__ == '__main__':

    if len(sys.argv)==4:
        datadir = sys.argv[1]
        pattern = sys.argv[2]
        configfile = sys.argv[3]
    elif len(sys.argv)==3:
        datadir = sys.argv[1]
        pattern = sys.argv[2]
        configfile = ''

    print(datadir)
    print(pattern)
    print(configfile)
    manual(datadir, pattern, configfile)
