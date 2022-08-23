#!/usr/bin/env python
## Last modified: Time-stamp: <2022-07-22 15:30:15 codar>
"""make_plist_config.py

Use this script to generate a config file for your specific use.
1. Reads qccodar_values dictionary from "tools/qccodar_default_configuration.py"
2. Specify the output location <outdir>
3. Run the script

Activate qccodar environment if have  one,
(base) % conda activate qccodar
(qccodar) % cd qccodar/tools
(qccodar) % python make_plist_config.py

For users, you can use the tool to generate a plist from scratch.
Then copy and edit the plist, rename it for specific conditions for
your site(s) and or time period.

For developers, you can edit the qccodar_default_configuration.py and
use this tool to generate a new default plist with any new values you
add to processing.  Then update the default qccodar.plist in
src/qccodar/config to be installed at software setup.

"""

import sys
import os
from pprint import PrettyPrinter
from plistlib import dump, FMT_XML

def write_plist(fullfn, qv):
    """ output dict of qccodar_values (qv) as plist """
    with open(ofn, 'wb') as fp:
        dump(qv,fp, fmt=FMT_XML, sort_keys=False) 

def get_config(name):
    """Usage Example >>>qccodar_values = get_config('qccodar_default_configuration.qccodar_values')"""
    components = name.split('.')
    mod = __import__(components[0])
    for comp in components[1:]:
        attr = getattr(mod, comp)
    return attr


indir = '.'
ifn = os.path.join(indir, 'qccodar_default_configuration.py')
cn = os.path.splitext(os.path.basename(ifn))[0]
# cn = 'qccodar_default_configuration'
print('getting values from ... %s' % ifn)

qv = get_config(cn+'.qccodar_values')
print('------ qccodar_values ------------')
pp = PrettyPrinter()
pp.pprint(qv)
print('-----------------------------------')

outdir = '.'
ofn = os.path.join(outdir, 'qccodar.plist')
print('saving plist to ... %s' % ofn)
write_plist(ofn, qv)
print('DONE')
