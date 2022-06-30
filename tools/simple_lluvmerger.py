#

import sys
import copy
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.codarutils import *
from qccodar.qcutils import *
from plistlib import load

configfile = '../src/qccodar/config/qccodar.plist'
plistfile = Path(configfile)
with open(plistfile, 'rb') as fp:
    qccodar_values = load(fp)

datadir = './2017_01'
pattern = 'IdealPattern'

#######
# test specifics
fn = 'RDLx_HATY_2017_01_01_0200.ruv'

ofn = run_LLUVMerger(datadir, fn, pattern, diag='4', **qccodar_values['merge'])
r = read_lluv_file(ofn)

for key in r._tables.keys():
    table = r._tables[key]
    if 'MRGS' in table['TableType']:
        filelist = table['data']['PATH']

shortpath = os.path.join(datadir, 'RadialShorts_qcd',pattern)
firstshortfn = os.path.join(shortpath, filelist[0])
rs = read_lluv_file(firstshortfn)

r = add_short_metadata(r, qccodar_values)

#############
# test general
fns = recursive_glob(os.path.join(datadir, 'RadialShorts_qcd', pattern), 'RDL'+ qccodar_values['merge']['shorts_minute_filter'] + '.ruv')
# run LLUVMerger for each
for fullfn in fns:
    print('###################')
    print('... input: %s' % fullfn)
    fn = os.path.basename(fullfn)
    ofn = do_merge(datadir, fn, pattern, qccodar_values)
    print('... output: %s' % ofn)

