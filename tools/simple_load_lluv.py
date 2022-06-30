#

import sys
import copy
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.codarutils import *
from qccodar.qcutils import *
from qccodar.app import load_configs

configfile = '../src/qccodar/config/qccodar.plist'
qccodar_values = load_configs(configfile)

# test default is same as distro src plist
# qv1 = load_configs(configfile)
# qv2 = load_configs('') # default if not specified
# qv1 == qv2

ifn = './2017_01/RadialMetric/IdealPattern/RDLv_HATY_2017_01_01_0000.ruv'

rs = read_lluv_file(ifn)

# read in other radial metric data to use in averaging over time
ixfns = find_files_to_concatenate(ifn, **qccodar_values['metric_concatenation'])

for xfn in ixfns:
    if xfn == ifn:
        continue

    rs1 = read_lluv_file(xfn)
    if len(rs.data.shape) == len(rs1.data.shape) == 2:
        if (rs.data.shape[1] == rs1.data.shape[1]):
            if all(rs.data.columns == rs1.data.columns):
                # if same number and order of columns as rs, then append the data rs
                if debug:
                    print('... ... include: %s' % xfn)
                rs.data = pd.concat([rs.data,rs1.data],ignore_index=True)

# (1) do threshold qc on radialmetric
rmqc = threshold_qc_all(rs, qccodar_values)
# In [3]: %timeit rmqc = threshold_qc_all(rs, qccodar_values)
# 24.4 ms ± 627 µs per loop (mean ± std. dev. of 7 runs, 10 loops each)


# (2) do weighted averaging of good, this is done when the radial short is created
rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
# In [4]: %timeit rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
# 48.3 s ± 598 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)


# (3) require a minimum numpoints used in to form cell average
rsx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])
# In [7]: %timeit rsxx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])
# 1.63 ms ± 73.1 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)


# test remove deepcopy in loop weighted_velocities
from test_change_weighted_velocities import *

rsx1 = generate_radialshort_SH(rmqc, **qccodar_values['weighted_shorts'])
# ***** some speedup (8 sec faster and not equal at Radial level but DataFrame is ******
# In [9]: %timeit rsx1 = generate_radialshort_SH(rmqc, **qccodar_values['weighted_shorts'])
# 40.7 s ± 138 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)

rsx1 = threshold_rsd_numpoints(rsx1, **qccodar_values['qc_radialshort_velocity_count'])

# rsx == rsx1
# False

# But
# In [7]: rsx.data == rsx1.data
# Out[7]: 
#       LOND  LATD  VELU  VELV  VFLG  ESPC  MAXV  MINV  EDVC  ERSC  XDST  YDST  RNGE  BEAR  VELO  HEAD  SPRC
# 9     True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True
# 10    True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True
# 11    True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True
# 12    True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True  True

numpy.where(rsx.data != rsx1.data)



# --------------------------------
r = read_lluv_file(ifn)

# formatfile = Path(__file__).parent.resolve() / 'file_formats' / 'radialshort_LLUV_RDL7.ruv'
formatfile = os.path.join('/Users/codar/Documents/qccodar_dev/src/qccodar', 'file_formats', 'radialshort_LLUV_RDL7.ruv')
rs = Radial(formatfile,empty_radial=True)

# different tables between RadialMetric (RM) read in and the "template" for RadialShorts (RS)
# tt = tabletype. e.g. 'LLUV RDM1'
# so want to retrain only the first of two words.  e.g. 'LLUV'
r_tts = [r._tables[key]['TableType'] for key in r._tables.keys()]
r_tts = [tt.split(' ')[0] for tt in r_tts]
print(r_tts)

rs_tts = [rs._tables[key]['TableType'] for key in rs._tables.keys()]
rs_tts = [tt.split(' ')[0] for tt in rs_tts]
print(rs_tts)


print(r_tts)
print(rs_tts)

#
r_tts = [r._tables[key]['TableType'].split(' ')[0] for key in r._tables.keys()]
rs_tts = [rs._tables[key]['TableType'].split(' ')[0] for key in rs._tables.keys()]
rs_keys= [key for key in rs._tables.keys()]
to_be_removed = [tt for tt in rs_tts if tt not in r_tts]

for tt in to_be_removed:
    idx = rs_tts.index(tt)
    rs._tables.pop(rs_keys[idx])

rs._tables.keys()

r = add_short_metadata(r,qccodar_values)


write_output(r,'/Users/codar/documents/qccodar_dev/tools/test_output.ruv')
# this is culprit for verbose ouput
# line 1108 self.to_lluv() in radials.py of hfradarpy needs to bee commented out
# get latest version of hfradarpy from github ?
# or comment out myself
