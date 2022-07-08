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

#################
# testing revisit_weighted_velocities
#################
from revisit_weighted_velocities import *

dev_path = '/Users/codar/Documents/qccodar_dev/src/qccodar'
files = os.path.join(dev_path, 'test', 'files')
ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')

r0 = read_lluv_file(ifn)
# weighted_velocities returns a pandas dataframe
# NOTE we are using the numpy array version
df1 = weighted_velocities_np(r0, numdegrees=1, weight_parameter='MP')
xd = df1.to_numpy()
xc = get_columns( ' '.join(df1.columns.to_list()) )

ifn2 = os.path.join(files, 'RadialShorts_mp_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
tr = read_lluv_file(ifn2)
td = tr.data.to_numpy()
tc = get_columns( ' '.join(tr.data.columns.to_list()) )

trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
trows, xrows = cell_intersect(trngbear, xrngbear)

subtd = td[trows, tc['VELO']]
subxd = xd[xrows, xc['VELO']]
# test that these are the same 
numpy.isclose(subxd, subtd, equal_nan=True).all()

#######################
# timing the two versions
#######################

ifn = '/Volumes/TRANSFER/reprocess_HATY/2017_01/RadialMetric/IdealPattern/RDLv_HATY_2017_01_01_0000.ruv'
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

# time the two versions
%timeit df1 = weighted_velocities_df(rmqc, numdegrees=3, weight_parameter='MP')
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
#  ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES
# 52.1 s � 223 ms per loop (mean � std. dev. of 7 runs, 1 loop each)


%timeit df2 = weighted_velocities_np(rmqc, numdegrees=3, weight_parameter='MP')
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
#  ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES
# 1.15 s � 14.6 ms per loop (mean � std. dev. of 7 runs, 1 loop each)


# (2) do weighted averaging of good, this is done when the radial short is created
# rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
# In [4]: %timeit rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
# 48.3 s ± 598 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)


# (3) require a minimum numpoints used in to form cell average
# rsx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])
# In [7]: %timeit rsxx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])
# 1.63 ms ± 73.1 µs per loop (mean ± std. dev. of 7 runs, 1000 loops each)


