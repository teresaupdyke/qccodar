#

import sys
import copy
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.codarutils import *
from qccodar.qcutils import *
from qccodar.app import load_configs

configfile = '../src/qccodar/config/qccodar.plist'
# configfile = '' # test default
qccodar_values = load_configs(configfile)

ifn = './2017_01/RadialMetric/IdealPattern/RDLv_HATY_2017_01_01_0000.ruv'

r = read_lluv_file(ifn)

def threshold_qc_doa_peak_power_Radial_copy(r, threshold=5.0):
#def threshold_qc_doa_peak_power(d, types_str, threshold=5.0):
    """Bad Flag any DOA peak power (dB) less than threshold value (default 5.0 dB).

    Flags any direction of arrival (DOA) peak power (dB) that falls
    below the input threshold value (default 5.0 dB).  Depending on
    the value of MSEL (1, 2, or 3), MSR1, MDR1, or MDR2 columns are
    evaluated.  Returns modified matrix with VFLG column the only
    changed values.

    """

    r1 = copy.deepcopy(r)
    havenan = numpy.isnan(r1.data['MSR1']) | numpy.isnan(r1.data['MDR1']) | numpy.isnan(r1.data['MDR2'])
    bad = (r1.data['MSEL']==1) & (r1.data['MSR1']<float(threshold))| \
          ((r1.data['MSEL']==2) & (r1.data['MDR1']<float(threshold))) | \
          ((r1.data['MSEL']==3) & (r1.data['MDR2']<float(threshold))) | havenan
    r1.data.loc[bad, 'VFLG'] = r1.data.loc[bad, 'VFLG'] + (1<<1)
    return r1


def threshold_qc_doa_peak_power_Pandas_copy(r, threshold=5.0):
#def threshold_qc_doa_peak_power(d, types_str, threshold=5.0):
    """Bad Flag any DOA peak power (dB) less than threshold value (default 5.0 dB).

    Flags any direction of arrival (DOA) peak power (dB) that falls
    below the input threshold value (default 5.0 dB).  Depending on
    the value of MSEL (1, 2, or 3), MSR1, MDR1, or MDR2 columns are
    evaluated.  Returns modified matrix with VFLG column the only
    changed values.

    """
    
    data = copy.deepcopy(r.data)
    havenan = numpy.isnan(data['MSR1']) | numpy.isnan(data['MDR1']) | numpy.isnan(data['MDR2'])
    bad = (data['MSEL']==1) & (data['MSR1']<float(threshold))| \
          ((data['MSEL']==2) & (data['MDR1']<float(threshold))) | \
          ((data['MSEL']==3) & (data['MDR2']<float(threshold))) | havenan
    data.loc[bad, 'VFLG'] = data.loc[bad, 'VFLG'] + (1<<1)
    r.data = data
    return r


rR = threshold_qc_doa_peak_power_Radial_copy(r, threshold=5.0)
rPD = threshold_qc_doa_peak_power_Pandas_copy(r, threshold=5.0)





