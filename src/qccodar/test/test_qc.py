#!/usr/bin/env python
#
# Last modified: Time-stamp: <2017-07-05 18:19:58 codar>

"""
Tests for qc thresholds and weighted averaging.
Using hfradarpy Radial Class
Using config qccodar_values

r0, d0 are the raw original data from using code
d1 is d0 processed by test

tr, td are the "truth" data from verified/saved files processed 
by each QC test used to confirm d1 is indeed the same

"""
import os
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')
# dev_path = '/Users/codar/Documents/qccodar_dev/src/qccodar'
# files = os.path.join(dev_path, 'test', 'files')

qccodar_values = dict(
  qc_doa_half_power_width=dict(doa_half_power_width_max=50.0),
  qc_doa_peak_power=dict(doa_peak_power_min=5.0),
  qc_monopole_snr=dict(monopole_snr_min=5.0),
  qc_loop_snr=dict(loop_snr_min=5.0),
  qc_radialshort_velocity_count=dict(radialshort_velocity_count_min=1.0),
  weighted_shorts=dict(numdegrees=3,weight_parameter='MP', table_type='LLUV RDL7'),
  merge=dict(css_interval_minutes=30.0,number_of_css=5.0)
  )


def test_no_qc():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    d0 = r0.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_test0', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    assert numpy.isclose(d0, td, equal_nan=True).all(), 'should be equal, including where NaN'

    
def test_threshold_qc_doa_peak_power():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # specify threshold in case default changes
    r1 = threshold_qc_doa_peak_power(r0, threshold=5.0)
    d1 = r1.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_test1', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_doa_half_power_width():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # specify threshold in case default changes
    r1 = threshold_qc_doa_half_power_width(r0, threshold=50.0)
    d1 = r1.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_test2', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_monopole_snr():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # specify threshold in case default changes
    r1 = threshold_qc_monopole_snr(r0, threshold=5.0)
    d1 = r1.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_test3', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_loop_snr():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # specify threshold in case default changes
    r1 = threshold_qc_loop_snr(r0, threshold=5.0)
    d1 = r1.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_test4', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'

def test_threshold_qc_all():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # specify threshold in case default changes -- previous version -- deprecated
    #_ r1 = threshold_qc_all(r0, thresholds=[5.0, 50.0, 5.0, 5.0])
    #
    # now using a structured dict() to hold / pass configuration info and test tolerances
    # specify thresholds here in case default changes

    r1 = threshold_qc_all(r0, qccodar_values)
    d1 = r1.data.to_numpy()
    #
    ifn2 = os.path.join(files, 'Radialmetric_testall', 'RDLv_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    #
    assert numpy.isclose(d1, td, equal_nan=True).all(), 'should be equal, including where NaN'


# Using early verified output RadialShorts as test data. These files
# were from earlier testing of weight function. They have fewer cells
# since cells were based on CODAR output RadialShorts.  For this test,
# we need to find where SPRC and BEAR are the same from td to xd.  So
# we are not using the whole array of xd.  We can add files later for
# whole array once we are comfortable with this function.


def test_weighted_average_mp_weight_angres1():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # weighted_velocities returns a pandas dataframe
    df1 = weighted_velocities(r0, numdegrees=1, weight_parameter='MP')
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
    assert numpy.isclose(subxd, subtd, equal_nan=True).all()

def _weighted_average_mp_weight_angres3():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # weighted_velocities returns a pandas dataframe
    df1 = weighted_velocities(r0, numdegrees=3, weight_parameter='MP')
    xd = df1.to_numpy()
    xc = get_columns( ' '.join(df1.columns.to_list()) )

    ifn2 = os.path.join(files, 'RadialShorts_mp_weight_angres3', 'RDLx_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    tc = get_columns( ' '.join(tr.data.columns.to_list()) )

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    assert numpy.isclose(subxd, subtd, rtol=1e-05, atol=1e-03, equal_nan=True).all()


def test_weighted_average_snr_weight_angres1():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # weighted_velocities returns a pandas dataframe
    df1 = weighted_velocities(r0, numdegrees=1, weight_parameter='SNR3')
    xd = df1.to_numpy()
    xc = get_columns( ' '.join(df1.columns.to_list()) )

    ifn2 = os.path.join(files, 'RadialShorts_snr_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    tc = get_columns( ' '.join(tr.data.columns.to_list()) )

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    assert numpy.isclose(subxd, subtd, equal_nan=True).all()

def _weighted_average_snr_weight_angres3():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # weighted_velocities returns a pandas dataframe
    df1 = weighted_velocities(r0, numdegrees=3, weight_parameter='SNR3')
    xd = df1.to_numpy()
    xc = get_columns( ' '.join(df1.columns.to_list()) )

    ifn2 = os.path.join(files, 'RadialShorts_snr_weight_angres3', 'RDLx_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    tc = get_columns( ' '.join(tr.data.columns.to_list()) )

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    assert numpy.isclose(subxd, subtd, rtol=1e-05, atol=1e-03, equal_nan=True).all()


def test_weighted_average_no_weight_angres1():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    r0 = read_lluv_file(ifn)
    # weighted_velocities returns a pandas dataframe
    df1 = weighted_velocities(r0, numdegrees=1, weight_parameter='NONE')
    xd = df1.to_numpy()
    xc = get_columns( ' '.join(df1.columns.to_list()) )

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    tr = read_lluv_file(ifn2)
    td = tr.data.to_numpy()
    tc = get_columns( ' '.join(tr.data.columns.to_list()) )

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)
    # trows, xrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], xd[:, [xc['SPRC'], xc['BEAR']]])

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    # subxd VELO is close to subtd VELO within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subxd, subtd, rtol=1e-05, atol=1e-03, equal_nan=True).all()


def _scratch():
    ofn = os.path.join(files, 'test1_output.txt')
    write_output(ofn, header, d1, footer)
    #
    idx = numpy.where(d1 != td)    
    assert numpy.isnan(d[idx]).all()
    assert numpy.isnan(d2[idx]).all()
    for i,j in numpy.array(idx).T:
        # if not numpy.isnan(d1[i,j]):
            print("(%4d, %4d) %5g %5g" % (i,j, d1[i,j], td[i,j]))    

def _generate_output():
    # used this subroutine to generate output when assured qc functions correct
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    # specify threshold in case default changes
    # e.g. test4
    d4 = threshold_qc_loop_snr(d, types_str, threshold=5.0)
    #
    ofn = os.path.join(files, 'Radialmetric_test4', 'RDLv_HATY_2013_11_05_0000.ruv')
    write_output(ofn, header, d4, footer)
