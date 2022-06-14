#!/usr/bin/env python
#
# Last modified: Time-stamp: <2017-07-05 13:24:32 codar>
"""
Test functions for generating radialshort data output.

"""
import os
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')

def test_generate_radialshort_array():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    
    # no qc thresholds and this no weighting and 1 deg angres (bearing_offset=0.0) is same as CODAR processing
    xd, xtypes_str = weighted_velocities(d, types_str, 1, 'NONE')
    ####################
    # TESTING
    rsd, rsdtypes_str = generate_radialshort_array(xd, xtypes_str, header)
    ####################
    rsc = get_columns(rsdtypes_str)

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    tc = get_columns(ttypes_str)

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    rsrngbear = rsd[:, [rsc['SPRC'], rsc['BEAR']]]
    trows, rsrows = cell_intersect(trngbear, rsrngbear)
    # trows, rsrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], rsd[:, [rsc['SPRC'], rsc['BEAR']]])
    # 
    subtd = td[trows, :]
    subrsd = rsd[rsrows, :]

    # 1st make sure simple Range cell and bearing columns okay, should be simplest
    tcol  = numpy.array([tc['BEAR'], tc['SPRC']])
    rscol = numpy.array([rsc['BEAR'], rsc['SPRC']])
    # subrsd is close to subtd  within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with BEAR or SPRC, not close to CODAR '

    # range is computed, now check this
    tcol  = numpy.array([tc['RNGE']])
    rscol = numpy.array([rsc['RNGE']])
    # subrsd is close to subtd  within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with RNGE, not close to CODAR '

    # LATD, LOND uses special geopy.distance.vincenty computation -- 
    # these are the columns generated from d that we want to compare with CODAR radialshorts
    tcol  = numpy.array([tc['LOND'], tc['LATD']])
    rscol = numpy.array([rsc['LOND'], rsc['LATD']])
    # subrsd is close to subtd  within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with LATD or LOND, not close to CODAR '
    
    # 2nd then test the columns that are filled 
    # these are the columns filled from xd VELO and that we want to compare with CODAR radialshorts
    tcol  = numpy.array([tc['VELO']])
    rscol = numpy.array([rsc['VELO']])
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with VELOs'

    tcol  = numpy.array([tc['XDST'], tc['YDST']])
    rscol = numpy.array([rsc['XDST'], rsc['YDST']])
    assert numpy.isclose(subrsd[:,rscol], subtd[:,tcol], rtol=1e-05, atol=1e-03, equal_nan=True).all(), \
        'something wrong with XDST, YDST'

    assert numpy.isclose(subrsd[:,rsc['VELU']], subtd[:,tc['VELU']], \
                         rtol=1e-01, atol=1e-01, equal_nan=True).all(), \
        'something wrong with VELU, is not close to CODAR VELU'

    ########################
    # I think that something is wrong with CODAR VELV computation.
    # Need to investigate further.  when BEAR is ~ 90 and HEAD is 270
    # is when these are not close or even approximate.  Might be true
    # for HATY orientation and proximity to GS. Might be different for
    # another site with different orientation and current with such
    # magnitude as in the GS.
    ########################

    # assert numpy.isclose(subrsd[:,rsc['VELV']], subtd[:,tc['VELV']], \
    #                      rtol=1e-01, atol=1e-01, equal_nan=True).all(), \
    #     'something wrong with VELV, is not close to CODAR VELV'


def test_generate_radialshort_header():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    # no qc thresholds and this no weighting and 1 deg angres (bearing_offset=0.0) is same as CODAR processing
    xd, xtypes_str = weighted_velocities(d, types_str, 1, 'NONE')
    rsd, rsdtypes_str = generate_radialshort_array(xd, xtypes_str, header) 
    ####################
    # TESTING    
    rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
    ####################
    assert rsdheader.split('\n')[0] == '%CTF: 1.00'
    assert re.search(r'%TableColumnTypes:.*\n', rsdheader).group() == \
        '%TableColumnTypes: LOND LATD VELU VELV VFLG ESPC MAXV MINV EDVC ERSC XDST YDST RNGE BEAR VELO HEAD SPRC\n'
    assert re.search(r'%TableType:.*\n', rsdheader).group() == '%TableType: LLUV RDL7\n'
    assert re.search(r'%TableStart:.*\n', rsdheader).group()== '%TableStart:\n'

    assert int(re.search(r'%TableRows:\s(\d*?)\n', rsdheader).groups()[0]) == rsd.shape[0]
    assert int(re.search(r'%TableColumns:\s(\d*?)\n', rsdheader).groups()[0]) == rsd.shape[1]

# tests for dealing with empty arrays or when there are no radials -- need empty radialshort data
def test_generate_radialshort_empty_array_when_no_radial_data():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_01', 'RDLv_HATY_2013_11_01_1830.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)

    # Make sure no data first
    assert d.size == 0

    ####################
    # TESTING
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str, header) 
    ####################

    assert rsd.size == 0
    assert rsdtypes_str == 'LOND LATD VELU VELV VFLG ESPC MAXV MINV EDVC ERSC XDST YDST RNGE BEAR VELO HEAD SPRC'

def test_generate_radialshort_header_for_empty_array_when_no_radial_data():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_01', 'RDLv_HATY_2013_11_01_1830.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    rsd, rsdtypes_str = generate_radialshort_array(d, types_str, header) 

    ####################
    # TESTING
    rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
    ####################

    rsc = get_columns(rsdtypes_str)

    assert rsdheader.split('\n')[0] == '%CTF: 1.00'
    assert re.search(r'%TableColumnTypes:.*\n', rsdheader).group() == \
        '%TableColumnTypes: LOND LATD VELU VELV VFLG ESPC MAXV MINV EDVC ERSC XDST YDST RNGE BEAR VELO HEAD SPRC\n'
    assert re.search(r'%TableType:.*\n', rsdheader).group() == '%TableType: LLUV RDL7\n'
    assert re.search(r'%TableStart:.*\n', rsdheader).group()== '%TableStart:\n'

    assert int(re.search(r'%TableRows:\s(\d*?)\n', rsdheader).groups()[0]) == 0
    assert int(re.search(r'%TableColumns:\s(\d*?)\n', rsdheader).groups()[0]) == len(rsc)

def _scratch():
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', 'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    xd, xtypes_str = weighted_velocities(d, types_str, numdegrees=1, weight_parameter='NONE')
    xc = get_columns(xtypes_str)

    ifn2 = os.path.join(files, 'RadialShorts_no_weight_angres1', 'RDLx_HATY_2013_11_05_0000.ruv')
    td, ttypes_str, theader, tfooter = read_lluv_file(ifn2)
    tc = get_columns(ttypes_str)

    trngbear = td[:, [tc['SPRC'], tc['BEAR']]]
    xrngbear = xd[:, [xc['SPRC'], xc['BEAR']]]
    trows, xrows = cell_intersect(trngbear, xrngbear)
    # trows, xrows = cell_intersect(td[:, [tc['SPRC'], tc['BEAR']]], xd[:, [xc['SPRC'], xc['BEAR']]])

    subtd = td[trows, tc['VELO']]
    subxd = xd[xrows, xc['VELO']]
    # subxd VELO is close to subtd VELO within 1/1000 th, since test data was output by CODAR
    assert numpy.isclose(subxd, subtd, rtol=1e-05, atol=1e-03, equal_nan=True).all()
