#!/usr/bin/env python
#
# Last modified: Time-stamp: <2017-07-05 13:25:48 codar>
"""
Tests for writing lluv files.

We may need to deal with different types of LLUV files eventually.

"""
import os
import numpy
numpy.set_printoptions(suppress=True)
from qccodar.qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')

def test_write_empty_output_by_readback():
    """ 
    Write empty LLUV file output, and test by comparing to readback.

    """
    # get the empty file example
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_01', \
                   'RDLv_HATY_2013_11_01_1830.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    ofn = os.path.join(files, 'test_output.txt')
    write_output(ofn, header, d, footer)
    ifn2 = ofn
    d2, types_str2, header2, footer2 = read_lluv_file(ifn2)

    assert header == header2
    assert footer == footer2 
    assert numpy.isclose(d, d2, equal_nan=True).all(), 'should be equal, including where NaN'
    assert d.size == 0, 'should be empty'
    assert d2.size == 0, 'should be emtpy'

def test_write_output_by_readback():
    """ 
    Write typical LLUV file output test by comparing to readback.

    """
    ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', \
                   'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    ofn = os.path.join(files, 'test_output.txt')
    write_output(ofn, header, d, footer)
    ifn2 = ofn
    d2, types_str2, header2, footer2 = read_lluv_file(ifn2)

    assert header == header2
    assert footer == footer2
    assert numpy.isclose(d, d2, equal_nan=True).all(), 'should be equal, including where NaN'

    # Asserting that the data arrays are equal, raised an interesting 
    # issue about NaN's.  It turns out where they are not equal is 
    # where d and d2 are nan's.
    # 
    # >>> idx = numpy.where(d!=d2)
    # They are all NaNs
    # >>> assert numpy.isnan(d[idx]).all()
    # >>> assert numpy.isnan(d2[idx]).all()
    #
    # but they do not equal 
    #
    # >>> assert numpy.array_equal(d, d2) 
    # False
    # >>> assert numpy.array_equiv(d, d2)
    # False
    #
    # >>> ix = idx[0]
    # >>> iy = idx[1]
    # >>> d[ix[0],iy[0]] == d2[ix[0],iy[0]]
    # False
    # >>> d[ix,iy] == d2[ix,iy]
    # array([False, False, False, ..., False, False, False], dtype=bool)
    #
    # Thought this might work but if any value is NaN it will return False
    # >>> assert numpy.allclose(d, d2)
    # False
    # 
    # Finally, found numpy.isclose where we can specify nan equal
    # help(numpy.isclose)
    # isclose(a, b, rtol=1e-05, atol=1e-08, equal_nan=False)
    #
    # >>> numpy.isclose(d, d2).all()
    # False
    # >>> numpy.isclose(d, d2, equal_nan=True).all()
    # True

