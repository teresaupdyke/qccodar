#!/usr/bin/env python
#
# Last modified: Time-stamp: <2017-07-05 13:25:32 codar>
"""
Tests for reading lluv files.

We may need to deal with different types of LLUV files eventually.

"""
import os
from qccodar.qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')
ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', \
                   'RDLv_HATY_2013_11_05_0000.ruv')

def test_load_data():
    """
    test_load_data -- Read the file as list of lines
    """
    lines = load_data(ifn)
    assert type(lines) is list
    assert type(lines[0]) is str
    assert len(lines) == 8136

def test_read_lluv_file():
    """
    test_read_lluv_file -- Read typical LLUV data that has data
    
    """
    d, types_str, header, footer = read_lluv_file(ifn)

    assert header.split('\n')[0] == '%CTF: 1.00'
    assert re.search(r'%TableType:.*\n', header).group() == '%TableType: LLUV RDM1\n'
    assert re.search(r'%TableStart:.*\n', header).group()== '%TableStart:\n'
    assert footer.split('\n')[0] == '%TableEnd:'

    assert types_str == 'LOND LATD VELU VELV VFLG RNGE BEAR VELO HEAD SPRC SPDC MSEL MSA1 MDA1 MDA2 MEGR MPKR MOFR MSP1 MDP1 MDP2 MSW1 MDW1 MDW2 MSR1 MDR1 MDR2 MA1S MA2S MA3S MEI1 MEI2 MEI3 MDRJ '

    # check some elements of return numpy.array
    assert d.shape == (8029, 34)
    # not ideal to pull specific values out, but how else to test?
    # first and last, lon, lat
    assert numpy.array_equal(d[0,0:2], [-75.3316214,  35.227066])
    assert numpy.array_equal(d[-1,0:2], [-75.820443,  36.972371])

def test_read_lluv_empty_file():
    """
    test_read_lluv_empty_file -- Read LLUV data that has NO radial data
    
    """
    empty_ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_01', \
                   'RDLv_HATY_2013_11_01_1830.ruv')

    d, types_str, header, footer = read_lluv_file(empty_ifn)

    # find header and type_str
    assert header.split('\n')[0] == '%CTF: 1.00'
    # the footer should be empty
    assert len(footer) == 0
    # should still have columntypes as types_str
    assert types_str == 'LOND LATD VELU VELV VFLG RNGE BEAR VELO HEAD SPRC SPDC MSEL MSA1 MDA1 MDA2 MSP1 MDP1 MDP2 MSW1 MDW1 MDW2 MSR1 MDR1 MDR2 MA1S MA2S MA3S MEI1 MEI2 MEI3 '
    
    # check size of data (no rows if no data)
    assert d.size == 0
    assert numpy.array_equal(d, [])
   
def test_get_columns():
    """
    test_get_columns -- Parse %TableColumnTypes into dict
    """
    d, types_str, header, footer = read_lluv_file(ifn)
    c = get_columns(types_str)
    # test the keys
    for i, k in enumerate(sorted(c, key=c.__getitem__)):
        # assert c['LOND'] == 0
        assert eval("c['%s'] == %d" % (k, i))
        
    # test the values
    # [0, 1, ... 34] == [0, 1, ... 34]
    assert sorted(c.values()) == list(range(34))

