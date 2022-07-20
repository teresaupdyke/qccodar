#!/usr/bin/env python
#
# Last modified: Time-stamp: <2022-06-29 18:26:03 codar>
"""
Tests for reading lluv files. 
Now using hfradarpy Radial Class

We may need to deal with different types of LLUV files eventually.

"""
import os

from qccodar.codarutils import *
from qccodar.qcutils import *

files = os.path.join(os.path.curdir, 'test', 'files')
ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_05', \
                   'RDLv_HATY_2013_11_05_0000.ruv')

def load_data(inFile):
    lines=None
    if os.path.exists(inFile):
        f = open(inFile, 'r')
        lines = f.readlines()
        f.close()
        if len(lines)<=0:
            print('Empty file: '+ inFile)
            raise EOFError('Empty File: %s' % inFile)
    else:
        print('File does not exist: '+ inFile)
        raise IOError('Error opening %s' % inFile)
    return lines

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
    test_read_lluv_file -- Read typical LLUV data that has data as a Radial object
    
    """
    r = read_lluv_file(ifn)
    d = r.data # now a pandas dataframe 

    assert r.metadata['CTF'] == '1.00'
    assert r.metadata['FileType'] == 'LLUV rdls "RadialMap"'

    # check first table (should be data in LLUV format)
    assert r._tables[1]['TableType'] == 'LLUV RDM1'
    assert r._tables[1]['TableColumns'] == '34'
    assert r._tables[1]['TableRows'] == '8029'
    
    # check column labels in dataframe as expected from RadialMetric data
    cols = d.columns.to_list()
    cols_str = 'LOND LATD VELU VELV VFLG RNGE BEAR VELO HEAD SPRC SPDC MSEL MSA1 MDA1 MDA2 MEGR MPKR MOFR MSP1 MDP1 MDP2 MSW1 MDW1 MDW2 MSR1 MDR1 MDR2 MA1S MA2S MA3S MEI1 MEI2 MEI3 MDRJ'
    assert ' '.join(cols) == cols_str

    # spot check some elements of data are read correctly
    assert d.shape == (8029, 34)
    # not ideal to pull specific values out, but how else to test?
    # first and last, lon, lat
    assert numpy.array_equal(d.iloc[0,0:2].to_numpy(), [-75.3316214,  35.227066])
    assert numpy.array_equal(d.iloc[-1,0:2].to_numpy(), [-75.820443,  36.972371])

def test_read_lluv_empty_file():
    """
    test_read_lluv_empty_file -- Read LLUV data that has NO radial data as a Radial object
    
    """
    empty_ifn = os.path.join(files, 'codar_raw', 'Radialmetric_HATY_2013_11_01', \
                   'RDLv_HATY_2013_11_01_1830.ruv')

    r = read_lluv_file(empty_ifn)
    d = r.data

    # check column labels in dataframe as expected from RadialMetric data
    cols = d.columns.to_list()
    cols_str = 'LOND LATD VELU VELV VFLG RNGE BEAR VELO HEAD SPRC SPDC MSEL MSA1 MDA1 MDA2 MSP1 MDP1 MDP2 MSW1 MDW1 MDW2 MSR1 MDR1 MDR2 MA1S MA2S MA3S MEI1 MEI2 MEI3'
    assert ' '.join(cols) == cols_str
    
    # check size of data (no rows if no data)
    assert d.size == 0
    shape = d.to_numpy().shape
    assert shape == (0, 30)
    assert numpy.array_equal(d.to_numpy(), numpy.empty(shape))
   
