#!/usr/bin/env python
#
# Last modified: Time-stamp: <2019-12-03 10:25:10 codar>

"""Quality control (QC) functions for CODAR SeaSonde Radialmetric data

QC categories:
A. Threshold Tests -- badflag any values that fall below or above a single threshold
B. Weighted Averaging -- average several values with weights based on signal quality parameters 

QC Threshold Tests:
1. DOA peak power (MSR1, MDR1, MDR2) < 5 dB default 
2. DOA 1/2 power width (3dB down) (MSW1, MDW1, MDW2) > 50 deg default
3. SNR on monopole (MA3S) < 5 dB default
4. SNR on both loops (MA1S and MA2S) < 5 dB

Weighted Averaging:
1. Weighting based on Music Power (MSP1, MDP1, MDP2)
2. Weighting based on SNR on monopole (MA3S)
3. No weight function (None) 

"""
import sys
import os
import re
import fnmatch
import datetime

import numpy
numpy.set_printoptions(suppress=True)

# 
from .codarutils import *

debug = 1

def _commonly_assigned_columns():
    """
    Commonly assigned CODAR RadialMetric columns
    """
    # for cut and pasting to other functions
    # help make some code more readable
    VFLG = c['VFLG']
    MSEL = c['MSEL']
    MSR1 = c['MSR1']
    MDR1 = c['MDR1']
    MDR2 = c['MDR2']
    MSW1 = c['MSW1'] 
    MDW1 = c['MDW1']
    MDW2 = c['MDW2']
    MA3S = c['MA3S']

def threshold_qc_doa_peak_power(d, types_str, threshold=5.0):
    """Bad Flag any DOA peak power (dB) less than threshold value (default 5.0 dB).

    Flags any direction of arrival (DOA) peak power (dB) that falls
    below the input threshold value (default 5.0 dB).  Depending on
    the value of MSEL (1, 2, or 3), MSR1, MDR1, or MDR2 columns are
    evaluated.  Returns modified matrix with VFLG column the only
    changed values.

    """
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MSEL = c['MSEL']
    MSR1 = c['MSR1']
    MDR1 = c['MDR1']
    MDR2 = c['MDR2']
    d1=numpy.copy(d)
    # 
    havenan = numpy.isnan(d[:,MSR1]) | numpy.isnan(d[:,MDR1]) | numpy.isnan(d[:,MDR2])
    bad = ((d[:,MSEL]==1) & (d[:,MSR1]<float(threshold))) | \
          ((d[:,MSEL]==2) & (d[:,MDR1]<float(threshold))) | \
          ((d[:,MSEL]==3) & (d[:,MDR2]<float(threshold))) | havenan
    d1[bad, VFLG] = d[bad,VFLG]+(1<<1)
    return d1

def threshold_qc_doa_half_power_width(d, types_str, threshold=50.0):
    """Bad Flag DOA 1/2 Power Width (degress) greater than threshold value (default 50.0 degrees).

    Flags any direction of arrival (DOA) 1/2 Power width (degress)
    that is wider than the input threshold value (default 50.0
    degrees).  Depending on the value of MSEL (1, 2, or 3), MSW1,
    MDW1, or MDW2 columns are evaluated.  Returns modified matrix with
    VFLG column the only changed values.

    """
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MSEL = c['MSEL']
    MSW1 = c['MSW1'] 
    MDW1 = c['MDW1']
    MDW2 = c['MDW2']
    d2=numpy.copy(d)
    # 
    havenan = numpy.isnan(d[:,MSW1]) | numpy.isnan(d[:,MDW1]) | numpy.isnan(d[:,MDW2])
    bad = ((d[:,MSEL]==1) & (d[:,MSW1]>float(threshold))) | \
          ((d[:,MSEL]==2) & (d[:,MDW1]>float(threshold))) | \
          ((d[:,MSEL]==3) & (d[:,MDW2]>float(threshold))) | havenan
    d2[bad, VFLG] = d[bad,VFLG]+(1<<2)
    return d2

def threshold_qc_monopole_snr(d, types_str, threshold=5.0):
    """Bad flag any SNR on monopole (dB)  less than threshold value (default 5.0 dB).

    Flags any signal-to-noise ratio (SNR) on monopole (dB) that falls
    below the input threshold value (default 5.0 dB).  No dependency on MSEL selections.

    """
    # Test 3 SNR on monopole (dB) for all selections
    # 
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MA3S = c['MA3S']
    d3=numpy.copy(d)
    # 
    bad = d[:,MA3S]<float(threshold)
    d3[bad, VFLG] = d[bad,VFLG]+(1<<3)
    return d3

def threshold_qc_loop_snr(d, types_str, threshold=5.0):
    """Bad flag if both loop SNR are less than threshold value (default 5.0 dB).

    Flags if signal-to-noise ratio (SNR) (dB) on loop1 AND on loop2 falls
    below the input threshold value (default 5.0 dB). No dependency on MSEL selections.

    """
    # Test 4 SNR on both loop antennas (dB) for all selections
    # 
    c = get_columns(types_str)
    VFLG = c['VFLG'] # help make the test more readable
    MA1S = c['MA1S']
    MA2S = c['MA2S']
    d4=numpy.copy(d)
    # 
    bad = (d[:,MA1S]<float(threshold)) & (d[:,MA2S]<float(threshold))
    d4[bad, VFLG] = d[bad,VFLG]+(1<<3)
    return d4

def threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0, 5.0]):
    """Combine all three threshold tests

    Returns modified matrix with VFLG column only changed values.

    """
    # TO DO: Use a dict or list to pass thresholds for all tests ??
    # if dict what key to use
    # if list how know the correct order for tests
    # 
    dall = threshold_qc_doa_peak_power(d, types_str, thresholds[0])
    dall = threshold_qc_doa_half_power_width(dall, types_str, thresholds[1])
    dall = threshold_qc_monopole_snr(dall, types_str, thresholds[2])
    dall = threshold_qc_loop_snr(dall, types_str, thresholds[3])
    
    return dall

def threshold_rsd_numpoints(rsd, rstypes_str, numpoints=1):
    """Bad flag any radialshort data with doppler velocity count (EDVC) less than "numpoints"

    Returns modified rsd matrix with VFLG column only changed if EDVC
    count is less than or equal to numpoints.  This threshold is
    checked after weighted_velocities()

    """
    if rsd.size == 0:
        return numpy.array([])

    rsc = get_columns(rstypes_str)
    VFLG = rsc['VFLG'] # help make the test more readable
    EDVC = rsc['EDVC']
    rsd1=numpy.copy(rsd)
    # 
    bad = rsd[:,EDVC]<int(numpoints)
    rsd1[bad, VFLG] = rsd[bad,VFLG]+(1<<12) # of dubious quality and should not be used or displayed
    return rsd1
   

def weighted_velocities(d, types_str, numdegrees=3, weight_parameter='MP'):
    """Calculates weighted average of radial velocities (VELO) at bearing and range.

    The weighted average of velocities found at given range and
    bearing based on weight_parameter.

    Paramters
    ---------
    d : ndarray
        The data from LLUV file(s). 
    types_str : string 
        The 'TalbleColumnTypes' string header of LLUV file(s) provide keys for each column.
    weight_parameter : string ('MP', 'SNR3', 'NONE'), optional 
        If 'MP' (default), uses MUSIC antenna peak power values for weighting function
           using MSEL to select one of (MSP1, MDP1, or MDP2).
        If 'SNR3', uses signal-to-noise ratio on monopole (MA3S).
        If 'NONE', just average with no weighting performed.
    numdegrees: int, optional (default 3 degree)
       The number of degrees of bearing from which to get velocities to spatially average over.
       For example, 
          If 1 deg, velocities from window of 1 deg will be averaged.
          If 3 deg, velocities from a window of 3 degrees will be averaged. This is the default.
          If 5 deg, velocities from a window of 5 degrees will be averaged.

    Returns
    -------
    xd : ndarray
       The averaged values with range and bearing.
       An array with averaged values, range, bearing, 
    xtypes_str : string 
        The order and key-labels for each column of xd array

    """
    # 
    # order of columns and labels for output data
    xtypes_str = 'VFLG SPRC BEAR VELO ESPC MAXV MINV EDVC ERSC'
    xc = get_columns(xtypes_str)
    #
    c = get_columns(types_str)
    offset = ((numdegrees-1)/2)
    # 
    ud = unique_rows(d[:,[c['SPRC'],c['BEAR'],c['VFLG']]].copy())
    # return only rows that have VFLG==0 (0 == good, >0 bad) so only get good data
    ud = ud[ud[:,2]==0]
    if ud.size == 0:
        return numpy.array([]), xtypes_str
    
    #
    allbearings = numpy.unique(ud[:,1])
    allranges = numpy.unique(ud[:,0])
    ud = numpy.array([[r,b] for r in allranges for b in allbearings])

    #
    nrows, _ = ud.shape
    ncols = len(xc)
    xd = numpy.ones(shape=(nrows,ncols))*numpy.nan
    #
    for irow, cell in enumerate(ud):
        rngcell, bearing = cell[0:2]
        # numpy.where() returns a tuple for condition so use numpy.where()[0]
        # also VFLG must equal 0 (0 == good, >0 bad) so only get good data
        xrow = numpy.where((d[:,c['SPRC']]==rngcell) & \
                           (d[:,c['BEAR']]>=bearing-offset) & \
                           (d[:,c['BEAR']]<=bearing+offset) & \
                           (d[:,c['VFLG']]==0))[0]
        # If no row matches rngcell AND bearing, then no VELO data, skip to next bearing
        if xrow.size == 0: 
            continue
        
        xcol = numpy.array([c['VELO'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2'], c['MA3S']])
        a = d[numpy.ix_(xrow, xcol)].copy()
        # if xrow.size == edvc:
        VELO = a[:,0] # all radial velocities found in cell
        SNR3 = a[:,5] # SNR on monopole for each velocity
        if weight_parameter.upper() == 'MP':
            # Create array to hold each Music Power (based on MSEL)
            MP = numpy.array(numpy.ones(VELO.shape)*numpy.nan) 
            # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
            for msel in [1, 2, 3]:
                which = a[:,1]==msel
                MP[which,] = a[which, msel+1]
            # convert MP from db to voltage for weighting
            MP = numpy.power(10, MP/10.)
            wts = MP/MP.sum()
            velo = numpy.dot(VELO,wts)
        elif weight_parameter.upper() == 'SNR3' or weight_parameter.upper() == 'SNR':
            wts = SNR3/SNR3.sum()
            velo = numpy.dot(VELO,wts)
        elif weight_parameter.upper() == 'NONE':
            # do no weighting and just compute the mean of all velo's
            velo = VELO.mean()
        # data
        xd[irow,xc['VFLG']] = 0
        xd[irow,xc['SPRC']] = rngcell
        xd[irow,xc['BEAR']] = bearing
        xd[irow,xc['VELO']] = velo
        # other stat output
        xd[irow,xc['ESPC']] = VELO.std() # ESPC
        xd[irow,xc['MAXV']] = VELO.max() # MAXV
        xd[irow,xc['MINV']] = VELO.min() # MINV
        # (EDVC and ERSC are the same in this subroutine's context)
        xd[irow,xc['EDVC']] = VELO.size # EDVC Velocity Count 
        xd[irow,xc['ERSC']] = VELO.size # ERSC Spatial Count
                
    # delete extra lines (nan) not filled above
    wherenan = numpy.where(numpy.isnan(xd[:,xc['VFLG']]))[0]
    xd = numpy.delete(xd, wherenan, axis=0)

    return xd, xtypes_str


def recursive_glob(treeroot, pattern):
    """ Glob-like search for filenames based on pattern by recursve walk 
    subdirectories starting in treeroot.

    Parameters
    ----------
    treeroot : string
       The top most directory path to begin search.
    pattern : string
       The pattern to match file in search.

    Return
    ------
    results : list of paths from treeroot
       The results of search.

    >>> files = os.path.join(os.path.curdir, 'test', 'files')
    >>> fns = recursive_glob(files, 'RDLx*.*')
    
    """
    # fnmatch gives you exactly the same patterns as glob, so this is
    # really an excellent replacement for glob.glob with very close
    # semantics.  Pasted from
    # <http://stackoverflow.com/questions/2186525/use-a-glob-to-find-files-recursively-in-python>
    results = [] 
    for base, dirs, files in os.walk(treeroot): 
        goodfiles = fnmatch.filter(files, pattern) 
        results.extend(os.path.join(base, f) for f in goodfiles)

    # sort the unsorted results
    results.sort()
    return results 

def filt_datetime(input_string, pattern=None):
    """Attempts to filter date and time from input string based on regex pattern.
    
    Default pattern follows the template, YYYY(-)MM(-)DD(-)(hh(:)(mm(:)(ss)))
    with minimum of YYYY MM and DD (date) supplied in descending order to 
    return its datetime object, otherwise returns None.

    Typical matches include, YYYYMMDD-hhmmss, YYYY-MM-DD-hh:mm:ss All
    the following will produce the corresponding datetime object. Any
    2-digit year will return None.

    Requires date with all three (year, month, day) in decreasing
    order as integers. Time is optional.

    >>> filt_datetime('RDLv_HATY_2013_11_05_000000.ruv')
    datetime.datetime(2013, 11, 5, 0, 0)

    >>> filt_datetime('RDLv_HATY_2013_11_05_0000.ruv')
    datetime.datetime(2013, 11, 5, 0, 0)
    
    >>> filt_datetime('RDLv_HATY_2013_11_05_00.ruv')
    datetime.datetime(2013, 11, 5, 0, 0)
    
    >>> filt_datetime('RDLv_HATY_2013_11_05.ruv')
    datetime.datetime(2013, 11, 5, 0, 0)
    
    # NOTE: returns None
    >>> filt_datetime('RDLv_HATY_13_11_05.ruv')
    
    >>> filt_datetime('RDLv_HATY_2013-11-05T00:00:00.ruv')
    datetime.datetime(2013, 11, 5, 0, 0)

    """

    #  the default pattern for a typical codar time stamp format
    if not pattern:
        pattern = r"""
        # YYYY(-)MM(-)DD(-)(hh(:)(mm(:)(ss)))
        (\d{4})           # 4-digit YEAR 
        \D?               # optional 1 character non-digit separator (e.g. ' ' or '-')
        (\d{2})           # 2-digit MONTH 
        \D?               # optional 1 character non-digit separator
        (\d{2})           # 2-digit DAY 
        \D?               # optional 1 character non-digit separator (e.g. ' ' or 'T')
        (\d{2})?          # optional 2-digit HOUR 
        \D?               # optional 1 character non-digit separator (e.g. ' ' or ':')
        (\d{2})?          # optional 2-digit MINUTE 
        \D?               # optional 1 character non-digit separator (e.g. ' ' or ':')
        (\d{2})?          # optional 2-digit SECOND
        """
    #         
    p = re.compile(pattern, re.VERBOSE)
    m = p.search(input_string) 
    # m.groups() # should be ('2013', '11', '05', '00', '00', None) for 'RDLv_HATY_2013_11_05_0000.ruv'
    if m:
        values = [int(yi) for yi in m.groups() if yi is not None] # [2013, 11, 5, 0, 0]
        # datetime.datetime(*v) requires mininum of year, month, day
        dt = datetime.datetime(*values) # datetime.datetime(2013, 11, 5, 0, 0)
    else:
        dt = None
    return dt

def find_files_to_merge(ifn, numfiles=3, sample_interval=30):
    """Finds the files that will be used in averaging in addition to ifn,
    based on sample_interval and numfiles to average over.
 
    Parameters:
    -----------
    ifn : string
       The complete path and filename of target date time to process.
    numfiles : int
       The number of files to average over, including input filename (ifn).
    sample_interval : int
       The sample interval in minutes. For radialmetric should be the
       same as CODAR's output rate.

    Return
    ------
    files : list of strings
       The files to process, including ifn.

    """

    indir = os.path.dirname(ifn)
    rdlstr = re.match(r'RDL[vw]', os.path.basename(ifn)).group()
    all_files = recursive_glob(os.path.join(indir), rdlstr+'*.ruv')

    delta_minutes = ((numfiles-1)/2)*sample_interval
    target_dt = filt_datetime(os.path.basename(ifn))
    dt_start = target_dt - datetime.timedelta(minutes=delta_minutes)
    dt_end = target_dt + datetime.timedelta(minutes=delta_minutes)
    
    files = []
    for fn in all_files:
        dt = filt_datetime(os.path.basename(fn))
        if dt is not None:
            if dt_start <= dt <= dt_end:
                files.append(fn)
    assert len(files) <= numfiles, \
        "Some duplicate files found since number found > numfiles needed "
    return files           

def do_qc(datadir, fn, patterntype):
    """ Do qc and then average over 3 sample_intervals (time), 3 degrees of bearing.
    """
    # read in the data
    rmfoldername = get_radialmetric_foldername(datadir)
    ifn = os.path.join(datadir, rmfoldername, patterntype, fn)
    d, types_str, header, footer = read_lluv_file(ifn)

    # test_str = 'testall_mp_weight_npts1'
    # test_str = 'testall_mp_weight_npts3'
    # test_str = 'qcd'

    # determine output directory and filename for radialshort data
    outdir = os.path.join(datadir, 'RadialShorts_qcd', patterntype)
    if patterntype=='IdealPattern':
        lluvtype = 'x'
    elif patterntype=='MeasPattern':
        lluvtype = 'y'
    else:
        print('Do not recognize patterntype='+patterntype+' -- must be IdealPattern or MeasPattern ') 
        return None
    # substitute RDLv(w) for RDLx(y) in filename
    rsdfn = re.sub(r'RDL[vw]', 'RDL'+lluvtype, fn)
    ofn = os.path.join(outdir, rsdfn)

    # handle empty radialmetric by outputting an empty radialshorts file
    # do not try to merge other files to fill empty radialmetric for this timestampe
    if d.size == 0:
        rsd, rsdtypes_str = generate_radialshort_array(d, types_str, header)
        rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
        rsdfooter = footer
        # 
        write_output(ofn, rsdheader, rsd, rsdfooter)
        return ofn

    # read in other data to use in averaging over time
    ixfns = find_files_to_merge(ifn, numfiles=3, sample_interval=30)
    for xfn in ixfns:
        if xfn == ifn:
            continue
        d1, types_str1, _, _ = read_lluv_file(xfn)
        if len(d.shape) == len(d1.shape) == 2:
            if (d.shape[1] == d1.shape[1]) & (types_str == types_str1):
                # if same number and order of columns as d, then append the data d
                if debug:
                    print('... ... include: %s' % xfn)
                d = numpy.vstack((d,d1))

    # (1) do threshold qc on radialmetric
    d = threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0, 5.0])
   
    # (2) do weighted averaging of good 
    xd, xtypes_str = weighted_velocities(d, types_str, numdegrees=3, weight_parameter='MP')

    # (3) require a minimum numpoints used in to form cell average
    xd = threshold_rsd_numpoints(xd, xtypes_str, numpoints=3)

    # create radialshort data, 
    rsd, rsdtypes_str = generate_radialshort_array(xd, xtypes_str, header)

    # create header from radialmetric, based on new radialshort data
    rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
    # not modifying the footer at this time
    rsdfooter = footer
    # 
    write_output(ofn, rsdheader, rsd, rsdfooter)
    return ofn

# for debugging
def _trial_qc():
    # read in the data
    ifn = os.path.join('.', 'test', 'files', 'codar_raw', \
                       'RadialMetric', 'IdealPattern', \
                       'RDLv_HATY_2013_11_05_0000.ruv')
    d, types_str, header, footer = read_lluv_file(ifn)
    # thresholding
    dall = threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0, 5.0])
    # weighting
    xd, xtypes_str = weighted_velocities(dall, types_str, numdegrees=3, weight_parameter='MP')
    
    # create radialshort data, first generate array then fill it
    rsd, rsdtypes_str = generate_radialshort_array(xd, xtypes_str, header)

    # badflag not enough numpoints
    rsd = threshold_rsd_numpoints(rsd, rsdtypes_str, numpoints=1)

    # modify header from radialmetric, based on new radialshort data
    rsdheader = generate_radialshort_header(rsd, rsdtypes_str, header)
    # not modifying the footer at this time
    rsdfooter = footer
    # output the radialshort data specified location
    ofn = os.path.join('.', 'test', 'files', 'test_output.txt')
    write_output(ofn, rsdheader, rsd, rsdfooter)
    
    #
    rsc = get_columns(rsdtypes_str)
    xc = get_columns(xtypes_str)

if __name__ == '__main__':
    # 
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    # datadir = '/Users/codar/Documents/reprocessing_2015/Reprocess_HATY_70_35/'
    # patterntype = 'MeasPattern' 
    # patterntype = 'IdealPattern' 
    # try:
    batch_qc(datadir, patterntype)
    # except:
    #     pass
    
