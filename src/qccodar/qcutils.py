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
import copy
import numpy
numpy.set_printoptions(suppress=True)
from qccodar3.codarutils import *

debug = 1


def threshold_qc_doa_peak_power(r, threshold=5.0):
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

def threshold_qc_doa_half_power_width(r, threshold=50.0):
#def threshold_qc_doa_half_power_width(d, types_str, threshold=50.0):
    """Bad Flag DOA 1/2 Power Width (degress) greater than threshold value (default 50.0 degrees).

    Flags any direction of arrival (DOA) 1/2 Power width (degress)
    that is wider than the input threshold value (default 50.0
    degrees).  Depending on the value of MSEL (1, 2, or 3), MSW1,
    MDW1, or MDW2 columns are evaluated.  Returns modified matrix with
    VFLG column the only changed values.

    """


    r1 = copy.deepcopy(r)
    havenan = numpy.isnan(r1.data['MSW1']) | numpy.isnan(r1.data['MDW1']) | numpy.isnan(r1.data['MDW2'])
    bad = (r1.data['MSEL']==1) & (r1.data['MSW1']>float(threshold))| \
          ((r1.data['MSEL']==2) & (r1.data['MDW1']>float(threshold))) | \
          ((r1.data['MSEL']==3) & (r1.data['MDW2']>float(threshold))) | havenan
    r1.data.loc[bad, 'VFLG'] = r1.data.loc[bad, 'VFLG'] + (1<<2)
    return r1

def threshold_qc_monopole_snr(r, threshold=5.0):
#def threshold_qc_monopole_snr(d, types_str, threshold=5.0):
    """Bad flag any SNR on monopole (dB)  less than threshold value (default 5.0 dB).

    Flags any signal-to-noise ratio (SNR) on monopole (dB) that falls
    below the input threshold value (default 5.0 dB).  No dependency on MSEL selections.

    """
    r1 = copy.deepcopy(r)
    bad = r1.data['MA3S'] < float(threshold)
    r1.data.loc[bad, 'VFLG'] = r1.data.loc[bad, 'VFLG'] + (1<<3)
    return r1

def threshold_qc_loop_snr(r, threshold=5.0):
#def threshold_qc_loop_snr(d, types_str, threshold=5.0):
    """Bad flag if both loop SNR are less than threshold value (default 5.0 dB).

    Flags if signal-to-noise ratio (SNR) (dB) on loop1 AND on loop2 falls
    below the input threshold value (default 5.0 dB). No dependency on MSEL selections.

    """
    r1 = copy.deepcopy(r)
    bad = (r1.data['MA1S']<float(threshold)) & (r1.data['MA2S']<float(threshold))
    r1.data.loc[bad, 'VFLG'] = r1.data.loc[bad, 'VFLG'] + (1<<3)
    return r1

def threshold_qc_all(r, qccodar_values = dict()):
#def threshold_qc_all(d, types_str, thresholds=[5.0, 50.0, 5.0, 5.0]):
    """Combine all three threshold tests

    Returns modified matrix with VFLG column only changed values.

    """

    r1 = copy.deepcopy(r)
    if r1.is_valid():

        qc_keys = qccodar_values.keys()

        # run high frequency radar qartod tests on open radial file
        if 'qc_doa_peak_power' in qc_keys:
            r1 = threshold_qc_doa_peak_power(r1, qccodar_values['qc_doa_peak_power']['doa_peak_power_min'])

        if 'qc_doa_half_power_width' in qc_keys:
            r1 = threshold_qc_doa_half_power_width(r1, qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max'] )

        if 'qc_monopole_snr' in qc_keys:
            r1 = threshold_qc_monopole_snr(r1, qccodar_values['qc_monopole_snr']['monopole_snr_min'])

        if 'qc_loop_snr' in qc_keys:
            r1 = threshold_qc_loop_snr(r1, qccodar_values['qc_loop_snr']['loop_snr_min'])

    return r1

def threshold_rsd_numpoints(rs, radialshort_velocity_count_min=1):
#def threshold_rsd_numpoints(rsd, rstypes_str, numpoints=1):
    """Bad flag any radialshort data with doppler velocity count (EDVC) less than "numpoints"

    Returns modified rsd matrix with VFLG column only changed if EDVC
    count is less than or equal to numpoints.  This threshold is
    checked after weighted_velocities()

    """

    rs1 = copy.deepcopy(rs)
    bad = rs.data['EDVC']<int(radialshort_velocity_count_min)
    rs1.data.loc[bad, 'VFLG'] = rs1.data.loc[bad, 'VFLG'] + (1<<12)
    return rs1


def weighted_velocities(r, numdegrees=3, weight_parameter='MP'):
    """Calculates weighted average of radial velocities (VELO) at bearing and range.

    The weighted average of velocities found at given range and
    bearing based on weight_parameter.

    Paramters
    ---------
    d : ndarray
        The data from LLUV file(s). 
    types_str : string 
        The 'TableColumnTypes' string header of LLUV file(s) provide keys for each column.
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
    xcols = ['VFLG', 'SPRC', 'BEAR', 'VELO', 'ESPC', 'MAXV', 'MINV', 'EDVC', 'ERSC']

    r1 = copy.deepcopy(r)
    d = copy.deepcopy(r.data)
    offset = ((numdegrees-1)/2)

    r1.data = r1.data.drop_duplicates(subset=['SPRC','BEAR','VFLG'])
    # return only rows that have VFLG==0 (0 == good, >0 bad) so only get good data
    good = r1.data['VFLG'] == 0
    ud = r1.data.loc[good, :]
    if ud.size == 0:
        r.data = numpy.array([])
        return r.data

    allbearings = numpy.unique(ud['BEAR'])
    allranges = numpy.unique(ud['SPRC'])
    ud = numpy.array([[r, b] for r in allranges for b in allbearings])
    nrows, _ = ud.shape
    xd = pd.DataFrame(columns = xcols ,index=numpy.arange(nrows))


    for irow, cell in enumerate(ud):
        rngcell, bearing = cell[0:2]
        # numpy.where() returns a tuple for condition so use numpy.where()[0]
        # also VFLG must equal 0 (0 == good, >0 bad) so only get good data

        xrow = numpy.where((d['SPRC']==rngcell) & \
                           (d['BEAR']>=bearing-offset) & \
                           (d['BEAR']<=bearing+offset) & \
                           (d['VFLG']==0))[0]

        # If no row matches rngcell AND bearing, then no VELO data, skip to next bearing
        if xrow.size == 0: 
            continue

        #xcol = numpy.array([['VELO'], ['MSEL'], ['MSP1'], ['MDP1'], ['MDP2'], ['MA3S']])

        #a = d[numpy.ix_(xrow, xcol)].copy()
        ao = copy.deepcopy(d)
        a = ao.loc[xrow, ['VELO', 'MSEL', 'MSP1', 'MDP1', 'MDP2', 'MA3S']]
        # if xrow.size == edvc:
        VELO = a['VELO']  # all radial velocities found in cell
        SNR3 = a['MA3S']  # SNR on monopole for each velocity
        if weight_parameter.upper() == 'MP':
            # Create array to hold each Music Power (based on MSEL)
            MP = numpy.array(numpy.ones(VELO.shape)*numpy.nan) 
            # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
            mselcol = ['','MSP1','MDP1','MDP2']
            for msel in [1, 2, 3]:
                which = a['MSEL'] == msel
                MP[which,] = a.loc[which, mselcol[msel]]
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
        xd.loc[irow,['VFLG']] = 0
        xd.loc[irow,['SPRC']] = rngcell
        xd.loc[irow,['BEAR']] = bearing
        xd.loc[irow,['VELO']] = velo
        # other stat output
        xd.loc[irow, ['ESPC']] = VELO.values.std()  # ESPC
        xd.loc[irow,['MAXV']] = VELO.max() # MAXV
        xd.loc[irow,['MINV']] = VELO.min() # MINV
        # (EDVC and ERSC are the same in this subroutine's context)
        xd.loc[irow,['EDVC']] = VELO.size # EDVC Velocity Count
        xd.loc[irow,['ERSC']] = VELO.size # ERSC Spatial Count
                
    # delete extra lines (nan) not filled above
    xd.dropna(axis=0, how='all', inplace=True)
    # ESPC had several nan's, replace these and other nan before returning?
    xd.fillna(999.000, inplace=True)

    return xd

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

def find_files_to_concatenate(ifn, numfiles=3, sample_interval=30):
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
    rdlstr = re.match(r'RDL[vwxy]', os.path.basename(ifn)).group()
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

def add_short_metadata(r,qccodar_values):
    # create a list of QCD metadata to write in the radial file header
    r.metadata['QCD'] = []

    r.metadata['QCD'].append('QCDReference: Radial Metric Quality Control Reference: doi:10.1175/JTECH-D-16-0203.1')

    # Not sure why but I had to define simple variables to put in {} for the formatted output.
    doa_peak_power_min = qccodar_values['qc_doa_peak_power']['doa_peak_power_min']
    r.metadata['QCD'].append((
        f'QCDTest: qc_doa_peak_power threshold ['
        f'doa_peak_power_min = {doa_peak_power_min}(dB)] '
    ))
    doa_half_power_width_max = qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max']
    r.metadata['QCD'].append((
        f'QCDTest: qc_doa_half_power_width threshold ['
        f'doa_half_power_width_max = {doa_half_power_width_max}(degrees)] '
    ))
    monopole_snr_min = qccodar_values['qc_monopole_snr']['monopole_snr_min']
    r.metadata['QCD'].append((
        f'QCDTest: qc_monopole_snr threshold ['
        f'monopole_snr_min = {monopole_snr_min}(dB)] '
    ))
    loop_snr_min = qccodar_values['qc_loop_snr']['loop_snr_min']
    r.metadata['QCD'].append((
        f'QCDTest: qc_loop_snr threshold ['
        f'loop_snr_min = {loop_snr_min}(dB)] '
    ))
    radialshort_velocity_count_min = qccodar_values['qc_radialshort_velocity_count']['radialshort_velocity_count_min']
    r.metadata['QCD'].append((
        f'QCDTest: qc_radialshort_velocity_count ['
        f'radialshort_velocity_count_min = {radialshort_velocity_count_min}(velocities)] '
    ))
    numfiles = qccodar_values['metric_concatenation']['numfiles']
    sample_interval = qccodar_values['metric_concatenation']['sample_interval']
    r.metadata['QCD'].append((
        f'QCDSettings: metric_concatenation ['
        f'numfiles = {numfiles}(files), '
        f'sample_interval = {sample_interval}(minutes)] '
    ))
    numdegrees = qccodar_values['weighted_shorts']['numdegrees']
    weight_parameter = qccodar_values['weighted_shorts']['weight_parameter']
    table_type = qccodar_values['weighted_shorts']['table_type']
    r.metadata['QCD'].append((
        f'QCDSettings: weighted_shorts ['
        f'numdegrees = {numdegrees}(degrees), '
        f'weight_parameter = {weight_parameter}, '
        f'table_type = {table_type}]'
    ))

    return r



def do_qc(datadir, fn, patterntype, qccodar_values = dict()):
    """ Do qc and then average over 3 sample_intervals (time), 3 degrees of bearing.
    """
    # read in the data
    rmfoldername = get_radialmetric_foldername(datadir)
    ifn = os.path.join(datadir, rmfoldername, patterntype, fn)
    rs = read_lluv_file(ifn)

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
    if rs.data.size == 0:
        rs = generate_radialshort(rs, **qccodar_values['weighted_shorts'])
        write_output(rs, ofn)
        return ofn

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

    # (2) do weighted averaging of good, this is done when the radial short is created
    rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])

    # (3) require a minimum numpoints used in to form cell average
    rsx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])

    # create a list of QCD metadata to write in the radial file header
    rsx = add_short_metadata(rsx,qccodar_values)

    write_output(rsx,ofn)
    return ofn

# for debugging
def _trial_qc():
    # read in the data
    ifn = os.path.join('.', 'test', 'files', 'codar_raw', \
                       'RadialMetric', 'IdealPattern', \
                       'RDLv_HATY_2013_11_05_0000.ruv')

    r = read_lluv_file(ifn)

    qccodar_values = dict(
        qc_doa_half_power_width=dict(doa_half_power_width_max=50.0),
        qc_doa_peak_power=dict(doa_peak_power_min=5.0),
        qc_monopole_snr=dict(monopole_snr_min=5.0),
        qc_loop_snr=dict(loop_snr_min=5.0),
        qc_radialshort_velocity_count=dict(radialshort_velocity_count_min=1.0),
        weighted_shorts=dict(numdegrees=3,weight_parameter='MP', table_type='LLUV RDL7'),
        merge=dict(css_interval_minutes=30.0,number_of_css=5.0)
    )

    # thresholding
    rall = threshold_qc_all(r, qccodar_values)

    # weighting, generating radial short
    rsd = generate_radialshort(rall, **qccodar_values['weighted_shorts'] )

    # badflag not enough numpoints
    rs = threshold_rsd_numpoints(rsd, qccodar_values['qc_radialshort_velocity_count']['qc_radialshort_velocity_count_min'])

    ofn = os.path.join('.', 'test', 'files', 'test_output.txt')
    write_output(rs, ofn)


if __name__ == '__main__':
    # 
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    # datadir = '/Users/codar/Documents/reprocessing_2015/Reprocess_HATY_70_35/'
    # patterntype = 'IdealPattern' 
    # try:
    #batch_qc(datadir, patterntype) #can't find batch_qc defined anywhere!
    # except:
    #     pass
