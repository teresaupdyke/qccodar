#!/usr/bin/env python
#
""" CODAR Utilities
"""
import sys
import os
import re
import copy
import fnmatch
import datetime
import geopy
import geopy.distance
import numpy
numpy.set_printoptions(suppress=True)
import pandas as pd
from hfradarpy.radials import Radial

debug = 1

def write_output(r,ofn,export_type='radial'):
    """Write LLUV file using HFRadarPy toolbox """
    r.export(ofn, export_type)

def read_lluv_file(ifn):
    """Reads LLUV file using HFRadarPy toolbox"""
    if os.path.exists(ifn):
        r = Radial(ifn, mask_over_land=False)
    else:
        print('File does not exist: '+ ifn)
        raise IOError('Error opening %s' % ifn)
    return r


def get_radialmetric_foldername(datadir, pattern='?adial*etric*'):
    """ Slightly different variances in the name of the folder for RadialMetric[s] data"""
    fns = os.listdir(datadir)
    mfns = fnmatch.filter(fns, pattern)
    if mfns:
        foldername = mfns[0]
    else:
        foldername = ''
    return foldername

def generate_radialshort(r, table_type='LLUV RDL7', numdegrees=3, weight_parameter='MP'):
    """Generates radialshort (rsd) data array.

    This function generates radialshort data (rsd) array based on data
    from weighted average data (xd).  It is the output CSV data table
    (middle) of CODAR LLUV format with header, middle, and footer.

    If xd is empty, return empty rsd for writing out empty LLUV files.

    Each LATD, LOND is computed using Vincenty's algorithm for
    destination along a BEAR (NCW) and a distance of RNGE (km) from
    point of origin (lat1,lon1).  Vincenty's GC is great circle
    distance on an WGS-84 ellipsoid model between two points. This
    requires having the CODAR Site location, and range resolution from
    the header data or config data.

    Parameters
    ----------
    xd : ndarray
       QC'd and weighted average radials for each range, bearing where data were found
    table_type : string

    Returns
    -------
    rs: Radial object
    """
    from qccodar3.qcutils import weighted_velocities

    if table_type == 'LLUV RDL7':

        rs = Radial('/Users/teresa/Desktop/Codar_Files/Software_Projects/git/qccodar3/src/qccodar3/file_formats/radialshort_LLUV_RDL7.ruv',empty_radial=True)
    else:
        print('generate_radial_array() : Unrecognized table_type "%s"' % (table_type,))
        return numpy.array([]), ''

    # copy over the file information, header, name, tables
    rs.metadata = r.metadata
    #  '% PatternMethod: 1 PatternVectors' should be added but I'm not sure how to
    #   insert at a specific position in Ordered dictionary
    rs.diagnostics_radial = r.diagnostics_radial
    rs.diagnostics_hardware = r.diagnostics_hardware
    # the range information table changes in the processing from radial metric to radial short
    # but I have simply copied over the metric version for now
    rs.range_information = r.range_information

    for key in rs._tables.keys():
        table = rs._tables[key]
        if 'LLUV' in table['TableType']:
            rs._tables[key]['data'] = rs.data
            rs._tables[key]['TableRows'] = str(rs.data.shape[0])
        elif 'rads' in table['TableType']:
            rs._tables[key]['data'] = rs.diagnostics_radial
            rs._tables[key]['TableRows'] = str(rs.diagnostics_radial.shape[0])
        elif 'rcvr' in table['TableType']:
            rs._tables[key]['data'] = rs.diagnostics_hardware
            rs._tables[key]['TableRows'] = str(rs.diagnostics_hardware.shape[0])
        elif 'RINF' in table['TableType']:
            rs._tables[key]['data'] = rs.range_information
            rs._tables[key]['TableRows'] = str(rs.range_information.shape[0])

    #rs.file_path = r.file_path
    fn = r.file_name.replace('RDLv', 'RDLx')
    fn = fn.replace('RDLw', 'RDLy')
    rs.file_name = fn
    #rs.full_file = r.full_file
    rs.is_wera = False
    rs._iscorrupt = False
    rs.time = datetime.datetime(*[int(s) for s in r.metadata['TimeStamp'].split()])

    origin = r.metadata['Origin']
    lat1, lon1 = [float(x) for x in origin.split()]
    range_resolution = float(r.metadata['RangeResolutionKMeters'])

    xd = weighted_velocities(r, numdegrees, weight_parameter)

    if xd.size == 0:
        return rs

    rs.data['VFLG'] = xd['VFLG']
    rs.data['SPRC'] = xd['SPRC']
    rs.data['BEAR'] = xd['BEAR']
    rs.data['VELO'] = xd['VELO']
    rs.data['ESPC'] = xd['ESPC']
    rs.data['MAXV'] = xd['MAXV']
    rs.data['MINV'] = xd['MINV']
    rs.data['EDVC'] = xd['EDVC']
    rs.data['ERSC'] = xd['ERSC']

    ############################
    # computations for filling in other columns of radialshort data
    ############################

    # create HEAD column based on BEAR+180
    rs.data['HEAD'] = numpy.mod(rs.data['BEAR'] + 180., 360.)
    # compute velocity components
    (rs.data['VELU'], rs.data['VELV']) = compass2uv(rs.data['VELO'], rs.data['HEAD'])

    rs.data['RNGE'] = range_resolution * rs.data['SPRC']
    #
    # Vincenty Great Circle destination point (LATD, LOND) based on rnge, bear from site origin
    origin = geopy.Point(lat1, lon1)
    pts = numpy.array([geopy.distance.geodesic(kilometers=r).destination(origin, b)[0:2] for (r, b) in zip(rs.data['RNGE'], rs.data['BEAR'])])
    rs.data['LATD'], rs.data['LOND'] = pts[:, 0], pts[:, 1]

    (rs.data['XDST'], rs.data['YDST']) = compass2uv(rs.data['RNGE'], rs.data['BEAR'])

    return rs

def cell_intersect(rngbear1, rngbear2):
    """ Return rows that match range and bearing data between the two nx2 matrices.

    Parameters
    ----------
    rngbear1 : nx2 array
       The first array of range and bearings that define each cell. 
       The first column is range cells.  The second column are the bearings.
    rngbear2 : nx2 array
       The second array of range and bearings that define each cell. 
       The first column is range cells.  The second column are the bearings.
    
    Return
    ------
    (rows1, rows2) : tuple of nx1 arrays 
       The rows where range and bearing cell are the same between two input matrices.
    """
    rows1 = []; rows2=[]
    for irow, cell in enumerate(rngbear1):
        rngcell, bearing = cell
        xrow = numpy.where( (rngbear2[:, 0] == rngcell) & \
                            (rngbear2[:, 1] == bearing) )[0]
        if xrow.size == 0:
            continue
        rows1.append(irow); rows2.append(xrow)
    rows1 = numpy.squeeze(rows1)
    rows2 = numpy.squeeze(rows2)
    return (rows1, rows2)


def compass2uv(wmag, wdir):
    """ Vector conversion from mag and direction (wmag,wdir) to x,y
    vector components (u,v)

    If inputs are lists, it is cast into an arrays.

    Parameters
    ----------
    wmag : array-like, same size as wdir
       The magnitude of the vector. 
    wdir : array-like, same size as wmag
       The compass direction, Clockwise from y-axis or North equals 0/360 deg

    Returns
    -------
    (u,v) : tuple of array-like u and v vectors the same size and shape as inputs.
       The x,y vector components. 

    >>> compass2uv(1.0, 0.0)
    (0.0, 1.0)
    >>> compass2uv([1., 1., 1., 1.], [0., 90., 180., 270.])
    (array([ 0.,  1.,  0., -1.]), array([ 1.,  0., -1., -0.]))
    
    """
    # calculate horizontal vector components (u,v) from magnitude and compass direction
    # cast the inputs into numpy.array
    wmag = numpy.array(wmag)
    wdir = numpy.array(wdir)
    assert wmag.shape == wdir.shape, 'wmag and wdir must be same size and shape'
    
    r = numpy.pi/180.
    u = wmag*numpy.sin(wdir*r)
    v = wmag*numpy.cos(wdir*r)
    return (u,v)

def run_LLUVMerger(datadir, fn, patterntype, css_interval_minutes=30.0, number_of_css=5.0, debug=2, diag = '4'):
    """ Run CODAR's LLUVMerger app in subprocess """

    import subprocess
    from .qcutils import filt_datetime

    ifn = os.path.join(datadir, 'RadialShorts_qcd', patterntype, fn)
    outdir = os.path.join(datadir, 'Radials_qcd', patterntype)

    if patterntype=='IdealPattern':
        lluvtype = 'i'
    elif patterntype=='MeasPattern':
        lluvtype = 'm'
    else:
        print('Do not recognize patterntype='+patterntype+' -- must be IdealPattern or MeasPattern ') 
        return

    rs_output_interval = datetime.timedelta(minutes=css_interval_minutes)
    rs_num = number_of_css
    # (merge average 5*30 min = 150 min or 2.5 hours)
    span_hrs = rs_num * (rs_output_interval.seconds/3600.) # hours, 2.5 hours
    span_hrs_str = '%f' % span_hrs # '2.5000'
    # if ifn (source file) is on the hour (00 min) expected time 
    expected_timedelta = datetime.timedelta(minutes=60)
    #
    #       22:30
    # 5\    23:00
    # 4 |   23:30
    # 3 |-- 00:00 <--expected time for merger of 5 files is 60 min behind source 
    # 2 |   00:30
    # 1 /   01:00 <-- source file time
    #       01:30
    
    # ordered list of args, order of some options is important,
    # e.g. -span and -startwith before -source
    args = ['/Codar/SeaSonde/Apps/Bin/LLUVMerger',
            '-span='+span_hrs_str,
            '-lluvtype='+lluvtype, 
            '-angres=5',
            '-angalign=2',
            '-angmethod=short',
            '-method=average',
            '-minvect=2',
            '-velcount',
            '-diag='+diag,
            '-source='+ifn,
            '-output='+outdir]

    if debug>=2:
        print(' '.join(args))

    # NEW WAY to capture both stderr and stdout
    p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    p.wait()
    (stdout_content, stderr_content) = p.communicate()

    # PIPE is a readable file-like object created when Popen attribute is defined.
    # communicate() sends stdin (if defined) and waits for exit. Next it gets stdout
    # and sterr, then it closes all three (stdin, stdout, stderr).
    # p.communicate() returns tuple (stdout_content, stderr_content). 

    # print error and return
    if stderr_content:
        print("Error running LLUVMerger()") 
        print(stderr_content)
        return

    # Check merged file has correct time-of-merge
    if stdout_content:
        # stdout_content should be
        # Running LLUVMerger 1.4.1
        # Merging 5 Sources...
        # Merging done.
        # MergedFile: "/Codar/SeaSonde/Data/Radials_qcd/IdealPattern/RDLi_HATY_2013_11_05_0000.ruv"
        if debug>=2:
            print(stdout_content.strip())

        lines = stdout_content.decode('utf-8').split('\n')
        # get line with MergedFile: path and filename from stdout_content
        line = [x for x in lines if 'MergedFile:' in x][0]
        if debug>=2:
            print(line)
        # mfn -- extract full path and file name of merged file 
        # /Codar/SeaSonde/Data/Radials_qcd/IdealPattern/RDLi_HATY_2013_11_05_0000.ruv
        m = re.match(r'^MergedFile:\s*\"(.*)\"$', line)
        if m: 
            mfn = m.groups()[0]
        else:
            mfn = ''
        ofn = mfn

        if not mfn:
            if debug>=2:
                print('No merged file found')
            return ofn

        # expected datetime
        dt_expected = filt_datetime(os.path.basename(ifn)) - expected_timedelta
        # actual datetime
        dt_actual = filt_datetime(os.path.basename(mfn))

        if dt_expected != dt_actual:
            # rename the file to dt_expected
            # >>> dt_actual.strftime('%Y_%m_%d_%H%M')
            # '2013_11_05_0115'
            outdir = os.path.dirname(mfn)
            (basename,suffix) = os.path.splitext(os.path.basename(mfn))
            p = re.match(r'^(.*)\d{4}_\d{2}_\d{2}_\d{4}', basename)
            if p:
                prefix = p.groups()[0]
            else:
                prefix = ''
            # outdir = '/Codar/SeaSonde/Data/Radials_qcd/IdealPattern'
            # prefix = 'RDLi_HATY_'
            # >>> dt_expected.strftime('%Y_%m_%d_%H%M')
            # '2013_11_05_0100'
            # suffix = '.ruv'
            newfn = os.path.join(outdir, prefix+dt_expected.strftime('%Y_%m_%d_%H%M')+suffix)
            if debug>=2:
                print('MergedFile renamed: "%s"' % newfn)
            os.rename(mfn, newfn)
            ofn = newfn
        return ofn


def check_headertime(fullfn):
    """ Loads merged file, checks if time in file name matches time in header
            and corrects if times do not match """
    from qccodar3.qcutils import filt_datetime

    r = read_lluv_file(fullfn)
    fn = os.path.split(fullfn)[1]
    fn_time = filt_datetime(r.file_name)
    if fn == r.file_name:
        header_time = filt_datetime(r.metadata['TimeStamp'])
        if fn_time != header_time:

        # if header timestamp is within 30 minutes of the file name time then update
        # header to match file name
            if header_time - fn_time < datetime.timedelta(seconds=1801):
                r.metadata['TimeStamp'] = fn_time.strftime("%Y %m %d %H %M %S")
                print('TimeStamp in header changed to match the time indicated in file name.')

        write_output(r, fullfn, export_type='radial')
        return

    else:
        print('File names do not match!')
        return

def add_diagnostic_tables(fullfn, shortpath):
    """ Uses diagnostic information in the radialshorts to build radial and hardware diagnostic tables
      for the merged radial file """

    from collections import OrderedDict

    r = read_lluv_file(fullfn)

    for key in r._tables.keys():
        table = r._tables[key]
        if 'MRGS' in table['TableType']:
            filelist = table['data']['PATH']

    firstshortfn = os.path.join(shortpath, filelist[0])
    rs = read_lluv_file(firstshortfn)
    rdtdata = rs.diagnostics_radial
    hdtdata = rs.diagnostics_hardware

    for shortfile in filelist:
        shortfullfile = os.path.join(shortpath, shortfile)
        print('... input: %s' % shortfullfile)

        rs= read_lluv_file(shortfullfile)
        for key in rs._tables.keys():
            table = rs._tables[key]
            if 'rads' in table['TableType']:
                rdt = table['data']
                rdtdata = pd.concat([rdtdata, rdt], ignore_index=True)
            if 'rcvr' in table['TableType']:
                hdt = table['data']
                hdtdata = pd.concat([hdtdata, hdt], ignore_index=True)

    # remove duplicates
    hdtdata = hdtdata.drop_duplicates(subset=['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'])
    rdtdata = rdtdata.drop_duplicates(subset=['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'])

    # adjust rdt for seconds from start
    # adjust hdt for minutes from start
    # create a new _tables dictionary for the merged file
    od = OrderedDict()

    for key in r._tables.keys():
        table = r._tables[key]
        if 'LLUV' in table['TableType']:
            od['1'] = r._tables[key]

    for key2 in rs._tables.keys():
        table2 = rs._tables[key2]
        if 'rads' in table2['TableType']:
            od['2'] = rs._tables[key2]
            od['2']['data'] = rdtdata
            od['2']['TableRows'] = rdtdata.shape[0]
        elif 'rcvr' in table2['TableType']:
            od['3'] = rs._tables[key2]
            od['3']['data'] = hdtdata
            od['3']['TableRows'] = hdtdata.shape[0]

    r._tables = od
    write_output(r, fullfn, export_type='radial')
    return

def do_merge(datadir, fn, pattern, qccodar_values):
    """ Calls LLUVMerger to create the merged file, checks that output filename includes the expected
        time, (if not the script corrects the file name and ensures that time in filename and header
        match), adds diagnostic tables to the end of the file """

    ofn = run_LLUVMerger(datadir, fn, pattern, diag='4', **qccodar_values['merge'])
    r = read_lluv_file(ofn)

    for key in r._tables.keys():
        table = r._tables[key]
        if 'MRGS' in table['TableType']:
            filelist = table['data']['PATH']

    if filelist.shape[0] < 2:
        # delete the file from the computer
        print('Not enough short files available for the merge.')
        return ofn
    else:
        shortpath = os.path.join(datadir, 'RadialShorts_qcd',pattern)
        check_headertime(ofn)
        add_diagnostic_tables(ofn, shortpath)
        return ofn