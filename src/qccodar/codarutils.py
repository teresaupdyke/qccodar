#!/usr/bin/env python
#
""" CODAR Utilities
"""
import os
import re
import fnmatch
import datetime
import geopy
import geopy.distance
import numpy
from pathlib import Path
numpy.set_printoptions(suppress=True)
import pandas as pd
from hfradarpy.radials import Radial

debug = 1

def write_output(r,ofn):
    """Write LLUV file using HFRadarPy toolbox """

    if not r._iscorrupt:
        r.to_ruv(ofn, validate=False, overwrite=True)
    else:
        print('Corrupt radial. File not written %s' % ofn)

def read_lluv_file(ifn):
    """Reads LLUV file using HFRadarPy toolbox"""
    if os.path.exists(ifn):
        r = Radial(ifn)
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
    from qccodar.qcutils import weighted_velocities

    if table_type == 'LLUV RDL7':
        formatfile = Path(__file__).parent.resolve() / 'file_formats' / 'radialshort_LLUV_RDL7.ruv'
        rs = Radial(formatfile,empty_radial=True)
    else:
        print('generate_radialshort() : Unrecognized table_type "%s"' % (table_type,))
        return numpy.array([]), ''

    # copy over the file information, header, name, tables
    rs.metadata = r.metadata
    #  '% PatternMethod: 1 PatternVectors' should be added but I'm not sure how to
    #   insert at a specific position in Ordered dictionary
    if hasattr(r,'diagnostics_radial'):
	    rs.diagnostics_radial = r.diagnostics_radial
    if hasattr(r,'diagnostics_hardware'):
	    rs.diagnostics_hardware = r.diagnostics_hardware
    # the range information table changes in the processing from radial metric to radial short
    # but I have simply copied over the metric version for now
    if hasattr(r,'range_information'):
	    rs.range_information = r.range_information

    # remove table in rs that is not in r
    r_tts = [r._tables[key]['TableType'].split(' ')[0] for key in r._tables.keys()]
    rs_tts = [rs._tables[key]['TableType'].split(' ')[0] for key in rs._tables.keys()]
    rs_keys= [key for key in rs._tables.keys()]
    to_be_removed = [tt for tt in rs_tts if tt not in r_tts]

    for tt in to_be_removed:
        idx = rs_tts.index(tt)
        rs._tables.pop(rs_keys[idx])
    
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

    if r.data.size == 0:
        return rs

    xd = weighted_velocities(r, numdegrees, weight_parameter)

    if xd.size == 0:
        return rs


    for key in rs._tables.keys():
        table = rs._tables[key]
        if 'LLUV' in table['TableType']:
            rs._tables[key]['TableRows'] = str(xd.shape[0])

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

def get_columns(types_str):
    # use dict to store column label and it's column number for indexing numpy array
    #c = col.defaultdict(int)
    c = {}
    column_labels = types_str.strip().split(' ')
    m = re.findall(r'\w{4}', types_str)
    for label in column_labels:
        c[label]=m.index(label) # c['VFLG']=4
    return c
    
def unique_rows(a):
    # http://stackoverflow.com/questions/8560440/removing-duplicate-columns-and-rows-from-a-numpy-2d-array?lq=1
    a = numpy.ascontiguousarray(a)
    unique_a = numpy.unique(a.view([('', a.dtype)]*a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

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

def run_LLUVMerger(datadir, fn, patterntype, qv_merge, debug=2):
# def run_LLUVMerger(datadir, fn, patterntype, css_interval_minutes=30.0, number_of_css=5.0, shorts_minute_filter='*00',debug=2, diag = '4'):
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

    # set defaults for required elements of LLUVMerger to generate Radials, in case missing from configfile
    css_interval_minutes=30.0
    number_of_css=5.0
    shorts_minute_filter='*00'
    method='average'
    minvect = '2'
    diag = '4'

    # items that are configurable in qccodar_values['merge'] passed as qv_merge
    if 'css_interval_minutes' in qv_merge.keys():
        css_interval_minutes=qv_merge['css_interval_minutes']
    if 'number_of_css' in qv_merge.keys():
        number_of_css=qv_merge['number_of_css']
    if 'shorts_minute_filter' in qv_merge.keys():
        shorts_minute_filter=qv_merge['shorts_minute_filter']
    if 'method' in qv_merge.keys():
        method=qv_merge['method']
    if 'minvect' in qv_merge.keys():
        minvect=qv_merge['minvect']
    if 'angalign' in qv_merge.keys():
        angalign = qv_merge['angalign']
    if 'diag' in qv_merge.keys():
        diag = qv_merge['diag']
    if 'reference' in qv_merge.keys():
        reference = qv_merge['reference']

    
    rs_output_interval = datetime.timedelta(minutes=css_interval_minutes)
    rs_num = number_of_css
    # (merge average 5*30 min = 150 min or 2.5 hours)
    span_hrs = rs_num * (rs_output_interval.seconds/3600.) # hours, 2.5 hours
    span_hrs_str = '%f' % span_hrs # '2.5000'

    minutes_behind_source = (number_of_css // 2) * css_interval_minutes
    expected_timedelta = datetime.timedelta(minutes=minutes_behind_source)

    # if ifn (source file) is on the hour (00 min) expected time
    #
    #       22:30
    # 5\    23:00
    # 4 |   23:30
    # 3 |-- 00:00 <--expected time for merge of 5 files with 30 minute intervals is 60 min behind source
    # 2 |   00:30
    # 1 /   01:00 <-- source file time
    #       01:30

    #       23:20
    # 7 \   23:30
    # 6 |   23:40
    # 5 |   23:50
    # 4 |-- 00:00 <--expected time for merge of 7 files with 10 minute intervals is 30 min behind source
    # 3 |   00:10
    # 2 |   00:20
    # 1 /   00:30 <-- source file time
    #       00:40

    # ordered list of args, order of some options is important,
    # -span and -startwith must be before -source
    
    # Must have -angres=5 and -angmethod=short to merge RadialShorts to Radials 
    # so keeping these static here and are not changeable in configfile
    args = ['/Codar/SeaSonde/Apps/Bin/LLUVMerger',
            '-span='+span_hrs_str,
            '-lluvtype='+lluvtype, 
            '-angres=5',
            '-angmethod=short',
            '-method='+method,
            '-minvect='+minvect,
            '-velcount',
            '-diag='+diag,
            '-source='+ifn,
            '-output='+outdir]

    # if any of these are set within this function as local variable then add it to args for running LLUVMerger
    if 'reference' in locals():
        args.append('-reference='+reference)
    if 'angalign' in locals():
        args.append('-angalign='+angalign)

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
        # stdout_content should look similar to
        # Running LLUVMerger 1.4.1
        # Merging 5 Sources...
        # Merging done.
        # MergedFile: "/Codar/SeaSonde/Data/Radials_qcd/IdealPattern/RDLi_HATY_2013_11_05_0000.ruv"
        if debug>=2:
            print(stdout_content.strip())

        lines = stdout_content.decode('utf-8').split('\n')
        # get line with MergedFile: path and filename from stdout_content
        if lines[0] == 'No source files found':
            line = lines[0]
        else:
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


def check_headertime(r,fullfn):
    """ Loads merged file, checks if time in file name matches time in header
            and corrects if times do not match """
    from qccodar.qcutils import filt_datetime

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

        return r

    else:
        print('File names do not match!')
        return

def add_diagnostic_tables(r, shortpath):
    """ Uses diagnostic information in the radialshorts to build radial and hardware diagnostic tables
      for the merged radial file """

    from collections import OrderedDict

    for key in r._tables.keys():
        table = r._tables[key]
        if 'MRGS' in table['TableType']:
            filelist = table['data']['PATH']

    firstshortfn = os.path.join(shortpath, filelist[0])
    rs = read_lluv_file(firstshortfn)
    if hasattr(rs,'diagnostics_radial'):
        rdtdata = rs.diagnostics_radial
    else:
        # use empty dataframe
        rdtdata = pd.DataFrame()
    if hasattr(rs,'diagnostics_hardware'):
        hdtdata = rs.diagnostics_hardware
    else:
    	hdtdata = pd.DataFrame()

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
    if not rdtdata.empty:
        rdtdata = rdtdata.drop_duplicates(subset=['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'])
    if not hdtdata.empty:
        hdtdata = hdtdata.drop_duplicates(subset=['TYRS', 'TMON', 'TDAY', 'THRS', 'TMIN', 'TSEC'])

    # adjust rdt for seconds from start
    # adjust hdt for minutes from start
    # create a new _tables dictionary for the merged file
    od = OrderedDict()

    for key in r._tables.keys():
        table = r._tables[key]
        if 'LLUV' in table['TableType']:
            od[1] = r._tables[key]

    for key2 in rs._tables.keys():
        table2 = rs._tables[key2]
        if 'rads' in table2['TableType']:
            od[2] = rs._tables[key2]
            od[2]['data'] = rdtdata
            od[2]['TableRows'] = rdtdata.shape[0]
            r.diagnostics_radial = rdtdata
        elif 'rcvr' in table2['TableType']:
            od[3] = rs._tables[key2]
            od[3]['data'] = hdtdata
            od[3]['TableRows'] = hdtdata.shape[0]
            r.diagnostics_hardware = hdtdata

    r._tables = od

    return r


def fix_empty_radial(r, table_type='LLUV RDL9'):
    """Fixes an empty radial (r) Radial object that be can be written out using 
    r.to_ruv() to avoid an error when there is missing table information.

    return r (Radial object) for writing out empty LLUV files.

    Parameters
    ----------
    r : Radial object read from LLUVMerger empty radial file 
    table_type = 'LLUV RDL9'  ensures correct file format used as template
       

    Returns
    -------
    r: Radial object, fixed by using re._tables[1]
    """

    if table_type == 'LLUV RDL9':
        formatfile = Path(__file__).parent.resolve() / 'file_formats' / 'radial_LLUV_RDL9.ruv'
        # formatfile = './test_empty_Radials_qcd/RDLi_HATY_2020_10_08_2000.ruv'
        re = Radial(formatfile,empty_radial=True)
    else:
        print('fix_empty_radial() : Unrecognized table_type "%s"' % (table_type,))
        return r

    ## it would be easy if we knew for sure they would always be key==1 for both
    ## r._tables[1] = re._tables[1]
    ## r.data = re._tables[1]['data']

    # all this to make sure index to r._tables and re._tables are getting 'LLUV' tables from both,
    # in case they are different keys

    # get all the table types
    r_tts = [r._tables[key]['TableType'].split(' ')[0] for key in r._tables.keys()]
    re_tts = [re._tables[key]['TableType'].split(' ')[0] for key in re._tables.keys()]

    # get list of keys -- r._tables is OrderedDict() and odict_keys are not subscriptable
    # so generate a list of keys
    r_keys= [key for key in r._tables.keys()]
    re_keys= [key for key in re._tables.keys()]

    idx = r_tts.index('LLUV')
    r_key = r_keys[idx]
    idx = re_tts.index('LLUV')
    re_key = re_keys[idx]
    
    # copy correctly formatted empty data table from re (Radial object) to r (Radial object)
    r._tables[r_key] = re._tables[re_key]
    r.data = re._tables[re_key]['data']

    return r

def do_merge(datadir, fn, pattern, qccodar_values):
    """ Calls LLUVMerger to create the merged file, checks that output filename includes the expected
        time, (if not the script corrects the file name and ensures that time in filename and header
        match), adds diagnostic tables to the end of the file """
    from qccodar.qcutils import add_short_metadata

    ofn = run_LLUVMerger(datadir, fn, pattern, qccodar_values['merge'])
    # ofn = run_LLUVMerger(datadir, fn, pattern, diag='4', **qccodar_values['merge'])
    if ofn:
        r = read_lluv_file(ofn)
        r = add_short_metadata(r, qccodar_values)

        css_interval_minutes = qccodar_values['merge']['css_interval_minutes']
        number_of_css = qccodar_values['merge']['number_of_css']
        if 'QCD' in r.metadata:
            r.metadata['QCD'].append((
                f'QCDSettings: merge ['
                f'css_interval_minutes = {css_interval_minutes}(minutes), '
                f'number_of_css = {number_of_css}(files)] '
            ))

        for key in r._tables.keys():
            table = r._tables[key]
            if 'MRGS' in table['TableType']:
                filelist = table['data']['PATH']

        if filelist.shape[0] < 2:
            # delete the file from the computer
            print('Not enough short files available for the merge.')
            os.remove(ofn)
            ofn=''
            return ofn
        else:
            shortpath = os.path.join(datadir, 'RadialShorts_qcd',pattern)
            r = check_headertime(r,ofn)
            r = add_diagnostic_tables(r, shortpath)
            # if there is offending data table header, fix the empty table
            if r._tables[1]['_TableHeader'] == [[''], ['']]:
                r = fix_empty_radial(r)
            write_output(r, ofn)
            return ofn
