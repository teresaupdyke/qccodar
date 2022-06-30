#
import sys
import copy

# # import os
# import re
# import fnmatch
# import datetime
# import geopy
# import geopy.distance
# import numpy
# from pathlib import Path
# numpy.set_printoptions(suppress=True)
# import pandas as pd
# from hfradarpy.radials import Radial

from qccodar.codarutils import *
from qccodar.qcutils import *


debug = 1

# test different deepcopy in weighted_velocities

def weighted_velocities_SH(r, numdegrees=3, weight_parameter='MP'):
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

        # just copy what is found in d this loop
        a = d.loc[xrow, ['VELO', 'MSEL', 'MSP1', 'MDP1', 'MDP2', 'MA3S']].copy(deep=True) # default is deep but say so for clarity

        # This copies all of d again and again. # I think we don't need this for every irow 
        # a0 = copy.deepcopy(d) 
        # a = a0.loc[xrow, ['VELO', 'MSEL', 'MSP1', 'MDP1', 'MDP2', 'MA3S']].copy(deep=True) # default is deep

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

def generate_radialshort_SH(r, table_type='LLUV RDL7', numdegrees=3, weight_parameter='MP'):
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
        formatfile = Path('..') / 'src' / 'qccodar' / 'file_formats' / 'radialshort_LLUV_RDL7.ruv'
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

    xd = weighted_velocities_SH(r, numdegrees, weight_parameter)

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

