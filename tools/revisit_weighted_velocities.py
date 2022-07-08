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

# test timing differences using pandas (dataframe) vs numpy array in weighted_velocities
#
# theory is that indexing dataframe is 100x slower than indexing numpy array
# current version with dataframe (weighted_velocities_df) takes ~40-50 sec for each RadialShort
# this version does weighted avg and indexes into a dataframe
# 
# so will test version using pure numpy array (weighted_velocities_np) to do weighted avg
# then return dataframe created from the numpy array
#
# time both by using 
# from revisit_weighted_velocities import weighted_velocities_df as weighted_velocities
# or
# from revisit_weighted_velocities import weighted_velocities_np as weighted_velocities

def weighted_velocities_df(r, numdegrees=3, weight_parameter='MP'):
    """Calculates weighted average of radial velocities (VELO) at bearing and range.

    The weighted average of velocities found at given range and
    bearing based on weight_parameter.

    Paramters
    ---------
    r : Radial object -- created by hfradarpy
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
    xd : pandas Dataframe with columns labeled
       The averaged values with range and bearing.

    """

    print(" ... DATAFRAME VERSION OF WEIGHTED_VELOCITIES")

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



def weighted_velocities_np(r, numdegrees=3, weight_parameter='MP'):
    """Calculates weighted average of radial velocities (VELO) at bearing and range.

    The weighted average of velocities found at given range and
    bearing based on weight_parameter.
    
    NOTE:  this version, like an earlier version, uses numpy array for building the RadialShort data array
    Indexing into the numpy array to fill data within the loop is much faster

    Paramters
    ---------
    r : Radial object -- created by hfradarpy with the data from LLUV file(s). 
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
    df : pandas Dataframe with columns labeled
       The averaged values with range and bearing.

    """

    print(" ... NUMPY ARRAY VERSION OF WEIGHTED_VELOCITIES")


    # 
    # order of columns and labels for output data
    xcol_labels = ['VFLG', 'SPRC', 'BEAR', 'VELO', 'ESPC', 'MAXV', 'MINV', 'EDVC', 'ERSC']
    xc = get_columns(' '.join(xcol_labels))
    
    # data and columns from input Radial object as numpy array
    d = copy.deepcopy(r.data.to_numpy()) # NOTE using np array, not the dataframe
    c = get_columns( ' '.join(r.data.columns.to_list()) ) # dict of column labels, their index
    offset = ((numdegrees-1)/2)

    ud = unique_rows(d[:,[c['SPRC'],c['BEAR'],c['VFLG']]].copy())
    # return only rows that have VFLG==0 (0 == good, >0 bad) so only get good data
    ud = ud[ud[:,2]==0]
    if ud.size == 0:
        r.data = numpy.array([])
        return r.data

    #
    allbearings = numpy.unique(ud[:,1])
    allranges = numpy.unique(ud[:,0])
    ud = numpy.array([[r, b] for r in allranges for b in allbearings])
    
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

    # create dataframe with the np array
    df = pd.DataFrame(xd)
    df.columns = xcol_labels
    # delete extra lines (nan) not filled above
    df.dropna(axis=0, how='all', inplace=True)
    # ESPC had several nan's, replace these and other nan before returning?
    df.fillna(999.000, inplace=True)

    return df

