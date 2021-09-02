#!/usr/bin/env pythonw
#
# Last modified: Time-stamp: <2017-10-13 16:29:02 codar>
""" Vizualization tool for QC process and settings

For a given bearing:
  First plot -- velocities with range cell for a given bearing
  Second plot -- SNR with range or other selected radialmetric parameter
  Third plot -- compass plot of current bearing direction 
Animate with bearing.

Other GUI:
  Change threshold values for each threhold test [5.0 50.0 5.0 5.0]
  Change numfiles 3 (odd int)
  Change numdegrees 3 (odd int)
  Select different weighting factor MP (MP, SNR, NONE)

Using IPython console, use magic to run code as if at unix prompt and provide datadir, patterntype, fn
e.g. %run qcviz.py [datadir] [patterntype] [fn]

In[]: cd workspace/qc-codar-radialmetric/
In[]: %run qcviz.py /Users/codar/Documents/reprocessing_2014_11/Reprocess_HATY_70_35/ IdealPattern RDLv_HATY_2013_11_05_0000.ruv
In[]: plt.show()

Using test dataset defaults
uses input file ./test/files/codar_raw/Radialmetric/IdealPattern/RDLv_HATY_2013_11_05_0000.ruv
In[]: cd workspace/qc-codar-radialmetric/
In[]: %run qcutil.py 
In[]: plt.show()

"""
from qccodar.qcutils import *
from qccodar.codarutils import *
from qccodar.qcviz.sliders import IndexedSlider

import sys

import matplotlib 
print('matplotlib backend: %s' % (matplotlib.get_backend(),))
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

# initialize the main parameters that can be changed on the plot
params = {'thresholds' : [5.0, 50.0, 5.0, 5.0],
          'numfiles' : 3,
          'numdegrees' :  3,
          'numpoints' :  2,
          'weight_parameter' : 'MP',
          'bearing' : 0,
}

pylab.rcParams['figure.figsize']= (11.0, 8.0)
fig, axs = plt.subplots(3,1)

axs[0].axhline(y=0, linewidth=1, color='k')
axs[0].set_ylim(-150, 150)
axs[0].set_ylabel('Radial Velocity (cm/s)')
ld_bad, = axs[0].plot([], [], 'ro', mec='red', mfc='None')
ld_good, = axs[0].plot([], [], 'go', mec='g', markersize=8)
lrs, = axs[0].plot([], [], 'bo-', mec='yellow')

axs[1].set_ylim(0, 45)
axs[1].set_xlabel('Range Cell')
axs[1].set_ylabel('MUSIC Power (dB)')
ls_bad, = axs[1].plot([], [], 'ro', mec='r', mfc='None')
ls_good, = axs[1].plot([], [], 'go', mec='g', markersize=8)

# Legend for upper plots
axleg = fig.legend((ld_good,ld_bad, lrs), ('good', 'badflagged', 'wtd averge'), 'upper right')

# change position of last plot to stay on left margin but make it square
bb2 = axs[2].get_position()
bb2.bounds = (bb2.bounds[0], bb2.bounds[1], bb2.height, bb2.height)
axs[2].set_position(bb2)

axs[2].set_ylim(-1,1)
axs[2].set_xlim(-1,1)
axs[2].axhline(y=0, linewidth=1, color='k')
axs[2].axvline(x=0, linewidth=1, color='k')
axs[2].set_aspect('equal')
axs[2].set_xticklabels('')
axs[2].set_yticklabels('')

lbear, = axs[2].plot([0,compass2uv(1,45)[0]], [0,compass2uv(1,45)[1]], 'b-')


# Widget functions
def sbear_change(val):
    global params
    params['bearing'] = sbear.seqvals[val]
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def wtdavg_change(label):
    global params
    global rsd
    params['weight_parameter'] = label

    if params['weight_parameter']=='MP':
        axs[1].set_ylim(-125, -75)
        axs[1].set_ylabel('MUSIC Power (dB)')
    elif params['weight_parameter']=='SNR':
        axs[1].set_ylim(0, 45)
        axs[1].set_ylabel('Monopole (A3) SNR (dB)')
    elif params['weight_parameter']=='NONE':
        axs[1].set_ylim(0, 1)
        axs[1].set_ylabel('No Weighting Param')
 
    rsd, rstypes_str = weighted_velocities(d, types_str, params['numdegrees'], params['weight_parameter'])
    rsd = threshold_rsd_numpoints(rsd, rstypes_str, params['numpoints'])
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def stest1_change(val):
    global params, d, types_str, rsd, rstypes_str
    params['thresholds'][0] = stest1.seqvals[val]
    d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def stest2_change(val):
    global params, d, types_str, rsd, rstypes_str
    params['thresholds'][1] = stest2.seqvals[val]
    d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def stest3_change(val):
    global params, d, types_str, rsd, rstypes_str
    params['thresholds'][2] = stest3.seqvals[val]
    d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def snumfiles_change(val):
    global params, d, types_str, rsd, rstypes_str
    numfiles = snumfiles.seqvals[val]
    params['numfiles'] = int(numfiles)
    d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def snumdegrees_change(val):
    global params, d, types_str, rsd, rstypes_str
    numdegrees = snumdegrees.seqvals[val]
    params['numdegrees'] = int(numdegrees)
    d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

def snumpoints_change(val):
    global params, rsd, rstypes_str
    numpoints = snumpoints.seqvals[val]
    params['numpoints'] = int(numpoints)
    rsd, rstypes_str = weighted_velocities(d, types_str, params['numdegrees'], params['weight_parameter'])
    rsd = threshold_rsd_numpoints(rsd, rstypes_str, params['numpoints'])
    plot_data(d, types_str, rsd, rstypes_str)
    fig.canvas.draw()

# Widgets
axbear = plt.axes([0.1, 0.05, 0.8, 0.03])
sbear = IndexedSlider(axbear, 'Bearing', seqvals=list(range(0,269,1)), valinit=0, valfmt='%03d (deg)')
sbear.on_changed(sbear_change)

axradio = plt.axes([0.4, 0.1, 0.15, 0.15], aspect='equal', title='Weighting Param')
rwtdavg = matplotlib.widgets.RadioButtons(axradio, ('MP', 'SNR', 'NONE'), active=1)
rwtdavg.on_clicked(wtdavg_change)

# outer grid to frame inner grid of sliders
ogs = matplotlib.gridspec.GridSpec(3,3)
# use lower-right ogs for sliders
igs = matplotlib.gridspec.GridSpecFromSubplotSpec(7,1,subplot_spec=ogs[-1,-1], hspace=0.0)

axtest1 = plt.subplot(igs[0], title='Thresholds')
stest1 = IndexedSlider(axtest1, 'DOA Power', seqvals=numpy.arange(0,25,0.1), \
                       valinit=params['thresholds'][0], valfmt='%3.1f (dB)')
stest1.on_changed(stest1_change)
axtest2 = plt.subplot(igs[1])
stest2 = IndexedSlider(axtest2, 'DOA Width', seqvals=list(range(100,0,-1)), \
                       valinit=params['thresholds'][1], valfmt='%3.1f (deg)')
stest2.on_changed(stest2_change)
axtest3 = plt.subplot(igs[2])
stest3 = IndexedSlider(axtest3, 'SNR Mono', seqvals=numpy.arange(0,25,0.1), \
                       valinit=params['thresholds'][2], valfmt='%3.1f (dB)')
stest3.on_changed(stest3_change)

axnf = plt.subplot(igs[4], title='Weighting Windows')
snumfiles = IndexedSlider(axnf, 'numfiles', seqvals=[1,3,5,7], valinit=params['numfiles'], valfmt='%d')
snumfiles.on_changed(snumfiles_change)

axnd = plt.subplot(igs[5])
snumdegrees = IndexedSlider(axnd, 'numdegrees', seqvals=[1,3,5,7], valinit=params['numdegrees'], valfmt='%d')
snumdegrees.on_changed(snumdegrees_change)

axnp = plt.subplot(igs[6])
snumpoints = IndexedSlider(axnp, 'numpoints', seqvals=list(range(1,11)), valinit=params['numpoints'], valfmt='%d')
snumpoints.on_changed(snumpoints_change)

def subset_rsdata(rsd, rsc, bearing):
    # get data from qc'd and averaged (now data for RadialShorts) array
    xrow = numpy.where( (rsd[:,rsc['BEAR']]==bearing) & (rsd[:,rsc['VFLG']]==0))[0]
    xcol = numpy.array([rsc['VELO'], rsc['SPRC'], rsc['BEAR']])
    a = rsd[numpy.ix_(xrow, xcol)]
    return a

def subset_data_good(d, c, bearing, numdegrees):
    # get GOOD data from RadialMetric array that is not badflagged 
    offset = ((numdegrees-1)/2)
    xrow = numpy.where( (d[:,c['BEAR']]>=bearing-offset) & \
                        (d[:,c['BEAR']]<=bearing+offset) & \
                        (d[:,c['VFLG']]==0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    # Create array to hold each Music Power (based on MSEL)
    MP = numpy.array(numpy.ones(a[:,0].shape)*numpy.nan) 
    # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
    for msel in [1, 2, 3]:
        which = a[:,4]==msel
        MP[which,] = a[which, msel+4]
    # append as last column
    a = numpy.hstack((a,MP.reshape(MP.size,1)))
    return a

def subset_data_bad(d, c, bearing, numdegrees):
    # get BAD data from RadialMetric array that is badflagged 
    offset = ((numdegrees-1)/2)
    xrow = numpy.where( (d[:,c['BEAR']]>=bearing-offset) & \
                        (d[:,c['BEAR']]<=bearing+offset) & \
                        (d[:,c['VFLG']]>0) )[0]
    xcol = numpy.array([c['VELO'], c['SPRC'], c['BEAR'], c['MA3S'], c['MSEL'], c['MSP1'], c['MDP1'], c['MDP2']])
    a = d[numpy.ix_(xrow, xcol)]
    # Create array to hold each Music Power (based on MSEL)
    MP = numpy.array(numpy.ones(a[:,0].shape)*numpy.nan) 
    # pluck the msel-based Music Power from MSP1, MDP1 or MPD2 column
    for msel in [1, 2, 3]:
        which = a[:,4]==msel
        MP[which,] = a[which, msel+4]
    # append as last column
    a = numpy.hstack((a,MP.reshape(MP.size,1)))
    return a

def compass2deg(az):
    """ Convert compass azimuth to cartesian angle in degrees

    https://en.wikipedia.org/wiki/Polar_coordinate_system#Converting_between_polar_and_Cartesian_coordinates

    Using arctan2 does this relative to (x0,y0)=(1,0)
    """
    x,y = compass2uv(1,az)
    # arctan2 does this relative to (x0,y0)=(1,0)
    return numpy.arctan2(y, x)*180./numpy.pi

def deg2compass(deg):
    """ Convert cartesian angle of degrees to compass azimuth"""
    compass = 90.-deg
    if compass < 0.:
        compass = compass + 360.
    return compass

def init_plot(d, types_str, rsd, rstypes_str):
    """ Set plot and slider limits
    """
    global params

    c = get_columns(types_str)
    rsc = get_columns(rstypes_str)
    xrow = numpy.where( (d[:,c['VFLG']]==0) )[0]
    allranges = numpy.unique(d[:,c['SPRC']] )
    allbearings = numpy.unique(d[xrow,c['BEAR']])
    thetas = numpy.array([compass2deg(b) for b in allbearings])
    lhs = int(deg2compass(thetas.max()))
    rhs = int(deg2compass(thetas.min()))
    if lhs > rhs:
        allbearings = numpy.arange(lhs, rhs+360.,1.) % 360
    else:
        allbearings = numpy.arange(lhs, rhs, 1.)
    params['bearing'] = allbearings[0]
    
    fig.suptitle(fn)
    axs[0].set_xlim(0, allranges.max()+2)
    axs[1].set_xlim(0, allranges.max()+2)

    # add radar patch limits
    axs[2].add_patch(matplotlib.patches.Wedge((0,0), 1, thetas.min(), thetas.max(), \
                                              zorder=-1, ec='None', fc=(.9,.9,.9)))

    # reset bearing slider since we now know what the bearings are from data
    sbear.seqvals = allbearings.tolist()
    sbear.valinit = sbear.seqvals.index(allbearings[0])
    sbear.minval = 0
    sbear.maxval = len(allbearings)
    sbear.set_val(sbear.valinit)

    # put the first bearing data into the plots
    plot_data(d, types_str, rsd, rstypes_str)


def plot_data(d, types_str, rsd, rstypes_str):
    # print bearing

    # update new bearing line
    lbear.set_xdata([0, compass2uv(1,params['bearing'])[0]])
    lbear.set_ydata([0, compass2uv(1,params['bearing'])[1]])

    c = get_columns(types_str)
    rsc = get_columns(rstypes_str)

    # update plots with new bearing
    rs = subset_rsdata(rsd, rsc, params['bearing'])
    gd = subset_data_good(d, c, params['bearing'], params['numdegrees'])
    bd = subset_data_bad(d, c, params['bearing'], params['numdegrees'])

    velo = 0
    rnge = 1
    mp = 8
    snr = 3
    if gd.size>0:
        ld_good.set_xdata(gd[:,rnge])
        ld_good.set_ydata(gd[:,velo])
        if params['weight_parameter']=='MP':
            ls_good.set_xdata(gd[:,rnge])
            ls_good.set_ydata(gd[:,mp])
        elif params['weight_parameter']=='SNR':
            ls_good.set_xdata(gd[:,rnge])
            ls_good.set_ydata(gd[:,snr])
        elif params['weight_parameter']=='NONE':
            ls_good.set_xdata([])
            ls_good.set_ydata([])

    if bd.size>0:
        ld_bad.set_xdata(bd[:,rnge])
        ld_bad.set_ydata(bd[:,velo])
        if params['weight_parameter']=='MP':
            ls_bad.set_xdata(bd[:,rnge])
            ls_bad.set_ydata(bd[:,mp])
        elif params['weight_parameter']=='SNR':
            ls_bad.set_xdata(bd[:,rnge])
            ls_bad.set_ydata(bd[:,snr])
        elif params['weight_parameter']=='NONE':
            ls_bad.set_xdata([])
            ls_bad.set_ydata([])

    if rs.size>0:
        lrs.set_xdata(rs[:,1])
        lrs.set_ydata(rs[:,0])

def debug_data(d, types_str, rsd, rstypes_str):
    # print bearing
    # update plots with new bearing
    rs = subset_rsdata(rsd, rsc, params['bearing'])
    gd = subset_data_good(d, c, params['bearing'], params['numdegrees'])
    bd = subset_data_bad(d, c, params['bearing'], params['numdegrees'])

def get_data(datadir, fn, patterntype):
    """
    """
    global params, d, types_str, rsd, rstypes_str

    # read in the data
    ifn = os.path.join(datadir, 'RadialMetric', patterntype, fn)
    d, types_str, header, footer = read_lluv_file(ifn)

    ixfns = find_files_to_merge(ifn, params['numfiles'], sample_interval=30)
    for xfn in ixfns:
        if xfn == ifn:
            continue
        d1, types_str1, _, _ = read_lluv_file(xfn)
        if len(d.shape) == len(d1.shape) == 2:
            if (d.shape[1] == d1.shape[1]) & (types_str == types_str1):
                # if same number and order of columns as d, then append the data d
                print('... ... merging: %s' % xfn)
                d = numpy.vstack((d,d1))

    # do any threshold qc
    d = threshold_qc_all(d, types_str, params['thresholds'])
   
    # do weighted averaging
    rsd, rstypes_str = weighted_velocities(d, types_str, params['numdegrees'], params['weight_parameter'])
        
    # check minimum numpoints 
    rsd = threshold_rsd_numpoints(rsd, rstypes_str, params['numpoints'])

    return d, types_str, rsd, rstypes_str

if len(sys.argv)==4:
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    fn = sys.argv[3]

elif len(sys.argv) < 2:
    datadir = os.path.join('.', 'test', 'files', 'codar_raw')
    # datadir = '/Users/codar/Documents/reprocessing_2015/Reprocess_HATY_70_35/'
    # patterntype = 'MeasPattern' 
    patterntype = 'IdealPattern' 
    fn = 'RDLv_HATY_2013_11_05_0000.ruv'

d, types_str, rsd, rstypes_str = get_data(datadir, fn, patterntype)
init_plot(d, types_str, rsd, rstypes_str)
plt.show()
