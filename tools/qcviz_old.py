## Last modified: Time-stamp: <2022-06-27 10:53:03 codar>
""" Vizualization tool for QC process and settings -- OLD VERSION --

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
e.g. %run qcviz_old.py [datadir] [patterntype] [fn]

In[]: cd ~/qccodar_dev/tools
In[]: %run qcviz_old.py /Users/codar/Documents/qccodar_dev/tools/2017_01 IdealPattern RDLv_HATY_2017_01_01_0400.ruv
In[]: plt.show()

Using test dataset defaults
uses input file ./2017_01/Radialmetric/IdealPattern/RDLv_HATY_2017_01_01_0400.ruv
In[]: cd ~/Documents/qccodar_dev/tools/
In[]: %run qcviz_old.py 
In[]: plt.show()

"""

from qccodar.qcutils import *
from qccodar.codarutils import *
from qccodar.app import load_configs

import sys

import matplotlib 
print('matplotlib backend: %s' % (matplotlib.get_backend(),))
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from matplotlib.widgets import Slider

# TODO eventaully use docopt for the commandline inputs argv's
if len(sys.argv)==5:
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    fn = sys.argv[3]
    configfile = sys.argv[4]
elif len(sys.argv)==4:
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    fn = sys.argv[3]
    configfile = ''
elif len(sys.argv) < 2:
    datadir = os.path.join('.', '2017_01')
    # datadir = '/Users/codar/Documents/qccodar_dev/tools/2017_01/'
    # patterntype = 'MeasPattern' 
    patterntype = 'IdealPattern' 
    fn = 'RDLv_HATY_2017_01_01_0400.ruv'
    configfile=''

# load the config file
qccodar_values = load_configs(configfile)
# print(qccodar_values)

# initialize the main parameters that can be changed on the plot
# params = {'thresholds' : [5.0, 50.0, 5.0, 5.0],
#           'numfiles' : 3,
#           'numdegrees' :  3,
#           'numpoints' :  2,
#           'weight_parameter' : 'MP',
#           'bearing' : 0,
# }

params = {'thresholds' : [qccodar_values['qc_doa_peak_power']['doa_peak_power_min'],
                          qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max'],
                          qccodar_values['qc_monopole_snr']['monopole_snr_min'],
                          qccodar_values['qc_loop_snr']['loop_snr_min']],
          'numfiles' : qccodar_values['metric_concatenation']['numfiles'],
          'numdegrees' :  qccodar_values['weighted_shorts']['numdegrees'],
          'numpoints' :  qccodar_values['qc_radialshort_velocity_count']['radialshort_velocity_count_min'],
          'weight_parameter' : qccodar_values['weighted_shorts']['weight_parameter'],
          'bearing' : 0,
}

def update_configs():
    """Update config settings from GUI params.  Need this can run qc tests """
    global params
    global qccodar_values
    qccodar_values['qc_doa_peak_power']['doa_peak_power_min']=params['thresholds'][0]
    qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max']=params['thresholds'][1]
    qccodar_values['qc_monopole_snr']['monopole_snr_min']=params['thresholds'][2]
    qccodar_values['qc_loop_snr']['loop_snr_min']=params['thresholds'][3]
    qccodar_values['metric_concatenation']['numfiles']=params['numfiles']
    qccodar_values['qc_radialshort_velocity_count']['radialshort_velocity_count_min']=params['numpoints']
    qccodar_values['weighted_shorts']['weight_parameter']=params['weight_parameter']
    return


pylab.rcParams['figure.figsize']= (11.0, 8.0)
fig, axs = plt.subplots(3,1)

axs[0].axhline(y=0, linewidth=1, color='k')
axs[0].set_ylim(-150, 150)
axs[0].set_ylabel('Radial Velocity (cm/s)')
ld_bad, = axs[0].plot([], [], 'ro', mec='red', mfc='None')
ld_good, = axs[0].plot([], [], 'go', mec='g', markersize=8)
lrs, = axs[0].plot([], [], 'bo-', mec='yellow')

axs[1].set_xlabel('Range Cell')
if params['weight_parameter']=='MP':
    axs[1].set_ylim(-125, -75)
    axs[1].set_ylabel('MUSIC Power (dB)')
elif params['weight_parameter']=='SNR':
    axs[1].set_ylim(0, 45)
    axs[1].set_ylabel('Monopole (A3) SNR (dB)')
elif params['weight_parameter']=='NONE':
    axs[1].set_ylim(0, 1)
    axs[1].set_ylabel('No Weighting Param')

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
    params['bearing'] = int(sbear.val)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def wtdavg_change(label):
    global params
    global rmqc, rsx
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

    update_configs()
    rsx.data = weighted_velocities(rmqc, params['numdegrees'], params['weight_parameter'])
    # rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
    rsx = threshold_rsd_numpoints(rsx, params['numpoints'])
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def stest0_change(val):
    global params, rmqc, rsx
    params['thresholds'][0] = stest0.val
    update_configs()
    
    rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def stest1_change(val):
    global params, rmqc, rsx
    params['thresholds'][1] = stest1.val
    update_configs()
    
    rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def stest2_change(val):
    global params, rmqc, rsx
    params['thresholds'][2] = stest2.val
    update_configs()
    
    rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def snumfiles_change(val):
    global params, rmqc, rsx
    numfiles = snumfiles.val
    params['numfiles'] = int(numfiles)
    update_configs()
    
    rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def snumdegrees_change(val):
    global params, rmqc, rsx
    numdegrees = snumdegrees.val
    params['numdegrees'] = int(numdegrees)
    update_configs()

    rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)
    fig.canvas.draw()

def snumpoints_change(val):
    global params, rmqc, rsx
    numpoints = snumpoints.val
    params['numpoints'] = int(numpoints)
    update_configs()
    
    rsx.data = weighted_velocities(rmqc, params['numdegrees'], params['weight_parameter'])
    # rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
    rsx = threshold_rsd_numpoints(rsx, params['numpoints'])
    plot_data(rmqc, rsx)
    fig.canvas.draw()

# Widgets
axbear = plt.axes([0.1, 0.05, 0.8, 0.03])
sbear = Slider(axbear, 'Bearing', 0, 270, valinit=0, valfmt='%03d (deg)')
sbear.on_changed(sbear_change)

axradio = plt.axes([0.4, 0.1, 0.15, 0.15], aspect='equal', title='Weighting Param')
rwtdavg = matplotlib.widgets.RadioButtons(axradio, ('MP', 'SNR', 'NONE'), active=0)
rwtdavg.on_clicked(wtdavg_change)

# outer grid to frame inner grid of sliders
ogs = matplotlib.gridspec.GridSpec(3,3)
# use lower-right ogs for sliders
igs = matplotlib.gridspec.GridSpecFromSubplotSpec(7,1,subplot_spec=ogs[-1,-1], hspace=0.0)

axtest1 = plt.subplot(igs[0], title='Thresholds')
stest0 = Slider(axtest1, 'DOA Power', 0, 25,  valstep=0.1, valinit=params['thresholds'][0], valfmt='%3.1f (dB)')
stest0.on_changed(stest0_change)
axtest2 = plt.subplot(igs[1])
stest1 = Slider(axtest2, 'DOA Width', 0, 180, valstep=10, valinit=params['thresholds'][1], valfmt='%3.1f (deg)')
stest1.on_changed(stest1_change)
axtest3 = plt.subplot(igs[2])
stest2 = Slider(axtest3, 'SNR Mono', 0, 25, valstep=0.1, valinit=params['thresholds'][2], valfmt='%3.1f (dB)')
stest2.on_changed(stest2_change)

axnf = plt.subplot(igs[4], title='Weighting Windows')
snumfiles = Slider(axnf, 'numfiles', 1, 7, valstep=2, valinit=params['numfiles'], valfmt='%d')
snumfiles.on_changed(snumfiles_change)

axnd = plt.subplot(igs[5])
snumdegrees = Slider(axnd, 'numdegrees', 1, 7, valstep=2, valinit=params['numdegrees'], valfmt='%d')
snumdegrees.on_changed(snumdegrees_change)

axnp = plt.subplot(igs[6])
snumpoints = Slider(axnp, 'numpoints', 1, 11, valinit=params['numpoints'], valfmt='%d')
snumpoints.on_changed(snumpoints_change)

def subset_rsdata(rsd, rsc, bearing):
    # get data from qc'd and averaged (now data for RadialShorts) array
    xrow = numpy.where( (rsd[:,rsc['BEAR']]==bearing) & (rsd[:,rsc['VFLG']]==0))[0]
    xcol = numpy.array([rsc['VELO'], rsc['SPRC'], rsc['BEAR']])
    a = rsd[numpy.ix_(xrow, xcol)]
    # print("rsd data at bearing %03d" % (bearing))
    # print(a)
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

def init_plot(rmqc, rsx):
    """ Set plot and slider limits
    """
    global params

    # data as numpy array
    d = rmqc.data.to_numpy()
    rsd = rsx.data.to_numpy()

    # get column labels
    c = get_columns( ' '.join(rmqc.data.columns.to_list()) )
    rsc = get_columns( ' '.join(rsx.data.columns.to_list()) )

    # 
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
    sbear.valinit = allbearings[0]
    sbear.minval = allbearings[0]
    sbear.maxval = allbearings[-1]
    sbear.set_val(sbear.valinit)

    # put the first bearing data into the plots
    plot_data(rmqc, rsx)


def plot_data(rmqc, rsx):
    # update new bearing line
    lbear.set_xdata([0, compass2uv(1,params['bearing'])[0]])
    lbear.set_ydata([0, compass2uv(1,params['bearing'])[1]])

    # get column labels
    c = get_columns( ' '.join(rmqc.data.columns.to_list()) )
    rsc = get_columns( ' '.join(rsx.data.columns.to_list()) )

    # get data for lines plotted in axs[0] and axs[1] along new bearing
    rsd_sub = subset_rsdata(rsx.data.to_numpy(), rsc, params['bearing'])
    gd_sub = subset_data_good(rmqc.data.to_numpy(), c, params['bearing'], params['numdegrees'])
    bd_sub = subset_data_bad(rmqc.data.to_numpy(), c, params['bearing'], params['numdegrees'])

    velo = 0
    rnge = 1
    mp = 8
    snr = 3

    ld_good.set_xdata(gd_sub[:,rnge])
    ld_good.set_ydata(gd_sub[:,velo])
    if params['weight_parameter']=='MP':
        ls_good.set_xdata(gd_sub[:,rnge])
        ls_good.set_ydata(gd_sub[:,mp])
    elif params['weight_parameter']=='SNR':
        ls_good.set_xdata(gd_sub[:,rnge])
        ls_good.set_ydata(gd_sub[:,snr])
    elif params['weight_parameter']=='NONE':
        ls_good.set_xdata([])
        ls_good.set_ydata([])

    ld_bad.set_xdata(bd_sub[:,rnge])
    ld_bad.set_ydata(bd_sub[:,velo])
    if params['weight_parameter']=='MP':
        ls_bad.set_xdata(bd_sub[:,rnge])
        ls_bad.set_ydata(bd_sub[:,mp])
    elif params['weight_parameter']=='SNR':
        ls_bad.set_xdata(bd_sub[:,rnge])
        ls_bad.set_ydata(bd_sub[:,snr])
    elif params['weight_parameter']=='NONE':
        ls_bad.set_xdata([])
        ls_bad.set_ydata([])

    lrs.set_xdata(rsd_sub[:,1])
    lrs.set_ydata(rsd_sub[:,0])


def get_data(datadir, fn, patterntype):
    """
    """
    global params, qccodar_values

    # read in the data
    ifn = os.path.join(datadir, 'RadialMetric', patterntype, fn)
    rs = read_lluv_file(ifn)

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

    return rmqc, rsx


rmqc, rsx = get_data(datadir, fn, patterntype)
init_plot(rmqc, rsx)
plt.draw()

