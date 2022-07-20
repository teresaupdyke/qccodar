## Last modified: Time-stamp: <2022-06-27 10:53:03 codar>
""" Vizualization tool for QC process and settings -- NEW VERSION --

For a given filename from list:
  First plot -- map of qc'd radialmetric radial velocities based qccodar_values OR which single and dual solutions
  Second plot -- map of resulting radialshort 
  Third plot -- histograms of individual tests and the single and dual solutions based on qccodar_values
Animate with each new filename or qccodar_values change.

Other GUI:
  Move forward and back with each file (based on datadir filenames)
  Change threshold values for each threhold test [5.0 50.0 5.0 5.0]
  Change numfiles 3 (odd int) [1-7]
  Change numdegrees 3 (odd int) [1-7]
  Select different weighting factor MP (MP, SNR, NONE)

Using IPython console, use magic to run code as if at unix prompt and provide datadir, patterntype, fn
e.g. %run qcviz.py [datadir] [patterntype] [configfile]

In[]: cd ~/qccodar_dev/tools
In[]: %run qcviz.py /Volumes/TRANSFER/reprocess_HATY/2017_01 IdealPattern 
In[]: plt.show()

Using test dataset defaults
uses input file ./2017_01/Radialmetric/IdealPattern/
In[]: cd /Volumes/TRANSFER/reprocess_HATY
In[]: %run qcviz.py 
In[]: plt.show()

"""

from qccodar.qcutils import *
from qccodar.codarutils import *
from qccodar.app import load_configs

import sys

import matplotlib 
print('matplotlib backend: %s' % (matplotlib.get_backend(),))
import matplotlib.pyplot as plt

# TODO eventaully use docopt for the commandline inputs argv's
if len(sys.argv)==4:
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    configfile = sys.argv[3]
elif len(sys.argv)==3:
    datadir = sys.argv[1]
    patterntype = sys.argv[2]
    configfile = ''
elif len(sys.argv) < 2:
    datadir = os.path.join('/Volumes/TRANSFER/reprocess_HATY', '2017_01')
    # datadir = '/Volumes/TRANSFER/reprocess_HATY/2017_01/'
    # patterntype = 'MeasPattern' 
    patterntype = 'IdealPattern' 
    configfile=''

# load the config file
qccodar_values = load_configs(configfile)
# print(qccodar_values)

params = {'thresholds' : [qccodar_values['qc_doa_peak_power']['doa_peak_power_min'],
                          qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max'],
                          qccodar_values['qc_monopole_snr']['monopole_snr_min'],
                          qccodar_values['qc_loop_snr']['loop_snr_min']],
          'numfiles' : qccodar_values['metric_concatenation']['numfiles'],
          'numdegrees' :  qccodar_values['weighted_shorts']['numdegrees'],
          'numpoints' :  qccodar_values['qc_radialshort_velocity_count']['radialshort_velocity_count_min'],
          'weight_parameter' : qccodar_values['weighted_shorts']['weight_parameter'],
          'fnidx' : 0,
}

# on this GUI we are keeping SNR the same for monopole and the loops for simplicity
params['thresholds'][3] = params['thresholds'][2]

def update_configs():
    """Update config settings from GUI params.  Need this can run qc tests """
    global params
    global qccodar_values
    qccodar_values['qc_doa_peak_power']['doa_peak_power_min']=params['thresholds'][0]
    qccodar_values['qc_doa_half_power_width']['doa_half_power_width_max']=params['thresholds'][1]
    qccodar_values['qc_monopole_snr']['monopole_snr_min']=params['thresholds'][2]
    qccodar_values['qc_loop_snr']['loop_snr_min']=params['thresholds'][3]
    qccodar_values['metric_concatenation']['numfiles']=params['numfiles']
    qccodar_values['weighted_shorts']['numdegrees']=params['numdegrees']
    qccodar_values['qc_radialshort_velocity_count']['radialshort_velocity_count_min']=params['numpoints']
    qccodar_values['weighted_shorts']['weight_parameter']=params['weight_parameter']
    return

def get_data(datadir, fn, patterntype):
    """
    Get data -- modified from qcutils.do_qc()
    """
    global params, qccodar_values

    # read in the data
    ifn = os.path.join(datadir, 'RadialMetric', patterntype, fn)
    rm = read_lluv_file(ifn)

    # read in other radial metric data to use in averaging over time
    ixfns = find_files_to_concatenate(ifn, **qccodar_values['metric_concatenation'])
    
    for xfn in ixfns:
        if xfn == ifn:
            continue
        
        rs1 = read_lluv_file(xfn)
        if len(rm.data.shape) == len(rs1.data.shape) == 2:
             if (rm.data.shape[1] == rs1.data.shape[1]):
                 if all(rm.data.columns == rs1.data.columns):
                     # if same number and order of columns as rs, then append the data rs
                     if debug:
                         print('... ... include: %s' % xfn)
                     rm.data = pd.concat([rm.data,rs1.data],ignore_index=True)
    
    # (0) rm is all the radialmetric data
    # remove any data from original radialmetric i.e. flagged by CODAR software (might speed things up too?)
    rm = remove_flagged_data(rm)

    # (1) do threshold qc on radialmetric
    rmqc = threshold_qc_all(rm, qccodar_values)

    # (2) do weighted averaging of good, this is done when the radial short is created
    rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])

    # (3) require a minimum numpoints used in to form cell average
    rsx = threshold_rsd_numpoints(rsx, **qccodar_values['qc_radialshort_velocity_count'])

    # for plotting we want to see which are good vectors, so remove badflagged based on qc tests
    rmqc = remove_flagged_data(rmqc)
    rsx = remove_flagged_data(rsx)

    return rm, rmqc, rsx

def remove_flagged_data(r):
    """ Remove data rows from Radial.data that has (bad) non-zero VFLG """
    r1 = copy.deepcopy(r)
    bad = numpy.where(r1.data['VFLG'] > 0)[0]
    bad = bad.tolist()
    r1.data.drop(r1.data.index[bad], axis=0, inplace=True)
    return r1

####################################
# get filenames and get data for first file
####################################
rmfoldername = get_radialmetric_foldername(datadir)
# get file listing of datadir
fns = recursive_glob(os.path.join(datadir, rmfoldername, patterntype), 'RDL*.ruv')
fn = os.path.basename(fns[0])
rm, rmqc, rsx = get_data(datadir, fn, patterntype)

####################################
# Setup figure layout for axes, and all the GUI
####################################
fig = plt.figure(figsize=(11.0, 8.0))

# outer grid to frame maps, inner grid of histogrames and inner grid of gui, 
# use the object handle of figure (fig) and method add_gridspec
ogs = fig.add_gridspec(4,3, left=0.05, right=0.95,  top=0.95, bottom=0.05, hspace=0.05, wspace=0.05)

# ogs[0,0] for rmqc map
# ogs[0,1] for rsx map
# ogs[0:,-1] far-right column for inner grid of GUI (dates/time, qccodar_values)
# ogs[1:,:-1] for inner grid of 3x3 of rmqc histograms

map_axs = [fig.add_subplot(ogs[0:2, 0]),
           fig.add_subplot(ogs[0:2, 1]),
           ]

# inner grid for 3x3 rmqc histograme
igs = matplotlib.gridspec.GridSpecFromSubplotSpec(3,3,subplot_spec=ogs[2:,:-1], hspace=0.15, wspace=0.15)
hist_axs = []
for i in range(9):
    hist_axs.append(fig.add_subplot(igs[i]))

# grids for GUI sliders and buttons on the right-most column
igs = matplotlib.gridspec.GridSpecFromSubplotSpec(25,1,subplot_spec=ogs[:,-1], hspace=0.1)

# igs[0] prev_fn, next_fn
# igs[1] sfn
# igs[2] ---- blank
# igs[3] stest0
# igs[4] stest1
# igs[5] stest2
# igs[6] ---- blank
# igs[7] snumfiles
# igs[8] snumdegrees
# igs[9] snumpoints
# igs[10] --- blank
# igs[11:] rwtavg radio buttons 

# another inner grid 2x1 for the prev and next buttons
igs2 = matplotlib.gridspec.GridSpecFromSubplotSpec(1,4,subplot_spec=igs[0,0], wspace=0.1)

# fn prev button
ax = fig.add_subplot(igs2[0])
bfnprev = matplotlib.widgets.Button(ax, '< prev')

# fn next button
ax = fig.add_subplot(igs2[1], title='xxx')
# and set it to current filename (fn) 
ax.title.set_text(fn)
bfnnext = matplotlib.widgets.Button(ax, 'next >')

ax = plt.subplot(igs[1])
sfn = matplotlib.widgets.Slider(ax, 'Index', 0, len(fns),  valstep=1, valinit=params['fnidx'], valfmt='%d')

# --------
ax = plt.subplot(igs[3], title='Thresholds')
stest0 = matplotlib.widgets.Slider(ax, 'DOA Power', 0, 25,  valstep=0.1, valinit=params['thresholds'][0], valfmt='%3.1f (dB)')

ax = plt.subplot(igs[4])
stest1 = matplotlib.widgets.Slider(ax, 'DOA Width', 0, 180, valstep=10, valinit=params['thresholds'][1], valfmt='%3.1f (deg)')

ax = plt.subplot(igs[5])
stest2 = matplotlib.widgets.Slider(ax, 'SNR Mono', 0, 25, valstep=0.1, valinit=params['thresholds'][2], valfmt='%3.1f (dB)')

# --------
ax = plt.subplot(igs[7], title='Weighting Windows')
snumfiles = matplotlib.widgets.Slider(ax, 'numfiles', 1, 7, valstep=2, valinit=params['numfiles'], valfmt='%d')

ax = plt.subplot(igs[8])
snumdegrees = matplotlib.widgets.Slider(ax, 'numdegrees', 1, 7, valstep=2, valinit=params['numdegrees'], valfmt='%d')

ax = plt.subplot(igs[9])
snumpoints = matplotlib.widgets.Slider(ax, 'numpoints', 1, 11, valinit=params['numpoints'], valfmt='%d')

ax = plt.subplot(igs[11:15], aspect='equal', title='Weighting Param')
rwtdavg = matplotlib.widgets.RadioButtons(ax, ('MP', 'SNR', 'NONE'), active=0)

# change where sliders are labeled
all_sliders = [sfn, stest0, stest1, stest2, snumfiles, snumdegrees, snumpoints]
for s in all_sliders:
    # s.label.set_horizontalalign('left')
    # x,y = s.label.get_position()
    s.label.set_position((1.00, 0.5)) # puts the label inside on the rhs


######################
# Initialize the maps
######################
# maps (vectors will change with GUI)

site_lat, site_lon = [float(x) for x in rmqc.metadata['Origin'].split('  ')]
map_axs[0].plot(site_lon, site_lat, 'o', markersize=10, markeredgecolor='black', color='red')
map_axs[0].set_title('RadialMetric')
rmuv = map_axs[0].quiver(rmqc.data.LOND, rmqc.data.LATD, rmqc.data.VELU, rmqc.data.VELV)
# map_axs[0].quiver(rmqc.data.LOND, rmqc.data.LATD, rmqc.data.VELU, rmqc.data.VELV, transform=ccrs.PlateCarree())
rmuv_key = map_axs[0].quiverkey(rmuv, 0.05, 0.9, 100, r'100 cm/s', labelpos='E', coordinates='axes')

map_axs[1].set_title('RadialShort')
map_axs[1].plot(site_lon, site_lat, 'o', markersize=10, markeredgecolor='black', color='red')
rsuv = map_axs[1].quiver(rsx.data.LOND, rsx.data.LATD, rsx.data.VELU, rsx.data.VELV)
rsuv_key = map_axs[1].quiverkey(rsuv, 0.05, 0.9, 100, r'100 cm/s', labelpos='E', coordinates='axes')

def plot_data(rmqc, rsx):
    """ Update plots """
    global rmuv, rsuv, rmuv_key, rsuv_key
    #
    # delete vectors and redraw new ones, since positions not always the same
    rmuv_key.remove()
    rsuv_key.remove()
    rmuv.remove()
    rsuv.remove()
    plt.draw()
    rmuv = map_axs[0].quiver(rmqc.data.LOND, rmqc.data.LATD, rmqc.data.VELU, rmqc.data.VELV)
    rmuv_key = map_axs[0].quiverkey(rmuv, 0.05, 0.9, 100, r'100 cm/s', labelpos='E', coordinates='axes')
    rsuv = map_axs[1].quiver(rsx.data.LOND, rsx.data.LATD, rsx.data.VELU, rsx.data.VELV)
    rsuv_key = map_axs[1].quiverkey(rsuv, 0.05, 0.9, 100, r'100 cm/s', labelpos='E', coordinates='axes')
    #
    plt.draw()


# histograms (bars and test tolerances that will change with GUI)
# 
# empty lists each 3 elements to hold object handles 
tols = bars = [None] * 9

###################
# Ready 3x3 histogram plots with labels and tolerances
###################

hist_axs[0].set_title('Single Soln' )
hist_axs[1].set_title('Dual Soln 1')
hist_axs[2].set_title('Dual Soln 2')

# Row of MSR1, MDR1, MDR2 histogram plots
hist_axs[0].annotate('DOA Pow\nMSR1', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[1].annotate('DOA Pow\nMDR1', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[2].annotate('DOA Pow\nMDR2', xy=(0.5,0.8), xycoords=u'axes fraction')

# Row of MSW1, MDW1, MDW2 histogram plots
hist_axs[3].annotate('DOA Width\nMSW1', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[4].annotate('DOA Width\nMDW1', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[5].annotate('DOA Width\nMDW2', xy=(0.5,0.8), xycoords=u'axes fraction')

# Row of MSW1, MDW1, MDW2 histogram plots
hist_axs[6].annotate('SNR Mono\nMA3S', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[7].annotate('SNR Mono\nMA3S', xy=(0.5,0.8), xycoords=u'axes fraction')
hist_axs[8].annotate('SNR Mono\nMA3S', xy=(0.5,0.8), xycoords=u'axes fraction')


def plot_hist_tolerances():
    global tols
    for tol in tols:
        if tol is not None:
            # remove tolerance polygons before we can add new updated ones
            tol.remove()
    
    tols[0] = hist_axs[0].axvspan(0, params['thresholds'][0], edgecolor=None, facecolor='0.5', zorder=0)
    tols[1] = hist_axs[1].axvspan(0, params['thresholds'][0], edgecolor=None, facecolor='0.5', zorder=0)
    tols[2] = hist_axs[2].axvspan(0, params['thresholds'][0], edgecolor=None, facecolor='0.5', zorder=0)
    
    tols[3] = hist_axs[3].axvspan(params['thresholds'][1], 180., edgecolor=None, facecolor='0.5', zorder=0)
    tols[4] = hist_axs[4].axvspan(params['thresholds'][1], 180., edgecolor=None, facecolor='0.5', zorder=0)
    tols[5] = hist_axs[5].axvspan(params['thresholds'][1], 180., edgecolor=None, facecolor='0.5', zorder=0)

    tols[6] = hist_axs[6].axvspan(0, params['thresholds'][2], edgecolor=None, facecolor='0.5', zorder=0)
    tols[7] = hist_axs[7].axvspan(0, params['thresholds'][2], edgecolor=None, facecolor='0.5', zorder=0)
    tols[8] = hist_axs[8].axvspan(0, params['thresholds'][2], edgecolor=None, facecolor='0.5', zorder=0)

    plt.draw()


def plot_hist_bars(rm):
    global bars
    for bar in bars:
        if bar is not None:
            # remove bar objects, before we can draw new updated ones
            bar.remove()

    bars = [None] * 9
    want = rm.data['MSEL'] == 1
    _,_,bars[0]= hist_axs[0].hist(rm.data.loc[want, 'MSR1'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")
    _,_,bars[3]= hist_axs[3].hist(rm.data.loc[want, 'MSW1'], bins=50, range=(0, 180), lw=1, ec="black", fc="blue")
    _,_,bars[6]= hist_axs[6].hist(rm.data.loc[want, 'MA3S'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")

    want = rm.data['MSEL'] == 2
    _,_,bars[1]= hist_axs[1].hist(rm.data.loc[want, 'MDR1'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")
    _,_,bars[4]= hist_axs[4].hist(rm.data.loc[want, 'MDW1'], bins=50, range=(0, 180), lw=1, ec="black", fc="blue")
    _,_,bars[7]= hist_axs[7].hist(rm.data.loc[want, 'MA3S'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")

    want = rm.data['MSEL'] == 3
    _,_,bars[2]= hist_axs[2].hist(rm.data.loc[want, 'MDR2'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")
    _,_,bars[5]= hist_axs[5].hist(rm.data.loc[want, 'MDW2'], bins=50, range=(0, 180), lw=1, ec="black", fc="blue")
    _,_,bars[8]= hist_axs[8].hist(rm.data.loc[want, 'MA3S'], bins=50, range=(0, 50), lw=1, ec="black", fc="blue")

    plt.draw()

##############
# initialize the bars and display tolerances
##############
plot_hist_bars(rm)
plot_hist_tolerances()

#####################

# GUI functions
def prev_fn(val):
    fnidx = int(sfn.val)
    sfn.set_val(fnidx-1)

def next_fn(val):
    fnidx = int(sfn.val)
    sfn.set_val(fnidx+1)

def sfn_change(val):
    # when fn slider changes
    global params, fn, rm, rmqc, rsx
    params['fnidx'] = int(sfn.val)
    idx = int(sfn.val)
    fn = os.path.basename(fns[idx])
    # set title to the current fn
    bfnnext.ax.set_title(fn)
    plt.draw()

    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_hist_bars(rm)
    plot_data(rmqc, rsx)


def wtdavg_change(label):
    global rm, rmqc, rsx
    params['weight_parameter'] = label
    update_configs()

    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_hist_bars(rm)
    plot_data(rmqc, rsx)

    # rsx.data = weighted_velocities(rmqc, params['numdegrees'], params['weight_parameter'])
    # rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
    # rsx = threshold_rsd_numpoints(rsx, params['numpoints'])
    # plot_data(rmqc, rsx)

def stest0_change(val):
    global params, rm, rmqc, rsx, tols
    params['thresholds'][0] = stest0.val
    update_configs()

    plot_hist_tolerances()
    
    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)

def stest1_change(val):
    global params, rm, mqc, rsx
    params['thresholds'][1] = stest1.val
    update_configs()

    plot_hist_tolerances()
    
    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)

def stest2_change(val):
    global params, rmqc, rsx
    # for simplicity keep SNR for monopole and loops the same
    params['thresholds'][2] = stest2.val
    params['thresholds'][3] = stest2.val
    update_configs()

    plot_hist_tolerances()

    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_data(rmqc, rsx)

def snumfiles_change(val):
    global params, rm, rmqc, rsx
    numfiles = snumfiles.val
    params['numfiles'] = int(numfiles)
    update_configs()

    plot_hist_tolerances()
    
    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_hist_bars(rm)
    plot_data(rmqc, rsx)

def snumdegrees_change(val):
    global params, rm, rmqc, rsx
    numdegrees = snumdegrees.val
    params['numdegrees'] = int(numdegrees)
    update_configs()

    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_hist_bars(rm)
    plot_data(rmqc, rsx)

def snumpoints_change(val):
    global params, rm, rmqc, rsx
    numpoints = snumpoints.val
    params['numpoints'] = int(numpoints)
    update_configs()
    
    rm, rmqc, rsx = get_data(datadir, fn, patterntype)
    plot_hist_bars(rm)
    plot_data(rmqc, rsx)
    
    # rsx.data = weighted_velocities(rmqc, params['numdegrees'], params['weight_parameter'])
    # rsx = generate_radialshort(rmqc, **qccodar_values['weighted_shorts'])
    # rsx = threshold_rsd_numpoints(rsx, params['numpoints'])
    # plot_data(rmqc, rsx)


# GUI actions
bfnprev.on_clicked(prev_fn)
bfnnext.on_clicked(next_fn)
sfn.on_changed(sfn_change)

stest0.on_changed(stest0_change)
stest1.on_changed(stest1_change)
stest2.on_changed(stest2_change)

snumfiles.on_changed(snumfiles_change)
snumdegrees.on_changed(snumdegrees_change)
snumpoints.on_changed(snumpoints_change)

rwtdavg.on_clicked(wtdavg_change)

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


plt.show()

