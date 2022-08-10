# Default configuration file for qccodar
# 
# 
qccodar_values = dict(
    qc_doa_half_power_width=dict(doa_half_power_width_max=50.0),
    qc_doa_peak_power=dict(doa_peak_power_min=5.0),
    qc_monopole_snr=dict(monopole_snr_min=5.0),
    qc_loop_snr=dict(loop_snr_min=5.0),
    qc_radialshort_velocity_count=dict(radialshort_velocity_count_min=2.0),
    metric_concatenation = dict(numfiles=3, sample_interval=30.0),
    weighted_shorts=dict(numdegrees=3,
                         weight_parameter='MP',
                         table_type='LLUV RDL7'),
    merge=dict(css_interval_minutes=30.0,
               number_of_css=5.0,
               shorts_minute_filter = '*00'
               )
 )


# NOTES on configurable items for merging RadialShorts_qcd to Radials_qcd
# See CODAR documentation /Codar/SeaSonde/Apps/Bin/BinDocs/Guide_LLUVMerger.pdf
#
# keys that can be added to merge step and used by qccodar to do the merge
# Not passing all possible parameters for LLUVMerger but enough to allow for differences in
# short- mid-, and long-range systems and some flexibility. 
# -------------------------
# method='average' # average/median/minimum/maximum/overlay/smallest/largest (Merge method)
#         Default is average.
# minvect='2'  # integer number as string (Minimum vectors for each source point) Default is 2.
# diag='4',    # integer number as string (output diagnostics 0=quiet, 1000=all) Default is 4.
# angalign='2' # Optional integer number as string (angular alignment, modulo of angres) 
# reference='/Users/codar/qccodar_files/HATY_lluv_grid_reference_for_merge.ruv'
#              # Optional lluv or grid file to use for output origin and reference location
#                vectors to merge. If not specified source file is origin.

# 
# -------------------------
# Settings for long-range system:
# -------------------------
#
# css_interval_minutes=30.0
# number_of_css=5.0
# shorts_minute_filter = '*00'
#
#       22:30
# 5\    23:00
# 4 |   23:30
# 3 |-- 00:00 <--expected time for merge of 5 files with 30 minute intervals is 60 min behind source
# 2 |   00:30
# 1 /   01:00 <-- source file time
#       01:30


# -------------------------
# Settings for mid-range system:
# -------------------------
#
# css_interval_minutes=10.0
# number_of_css=7.0
# shorts_minute_filter = '*30'
#
#       23:20
# 7 \   23:30
# 6 |   23:40
# 5 |   23:50
# 4 |-- 00:00 <--expected time for merge of 7 files with 10 minute intervals is 30 min behind source
# 3 |   00:10
# 2 |   00:20
# 1 /   00:30 <-- source file time
#       00:40


