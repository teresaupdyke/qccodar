#!/bin/bash

# ============================================================================
# Run qccodar script
# Last updated: 03/06/2025 (Teresa Updyke)
# Features:
# 1) Writes output to a log file in /Users/codar/qccodar_files/logs
# 2) Uses a lock file to prevent multiple instances of the script from running 
#    at the same time, which could cause problems with data processing.
# 3) Handles miniforge3, mambaforge3 or miniconda3 environments
# ============================================================================

# Set up environment
PATH=/bin:/usr/bin;
TMPDIR=/tmp;
LOGDIR='/Users/codar/qccodar_files/logs';


# Script name with path stripped off
script=$(basename $0);

# Look for lock file.  If it exists, this script is running already so we need
# to exit
lockFile="${TMPDIR}/${script}.lock";
staleHours=12; # Number of hours after which a lock file is considered stale

# Delete the lock file if the scripts exits due to kill, etc.
trap "rm $lockFile; exit 255;" SIGINT SIGHUP SIGKILL SIGTERM SIGQUIT;

# Script logging
if [ -f "$lockFile" ]
then
    echo "Script currently locked: $lockFile";
    echo "Another instance of the script is already running. Exiting."    
    exit 1
fi

# Create the lock file to lock the script from executing again
echo -n "Locking script with lock file: $lockFile...";
touch $lockFile;
if [ "$?" -ne 0 ]
then
    echo "$0: Failed to lock script...EXITING" >&2;
    exit 1;
fi
echo "LOCKED.";

# Setup logfile
LOGFILE="${LOGDIR}/${script}_$(date +%Y-%m-%d).log";
echo "LOG FILE: $LOGFILE";

# Start script
date >> $LOGFILE

# Run qccodar
if [ -f "/Users/codar/miniforge3/envs/qccodar/bin/qccodar" ]; then
    /Users/codar/miniforge3/envs/qccodar/bin/qccodar auto -p IdealPattern  >> $LOGFILE
elif [ -f "/Users/codar/mambaforge3/envs/qccodar/bin/qccodar" ]; then
    /Users/codar/mambaforge3/envs/qccodar/bin/qccodar auto -p IdealPattern  >> $LOGFILE
elif [ -f "/Users/codar/miniconda3/envs/qccodar/bin/qccodar" ]; then
    /Users/codar/miniconda3/envs/qccodar/bin/qccodar auto -p IdealPattern  >> $LOGFILE
else
    echo "qccodar not found in miniforge3, mambaforge3 or miniconda3 environment.  Please check your installation." >> $LOGFILE
fi

date >> $LOGFILE

echo -n "Removing lock file: $lockFile...";
rm $lockFile;
if [ "$?" -ne 0 ]
then
    echo "Lock file ERROR: $lockFile" >&2;
fi

