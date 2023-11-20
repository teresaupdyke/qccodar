# qccodar

This python code applies several quality (QC) functions based
on CODAR SeaSonde (COS) Radial Metric data currently output in COS
RadialSuite version 7.x (and versions 8.x, 21.x and 22.x upon request). There are two 
modes to this code: an auto-and manual-mode.  Auto-mode is for realtime
processing. Manual-mode is for processing all RadialMetric files in a
specified folder in post-processing. qccodar is intended to run beside 
(in parallel) to SeaSonde Analysis and not supplant any processing.  
In fact, qccodar uses the LLUVMerger.app provided by SeaSonde to merge 
the data back into standard SeaSonde processing methodology.

## Installation

qccodar is a python package that runs under Python 3. 
The following instructions describe how to install qccodar using [Miniconda3](https://conda.io/miniconda.html) and ensure that system requirements are met.
At this point, please check the CODAR software requirements at the end of this file and request any key files or add-on software if necessary.

### 1. Download
Use one of the two instructions below depending the type of Mac computer chip.  If unsure about which to use, select "About this Mac" from the Apple icon menu at top left of the screen, and a display will appear that identifies the type of computer chip. 
Open a Terminal window to execute the following commands.

#### For M1 chip:
```zsh
   cd ~/Downloads
   curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o "miniconda3.sh"
 ```
#### For Intel chip:
```bash
   cd ~/Downloads
   curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -o "miniconda3.sh"
```



Now to install this code, update the permissions on the file and run it with the following commands:

```bash
   chmod 755 miniconda3.sh
   ./miniconda3.sh
```

You must press Enter several times to read through the license agreement.  When prompted, type 'yes' to accept the license terms. 
Press ENTER to confirm the install location /Users/codar/miniconda3. 
You will be asked if you wish to update your shell profile to automatically initialize conda. 
Type 'yes'. After installation type 'exit' to close the shell and then open another Terminal window. 

### 2. Shell Setup Instructions

Type the following:
```bash
   conda -V
```
If you recieve a response that includes a conda version number (for example, conda 23.10.0), 
skip the rest of this section and proceed to Step 3 Environment Setup. 
If you receive a response that includes "command not found", then continue with the next instruction.

Use one of the two sets of instructions below depending on the default Terminal shell (zsh or bash).
Open a Terminal window and the name of the shell will be shown on the top bar of the window.

#### For zsh:
As user `codar`, edit or create the file /Users/codar/.zshrc

Note: If you navigate to /Users/codar in Finder, you can show this hidden file 
by pressing Command + Shift + . (the period key).
When finished with this step, you can repeat the keystroke to restore normal file viewing.

Make sure the following lines are included in the file:

```zsh
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/codar/miniconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/codar/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/codar/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/codar/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
# <<< conda initialize <<< 
 ```
Save the file. Type 'exit' to close the shell and then open another Terminal window. 
Try the conda -V command again.
If it doesn't work, try sourcing the profile from the command prompt with this command.
```zsh
   source ~/.zshrc
```

#### For bash:

As user `codar`, edit or create the file `/Users/codar/.bash_profile`

Note: If you navigate to /Users/codar in Finder, you can show this hidden file 
by pressing Command + Shift + . (the period key).  
When finished with this step, you can repeat the keystroke to restore normal file viewing.

Make sure the following lines are included in the file:

```bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/codar/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/codar/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/codar/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/codar/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<< 
 ```

Save the file. Type 'exit' to close the shell and then open another Terminal window. 
Try the conda -V command again.
If it doesn't work, try sourcing the profile from the command prompt with this command.
```bash
   source ~/.bash_profile
```


### 3. Environment Setup Instructions
 
Create a directory to store the qccodar files and a subdirectory to store the log files.
```bash
   mkdir /Users/codar/qccodar_files
   mkdir /Users/codar/qccodar_files/logs
```

Download the qccodar code by going to https://github.com/teresaupdyke/qccodar
On this website, click the "Code button" and choose "Download ZIP" from the menu.
In Finder, open the Downloads folder (/Users/codar/Downloads).
Look for the downloaded code as a zip file. Doubleclick on qccodar-main.zip to unzip the file.  
(The computer may automatically unzip the file for you, in which case you'll simply see the qccodar-main folder.) 
Move the entire qccodar-main folder into /Users/codar/qccodar_files/

Back in the Terminal window, use the following instructions to create a conda environment that allows the qccodar module and its
dependencies to run isolated from the installed Python avoiding
potential module version conflicts. 

```bash
   cd /Users/codar/qccodar_files/qccodar-main
   conda env create -f environment.yml
```
```bash
   conda activate qccodar
```
Now the command prompt should be preceded by the label (qccodar).
Run this command:
```bash
   which python
```
and the result should be  
```Users/codar/miniconda3/envs/qccodar/bin/python```

While still in the same directory, install qccodar to the environment: 
```bash
   pip install .
```
After the installation is complete, you can deactivate the environment:
```bash
   conda deactivate
```

## 4. Configuration for Realtime QC

First, enable RadialMetric output:

1. Edit line 21 of `AnalysisOptions.txt` in /Codar/SeaSonde/Configs/RadialConfigs.
1. Restart AnalyzeSpectra

```
   1           !21 Enable Radial Metric Output: 0(Off), 1(Enable), 2(Enable MaxVel.)
```


Copy the default qccodar configuration file into a standardized location with standardized name.  
For a 5 MHz system:
```bash
   cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_5MHz.plist /Users/codar/qccodar_files/qccodar.plist
```
For a 13 MHz system:
```bash
   cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_13MHz.plist /Users/codar/qccodar_files/qccodar.plist
```
For a 25 MHz system:
```bash
   cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_25MHz.plist /Users/codar/qccodar_files/qccodar.plist
```
If a customized configuration is needed, edit /Users/codar/qccodar_files/qccodar.plist

For example if you have 30 minute output, use a different shorts_minute_filter:
		
                <key>shorts_minute_filter</key>
                <string>*[0,3]0</string>

An analysis with the qcviz.py program in this software package can inform decisions on which thresholds to use for a particular station.  
Instructions for qcviz will be included in a future release.


## 5. Crontab Entry for Realtime QC

Finally, set crontab entry to run qccodar:

2. Place entry in crontab to run every 15 minutes and log the output.
If you are familiar with how to edit the crontab
with ```crontab -e```, then go ahead and add one or both of the following line to the crontab file 
```
00,15,30,45 * * * * /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_ideal.sh
00,15,30,45 * * * * /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_meas.sh
```
If you do not know how to edit the crontab this way, then use the following commands:
```bash
crontab -l > /Users/codar/crontab_backup_copy.txt
crontab -l > /Users/codar/mycron.txt
```
Open the mycron.txt file in whatever editor you choose. It may be an empty file if there were no scheduled jobs.
Add one or both of these lines 
```
00,15,30,45 * * * * /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_ideal.sh
00,15,30,45 * * * * /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_meas.sh
```
and save the file.
 
Then install the updated file with this command:
```bash
crontab /Users/codar/mycron.txt
```
If you encountered any errors with crontab and need to start over, you can always restore your original file with this command:
```bash
crontab /Users/codar/crontab_backup_copy.txt
```

After the task runs, you should see new radial files being generated in /Codar/SeaSonde/Data/RadialShorts_qcd and
/Codar/SeaSonde/Data/Radials_qcd




## 6. Archivalist Setup (very important!)

You must set up archiving for all of these new files, otherwise you risk filling up your computer's hard disk.

The following command will set up a radial metric task list for Archivalist. (Please note that if you already have the Archivalist_RadialMetric.plist file in RadialConfigs, this will overwrite your existing file.):
```bash
   cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/Archivalist_RadialMetric.plist /Codar/SeaSonde/Configs/RadialConfigs/Archivalist_RadialMetric.plist
```

In the SeaSonde Archivalist application, check that you have a task list called RadialMetric 
with tasks for the following files and make any adjustments for your preferences as needed:

See the CODAR Archivalist application guide if you are unfamiliar with the Archivalist program and would like further instructions on how to edit the archives.

/Codar/SeaSonde/Data/RadialMetrics/IdealPattern                   
/Codar/SeaSonde/Data/RadialMetrics/MeasPattern                   
/Codar/SeaSonde/Data/Radials_qcd/IdealPattern                   
/Codar/SeaSonde/Data/Radials_qcd/MeasPattern                   
/Codar/SeaSonde/Data/RadialShorts_qcd/IdealPattern                   
/Codar/SeaSonde/Data/RadialShorts_qcd/MeasPattern                   
/Codar/SeaSonde/Data/RadialResponses/IdealPattern   (see note below)            
/Codar/SeaSonde/Data/RadialResponses/MeasPattern    (see note below) 

Note: The RadialResponses tasks may not be necessary as later versions of the SeaSonde radial software do 
not include these files.  If the RadialResponses folders contain data, it is *extremely important* to 
set up archiving tasks for those files, because they take up alot of space and will fill up your disk quickly!



## Background

For information about radial metric QC processing, the tests, and the test thresholds:

Haines, Sara, Harvey Seim, and Mike Muglia. "Implementing quality control of high-frequency radar estimates and application to Gulf Stream surface currents." Journal of Atmospheric and Oceanic Technology 34.6 (2017): 1207-1224.

### Notes

| QC Function        | Description |
| -----------        | ----------- |
| Threshold Tests    | badflag any values that fall below or above a single threshold |
| Weighted Averaging | average several values with weights based on signal quality parameters |

#### QC Threshold Tests:
1. DOA peak power (MSR1, MDR1, MDR2) < 5 dB default 
1. DOA 1/2 power width (3dB down) (MSW1, MDW1, MDW2) > 50 deg default
1. SNR on monopole (MA3S) < 5 dB default
1. SNR on both loops (MA1S and MA2S) < 5 dB

#### Weighted Averaging:
1. Weighting based on Music Power (MSP1, MDP1, MDP2)
1. Weighting based on SNR on monopole (MA3S)
1. No weight function (None) 
 
### System Requirements
Python 3 is required. However, earlier versions of some of the other packages may be okay.

- Python >=3.7, <3.11
- numpy
    - https://pypi.python.org/pypi/numpy
    - Data read into memory are stored in the N-dimensional array datatype (ndarray) for indexing and computation.
- geopy
    - https://pypi.python.org/pypi/geopy
    - geopy.distance.geodesic()
    - Used to compute (LAT, LON) based on range and bearing from site origin in generating RadialShorts file
- watchdog
    - Used to monitor a directory for new files and trigger qc and merge process when new RadialMetric file is created (NOT YET IMPLEMENTED)
- docopt
- hfradarpy

### CODAR Software Requirements
- CODAR SeaSonde RadialSuite 7.x -- supports RadialMetric output out of the box
- CODAR SeaSonde RadialSuite 8.x -- does not support RadialMetric output unless requested
   - Requires special key file specifically to enable RadialMetric output (contact CODAR to get)
   - Requires RadialMetric R2 Addon (contact CODAR to get)
- CODAR SeaSonde RadialSuite 21 and 22 -- does not support RadialMetric output unless requested
   - Requires special key file specifically to enable RadialMetric output (contact CODAR to get)
   - DOES NOT require a RadialMetric Addon
-
    - /Codar/SeaSonde/Apps/Bin/LLUVMerger.app
    - Used to merge spatial and temporal RadialShorts data to final Radial
