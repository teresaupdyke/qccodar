#!/bin/bash

cd /Users/codar/Downloads

timestamp=`date "+%Y%m%d_%H%M%S"`
if [ $? -ne 0 ]; then
    backup_dir="/Users/codar/qccodar_backup"
else
    backup_dir="/Users/codar/qccodar_backup/${timestamp}"
fi
echo "Backup directory is ${backup_dir}"

# Create directories for qccodar files
echo "Creating directories for qccodar files..."
mkdir -p /Users/codar/qccodar_files/logs
mkdir -p ${backup_dir}

# Detect architecture (Apple Silicon or Intel)
ARCH=$(uname -m)

# Detect the default shell (zsh or bash)
SHELL_TYPE=$(basename "$SHELL")
if [[ "$SHELL_TYPE" == "zsh" ]]; then
      echo "Detected zsh shell."
      SHELL_FILE=".zshrc"
elif [[ "$SHELL_TYPE" == "bash" ]]; then
      echo "Detected bash shell."
      SHELL_FILE=".bash_profile"
fi
BACKUP_CONFIG_FILE=${backup_dir}/backup${SHELL_FILE}
CONFIG_FILE=/Users/codar/${SHELL_FILE}
if [ -f "${CONFIG_FILE}" ]; then
   cp "${CONFIG_FILE}" "${BACKUP_CONFIG_FILE}" 
else
   touch "${BACKUP_CONFIG_FILE}" 
fi

condaversion="none"
if [ -d "/Users/codar/miniforge3" ]; then
  condaversion="miniforge3"
elif [ -d "/Users/codar/mambaforge3" ]; then 
  condaversion="mambaforge3"
elif [ -d "/Users/codar/miniconda3" ]; then 
  condaversion="miniconda3"  
fi

# Check to see if miniforge is already installed. We don't want to automatically overwrite an existing installation.
if [[ ${condaversion} != "none" ]]; then
    echo "A minimal installation of conda called ${condaversion} is already installed and will be used by qccodar."
else
    echo "Installing Miniforge..."

    if [[ "$ARCH" == "arm64" ]]; then
        echo "Detected Apple Silicon (M1/M2). Downloading Miniforge for arm64..."
        curl -o /Users/codar/Downloads/miniforge3.sh -L -H "User-Agent: Safari/537.36" https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-24.11.3-0-MacOSX-arm64.sh
    else
        echo "Detected Intel architecture. Downloading Miniforge for x86_64..."
        curl -o /Users/codar/Downloads/miniforge3.sh -L -H "User-Agent: Safari/537.36" https://github.com/conda-forge/miniforge/releases/download/24.11.3-0/Miniforge3-24.11.3-0-MacOSX-x86_64.sh
    fi

    # Prompt user to accept license agreement and install location manually
    echo " "
    echo "*******************************************"
    echo " "
    echo "PLEASE READ"
    echo "Miniforge installer will now run." 
    echo "Follow the instructions, accept the license, and press ENTER to confirm installation location as /Users/codar/miniforge3"
    echo "Afterwards, you may be asked if you wish to update your shell profile to automatically initialize conda.  Answering either yes or no should be okay"
    echo "Press any key to continue..."
    read -n 1 -s


    # Install Miniforge
    chmod 755 miniforge3.sh
    ./miniforge3.sh


    # Detect the default shell (zsh or bash)
    SHELL_TYPE=$(basename "$SHELL")
    if [[ "$SHELL_TYPE" == "zsh" ]]; then
      echo "Configuring .zshrc for conda."
      cp ${BACKUP_CONFIG_FILE} "${CONFIG_FILE}" 
      INIT_COMMAND='source /Users/codar/.zshrc'

      # Add conda initialization to the shell configuration file
      echo "Updating $CONFIG_FILE to initialize conda..."
cat <<EOT >> "$CONFIG_FILE"    
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/codar/miniforge3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/codar/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/Users/codar/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/codar/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/Users/codar/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/Users/codar/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
EOT
    
    else
      echo "Configuring .bash_profile for conda."
      cp ${BACKUP_CONFIG_FILE} "${CONFIG_FILE}"
      INIT_COMMAND='source /Users/codar/.bash_profile'

      # Add conda initialization to the shell configuration file
      echo "Updating $CONFIG_FILE to initialize conda..."
cat <<EOT >> "$CONFIG_FILE"    
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/codar/miniforge3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/codar/miniforge3/etc/profile.d/conda.sh" ]; then
        . "/Users/codar/miniforge3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/codar/miniforge3/bin:$PATH"
    fi
fi
unset __conda_setup

if [ -f "/Users/codar/miniforge3/etc/profile.d/mamba.sh" ]; then
    . "/Users/codar/miniforge3/etc/profile.d/mamba.sh"
fi
# <<< conda initialize <<<
EOT

    fi


    # Apply the changes
echo " "
echo "*******************************************"
echo " "
#echo "PLEASE READ"
#    echo "This shouldn't happen BUT if you encounter an unexpected error in this next step (an error related to conda that might involve a plug-in) then WAIT 1 MINUTE for the script to continue.  After the install is over, you must edit your shell profile following step 2 of the manual installation instructions in the README file, otherwise you may have problems with CODAR radial processing."
#    echo "Press any key to continue..."
#    read -n 1 -s
    echo "Sourcing $CONFIG_FILE..."
    $INIT_COMMAND

fi

if [ -d "/Users/codar/qccodar_files/qccodar-main" ]; then
  echo " "
  echo "*******************************************"
  read -p "/Users/codar/qccodar_files/qccodar-main already exists. Do you want to overwrite your existing copy of qccodar (yes or no)? " ans_overwrite

     # Convert the input to lowercase for simplicity
     ans_overwrite=$(echo "$choice2" | tr '[:upper:]' '[:lower:]')


     if [[ "$ans_overwrite" == "y" || "$ans_overwrite" == "yes" ]]; then
         rm -rf /Users/codar/qccodar_files/qccodar-main
         # Download qccodar from GitHub
         echo "Downloading qccodar code from GitHub..."
         cd /Users/codar/qccodar_files
         curl -L https://github.com/rowg/qccodar/archive/refs/heads/main.zip -o /Users/codar/qccodar_files/qccodar-main.zip
         unzip qccodar-main.zip
         rm -rf qccodar-main.zip
         break
     else
         echo "Keeping previously installed copy of qccodar."
     fi
else

         # Download qccodar from GitHub
         echo "Downloading qccodar code from GitHub..."
         cd /Users/codar/qccodar_files
         curl -L https://github.com/rowg/qccodar/archive/refs/heads/main.zip -o /Users/codar/qccodar_files/qccodar-main.zip
         unzip qccodar-main.zip
         rm -rf qccodar-main.zip
   
fi




# Set up the Conda environment for qccodar
echo "Setting up Conda environment for qccodar..."
if [[ "$SHELL_TYPE" == "zsh" ]]; then
    conda init zsh
elif [[ "$SHELL_TYPE" == "bash" ]]; then
    conda init bash
fi
source ~/${condaversion}/etc/profile.d/conda.sh

cd /Users/codar/qccodar_files/qccodar-main
conda env create -f environment.yml
conda activate qccodar

if [ -d "/Users/codar/${condaversion}/envs/qccodar" ]; then
    echo "qccodar environment has been created."
else
    echo "Failed to create qccodar environment."
    exit 1
fi

# Install qccodar in the environment
echo "Installing qccodar..."
pip install .
conda deactivate

if [ -f "/Users/codar/${condaversion}/envs/qccodar/bin/qccodar" ]; then
    echo "Installation successful"
else
    echo "Installation failed. Aborting..."
    exit 1
fi

# Radial Metric Output Configuration
echo "Configuring Radial Metric Output..."
# Check Header.txt file to get the frequency to set system type (5 MHz, 13 MHz, or 25 MHz)
FREQ_OPTIONS=(5 13 25)
FREQ=`grep "! 7" /Codar/SeaSonde/Configs/RadialConfigs/Header.txt`
HEADER_FREQ=${FREQ:0:5}
HEADER_FREQ=$(printf "%.0f" "$HEADER_FREQ")
echo "Rounded frequency from Header.txt file is $HEADER_FREQ"

# Initialize min_value to a high value and min_index to 0
min_value=999
min_index=0

# Loop through the array to find the minimum value and its index
for i in "${!FREQ_OPTIONS[@]}"; do
  num=$(($HEADER_FREQ - FREQ_OPTIONS[$i]))
  # Check if number is negative
  if (( num < 0 )); then
    abs=$(( -num ))
  else
    abs=$num
  fi
  if (( abs < min_value )); then
    #echo "Index $i: Absolute Difference $abs is less than current min of $min_value"
    min_value=${abs}
    #echo "Min has been updated to ${min_value}"
    min_index=$i
  fi
done
#echo $min_index
SYSTEM_FREQ=${FREQ_OPTIONS[$min_index]}
echo "The closest nominal operating frequency is $SYSTEM_FREQ"


if [ "$SYSTEM_FREQ" == "5" ]; then
    cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_5MHz.plist /Users/codar/qccodar_files/qccodar.plist
elif [ "$SYSTEM_FREQ" == "13" ]; then
    cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_13MHz.plist /Users/codar/qccodar_files/qccodar.plist
elif [ "$SYSTEM_FREQ" == "25" ]; then
    cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/qccodar_25MHz.plist /Users/codar/qccodar_files/qccodar.plist
else
    echo "Error: Invalid system frequency selected!"
fi

# Edit AnalysisOptions.txt to enable Radial Metric Output
ANALYSIS_OPTIONS_FILE="/Codar/SeaSonde/Configs/RadialConfigs/AnalysisOptions.txt"
cp ${ANALYSIS_OPTIONS_FILE} ${backup_dir}/.
if [ -f "$ANALYSIS_OPTIONS_FILE" ]; then
    echo "Editing line 21 of AnalysisOptions.txt to enable Radial Metric Output..."
    sed -i '' '21s/.*/1           !21 Enable Radial Metric Output: 0(Off), 1(Enable), 2(Enable MaxVel.)/' "$ANALYSIS_OPTIONS_FILE"
else
    echo "Error: AnalysisOptions.txt not found at $ANALYSIS_OPTIONS_FILE"
fi

# Crontab setup for Realtime QC
echo " "
echo "*******************************************"
echo " "
echo "Setting up crontab for Realtime QC..."
crontab -l > ${backup_dir}/crontab_backup.txt
CRON_FILE="/Users/codar/mycron.txt"

crontab -l > "$CRON_FILE" 2>/dev/null || touch "$CRON_FILE"

#Remove any previous qccodar entries to avoid duplicates
sed -i '' '/qccodar/d' "$CRON_FILE"
crontab "$CRON_FILE"

#Add the new qccodar entries
cat <<EOT >> $CRON_FILE
00,15,30,45 * * * * source /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_ideal.sh
00,15,30,45 * * * * source /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar_meas.sh
EOT

crontab $CRON_FILE
rm -rf $CRON_FILE
echo "Crontab entries added. Contents of crontab: "
crontab -l

# Ensure that the scripts have permissions to execute
chmod 755 /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/run_qccodar*.sh

# Archivalist setup
echo " "
echo "Setting up Archivalist for Radial Metric tasks..."
if [ -f "/Codar/SeaSonde/Configs/RadialConfigs/Archivalist_RadialMetric.plist" ]; then
     echo " "
     echo "*******************************************"
     read -p "Do you want use qccodar recommended settings and overwrite your existing Archivalist_RadialMetric.plist (yes or no)? " choice2

     # Convert the input to lowercase for simplicity
     choice2=$(echo "$choice2" | tr '[:upper:]' '[:lower:]')

     if [[ "$choice2" == "y" || "$choice2" == "yes" ]]; then
         cp /Codar/SeaSonde/Configs/RadialConfigs/Archivalist_RadialMetric.plist ${backup_dir}/.
         cp /Users/codar/qccodar_files/qccodar-main/src/qccodar/config/Archivalist_RadialMetric.plist /Codar/SeaSonde/Configs/RadialConfigs/Archivalist_RadialMetric.plist
         echo "Archivalist settings installed. After install, please verify the Archivalist tasks in the SeaSonde Archivalist application!"

     else
         echo " "
         echo "PLEASE READ"
         echo "Keeping previously installed Archivalist_RadialMetric.plist.  Be sure these settings will not fill up your hard disk! Please pay careful attention to how you handle antenna response files!  Also, do not keep more metric files in the RadialMetric source directories than short files in the RadialShorts_qcd source directories!"
     fi
   
else
    echo "Archivalist_RadialMetric.plist has NOT been modified.  Be sure the settings will not fill up your hard disk! Please pay careful attention to how you handle antenna response files!  Also, do not keep more metric files in the RadialMetric source directories than short files in the RadialShorts_qcd source directories!"
fi

echo " "
echo " "
    echo "Press any key to continue..."
    read -n 1 -s

echo " "
echo "*******************************************"
echo " "
echo "PLEASE READ"
echo "Backup copies of the shell profile, the crontab, AnalysisOptions.txt, and Archivalist_RadialMetric.plist (if available) can be found in ${backup_dir}  Use if you need to restore a version of a file that existed before qccodar was installed!"
echo " "
    echo "Press any key to continue..."
    read -n 1 -s

# Final message
echo " "
echo "*******************************************"
echo " "
echo "qccodar has been successfully installed and configured for Realtime QC."
echo "Log files are located at /Users/codar/qccodar_files/logs/"
echo "Edit /Users/codar/qccodar_files/qccodar.plist if you need to use custom settings/thresholds."  
