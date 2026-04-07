#!/bin/bash

# Prepare run directories, submission scripts, and launch adjoint experiments

# Variables defined in clear_run.batch.sh
export CMASK="___CMASK___"
export CMASKDIR="___CMASKDIR___"
export TMASK="___TMASK___"
export PERTLABEL="___PERTLABEL___"
export LABEL="___LABEL___"

# Variables held constant across all experiments
export BASEDIR="/work/n01/n01/afstyles/MYECCO_fork/MITgcm/ECCOV4/release4/"
export SCRATCH="/mnt/lustre/a2fs-nvme/work/n01/n01/afstyles/"
export ADJOINT=true
export DATA="data"
export DATA_ECCO="data.ecco.highfreq.template"
export DATA_CTRL="data.ctrl.highfreq.noinit"
export DATA_CAL="data.cal.2004"
export DATA_EXCH2="data.exch2.360"
export DATA_PKG="data.pkg.nodiag"
export DATA_LAYERS="data.layers"
export DATA_DIAGNOSTICS="data.diagnostics.lowfreqperts"
export DAILY_CTRL=true

# Mask configuration
export MASKTYPE='C'  # Specify masktype. 'C' for boxmeans or 'WS' for m_horflux

# Pickup files (make sure to update data file, DON'T USE PICKUPS IN ADJOINT MODE)
export PICKUPDIR="${BASEDIR}/pickup/ecco_layers_26y/"
export PICKUPNO="none"  #0000035052  #Set to none if you don't want to use a pickup
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

module load cray-python

#if [ "${PERTLABEL}" = none ] ; then
#   export LABEL=${CMASK}.${TMASK}
#else
#   export LABEL=${PERTLABEL}.${CMASK}.${TMASK}
#fi

if [ "${ADJOINT}" = true ] ; then
   export RUNDIR="run_ad.${LABEL}"
   export EXECUTABLE="mitgcmuv_ad"
else
   export RUNDIR="run.${LABEL}"
   export EXECUTABLE="mitgcmuv"
fi

echo ${RUNDIR}

if [ -d "${SCRATCH}/${RUNDIR}" ]; then #If run directory exists clear it
   #rm -rv ${SCRATCH}/${RUNDIR}/*
   echo "${SCRATCH}/${RUNDIR}/* would be deleted normally"
else
   mkdir -v ${SCRATCH}/${RUNDIR} #Otherwise create new directory
fi

cd ${SCRATCH}/${RUNDIR}


# Link the data files
ln -s ${BASEDIR}/input_init/NAMELIST/* .
ln -s ${BASEDIR}/input_init/error_weight/ctrl_weight/* .
ln -s ${BASEDIR}/input_init/error_weight/data_error/* .
ln -s ${BASEDIR}/input_init/* .
ln -s ${BASEDIR}/input_init/tools/* .
ln -s ${BASEDIR}/input_ecco/*/* .
ln -s ${BASEDIR}/input_forcing/eccov4r4* .

# Replace specific files with custom versions
rm data
cp -v ${BASEDIR}/namelist/${DATA} data

rm data.ecco 
cp -v ${BASEDIR}/namelist/${DATA_ECCO} data.ecco
sed -i "s/___MASKNAME___/${CMASK}/g" data.ecco

rm data.ctrl
cp -v ${BASEDIR}/namelist/${DATA_CTRL} data.ctrl

rm data.cal
cp -v ${BASEDIR}/namelist/${DATA_CAL} data.cal

rm data.exch2
cp -v ${BASEDIR}/namelist/${DATA_EXCH2} data.exch2

rm data.pkg
cp -v ${BASEDIR}/namelist/${DATA_PKG} data.pkg

rm data.layers
cp -v ${BASEDIR}/namelist/${DATA_LAYERS} data.layers 

rm data.diagnostics
cp -v ${BASEDIR}/namelist/${DATA_DIAGNOSTICS} data.diagnostics

# Copy over mask files
cp -v ${BASEDIR}/my_masks/${TMASK}T .
ln -sv ${TMASK}T ${CMASK}T

if [ "${MASKTYPE}" == C ] ; then
   cp -v ${BASEDIR}/my_masks/${CMASKDIR}/${CMASK}C .
else
   cp -v ${BASEDIR}/my_masks/${CMASKDIR}/${CMASK}W .
   cp -v ${BASEDIR}/my_masks/${CMASKDIR}/${CMASK}S .
fi

# Add additional xx_* files which are just zeros
if [ "${DAILY_CTRL}" = false ] ; then
   cp -v ${BASEDIR}/xx_zero.0000000129.data .
   cp -v ${BASEDIR}/xx_zero.0000000129.meta .
else
   ln -sv ${BASEDIR}/xx_zero2ddaily.0000000129.data xx_zero.0000000129.data
   ln -sv ${BASEDIR}/xx_zero2ddaily.0000000129.meta xx_zero.0000000129.meta
fi

cp -v ${BASEDIR}/xx_zero2dparam.0000000129.data .
cp -v ${BASEDIR}/xx_zero2dparam.0000000129.meta .

cp -v ${BASEDIR}/xx_zero3dparam.0000000129.data .
cp -v ${BASEDIR}/xx_zero3dparam.0000000129.meta .

# Copy over pickup file (if necessary)
if [ "${PICKUPNO}" != none ] ; then
   cp -v ${PICKUPDIR}/pickup*.${PICKUPNO}.data .
   cp -v ${PICKUPDIR}/pickup*.${PICKUPNO}.meta .
fi

ln -sv xx_zero.0000000129.data xx_qnet.0000000129.data
ln -sv xx_zero.0000000129.data xx_empmr.0000000129.data
ln -sv xx_zero.0000000129.data xx_tauu.0000000129.data
ln -sv xx_zero.0000000129.data xx_tauv.0000000129.data
ln -sv xx_zero2dparam.0000000129.data xx_bottomdrag.0000000129.data
ln -sv xx_zero2dparam.0000000129.data xx_nlbottomdrag.0000000129.data

ln -sv xx_zero.0000000129.meta xx_qnet.0000000129.meta
ln -sv xx_zero.0000000129.meta xx_empmr.0000000129.meta
ln -sv xx_zero.0000000129.meta xx_tauu.0000000129.meta
ln -sv xx_zero.0000000129.meta xx_tauv.0000000129.meta
ln -sv xx_zero2dparam.0000000129.meta xx_bottomdrag.0000000129.meta
ln -sv xx_zero2dparam.0000000129.meta xx_nlbottomdrag.0000000129.meta

# Replace zero xx_* files if perturbing fwd run
if [ "${PERTLABEL}" != none ] ; then
   for FILE in ${BASEDIR}/perts/xx_*.data.${PERTLABEL}; do
      echo 'FILENAME: ' ${FILE}
      LINKNAME=$(basename ${FILE})
      LINKNAME=${LINKNAME%.*}
      rm -v ${LINKNAME}
      echo 'LINKNAME: ' ${LINKNAME}
      cp -v ${FILE} ${LINKNAME}
   done
fi

# Add weights_one files
ln -sv ${BASEDIR}/weights_ones.data .

# Add additional climatology files
ln -sv ${BASEDIR}/S_OWPv1_M_eccollc_90x50.bin .
ln -sv ${BASEDIR}/T_OWPv1_M_eccollc_90x50.bin .

# Copy executable file
if [ "${ADJOINT}" = true ] ; then
   cp -pv ${BASEDIR}/build_ad/${EXECUTABLE} .
else
   cp -pv ${BASEDIR}/build/${EXECUTABLE} .
fi

#Replace and run mkdir_subdir_diags.py
rm mkdir_subdir_diags.py
ln -s ${BASEDIR}/mkdir_subdir_diags.py .
python mkdir_subdir_diags.py

module unload cray-python

# Return to BASE directory
cd ${BASEDIR}

#Prepare submission scripts
cp -v runCommand.template.sh runCommand.${LABEL}.sh
cp -v run_submit.template.slurm run_submit.${LABEL}.slurm

sed -i "s/___RUNDIR___/${RUNDIR}/g" runCommand.${LABEL}.sh
sed -i "s/___RUNSUBMIT___/run_submit.${LABEL}.slurm/g" runCommand.${LABEL}.sh
sed -i "s/___RUNNAME___/e4_${LABEL}/g" runCommand.${LABEL}.sh

sed -i "s/___RUNDIR___/${RUNDIR}/g" run_submit.${LABEL}.slurm
sed -i "s/___EXECUTABLE___/${EXECUTABLE}/g" run_submit.${LABEL}.slurm

