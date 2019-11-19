
#!/usr/bin/env bash

# set -x

# RECON VEL VOLUME

# e.g., bash recon_vel_vol.bash ~/path/to/top/level/recon/directory/ vel_vol


# Input 

RECONDIR=$1
VOLDESC=$2


# Check that Recon Directory Exists

if [ ! -d "$RECONDIR" ]; then
  echo directory $RECONDIR does not exist
  exit 1
else


# Manage Paths to Allow Queueing of Jobs using TaskSpooler

ORIG_PATH=$(pwd)
SCRIPT_PATH=$(dirname `which $0`)

RECONVOLDIR=$RECONDIR/$VOLDESC
mkdir -p $RECONVOLDIR
cd $RECONVOLDIR

echo RECON CINE VOLUME
echo $RECONVOLDIR


# Variables 

RECON=$VOLDESC.nii.gz
STACKS="../data/s*_rlt_ph_corr_uterus.nii.gz"
GRADMOMVALSFILE="../data/grad_moment_vals.txt"
GRADMOMDIRSFILE="../data/grad_moment_dirs.txt"
THICKNESSFILE="../data/slice_thickness.txt"
EXCLUDESLICEFILE="../data/force_exclude_slice.txt"
EXCLUDECINEVOLFILE="../data/force_exclude_cine_vol.txt"
RESOLUTION=1.25
NMC=0
NSR=40
NSRLAST=20
NUMCARDPHASE=25
ALPHA=3
STACKDOFDIR="../dc_vol/stack_transformations"
CINETRANSDIR="../cine_vol/transformations"
MASKCINEVOL="../cine_vol/mask_cine_vol.nii.gz"
MEANRRFILE="../cardsync/mean_rrinterval.txt"
RRINTERVALSFILE="../cardsync/rrintervals.txt"
CARDPHASESFILE="../cardsync/cardphases_interslice_cardsync.txt"


# Setup

ITER=$(($NMC+1))
THICKNESS=$(cat $THICKNESSFILE)
MEANRR=$(cat $MEANRRFILE)
RRINTERVALS=$(cat $RRINTERVALSFILE)
CARDPHASES=$(cat $CARDPHASESFILE)
NUMSTACK=$(ls -1 ../data/s*_dc_ab.nii.gz | wc -l);
NUMSLICE=$(eval "wc -w $RRINTERVALSFILE | awk -F' ' '{print \$1}'" )
NUMFRAME=$(eval "wc -w $CARDPHASESFILE | awk -F' ' '{print \$1}'" )

EXCLUDESLICE=$(cat $EXCLUDESLICEFILE)
NUMEXCLUDESLICE=$(eval "wc -w $EXCLUDESLICEFILE | awk -F' ' '{print \$1}'" )
EXCLUDECINEVOL=$(cat $EXCLUDECINEVOLFILE)
NUMEXCLUDECINEVOL=$(eval "wc -w $EXCLUDECINEVOLFILE | awk -F' ' '{print \$1}'" )


# Recon Velocity Cine Volume

echo reconstructing velocity cine volume:
CMD="mirtk reconstructCardiacVelocity $NUMSTACK $STACKS $GRADMOMVALSFILE $GRADMOMDIRSFILE -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -transformations $CINETRANSDIR -mask $MASKCINEVOL -alpha $ALPHA -limit_intensities -rec_iterations $NSR -resolution $RESOLUTION -force_exclude_stack 0 -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -force_exclude $NUMEXCLUDECINEVOL $EXCLUDECINEVOL -numcardphase $NUMCARDPHASE -rrinterval $MEANRR -rrintervals $NUMSLICE $RRINTERVALS -cardphase $NUMFRAME $CARDPHASES -debug > log-main.txt"
echo $CMD > recon.bash
eval $CMD


# Clean Up
mkdir sr_iterations
mv addon-velocity-*.nii.gz sr_iterations/
mv confidence-map-velocity*.nii.gz sr_iterations/
mv sum-velocity-*.nii.gz sr_iterations/
mv sr_iterations/sum-velocity-final.nii.gz .
mv velocity-0*.nii.gz sr_iterations/
mv velocity-1*.nii.gz sr_iterations/
mv velocity-2*.nii.gz sr_iterations/


# Finish

echo "velocity volume reconstruction complete"

cd $ORIG_PATH


fi

