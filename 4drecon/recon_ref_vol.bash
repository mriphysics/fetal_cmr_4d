#!/usr/bin/env bash


# RECON REF VOLUME

# e.g., bash recon_ref_vol.bash ~/path/to/top/level/recon/directory/ ref_vol


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

echo RECON DC VOLUME
echo $RECONVOLDIR


# Variables 

RECON=$VOLDESC.nii.gz
STACKS="../data/s*_dc_ab.nii.gz"
THICKNESS=$(cat ../data/slice_thickness.txt)
MASKSTACKRREG="mask_stack_rreg.nii.gz"
MASKSLICERREG="mask_slice_rreg.nii.gz"
MASKDCVOL="mask_vol.nii.gz"
TGTSTACKNO=$(cat ../data/tgt_stack_no.txt)
EXCLUDESTACKFILE="../data/force_exclude_stack.txt"
EXCLUDESLICEFILE="../data/force_exclude_slice.txt"
RESOLUTION=1.25
NMC=3
NSR=10
NSRLAST=20
NUMCARDPHASE=1
STACKDOFDIR="stack_transformations"
ROIBLURSTACKRREG=35 # mm
ROIBLURDCVOL=80  # mm
ROIBLURTHRESH=0.002699796063  # 3*sigma
ROIBLURSTACKRREGSIGMA=$(echo "$ROIBLURSTACKRREG/3" | bc -l)
ROIBLURDCVOLSIGMA=$(echo "$ROIBLURDCVOL/3" | bc -l)
NUMSLICEENLARGE=40
NUMROICLOSINGITER=9

echo reconstructing DC volume: $RECONVOLDIR/$RECON


# Setup

ITER=$(($NMC+1))
NUMSTACK=$(ls -1 ../data/s*_dc_ab.nii.gz | wc -l);
STACKMASKALL="../mask/s*_mask_heart.nii.gz"
STACKMASKFILES=(../mask/s*_mask_heart.nii.gz)
STACKMASKTGT=${STACKMASKFILES[$TGTSTACKNO-1]}
EXCLUDESTACK=$(cat $EXCLUDESTACKFILE)
NUMEXCLUDESTACK=$(eval "wc -w $EXCLUDESTACKFILE | awk -F'[ ]' '{print \$1}'" )
EXCLUDESLICE=$(cat $EXCLUDESLICEFILE)
NUMEXCLUDESLICE=$(eval "wc -w $EXCLUDESLICEFILE | awk -F'[ ]' '{print \$1}'" )
echo "   target stack no.: "$TGTSTACKNO
echo "   target stack mask: "$STACKMASKTGT


# Generate Stack-Stack RReg Mask

echo generating stack-stack rreg mask: $MASKSTACKRREG
enlarge_image $STACKMASKTGT $MASKSTACKRREG -z $NUMSLICEENLARGE > /dev/null
resample $MASKSTACKRREG $MASKSTACKRREG -size $RESOLUTION $RESOLUTION $RESOLUTION > /dev/null
combineimages $MASKSTACKRREG $STACKMASKALL $MASKSTACKRREG > /dev/null
threshold $MASKSTACKRREG $MASKSTACKRREG 0 > /dev/null
echo "   blurring ROI with radius =" $ROIBLURSTACKRREG "mm (sigma =" $ROIBLURSTACKRREGSIGMA "mm)"
SIGMA=$(echo "$ROIBLURSTACKRREGSIGMA+2*$RESOLUTION" | bc -l)  # NOTE: additional blurring+erosion to smooth mask
blur $MASKSTACKRREG $MASKSTACKRREG $SIGMA -3D > /dev/null
threshold $MASKSTACKRREG $MASKSTACKRREG $ROIBLURTHRESH > /dev/null 
erosion $MASKSTACKRREG $MASKSTACKRREG -iterations 4 > /dev/null 
closing $MASKSTACKRREG $MASKSTACKRREG -iterations $NUMROICLOSINGITER > /dev/null

# Estimate Stack-Stack Transformations

CMD="reconstructionCardiac $RECON $NUMSTACK $STACKS -thickness $THICKNESS -stack_registration -target_stack $TGTSTACKNO -mask $MASKSTACKRREG -iterations 1 -rec_iterations_last 8 -resolution $RESOLUTION -force_exclude_stack $NUMEXCLUDESTACK $EXCLUDESTACK -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -force_exclude $NUMEXCLUDEFRAME $EXCLUDEFRAME -no_robust_statistics -numcardphase $NUMCARDPHASE -debug > log-main.txt"
echo estimating stack-stack transformations: $CMD
eval $CMD


# Clean Up

mkdir -p $STACKDOFDIR;
mv $RECON $STACKDOFDIR;
mv $MASKSTACKRREG $STACKDOFDIR;
mv mask.nii.gz $STACKDOFDIR;
mv log*.txt $STACKDOFDIR;
mv stack-transformation0*.dof $STACKDOFDIR;
mv info.tsv $STACKDOFDIR;
rm *.*


# Generate Volume Recon Mask

echo volume reconstruction mask: $MASKDCVOL
enlarge_image $STACKMASKTGT $MASKDCVOL -z $NUMSLICEENLARGE > /dev/null
resample $MASKDCVOL $MASKDCVOL -size $RESOLUTION $RESOLUTION $RESOLUTION > /dev/null
combineimages $MASKDCVOL $STACKMASKALL $MASKDCVOL > /dev/null
threshold $MASKDCVOL $MASKDCVOL 0 > /dev/null
echo "   blurring ROI with radius =" $ROIBLURDCVOL "mm (sigma =" $ROIBLURDCVOLSIGMA "mm)"
SIGMA=$(echo "$ROIBLURDCVOLSIGMA+2*$RESOLUTION" | bc -l)  # NOTE: additional blurring+erosion to smooth mask
blur $MASKDCVOL $MASKDCVOL $SIGMA -3D > /dev/null
threshold $MASKDCVOL $MASKDCVOL $ROIBLURTHRESH > /dev/null 
erosion $MASKDCVOL $MASKDCVOL -iterations 4 > /dev/null 
closing $MASKDCVOL $MASKDCVOL -iterations $NUMROICLOSINGITER > /dev/null


# Reconstruct Reference Volume

CMD="reconstructionCardiac $RECON $NUMSTACK $STACKS -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -mask $MASKDCVOL -iterations 1 -rec_iterations_last 8 -resolution $RESOLUTION -force_exclude_stack $NUMEXCLUDESTACK $EXCLUDESTACK -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -no_robust_statistics -numcardphase $NUMCARDPHASE -debug > log-main.txt"
echo reconstructing volume: $CMD
eval $CMD


# Clean Up

mkdir -p tmp;
mv $RECON tmp;
mv $MASKDCVOL tmp;
mv mask.nii.gz tmp;
mv log*.txt tmp;
mv info.tsv tmp;
rm *.*;
mv tmp/* .;
rm -r tmp;


# Finish

echo "volume reconstruction complete"

cd $ORIG_PATH


fi




