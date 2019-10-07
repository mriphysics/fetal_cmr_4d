#!/usr/bin/env bash


# RECON CINE VOLUME

# e.g., bash recon_cine_vol.bash ~/path/to/top/level/recon/directory/ cine_vol


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
STACKS="../data/s*_rlt_ab.nii.gz"
REFVOL="../mask/mask_chest.nii.gz"
MASKCINEVOL="mask_cine_vol.nii.gz"
MASKCHEST="../mask/mask_chest.nii.gz"
THICKNESSFILE="../data/slice_thickness.txt"
EXCLUDESTACKFILE="../data/force_exclude_stack.txt"
EXCLUDESLICEFILE="../data/force_exclude_slice.txt"
EXCLUDEFRAMEFILE="../data/force_exclude_frame.txt"
RESOLUTION=1.25
NMC=2 # TAR --- adjusted to 2 so that ITER = 3. ITER > 3 crashes MIRTK
NSR=10
NSRLAST=20
NUMCARDPHASE=25
STACKDOFDIR="../dc_vol/stack_transformations"
SLICEDOFDIR="../dc_vol/slice_transformations"
NUMSLICEENLARGE=20  # num. slices
ROIBLURCINEVOL=10  # mm
ROIBLURTHRESH=0.002699796063 # 3*sigma
ROIBLURCINEVOLSIGMA=$(echo "$ROIBLURCINEVOL/3" | bc -l)
NUMROICLOSINGITER=9
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
STACKFILES=($STACKS)
STACKMASKALL="../mask/s*_mask_heart.nii.gz"
STACKMASKFILES=(../mask/s*_mask_heart.nii.gz)
EXCLUDESTACK=$(cat $EXCLUDESTACKFILE)
NUMEXCLUDESTACK=$(eval "wc -w $EXCLUDESTACKFILE | awk -F' ' '{print \$1}'" )
EXCLUDESLICE=$(cat $EXCLUDESLICEFILE)
NUMEXCLUDESLICE=$(eval "wc -w $EXCLUDESLICEFILE | awk -F' ' '{print \$1}'" )
EXCLUDEFRAME=$(cat $EXCLUDEFRAMEFILE)
NUMEXCLUDEFRAME=$(eval "wc -w $EXCLUDEFRAMEFILE | awk -F' ' '{print \$1}'" )


# Identify Number of Slices in Each Stack

declare -a ARRAYNUMSLICEINSTACK
STACKINDEX=0
for STACK in ${STACKFILES[@]}
do
  CMD="mirtk info $STACK | grep \"Image dimensions\" | awk -F' ' '{print \$6}'"
  NUMSLICEINSTACK=$(eval $CMD)
  ARRAYNUMSLICEINSTACK[$STACKINDEX]=$NUMSLICEINSTACK
  ((STACKINDEX++))
done


# Generate Mask

echo generating mask: $MASKCINEVOL
STACKINDEX=0
SLICEINDEX=0
for STACKMASK in ${STACKMASKFILES[@]}
do
  for STACKSLICEINDEX in $(seq 1 ${ARRAYNUMSLICEINSTACK[$STACKINDEX]})
  do
    STACKMASKSUFFIX='_mask_heart'
    SLICEMASKSUFFIX=$(printf '%s_slice%02i' $STACKMASKSUFFIX $STACKSLICEINDEX)
    SLICEMASK=${STACKMASK/$STACKMASKSUFFIX/$SLICEMASKSUFFIX}
    CMD=$(printf 'mirtk extract-image-region %s %s -Rz1 %i -Rz2 %i' $STACKMASK $SLICEMASK $(($STACKSLICEINDEX-1)) $(($STACKSLICEINDEX-1)))
    eval $CMD
    SLICEDOF=$(printf '%s/transformation%05i.dof' $SLICEDOFDIR $SLICEINDEX)
    CMD="mirtk edit-image $SLICEMASK $SLICEMASK -dofin $SLICEDOF"
    eval $CMD > /dev/null
    ((SLICEINDEX++))
  done
  ((STACKINDEX++))
done

echo "   combining" $STACKMASKALL "using" $SLICEDOFDIR"/transformation*.dof"

mirtk combine_masks  $REFVOL ../mask/s*_mask_heart_slice*.nii.gz $MASKCINEVOL > /dev/null
mirtk threshold_image $MASKCINEVOL $MASKCINEVOL 0 > /dev/null

echo "   blurring ROI with radius =" $ROIBLURCINEVOL "mm (sigma =" $ROIBLURCINEVOLSIGMA "mm)"
SIGMA=$(echo "$ROIBLURCINEVOLSIGMA+2*$RESOLUTION" | bc -l)  # NOTE: additional blurring+erosion to smooth mas

mirtk smooth-image $MASKCINEVOL $MASKCINEVOL $SIGMA -3D > /dev/null
mirtk threshold_image $MASKCINEVOL $MASKCINEVOL $ROIBLURTHRESH > /dev/null
mirtk erode-image $MASKCINEVOL $MASKCINEVOL -iterations 2 > /dev/null
mirtk close-image $MASKCINEVOL $MASKCINEVOL -iterations $NUMROICLOSINGITER > /dev/null

rm ../mask/s*_mask_heart_slice*.nii.gz

mirtk transform-image $MASKCHEST mask_chest_vol.nii.gz -target $MASKCINEVOL -interp NN > /dev/null
mirtk pad_image $MASKCINEVOL mask_chest_vol.nii.gz $MASKCINEVOL 0 0 > /dev/null

# Recon Cine Volume

echo reconstructing cine volume: $RECON
CMD="mirtk reconstructCardiac $RECON $NUMSTACK $STACKS -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -slice_transformations $SLICEDOFDIR -mask $MASKCINEVOL -iterations $ITER -rec_iterations $NSR -rec_iterations_last $NSRLAST -resolution $RESOLUTION -force_exclude_stack $NUMEXCLUDESTACK $EXCLUDESTACK -force_exclude_sliceloc $NUMEXCLUDESLICE $EXCLUDESLICE -force_exclude $NUMEXCLUDEFRAME $EXCLUDEFRAME -numcardphase $NUMCARDPHASE -rrinterval $MEANRR -rrintervals $NUMSLICE $RRINTERVALS -cardphase $NUMFRAME $CARDPHASES -debug > log-main.txt"
echo $CMD > recon.bash
eval $CMD


# > log-main.txt

# Clean Up

echo "" >> recon.bash
CMD="mkdir transformations; mv transformation*.dof transformations;"
echo $CMD >> recon.bash
eval $CMD
echo "" >> recon.bash
CMD="mkdir sr_iterations; mv *_mc*sr* sr_iterations;"
echo $CMD >> recon.bash
eval $CMD


# Finish

echo "volume reconstruction complete"

cd $ORIG_PATH





fi

