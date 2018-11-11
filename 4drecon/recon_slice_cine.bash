#!/usr/bin/env bash


# RECON SINGLE SLICE CINE VOLUMES

# e.g., bash recon_slice_vol.bash ~/path/to/top/level/recon/directory/

# TODO: add FORCEEXCLUDEFRAME as used in other recon scripts

# Input 

RECONDIR=$1


# Check that Recon Directory Exists

if [ ! -d "$RECONDIR" ]; then
  echo directory $RECONDIR does not exist
  exit 1
else


# Manage Paths to Allow Queueing of Jobs using TaskSpooler

ORIG_PATH=$(pwd)
SCRIPT_PATH=$(dirname `which $0`)

RECONSLICECINEDIR=$RECONDIR/slice_cine_vol
MASKDIR=$RECONSLICECINEDIR/mask
mkdir -p $MASKDIR
cd $MASKDIR

echo RECON SINGLE SLICE CINE VOLUMES
echo $RECONSLICECINEDIR


# Variables 

RECON="cine.nii.gz"
STACKS="../../data/s*_rlt_ab.nii.gz"
REFVOL="../../dc_vol/mask.nii.gz"
MASKCINEVOL="../mask/mask_vol.nii.gz"
MASKCARDSYNC="../mask/mask_cardsync.nii.gz"
THICKNESSFILE="../../data/slice_thickness.txt"
RESOLUTION=1.25
NMC=0
NSRLAST=10
NUMCARDPHASE=25
STACKDOFDIR="../../dc_vol/stack_transformations"
SLICEDOFDIR="../../dc_vol/slice_transformations"
NUMSLICEENLARGE=20  # num. slices
ROIBLURCINEVOL=4  # mm
ROIBLURTHRESH=0.002699796063  # 3*sigma
ROIBLURCINEVOLSIGMA=$(echo "$ROIBLURCINEVOL/3" | bc -l)
NUMROICLOSINGITER=9
MEANRRFILE="../../cardsync/mean_rrinterval.txt"
RRINTERVALSFILE="../../cardsync/rrintervals.txt"
CARDPHASESFILE="../../cardsync/cardphases_intraslice_cardsync.txt"


# Setup

ITER=$(($NMC+1))
THICKNESS=$(cat $THICKNESSFILE)
MEANRR=$(cat $MEANRRFILE)
RRINTERVALS=$(cat $RRINTERVALSFILE)
CARDPHASES=$(cat $CARDPHASESFILE)
NUMSTACK=$(ls -1 ../../data/s*_dc_ab.nii.gz | wc -l);
NUMSLICE=$(eval "wc -w $RRINTERVALSFILE | awk -F'[ ]' '{print \$1}'" )
NUMFRAME=$(eval "wc -w $CARDPHASESFILE | awk -F'[ ]' '{print \$1}'" )
STACKFILES=($STACKS)
STACKMASKALL="../../mask/s*_mask_heart.nii.gz"
STACKMASKFILES=(../../mask/s*_mask_heart.nii.gz)


# Identify Number of Slices in Each Stack

declare -a ARRAYNUMSLICEINSTACK
STACKINDEX=0
for STACK in ${STACKFILES[@]}
do
  CMD="info $STACK | grep \"Image dimensions\" | awk -F'[ ]' '{print \$6}'"
  NUMSLICEINSTACK=$(eval $CMD)
  ARRAYNUMSLICEINSTACK[$STACKINDEX]=$NUMSLICEINSTACK
  ((STACKINDEX++))
done


# Generate Masks

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
    CMD=$(printf 'region %s %s -Rz1 %i -Rz2 %i' $STACKMASK $SLICEMASK $(($STACKSLICEINDEX-1)) $STACKSLICEINDEX)
    eval $CMD
    SLICEDOF=$(printf '%s/transformation%05i.dof' $SLICEDOFDIR $SLICEINDEX)
    CMD="headertool $SLICEMASK $SLICEMASK -dofin $SLICEDOF"
    eval $CMD > /dev/null
    ((SLICEINDEX++))
  done
  ((STACKINDEX++))
done
echo "   combining" $STACKMASKALL "using" $SLICEDOFDIR"/transformation*.dof"
enlarge_image $REFVOL tgt.nii.gz -x 20 -y 20 -z 20
combineimages tgt.nii.gz ../../mask/s*_mask_heart_slice*.nii.gz $MASKCINEVOL > /dev/null
threshold $MASKCINEVOL $MASKCINEVOL 0 > /dev/null
echo "   blurring ROI with radius =" $ROIBLURCINEVOL "mm (sigma =" $ROIBLURCINEVOLSIGMA "mm)"
SIGMA=$(echo "$ROIBLURCINEVOLSIGMA+2*$RESOLUTION" | bc -l)  # NOTE: additional blurring+erosion to smooth mask
blur $MASKCINEVOL $MASKCINEVOL $SIGMA -3D > /dev/null
threshold $MASKCINEVOL $MASKCINEVOL $ROIBLURTHRESH > /dev/null 
erosion $MASKCINEVOL $MASKCINEVOL -iterations 4 > /dev/null 
closing $MASKCINEVOL $MASKCINEVOL -iterations $NUMROICLOSINGITER > /dev/null
rm ../../mask/s*_mask_heart_slice*.nii.gz
rm tgt.nii.gz
cp $MASKCINEVOL $MASKCARDSYNC
dilation $MASKCINEVOL $MASKCINEVOL -iterations 4


# Recon Slice Cines

STACKINDEX=0
SLICEINDEX=0
for STACKMASK in ${STACKMASKFILES[@]}
do
  CMD="info $STACKMASK | grep \"Voxel dimensions are\" | awk -F'[ ]' '{print \$4}'"
  DX=$(eval $CMD)
  CMD="info $STACKMASK | grep \"Voxel dimensions are\" | awk -F'[ ]' '{print \$5}'"
  DY=$(eval $CMD)
  DZ=18  # TODO: fix this to read from THICKNESSFILE, .e.g, CMD=$(printf 'cat %s | awk -F'[ ]' '{print \$%i}'' $THICKNESSFILE $STACKINDEX); DZ=$(eval $CMD);
  if [ "$STACKINDEX" = "0" ]; then
    echo "NOTE: using hard-coded dz = $DZ mm for all slice masks volumes"
  fi
  for STACKSLICEINDEX in $(seq 1 ${ARRAYNUMSLICEINSTACK[$STACKINDEX]})
  do
    SLICEDIR=$(printf '%s/slice%05i' $RECONSLICECINEDIR $SLICEINDEX)
    mkdir -p $SLICEDIR
    cd $SLICEDIR
    echo $(printf 'slice%05i' $SLICEINDEX)
    # Identify Excluded Slices
    declare -a EXCLUDEDSLICES
    NUMEXCLUDEDSLICES=0
    ISLC=0
    for ISTK in $(seq 1 $NUMSTACK)
    do
      for ISTKSLC in $(seq 1 ${ARRAYNUMSLICEINSTACK[$ISTK-1]})
      do
        if [ ! "$ISLC" -eq "$SLICEINDEX" ]; then
          EXCLUDEDSLICES[$NUMEXCLUDEDSLICES]=$ISLC
          ((NUMEXCLUDEDSLICES++))
        fi
        ((ISLC++))
      done
    done
    # Recon Cine
    echo "  reconstructing $RECON"
    CMD="reconstructionCardiac $RECON $NUMSTACK $STACKS -thickness $THICKNESS -dofin $STACKDOFDIR/stack-transformation*.dof -slice_transformations $SLICEDOFDIR -mask $MASKCINEVOL -iterations $ITER -rec_iterations_last $NSRLAST -resolution $RESOLUTION -force_exclude_sliceloc $NUMEXCLUDEDSLICES ${EXCLUDEDSLICES[@]} -exclude_slices_only -numcardphase $NUMCARDPHASE -rrinterval $MEANRR -rrintervals $NUMSLICE $RRINTERVALS -cardphase $NUMFRAME $CARDPHASES -debug > log-main.txt"
    echo $CMD > recon.bash
    eval $CMD
    # Extract Heart Mask
    echo "  generating mask_heart_slice.nii.gz"
    CMD=$(printf 'region %s mask_heart_slice.nii.gz -Rz1 %i -Rz2 %i' $STACKMASK $(($STACKSLICEINDEX-1)) $STACKSLICEINDEX)
    echo $CMD >> recon.bash
    eval $CMD
    CMD=$(printf 'headertool mask_heart_slice.nii.gz mask_heart_slice.nii.gz -size %.3f %.3f %.3f -dofin %s/transformation%05i.dof' $DX $DY $DZ $SLICEDOFDIR $SLICEINDEX)
    echo $CMD >> recon.bash
    eval $CMD
    CMD="transformation mask_heart_slice.nii.gz mask_heart_slice.nii.gz -target cine.nii.gz -nn -matchInputType"
    echo $CMD >> recon.bash
    eval $CMD
    # Transform Cardsync Mask to Cine Volume
    transformation $MASKCARDSYNC mask_cardsync.nii.gz -target $RECON -linear -matchInputType
    # Clean Up
    CMD="rm transformation*.dof; rm info_mc*.tsv; rm log-registration*.txt; rm *_mc*sr*; rm average*.nii.gz; rm bias*.nii.gz; rm cropped*.nii.gz; rm init_mc*.nii.gz; rm maskformatching*.nii.gz; rm reconstructed_mc*.nii.gz; rm rescaledstack*.nii.gz; rm simstack*.nii.gz; rm stack*.nii.gz; rm weight*.nii.gz; rm mask0*.nii.gz;"
    echo $CMD >> recon.bash
    eval $CMD
    # Iterate Slice Index
    ((SLICEINDEX++))
  done
  ((STACKINDEX++))
done


# Finish

echo "volume reconstruction complete"

cd $ORIG_PATH


fi

