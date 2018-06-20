#!/bin/bash
set -euo pipefail
set -x

#Set local parallelism inherited from QBATCH
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS_PER_COMMAND:-$(nproc)}

#Locate priors for processing
#TODO Allow external specification
BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONBRAINMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"

#Mangle inputs
#Order is output file, then T1, then any others
output=$1
originput=$2
shift 2
multispectral_inputs=("$@")

#Add setting to allow override
tmpdir=$(mktemp -d)

#Internal resampled input used for processing
input=${tmpdir}/input.mnc

#If lesion mask exists, negate it to produce a multiplicative exlcusion mask
if [[ -s $(dirname ${originput})/$(basename ${originput} _t1.mnc)_lesion.mnc ]]; then
  ImageMath 3 ${tmpdir}/lesion.mnc Neg $(dirname ${originput})/$(basename ${originput} _t1.mnc)_lesion.mnc
  lesionmask=${tmpdir}/lesion.mnc
else
  lesionmask=""
fi

#Function used to do bias field correction
function do_N4_correct {
  #input fov mask weight maxval output bias shrink
  #Do first round of masked bias field correction, use brain mask as weight
  N4BiasFieldCorrection ${N4_VERBOSE:+--verbose} -d 3 -s $8 -w $4 -x $2 \
    -b [200] -c [300x300x300x300,1e-5] --histogram-sharpening [0.05,0.01,200] -i $1 -o [$6,$7] -r 0

  #Normalize and rescale intensity
  ImageMath 3 $6 TruncateImageIntensity $6 0.0005 0.9995 1024 $3
  ImageMath 3 $6 RescaleImage $6 0 $5
}

#Convert classify image into a mask
function classify_to_mask {
  #Breakup classification and drop try to exclude small misclassified chunks
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class3.mnc 3 3 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class2.mnc 2 2 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class1.mnc 1 1 1 0
  ImageMath 3 ${tmpdir}/${n}/class3.mnc GetLargestComponent ${tmpdir}/${n}/class3.mnc
  ImageMath 3 ${tmpdir}/${n}/class2.mnc GetLargestComponent ${tmpdir}/${n}/class2.mnc
  iMath 3 ${tmpdir}/${n}/class1.mnc ME ${tmpdir}/${n}/class1.mnc 1 1 box 1
  ImageMath 3 ${tmpdir}/${n}/class1.mnc GetLargestComponent ${tmpdir}/${n}/class1.mnc
  iMath 3 ${tmpdir}/${n}/class1.mnc MD ${tmpdir}/${n}/class1.mnc 1 1 ball 1

  #Reconstruct a closed mask from the classification
  ImageMath 3 ${tmpdir}/${n}/classifymask.mnc addtozero ${tmpdir}/${n}/class2.mnc ${tmpdir}/${n}/class3.mnc
  ImageMath 3 ${tmpdir}/${n}/classifymask.mnc addtozero ${tmpdir}/${n}/classifymask.mnc ${tmpdir}/${n}/class1.mnc
  iMath 3 ${tmpdir}/${n}/classifymask.mnc ME ${tmpdir}/${n}/classifymask.mnc 2 1 ball 1
  ImageMath 3 ${tmpdir}/${n}/classifymask.mnc GetLargestComponent ${tmpdir}/${n}/classifymask.mnc
  iMath 3 ${tmpdir}/${n}/classifymask.mnc MD ${tmpdir}/${n}/classifymask.mnc 3 1 ball 1
  iMath 3 ${tmpdir}/${n}/classifymask.mnc MC ${tmpdir}/${n}/classifymask.mnc 5 1 ball 1
  ImageMath 3  ${tmpdir}/${n}/classifymask.mnc m ${tmpdir}/${n}/classifymask.mnc ${tmpdir}/${n}/hotmask.mnc
  ImageMath 3 ${tmpdir}/${n}/classifymask.mnc FillHoles ${tmpdir}/${n}/classifymask.mnc 2

}

function cleanup_posteriors {
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class3.mnc 3 3 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class2.mnc 2 2 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class1.mnc 1 1 1 0
  for num in 1 2 3; do
    ImageMath 3 ${tmpdir}/${n}/class${num}.mnc GetLargestComponent ${tmpdir}/${n}/class${num}.mnc
    ImageMath 3 ${tmpdir}/${n}/SegmentationPosteriors${num}.mnc m ${tmpdir}/${n}/SegmentationPosteriors${num}.mnc ${tmpdir}/${n}/class${num}.mnc
    ImageMath 3 ${tmpdir}/${n}/SegmentationPosteriors${num}.mnc m ${tmpdir}/${n}/SegmentationPosteriors${num}.mnc ${tmpdir}/${n}/classifymask.mnc
  done
}


#Generate list of extra Atropos inputs for multispectral segmentation
multispectral_atropos_inputs=""
if (( ${#multispectral_inputs[@]} > 0 )); then
  for file in "${multispectral_inputs[@]}"; do
    multispectral_atropos_inputs+="-a ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc "
  done
fi

#Find maximum value of scan to rescale to for final output
maxval=$(mincstats -max -quiet ${originput})

#Isotropize, crop, and pad input volume
isostep=1
ResampleImage 3 ${originput} ${input} ${isostep}x${isostep}x${isostep} 0 4
mincmath -quiet -clamp -const2 0 ${maxval} ${input} ${tmpdir}/input.clamp.mnc
mv -f ${tmpdir}/input.clamp.mnc ${input}
ImageMath 3 ${tmpdir}/cropmask.mnc ThresholdAtMean ${input} 1
ExtractRegionFromImageByMask 3 ${input} ${tmpdir}/input.crop.mnc ${tmpdir}/cropmask.mnc 1 10
mv -f ${tmpdir}/input.crop.mnc ${input}

if [[ ! -z ${lesionmask} ]]; then
  antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${lesionmask} -r ${input} -n GenericLabel -o ${lesionmask}
fi

#Compute final round shink factor for N4
dx=$(mincinfo -attvalue xspace:step ${originput})
dy=$(mincinfo -attvalue yspace:step ${originput})
dz=$(mincinfo -attvalue zspace:step ${originput})
shrink_final_round=$(python -c "import math; print(max(2,int(math.ceil(2.0 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0)))))")

#Generate a whole-image mask to force N4 to always do correction over whole image
minccalc -quiet -unsigned -byte -expression 'A[0]?1:1' ${input} ${tmpdir}/initmask.mnc
minccalc -quiet -unsigned -byte -expression 'A[0]>0?1:0' ${input} ${tmpdir}/initweight.mnc
################################################################################
#Round 0, N4 across areas greater than 1% of mean
################################################################################
n=0
mkdir -p ${tmpdir}/${n}

#Correct entire image domain
ImageMath 3 ${tmpdir}/${n}/weight.mnc ThresholdAtMean ${input} 1
iMath 3 ${tmpdir}/${n}/weight.mnc ME ${tmpdir}/${n}/weight.mnc 3 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/weight.mnc 4 1 ball 1

do_N4_correct ${input} ${tmpdir}/initmask.mnc ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc ${maxval} ${tmpdir}/${n}/corrected.mnc ${tmpdir}/${n}/bias.mnc 8

################################################################################
#Round 1, N4 across areas greater than 1% of mean, repeat
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc

#Correct entire image domain
ImageMath 3 ${tmpdir}/${n}/weight.mnc ThresholdAtMean ${tmpdir}/${n}/input.mnc 1
iMath 3 ${tmpdir}/${n}/weight.mnc ME ${tmpdir}/${n}/weight.mnc 3 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/weight.mnc 4 1 ball 1

do_N4_correct ${input} ${tmpdir}/initmask.mnc ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc ${maxval} ${tmpdir}/${n}/corrected.mnc ${tmpdir}/${n}/bias.mnc 8

minccalc -zero -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print(float(\"$(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc)\") / float(\"$(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)\"))" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Round 2, N4 with brain mask intersected with Otsu mask
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
minc_anlm --mt ${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS} ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

#Register to MNI space
antsRegistration ${N4_VERBOSE:+--verbose} -d 3 --float 1 --collapse-output-transforms 1 --minc  \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform [${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1] \
  --transform Rigid[0.5] \
  --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,32,Regular,0.95] \
  --convergence [409600x51200x6400,1e-6,10] \
  --shrink-factors 16x8x4 \
  --smoothing-sigmas 6.7945744023041525x3.3972872011520763x1.6986436005760381mm \
  --masks [NULL,NULL] \
  --transform Similarity[0.1] \
  --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,64,Regular,0.95] \
  --convergence [51200x6400x800,1e-6,10] \
  --shrink-factors 8x4x2 \
  --smoothing-sigmas 3.3972872011520763x1.6986436005760381x0.8493218002880191mm \
  --masks [NULL,NULL] \
  --transform Affine[0.1] \
  --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,128,Regular,0.95] \
  --convergence [6400x800x100,1e-6,10] \
  --shrink-factors 4x2x1 \
  --smoothing-sigmas 1.6986436005760381x0.8493218002880191x0.42466090014400953mm \
  --masks [${REGISTRATIONBRAINMASK},NULL]

#Repeat with nuyl matched registration
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${REGISTRATIONBRAINMASK} -o ${tmpdir}/${n}/mnimask.mnc -n GenericLabel
iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${tmpdir}/${n}/mnimask.mnc 4 1 ball 1
minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/${n}/shrinkmask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 0 --fix_zero_padding --steps 1024

antsRegistration ${N4_VERBOSE:+--verbose} -d 3 --float 1 --collapse-output-transforms 1 --minc \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform ${tmpdir}/${n}/mni0_GenericAffine.xfm \
  --transform Affine[0.1] \
  --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,256,Regular,0.95] \
  --convergence [800x100x100x100,1e-6,10] \
  --shrink-factors 2x1x1x1 \
  --smoothing-sigmas 0.8493218002880191x0.42466090014400953x0.21233045007200477x0mm \
  --masks [${REGISTRATIONBRAINMASK},NULL]

antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/input.mnc -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc -r ${REGISTRATIONMODEL}

#BSpline[5] does weird things to intensity, clip back to positive range
mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/${n}/mni.mnc) ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

#Shrink the MNI mask for the first intensity matching
rm -f ${tmpdir}/${n}/shrinkmask.mnc
iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${REGISTRATIONBRAINMASK} 4 1 ball 1

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/${n}/shrinkmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

#Run a quick beast to get a brain mask
mincbeast ${N4_VERBOSE:+-verbose} -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.2mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

#Resample beast mask and MNI mask to native space
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${REGISTRATIONBRAINMASK} -o ${tmpdir}/${n}/mnimask.mnc -n GenericLabel
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${tmpdir}/${n}/beastmask.mnc -o ${tmpdir}/${n}/bmask.mnc -n GenericLabel

#Combine the masks because sometimes beast misses badly biased cerebellum
ImageMath 3 ${tmpdir}/${n}/mask.mnc addtozero ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc

#Expand the mask a bit
iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 3 1 ball 1

#Create hotmask and exlcude hot voxels from weight
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  0 $(mincstats -quiet -pctT 99.95 -mask ${tmpdir}/${n}/mask_D.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 1 0

#Exclude hotspots to avoid kmeans classifying skull as brain
ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/${n}/hotmask.mnc

ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/kmeans.mnc Kmeans 2 ${tmpdir}/${n}/mask_D.mnc
ThresholdImage 3 ${tmpdir}/${n}/kmeans.mnc ${tmpdir}/${n}/weight.mnc 2 3 1 0
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc

ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

do_N4_correct ${input} ${tmpdir}/initmask.mnc  ${tmpdir}/${n}/mask.mnc ${tmpdir}/${n}/primary_weight.mnc ${maxval} ${tmpdir}/${n}/corrected.mnc ${tmpdir}/${n}/bias.mnc 6

#Align multispectral inputs into T1 space
if (( ${#multispectral_inputs[@]} > 0 )); then
    echo "Needs to be reimplemented"
else
  cp -f ${tmpdir}/initweight.mnc ${tmpdir}/global_exclude.mnc
fi

minccalc -zero -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print(float(\"$(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc)\") / float(\"$(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)\"))" >>${tmpdir}/convergence.txt
#python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Round 3, N4 with nonlinearly MNI-bootstrapped WM/GM segmentation proabilities
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
minc_anlm --mt ${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}  ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/$((n - 1))/mask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 0 --fix_zero_padding --steps 1024

#Affine register to MNI space, tweak registration
antsRegistration ${N4_VERBOSE:+--verbose} -d 3 --float 1 --collapse-output-transforms 1 --minc \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform ${tmpdir}/$((n - 1))/mni0_GenericAffine.xfm \
  --transform Affine[0.1] \
  --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,256,Regular,0.95] \
  --convergence [800x100x100x100,1e-6,10] \
  --shrink-factors 2x1x1x1 \
  --smoothing-sigmas 0.8493218002880191x0.42466090014400953x0.21233045007200477x0mm \
  --masks [${REGISTRATIONBRAINMASK},${tmpdir}/$((n - 1))/mask_D.mnc]

antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/input.mnc -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc -r ${REGISTRATIONMODEL}

mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/${n}/mni.mnc) ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

#Shrink last round's beastmask for normalization
iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${tmpdir}/$((n - 1))/beastmask.mnc 5 1 ball 1

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/${n}/shrinkmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

#Run a quick beast to get a brain mask
mincbeast ${N4_VERBOSE:+-verbose} -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.1mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

antsApplyTransforms -i ${tmpdir}/${n}/beastmask.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/bmask.mnc ${N4_VERBOSE:+--verbose} -d 3 -n GenericLabel
antsApplyTransforms -i ${REGISTRATIONBRAINMASK} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/mniaffinemask.mnc ${N4_VERBOSE:+--verbose} -d 3 -n GenericLabel
cp -f ${tmpdir}/${n}/bmask.mnc ${tmpdir}/bmask.mnc

iMath 3 ${tmpdir}/${n}/nonlinregmask.mnc MD ${tmpdir}/${n}/bmask.mnc 1 1 ball 1

#Non linearly register priors
antsRegistration ${N4_VERBOSE:+--verbose} -d 3 --float 1 --minc \
  --output [${tmpdir}/${n}/nonlin] \
  --initial-moving-transform ${tmpdir}/${n}/mni0_GenericAffine.xfm \
  --use-histogram-matching 0 \
  --transform SyN[0.1,3,0] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,256,Regular,0.95] \
  --convergence [102400x68600x43200x25000x12800x5400x1600x200x50x50,1e-6,10] \
  --shrink-factors 16x14x12x10x8x6x4x2x1x1 \
  --smoothing-sigmas 6.7945744023041525x5.945252602016134x5.095930801728114x4.2466090014400955x3.3972872011520763x2.547965400864057x1.6986436005760381x0.8493218002880191x0.42466090014400953x0mm \
  --masks [${REGISTRATIONBRAINMASK},${tmpdir}/${n}/nonlinregmask.mnc]

#Resample MNI Priors to Native space for classification
antsApplyTransforms -i ${WMPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior3.mnc ${N4_VERBOSE:+--verbose} -d 3 -n Linear
antsApplyTransforms -i ${GMPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior2.mnc ${N4_VERBOSE:+--verbose} -d 3 -n Linear
antsApplyTransforms -i ${CSFPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior1.mnc ${N4_VERBOSE:+--verbose} -d 3 -n Linear
#Masks
antsApplyTransforms -i ${REGISTRATIONBRAINMASK} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/mnimask.mnc ${N4_VERBOSE:+--verbose} -d 3 -n GenericLabel

cp -f ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/mnimask.mnc

#Combine the masks because sometimes beast misses badly biased cerebellum
ImageMath 3 ${tmpdir}/${n}/mask.mnc addtozero ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc

#Expand the mask a bit
iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 2 1 ball 1

#Create hotmask
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  0 $(mincstats -quiet -pctT 99.95 -mask ${tmpdir}/${n}/mask.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 1 0

#Exclude lesions and such
ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/global_exclude.mnc ${lesionmask}
#ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m  ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/${n}/hotmask.mnc

#Do an initial classification using the MNI priors
Atropos ${N4_VERBOSE:+--verbose} -d 3 -x ${tmpdir}/${n}/mask_D.mnc -c [5,0] -a ${tmpdir}/${n}/input.nuyl.mnc ${multispectral_atropos_inputs}  \
  -i PriorProbabilityImages[3,${tmpdir}/${n}/SegmentationPrior%d.mnc,0.25] -k HistogramParzenWindows -m [0.1,1x1x1] \
  -o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[0] 1 --winsorize-outliers BoxPlot

classify_to_mask

#Combine GM and WM probably images into a N4 mask,
ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/classifymask.mnc #${lesionmask} ${tmpdir}/global_exclude.mnc
#SmoothImage 3 ${tmpdir}/${n}/weight.mnc 0.5 ${tmpdir}/${n}/weight.mnc 0 0
ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc
ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc RescaleImage ${tmpdir}/${n}/primary_weight.mnc 0 1

cleanup_posteriors
do_N4_correct ${input} ${tmpdir}/initmask.mnc  ${tmpdir}/${n}/mask.mnc ${tmpdir}/${n}/primary_weight.mnc ${maxval} ${tmpdir}/${n}/corrected.mnc ${tmpdir}/${n}/bias.mnc 4

if (( ${#multispectral_inputs[@]} > 0 )); then
    echo "Need to implement"
fi

minccalc -zero -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print(float(\"$(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc)\") / float(\"$(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)\"))" >>${tmpdir}/convergence.txt
#python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Remaining rounds, N4 with segmentations bootstrapped from prior run
################################################################################
atropos_prior_weight=0.25
while true; do
  ((++n))
  mkdir -p ${tmpdir}/${n}

  cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
  minc_anlm --mt ${ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS}  ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
  mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

  minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/$((n - 1))/mask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 0 --fix_zero_padding --steps 1024

  #Resample beast mask and MNI mask to native space
  cp -f ${tmpdir}/mnimask.mnc ${tmpdir}/${n}/mnimask.mnc

  #Combine the masks because sometimes beast misses badly biased cerebellum
  cp -f ${tmpdir}/$((n - 1))/classifymask.mnc ${tmpdir}/${n}/mask.mnc
  iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 1 1 ball 1

  #Create hotmask
  ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
    0 $(mincstats -quiet -pctT 99.95 -mask ${tmpdir}/${n}/mask_D.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 1 0

  ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/global_exclude.mnc ${lesionmask}
  #ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/${n}/hotmask.mnc

  #Do an initial classification using the last round posteriors, remove outliers
  Atropos ${N4_VERBOSE:+--verbose} -d 3 -x ${tmpdir}/${n}/mask_D.mnc -c [5,0.0] -a ${tmpdir}/${n}/input.nuyl.mnc ${multispectral_atropos_inputs} -s 1x2 -s 2x3 \
    -i PriorProbabilityImages[3,${tmpdir}/$((n - 1))/SegmentationPosteriors%d.mnc,${atropos_prior_weight}] -k HistogramParzenWindows -m [0.1,1x1x1] \
    -o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[1] --winsorize-outliers BoxPlot

  classify_to_mask

  #Combine GM and WM probably images into a N4 mask,
  ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc
  ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/classifymask.mnc #${lesionmask} ${tmpdir}/global_exclude.mnc
  #SmoothImage 3 ${tmpdir}/${n}/weight.mnc 0.5 ${tmpdir}/${n}/weight.mnc 0 0
  ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc
  ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc RescaleImage ${tmpdir}/${n}/primary_weight.mnc 0 1

  cleanup_posteriors
  do_N4_correct ${input} ${tmpdir}/initmask.mnc  ${tmpdir}/${n}/mask.mnc ${tmpdir}/${n}/primary_weight.mnc ${maxval} ${tmpdir}/${n}/corrected.mnc ${tmpdir}/${n}/bias.mnc 4

  #Loop over multispectral inputs
  if (( ${#multispectral_inputs[@]} > 0 )); then
    echo "Need to implement"
  fi

  #Compute coeffcient of variation
  minccalc -zero -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
  python -c "print(float(\"$(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc)\") / float(\"$(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)\"))" >>${tmpdir}/convergence.txt
  #python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt

  rm -rf ${tmpdir}/$((n - 1))

  # Maximum number of iterations is 10, CV difference less than 0.005
  [[ (${n} -lt 10) && ($(python -c "print $(tail -1 ${tmpdir}/convergence.txt) > 0.005") == "True") ]] || break

done

echo "Convergence results:"
cat ${tmpdir}/convergence.txt

#Transform all the working files into the original input space
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/primary_weight.mnc -o ${tmpdir}/finalweight.mnc -r ${originput} -n Linear
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/mask.mnc -o ${tmpdir}/finalmask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/bmask.mnc -o ${tmpdir}/finalbmask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/mnimask.mnc -o ${tmpdir}/finalmnimask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/classifymask.mnc -o ${tmpdir}/finalclassifymask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms ${N4_VERBOSE:+--verbose} -d 3 --float 1 -i ${tmpdir}/${n}/classify.mnc -o ${tmpdir}/finalclassify.mnc -r ${originput} -n GenericLabel

#Create a FOV mask for the original input
minccalc -quiet -unsigned -byte -expression 'A[0]?1:1' ${originput} ${tmpdir}/originitmask.mnc

do_N4_correct ${originput} ${tmpdir}/originitmask.mnc ${tmpdir}/finalmask.mnc ${tmpdir}/finalweight.mnc ${maxval} ${tmpdir}/corrected.mnc ${tmpdir}/bias.mnc ${shrink_final_round}

if (( ${#multispectral_inputs[@]} > 0 )); then
    echo "Need to implement"
fi

cp -f ${tmpdir}/corrected.mnc ${output}
if [[ -z "${N4_STANDALONE-}" ]]; then
  cp -f ${tmpdir}/finalbmask.mnc $(dirname $output)/$(basename $output .mnc).beastmask.mnc
  cp -f ${tmpdir}/finalmnimask.mnc $(dirname $output)/$(basename $output .mnc).mnimask.mnc
  cp -f ${tmpdir}/finalclassifymask.mnc $(dirname $output)/$(basename $output .mnc).classifymask.mnc
  cp -f ${tmpdir}/finalclassify.mnc $(dirname $output)/$(basename $output .mnc).classify.mnc
  cp -f ${tmpdir}/finalmask.mnc $(dirname $output)/$(basename $output .mnc).mask.mnc
fi

rm -rf ${tmpdir}
