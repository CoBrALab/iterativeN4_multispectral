#!/bin/bash
set -euo pipefail
set -x

export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=${THREADS_PER_COMMAND:-$(nproc)}

BEASTLIBRARY_DIR="${QUARANTINE_PATH}/resources/BEaST_libraries/combined"
REGISTRATIONMODEL="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c.mnc"
REGISTRATIONBRAINMASK="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_t1_tal_nlin_sym_09c_mask.mnc"
WMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_wm_tal_nlin_sym_09c.mnc"
GMPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_gm_tal_nlin_sym_09c.mnc"
CSFPRIOR="${QUARANTINE_PATH}/resources/mni_icbm152_nlin_sym_09c_minc2/mni_icbm152_csf_tal_nlin_sym_09c.mnc"

#Converting FWHM to Sigma
blur16mm=6.79457440230415234177
blur14mm=5.94525260201613329905
blur12mm=5.09593080172811425633
blur10mm=4.24660900144009521361
blur8mm=3.39728720115207617089
blur6mm=2.54796540086405712816
blur4mm=1.69864360057603808544
blur3mm=1.27398270043202856408
blur2mm=0.84932180028801904272
blur1mm=0.42466090014400952136
blur05mm=0.21233045007200476068

output=$1
shift
input=$1
shift
inputs=("$@")

tmpdir=$(mktemp -d)
#tmpdir=temp

#If lesion mask exists, negate it to produce a multiplicative exlcusion mask
if [[ -s $(dirname $input)/$(basename $input _t1.mnc)_lesion.mnc ]]; then
  ImageMath 3 $tmpdir/lesion.mnc Neg $(dirname $input)/$(basename $input _t1.mnc)_lesion.mnc
  lesionmask=$tmpdir/lesion.mnc
else
  lesionmask=""
fi

#Generate list of extra Atropos inputs for multispectral segmentation
extra_atropos_inputs=""
if (( ${#inputs[@]} > 0 )); then
  mkdir -p ${tmpdir}/multispectral
  for file in "${inputs[@]}"; do
    extra_atropos_inputs+="-a ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc "
  done
fi

#Find maximum value of scan to rescale to
maxval=$(mincstats -max -quiet ${input})

#Calculate shrink values for N4
#Target is 4mm steps for first round, 2mm for all subsequent
dx=$(mincinfo -attvalue xspace:step ${input})
dy=$(mincinfo -attvalue yspace:step ${input})
dz=$(mincinfo -attvalue zspace:step ${input})

#shrinkround1=$(python -c "import math; print max(4,int(math.ceil(4 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")
#shrinkround2=$(python -c "import math; print max(2,int(math.ceil(2 / ( ( abs($dx) + abs($dy) + abs($dz) ) / 3.0))))")
shrinkround1=$(python -c "import math; print max(4,int(math.ceil(4 / ( min([abs($dx),abs($dy),abs($dz)]) ) )))")
shrinkround2=$(python -c "import math; print max(2,int(math.ceil(2 / ( min([abs($dx),abs($dy),abs($dz)]) ) )))")

#Isotropize and pad input volume
isostep=$(python -c "print min([abs($dx),abs($dy),abs($dz)])")
ResampleImage 3 $input ${tmpdir}/input.mnc ${isostep}x${isostep}x${isostep} 0 4
ImageMath 3 ${tmpdir}/input.mnc PadImage ${tmpdir}/input.mnc 10
mincmath -quiet -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/input.mnc) ${tmpdir}/input.mnc ${tmpdir}/input.clamp.mnc
mv -f ${tmpdir}/input.clamp.mnc ${tmpdir}/input.mnc

originput=${input}
input=${tmpdir}/input.mnc

#Generate a whole-image mask
minccalc -quiet -unsigned -byte -expression 'A[0]?1:1' ${input} ${tmpdir}/initmask.mnc

################################################################################
#Round 0, N4 across entire FOV
################################################################################
n=0
mkdir -p ${tmpdir}/${n}

#Correct entire image domain
minccalc -quiet -unsigned -byte -expression 'A[0]>0.5?1:0' ${input} ${tmpdir}/${n}/weight.mnc
cp -f ${tmpdir}/${n}/weight.mnc ${tmpdir}/initweight.mnc

N4BiasFieldCorrection -d 3 -s $shrinkround1 -w ${tmpdir}/${n}/weight.mnc -x ${tmpdir}/initmask.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

#Clip the intensity range
ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/weight.mnc

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

################################################################################
#Round 1, N4 with Otsu Mask
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc

#Threshold Image and grab largest component (hopefully brain)
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/headmask.mnc Otsu 1 ${tmpdir}/$((n - 1))/weight.mnc
ThresholdImage 3 ${tmpdir}/${n}/headmask.mnc ${tmpdir}/${n}/headmask.mnc 2 2 1 0
iMath 3 ${tmpdir}/${n}/headmask.mnc ME ${tmpdir}/${n}/headmask.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/headmask.mnc GetLargestComponent ${tmpdir}/${n}/headmask.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/headmask.mnc 2 1 ball 1

#Find hotspots to exclude
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) \
  $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 0 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

#Correct entire image domain
N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/${n}/weight.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

#Clip the intensity range
ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/weight.mnc

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

#Compute coeffcient of variation
minccalc -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Round 2, N4 with Otsu Mask repeated
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc

#Threshold Image and grab largest component (hopefully brain)
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/headmask.mnc Otsu 1 ${tmpdir}/initweight.mnc
ThresholdImage 3 ${tmpdir}/${n}/headmask.mnc ${tmpdir}/${n}/headmask.mnc 2 2 1 0
iMath 3 ${tmpdir}/${n}/headmask.mnc ME ${tmpdir}/${n}/headmask.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/headmask.mnc GetLargestComponent ${tmpdir}/${n}/headmask.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/headmask.mnc 2 1 ball 1

#Find hotspots to exclude
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) \
  $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 0 1

ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

#Correct entire image domain
N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/${n}/weight.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

#Compute coeffcient of variation
minccalc -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Round 3, N4 with brain mask intersected with Otsu mask
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
minc_anlm ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

#Register to MNI space
antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc -n BSpline[5] \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform [${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1] \
  --transform Rigid[0.1] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,32,Regular,0.25] --convergence [1000x500x250,1e-6,10] --shrink-factors 8x4x4 --smoothing-sigmas ${blur16mm}x${blur8mm}x${blur4mm}mm --masks [NULL,NULL] \
  --transform Similarity[0.1] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,32,Regular,0.5] --convergence [500x250x100,1e-6,10] --shrink-factors 4x4x2 --smoothing-sigmas ${blur8mm}x${blur4mm}x${blur2mm}mm --masks [NULL,NULL] \
  --transform Affine[0.1] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,32,Regular,0.75] --convergence [250x100x50,1e-6,10] --shrink-factors 4x2x2 --smoothing-sigmas ${blur4mm}x${blur2mm}x${blur1mm}mm --masks [NULL,NULL] \
  --transform Affine[0.1] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.mnc,1,32,Regular,1] --convergence [100x50x25x10,1e-6,10] --shrink-factors 2x2x1x1 --smoothing-sigmas ${blur2mm}x${blur1mm}x${blur05mm}x0mm \
  --masks [${REGISTRATIONBRAINMASK},NULL]

#Repeat with nuyl matched registration
antsApplyTransforms -d 3 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${REGISTRATIONBRAINMASK} -o ${tmpdir}/${n}/mnimask.mnc -n GenericLabel
minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/${n}/mnimask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 99.8 --fix_zero_padding --steps 1024
antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc -n BSpline[5] \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform ${tmpdir}/${n}/mni0_GenericAffine.xfm \
  --transform Affine[0.05] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,32,Regular,1] --convergence [500x250x100x50x25,1e-6,10] --shrink-factors 4x2x2x1x1 --smoothing-sigmas ${blur4mm}x${blur2mm}x${blur1mm}x${blur05mm}x0mm \
  --masks [${REGISTRATIONBRAINMASK},${tmpdir}/${n}/mnimask.mnc]
antsApplyTransforms -d 3 -i ${tmpdir}/${n}/input.mnc -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc -r ${REGISTRATIONMODEL}

#BSpline[5] does weird things to intensity, clip back to positive range
mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/${n}/mni.mnc) ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

#Shrink the MNI mask for the first intensity matching
iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${REGISTRATIONBRAINMASK} 5 1 ball 1

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/${n}/shrinkmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

#Run a quick beast to get a brain mask
mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.4mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

#Resample beast mask and MNI mask to native space
antsApplyTransforms -d 3 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${REGISTRATIONBRAINMASK} -o ${tmpdir}/${n}/mnimask.mnc -n GenericLabel
antsApplyTransforms -d 3 -r ${tmpdir}/${n}/input.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -i ${tmpdir}/${n}/beastmask.mnc -o ${tmpdir}/${n}/bmask.mnc -n GenericLabel

#Combine the masks because sometimes beast misses badly biased cerebellum
ImageMath 3 ${tmpdir}/${n}/mask.mnc addtozero ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc

#Expand the mask a bit
iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 3 1 ball 1

#Extract brain, threshold with otsu, create weight mask
ImageMath 3 ${tmpdir}/${n}/extracted.mnc m ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/mask_D.mnc
ThresholdImage 3 ${tmpdir}/${n}/extracted.mnc ${tmpdir}/${n}/weight.mnc Otsu 1 ${tmpdir}/${n}/mask_D.mnc
ThresholdImage 3 ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/weight.mnc 2 2 1 0
iMath 3 ${tmpdir}/${n}/weight.mnc ME ${tmpdir}/${n}/weight.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/weight.mnc GetLargestComponent ${tmpdir}/${n}/weight.mnc
iMath 3 ${tmpdir}/${n}/weight.mnc MD ${tmpdir}/${n}/weight.mnc 2 1 ball 1
#ImageMath 3 ${tmpdir}/${n}/weight.mnc FillHoles ${tmpdir}/${n}/weight.mnc 2

#Create hotmask and exlcude hot voxels from weight
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) \
  $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 0 1

ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

#Do first round of masked bias field correction, use brain mask as weight
N4BiasFieldCorrection -d 3 -s $shrinkround1 -w ${tmpdir}/${n}/primary_weight.mnc -x ${tmpdir}/initmask.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

#Align multispectral inputs into T1 space
if (( ${#inputs[@]} > 0 )); then
  for file in "${inputs[@]}"; do
    #Round 1, do whole-region N4
    minccalc -quiet -expression 'A[0]?1:1' ${file} ${tmpdir}/multispectral/$(basename $file .mnc).initmask.mnc
    minccalc -quiet -expression 'A[0]>0.5?1:0' ${file} ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc

    ThresholdImage 3 ${file} ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc \
      $(mincstats -quiet -pctT 99.5 ${file}) \
      $(mincstats -quiet -max ${file}) 0 1

    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc m ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc

    N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/multispectral/$(basename $file .mnc).initmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc \
      -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${file} -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0

    ThresholdImage 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).otsu.mnc Otsu 1

    ThresholdImage 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc \
      $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/multispectral/$(basename $file .mnc).otsu.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) \
      $(mincstats -quiet -max -mask ${tmpdir}/multispectral/$(basename $file .mnc).otsu.mnc -mask_binvalue 1  ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) 0 1

    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc m ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc ${tmpdir}/multispectral/$(basename $file .mnc).otsu.mnc

    N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/multispectral/$(basename $file .mnc).initmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).otsu.mnc \
      -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${file} -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0
    #Register corrected brains to T1, do affine because distorations of different scan types
    antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc \
      --output ${tmpdir}/multispectral/$(basename $file .mnc) \
      --use-histogram-matching 0 \
      --transform Affine[0.1] --metric Mattes[${tmpdir}/${n}/corrected.mnc,${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc,1,32,Regular,1] --convergence [2000x2000x2000x2000,1e-6,10] \
      --shrink-factors 4x2x1x1 --smoothing-sigmas ${blur2mm}x${blur1mm}x${blur05mm}x0mm --masks [${tmpdir}/${n}/mask.mnc,NULL]

    antsApplyTransforms -d 3 -i $file -r ${tmpdir}/${n}/corrected.mnc -t ${tmpdir}/multispectral/$(basename $file .mnc)0_GenericAffine.xfm -o ${tmpdir}/multispectral/$(basename $file) -n BSpline[5]
    antsApplyTransforms -d 3 -i ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc -r ${tmpdir}/initmask.mnc -t ${tmpdir}/multispectral/$(basename $file .mnc)0_GenericAffine.xfm \
      -o ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc -n GenericLabel

    mincmath -clamp -const2 0 $(mincstats -max -quiet ${tmpdir}/multispectral/$(basename $file)) ${tmpdir}/multispectral/$(basename $file) ${tmpdir}/multispectral/$(basename $file).clamp.mnc
    mv ${tmpdir}/multispectral/$(basename $file).clamp.mnc ${tmpdir}/multispectral/$(basename $file)
  done

  ImageMath 3 ${tmpdir}/global_exclude.mnc m ${tmpdir}/initweight.mnc ${tmpdir}/multispectral/*.initweight.mnc
  for file in "${inputs[@]}"; do
    ThresholdImage 3 ${tmpdir}/multispectral/$(basename $file) ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc \
      $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file)) \
      $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file)) 0 1
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc m ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc ${tmpdir}/${n}/weight.mnc \
      ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc
    N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc \
      -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${tmpdir}/multispectral/$(basename $file) -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0
    multimaxval=$(mincstats -max -quiet $file)
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc TruncateImageIntensity ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc RescaleImage ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0 ${multimaxval}
    minc_anlm ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc
    mv -f ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc
  done
else
  cp -f ${tmpdir}/initweight.mnc ${tmpdir}/global_exclude.mnc
fi

minccalc -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Round 4, N4 with nonlinearly MNI-bootstrapped WM/GM segmentation proabilities
################################################################################
((++n))
mkdir -p ${tmpdir}/${n}

cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
minc_anlm ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/$((n - 1))/mask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 99.8 --fix_zero_padding --steps 1024

#Affine register to MNI space, tweak registration
antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc -n BSpline[5] \
  --output [${tmpdir}/${n}/mni] \
  --use-histogram-matching 0 \
  --initial-moving-transform ${tmpdir}/$((n - 1))/mni0_GenericAffine.xfm \
  --transform Affine[0.05] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,32,Regular,1] --convergence [500x250x100x50x25,1e-6,10] --shrink-factors 4x2x2x1x1 --smoothing-sigmas ${blur4mm}x${blur2mm}x${blur1mm}x${blur05mm}x0mm \
  --masks [${REGISTRATIONBRAINMASK},${tmpdir}/$((n - 1))/mask.mnc]

antsApplyTransforms -d 3 -i ${tmpdir}/${n}/input.mnc -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc -r ${REGISTRATIONMODEL}

mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/${n}/mni.mnc) ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

#Shrink last round's beastmask for normalization
iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${tmpdir}/$((n - 1))/beastmask.mnc 5 1 ball 1

#Intensity normalize
volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/${n}/shrinkmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

#Run a quick beast to get a brain mask
mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.2mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

#Non linearly register priors
antsRegistration --dimensionality 3 --float 1 --collapse-output-transforms 1 -a 0 --minc -u 0 -n BSpline[5] \
  --output [${tmpdir}/${n}/nonlin,${tmpdir}/${n}/mninonlin.mnc] \
  --initial-moving-transform ${tmpdir}/${n}/mni0_GenericAffine.xfm \
  --use-histogram-matching 0 \
  --transform SyN[0.1,3,0] --metric CC[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,4] \
  --convergence [500x250x100x50x25x10,1e-6,10] \
  --shrink-factors 8x4x2x2x2x1 \
  --smoothing-sigmas ${blur16mm}x${blur8mm}x${blur4mm}x${blur2mm}x${blur1mm}x0mm \
  --masks [${REGISTRATIONBRAINMASK},${tmpdir}/$((n - 1))/mask.mnc]

#Resample MNI Priors to Native space for classification
antsApplyTransforms -i ${WMPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior3.mnc -d 3 -n Linear
antsApplyTransforms -i ${GMPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior2.mnc -d 3 -n Linear
antsApplyTransforms -i ${CSFPRIOR} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/SegmentationPrior1.mnc -d 3 -n Linear
#Masks
antsApplyTransforms -i ${REGISTRATIONBRAINMASK} -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -t ${tmpdir}/${n}/nonlin1_inverse_NL.xfm -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/mnimask.mnc -d 3 -n GenericLabel
antsApplyTransforms -i ${tmpdir}/${n}/beastmask.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/bmask.mnc -d 3 -n GenericLabel

cp -f ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/mnimask.mnc

#Combine the masks because sometimes beast misses badly biased cerebellum
ImageMath 3 ${tmpdir}/${n}/mask.mnc addtozero ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc

#Expand the mask a bit
iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 2 1 ball 1

#Create hotmask
ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
  $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/mask.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) \
  $(mincstats -quiet -max -mask ${tmpdir}/${n}/mask.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) 0 1

#Exclude lesions and such
ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/global_exclude.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

#Do an initial classification using the MNI priors
Atropos -d 3 -x ${tmpdir}/${n}/mask_D.mnc -c [15,0.0] -a ${tmpdir}/${n}/input.mnc ${extra_atropos_inputs} \
  -i PriorProbabilityImages[3,${tmpdir}/${n}/SegmentationPrior%d.mnc,0.25] -k HistogramParzenWindows -m [0.1,1x1x1] \
  -o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[0] 1 --winsorize-outliers BoxPlot

#Exlcude small bits of GM and WM classification from eventual weight mask
ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class3.mnc 3 3 1 0
ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class2.mnc 2 2 1 0
ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class1.mnc 1 1 1 0
ImageMath 3 ${tmpdir}/${n}/class3.mnc GetLargestComponent ${tmpdir}/${n}/class3.mnc
ImageMath 3 ${tmpdir}/${n}/class2.mnc GetLargestComponent ${tmpdir}/${n}/class2.mnc
ImageMath 3 ${tmpdir}/${n}/class2.mnc FillHoles ${tmpdir}/${n}/class2.mnc 2

iMath 3 ${tmpdir}/${n}/class1.mnc ME ${tmpdir}/${n}/class1.mnc 10 1 ball 1

ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc addtozero ${tmpdir}/${n}/class1.mnc ${tmpdir}/${n}/class2.mnc
ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc addtozero ${tmpdir}/${n}/tissuemask.mnc ${tmpdir}/${n}/class3.mnc
iMath 3 ${tmpdir}/${n}/tissuemask.mnc ME ${tmpdir}/${n}/tissuemask.mnc 2 1 ball 1
ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc GetLargestComponent ${tmpdir}/${n}/tissuemask.mnc
iMath 3 ${tmpdir}/${n}/tissuemask.mnc MD ${tmpdir}/${n}/tissuemask.mnc 4 1 ball 1
ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc FillHoles ${tmpdir}/${n}/tissuemask.mnc 2

iMath 3 ${tmpdir}/${n}/mask_E.mnc ME ${tmpdir}/${n}/mask.mnc 5 1 ball 1
ImageMath 3 ${tmpdir}/${n}/classifymask.mnc addtozero ${tmpdir}/${n}/tissuemask.mnc ${tmpdir}/${n}/mask_E.mnc

#Combine GM and WM probably images into a N4 mask,
ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc
ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/tissuemask.mnc ${lesionmask} ${tmpdir}/global_exclude.mnc
ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc

#Perform bias field correction with weight mask
N4BiasFieldCorrection -d 3 -s $shrinkround1 -w ${tmpdir}/${n}/primary_weight.mnc -x ${tmpdir}/initmask.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

if (( ${#inputs[@]} > 0 )); then
  for file in "${inputs[@]}"; do
    ThresholdImage 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc \
      $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) \
      $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) 0 1
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc m ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc ${tmpdir}/${n}/weight.mnc \
      ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc
    N4BiasFieldCorrection -d 3 -s $shrinkround1 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc \
      -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${tmpdir}/multispectral/$(basename $file) -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0
    multimaxval=$(mincstats -max -quiet $file)
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc TruncateImageIntensity ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc RescaleImage ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0 ${multimaxval}
    minc_anlm ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc
    mv -f ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc
  done
fi

#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

minccalc -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt
rm -rf ${tmpdir}/$((n - 1))

################################################################################
#Remaining rounds, N4 with segmentations bootstrapped from prior run
################################################################################
while true; do
  ((++n))
  mkdir -p ${tmpdir}/${n}

  cp -f ${tmpdir}/$((n - 1))/corrected.mnc ${tmpdir}/${n}/input.mnc
  minc_anlm ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/denoise.mnc
  mv -f ${tmpdir}/${n}/denoise.mnc ${tmpdir}/${n}/input.mnc

  minc_nuyl ${tmpdir}/${n}/input.mnc ${REGISTRATIONMODEL} ${tmpdir}/${n}/input.nuyl.mnc --source-mask ${tmpdir}/$((n - 1))/mask.mnc --target-mask ${REGISTRATIONBRAINMASK} --cut-off 99.8 --fix_zero_padding --steps 1024

  #Register to MNI space
  antsRegistration --dimensionality 3 --float 0 --collapse-output-transforms 1 --minc -n BSpline[5] \
    --output [${tmpdir}/${n}/mni] \
    --use-histogram-matching 0 \
    --initial-moving-transform ${tmpdir}/$((n - 1))/mni0_GenericAffine.xfm \
    --transform Affine[0.05] --metric Mattes[${REGISTRATIONMODEL},${tmpdir}/${n}/input.nuyl.mnc,1,32,Regular,1] --convergence [500x250x100x50x25,1e-6,10] --shrink-factors 4x2x2x1x1 --smoothing-sigmas ${blur4mm}x${blur2mm}x${blur1mm}x${blur05mm}x0mm \
    --masks [${REGISTRATIONBRAINMASK},${tmpdir}/$((n - 1))/mask.mnc]

  antsApplyTransforms -d 3 -i ${tmpdir}/${n}/input.mnc -t ${tmpdir}/${n}/mni0_GenericAffine.xfm -n BSpline[5] -o ${tmpdir}/${n}/mni.mnc -r ${REGISTRATIONMODEL}

  mincmath -clamp -const2 0 $(mincstats -quiet -max ${tmpdir}/${n}/mni.mnc) ${tmpdir}/${n}/mni.mnc ${tmpdir}/${n}/mni.clamp.mnc
  mv -f ${tmpdir}/${n}/mni.clamp.mnc ${tmpdir}/${n}/mni.mnc

  #Shrink last round's beastmask for normalization
  iMath 3 ${tmpdir}/${n}/shrinkmask.mnc ME ${tmpdir}/$((n - 1))/beastmask.mnc 5 1 ball 1

  #Intensity normalize
  volume_pol --order 1 --min 0 --max 100 --noclamp ${tmpdir}/${n}/mni.mnc ${REGISTRATIONMODEL} --source_mask ${tmpdir}/${n}/shrinkmask.mnc --target_mask ${REGISTRATIONBRAINMASK} ${tmpdir}/${n}/mni.norm.mnc

  #Run a quick beast to get a brain mask
  mincbeast -clobber -fill -median -same_res -flip -conf ${BEASTLIBRARY_DIR}/default.2mm.conf ${BEASTLIBRARY_DIR} ${tmpdir}/${n}/mni.norm.mnc ${tmpdir}/${n}/beastmask.mnc

  #Resample beast mask and MNI mask to native space
  cp -f ${tmpdir}/mnimask.mnc ${tmpdir}/${n}/mnimask.mnc
  antsApplyTransforms -i ${tmpdir}/${n}/beastmask.mnc -t [${tmpdir}/${n}/mni0_GenericAffine.xfm,1] -r ${tmpdir}/${n}/input.mnc -o ${tmpdir}/${n}/bmask.mnc -d 3 -n GenericLabel

  #Combine the masks because sometimes beast misses badly biased cerebellum
  ImageMath 3 ${tmpdir}/${n}/mask.mnc addtozero ${tmpdir}/${n}/mnimask.mnc ${tmpdir}/${n}/bmask.mnc
  iMath 3 ${tmpdir}/${n}/mask_D.mnc MD ${tmpdir}/${n}/mask.mnc 1 1 ball 1

  #Create hotmask
  ThresholdImage 3 ${tmpdir}/${n}/input.mnc ${tmpdir}/${n}/hotmask.mnc \
    $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/mask.mnc -mask_binvalue 1 ${tmpdir}/${n}/input.mnc) \
    $(mincstats -quiet -max -mask ${tmpdir}/${n}/mask.mnc -mask_binvalue 1  ${tmpdir}/${n}/input.mnc) 0 1

  ImageMath 3 ${tmpdir}/${n}/mask_D.mnc m ${tmpdir}/${n}/mask_D.mnc ${tmpdir}/global_exclude.mnc ${tmpdir}/${n}/hotmask.mnc ${lesionmask}

  #Do an initial classification using the MNI priors, remove outliers
  Atropos -d 3 -x ${tmpdir}/${n}/mask_D.mnc -c [15,0.0] -a ${tmpdir}/${n}/input.mnc ${extra_atropos_inputs} -s 1x2 -s 2x3 \
    -i PriorProbabilityImages[3,${tmpdir}/$((n - 1))/SegmentationPosteriors%d.mnc,0] -k HistogramParzenWindows -m [0.1,1x1x1] \
    -o [${tmpdir}/${n}/classify.mnc,${tmpdir}/${n}/SegmentationPosteriors%d.mnc] -r 1 -p Socrates[1]

  #Exlcude small bits of GM and WM classification from eventual weight mask
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class3.mnc 3 3 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class2.mnc 2 2 1 0
  ThresholdImage 3 ${tmpdir}/${n}/classify.mnc ${tmpdir}/${n}/class1.mnc 1 1 1 0
  ImageMath 3 ${tmpdir}/${n}/class3.mnc GetLargestComponent ${tmpdir}/${n}/class3.mnc
  ImageMath 3 ${tmpdir}/${n}/class2.mnc GetLargestComponent ${tmpdir}/${n}/class2.mnc
  ImageMath 3 ${tmpdir}/${n}/class2.mnc FillHoles ${tmpdir}/${n}/class2.mnc 2

  iMath 3 ${tmpdir}/${n}/class1.mnc ME ${tmpdir}/${n}/class1.mnc 10 1 ball 1

  ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc addtozero ${tmpdir}/${n}/class1.mnc ${tmpdir}/${n}/class2.mnc
  ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc addtozero ${tmpdir}/${n}/tissuemask.mnc ${tmpdir}/${n}/class3.mnc
  iMath 3 ${tmpdir}/${n}/tissuemask.mnc ME ${tmpdir}/${n}/tissuemask.mnc 2 1 ball 1
  ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc GetLargestComponent ${tmpdir}/${n}/tissuemask.mnc
  iMath 3 ${tmpdir}/${n}/tissuemask.mnc MD ${tmpdir}/${n}/tissuemask.mnc 4 1 ball 1
  ImageMath 3 ${tmpdir}/${n}/tissuemask.mnc FillHoles ${tmpdir}/${n}/tissuemask.mnc 2

  iMath 3 ${tmpdir}/${n}/mask_E.mnc ME ${tmpdir}/${n}/mask.mnc 5 1 ball 1
  ImageMath 3 ${tmpdir}/${n}/classifymask.mnc addtozero ${tmpdir}/${n}/tissuemask.mnc ${tmpdir}/${n}/mask_E.mnc

  #Combine GM and WM probably images into a N4 mask,
  ImageMath 3 ${tmpdir}/${n}/weight.mnc PureTissueN4WeightMask ${tmpdir}/${n}/SegmentationPosteriors2.mnc ${tmpdir}/${n}/SegmentationPosteriors3.mnc
  ImageMath 3 ${tmpdir}/${n}/weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/tissuemask.mnc ${lesionmask} ${tmpdir}/global_exclude.mnc
  ImageMath 3 ${tmpdir}/${n}/primary_weight.mnc m ${tmpdir}/${n}/weight.mnc ${tmpdir}/${n}/hotmask.mnc

  #Perform bias field correction with weight mask
  N4BiasFieldCorrection -d 3 -s $shrinkround2 -w ${tmpdir}/${n}/weight.mnc -x ${tmpdir}/initmask.mnc \
    -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${input} -o [${tmpdir}/${n}/corrected.mnc,${tmpdir}/${n}/bias.mnc] -r 0

  if (( ${#inputs[@]} > 0 )); then
    for file in "${inputs[@]}"; do
      ThresholdImage 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc \
        $(mincstats -quiet -pctT 99.5 -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) \
        $(mincstats -quiet -max -mask ${tmpdir}/${n}/weight.mnc -mask_binvalue 1 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc) 0 1
      ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc m ${tmpdir}/multispectral/$(basename $file .mnc).hotmask.mnc ${tmpdir}/${n}/weight.mnc \
        ${tmpdir}/multispectral/$(basename $file .mnc).initweight.mnc
      N4BiasFieldCorrection -d 3 -s $shrinkround2 -x ${tmpdir}/initmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc \
        -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${tmpdir}/multispectral/$(basename $file) -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0
      multimaxval=$(mincstats -max -quiet $file)
      ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc TruncateImageIntensity ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
      ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc RescaleImage ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0 ${multimaxval}
      minc_anlm ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc
      mv -f ${tmpdir}/multispectral/$(basename $file .mnc).N4.denoise.mnc ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc
    done
  fi

  #Normalize and rescale intensity
  ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/${n}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/${n}/mask.mnc
  ImageMath 3 ${tmpdir}/${n}/corrected.mnc RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

  #Compute coeffcient of variation
  minccalc -quiet -clobber -expression 'A[0]/A[1]' ${tmpdir}/$((n - 1))/bias.mnc ${tmpdir}/${n}/bias.mnc ${tmpdir}/${n}/ratio.mnc
  python -c "print $(mincstats -quiet -stddev ${tmpdir}/${n}/ratio.mnc) / $(mincstats -quiet -mean ${tmpdir}/${n}/ratio.mnc)" >>${tmpdir}/convergence.txt

  rm -rf ${tmpdir}/$((n - 1))

  [[ (${n} -lt 10) && ($(python -c "print $(tail -1 ${tmpdir}/convergence.txt) > 0.005") == "True") ]] || break

done

echo "Convergence results:"
cat ${tmpdir}/convergence.txt

antsApplyTransforms -d 3 -i ${tmpdir}/${n}/weight.mnc -o ${tmpdir}/finalweight.mnc -r ${originput} -n Linear
antsApplyTransforms -d 3 -i ${tmpdir}/${n}/mask.mnc -o ${tmpdir}/finalmask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms -d 3 -i ${tmpdir}/${n}/bmask.mnc -o ${tmpdir}/finalbmask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms -d 3 -i ${tmpdir}/mnimask.mnc -o ${tmpdir}/finalmnimask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms -d 3 -i ${tmpdir}/${n}/classifymask.mnc -o ${tmpdir}/finalclassifymask.mnc -r ${originput} -n GenericLabel
antsApplyTransforms -d 3 -i ${tmpdir}/${n}/classify.mnc -o ${tmpdir}/finalclassify.mnc -r ${originput} -n GenericLabel

minccalc -quiet -unsigned -byte -expression 'A[0]?1:1' ${originput} ${tmpdir}/originitmask.mnc
N4BiasFieldCorrection -d 3 -s 2 -w ${tmpdir}/finalweight.mnc -x ${tmpdir}/originitmask.mnc \
  -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${originput} -o ${tmpdir}/corrected.mnc -r 0
#Normalize and rescale intensity
ImageMath 3 ${tmpdir}/${n}/corrected.mnc TruncateImageIntensity ${tmpdir}/corrected.mnc 0.0005 0.9995 1024 ${tmpdir}/finalmask.mnc
ImageMath 3 ${output} RescaleImage ${tmpdir}/${n}/corrected.mnc 0 ${maxval}

if (( ${#inputs[@]} > 0 )); then
  for file in "${inputs[@]}"; do
    antsApplyTransforms -d 3 -i $file -r ${output} -t ${tmpdir}/multispectral/$(basename $file .mnc)0_GenericAffine.xfm -o ${tmpdir}/multispectral/$(basename $file) -n BSpline[5]
    antsApplyTransforms -d 3 -i ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc -o ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc -r ${output}
    mincmath -clamp -const2 0 $(mincstats -max -quiet ${tmpdir}/multispectral/$(basename $file)) ${tmpdir}/multispectral/$(basename $file) ${tmpdir}/multispectral/$(basename $file).clamp.mnc
    mv ${tmpdir}/multispectral/$(basename $file).clamp.mnc ${tmpdir}/multispectral/$(basename $file)
    N4BiasFieldCorrection -d 3 -s 2 -x ${tmpdir}/originitmask.mnc -w ${tmpdir}/multispectral/$(basename $file .mnc).weight.mnc \
      -b [200] -c [300x300x300x300,1e-7] --histogram-sharpening [0.05,0.01,200] -i ${tmpdir}/multispectral/$(basename $file) -o ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc -r 0
    multimaxval=$(mincstats -max -quiet $file)
    ImageMath 3 ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc TruncateImageIntensity ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0.0005 0.9995 1024 ${tmpdir}/finalmask.mnc
    ImageMath 3 $(dirname $output)/$(basename $output .mnc).$(basename $file) RescaleImage ${tmpdir}/multispectral/$(basename $file .mnc).N4.mnc 0 ${multimaxval}
  done
fi

cp -f ${tmpdir}/finalbmask.mnc $(dirname $output)/$(basename $output .mnc).beastmask.mnc
cp -f ${tmpdir}/finalmnimask.mnc $(dirname $output)/$(basename $output .mnc).mnimask.mnc
cp -f ${tmpdir}/finalclassifymask.mnc $(dirname $output)/$(basename $output .mnc).classifymask.mnc
cp -f ${tmpdir}/finalclassify.mnc $(dirname $output)/$(basename $output .mnc).classify.mnc


rm -rf ${tmpdir}
