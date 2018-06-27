iterativeN4_multispectral.sh
------------------------------

A preprocessing pipeline which performs iterative N4 bias field correction
and tissue classification until stability of the bias field is achieved.

Utilizes the MNI ICBM priors (brainmask, CSF, GM, WM probabilities) and the
BEaST patch based segmentation tool.

```sh
> iterativeN4_multispectral.sh <output.mnc> <T1.mnc> [<other input spectra>]
```

Note: multispectral support is currently in flux

## Configuration

Currently the functionality of `iterativeN4_multispectral.sh` is mostly automatic.

If the environment variable `N4_VERBOSE` is set, all bash commands will be printed
and any programs which can have progress reported will executed as such.

If the environment variable `N4_STANDALONE` is set, `iterativeN4_multispectral.sh`
will save the final outputs from the internal posteriors in addition to the corrected
file (masks and classification)
