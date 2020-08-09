iterativeN4_multispectral.sh
------------------------------

A preprocessing pipeline which performs iterative N4 bias field correction
and tissue classification until stability of the bias field is achieved.

Utilizes the MNI ICBM priors (brainmask, CSF, GM, WM probabilities) and the
BEaST patch based segmentation tool.

Insprired by the https://github.com/ANTsX/ANTs ``antsAtroposN4.sh`` tool.

```
> iterativeN4_multispectral.sh -h 
iterativeN4_multispectral.sh is script which performs iterative inhomogeneity (bias field) correction and classification on T1w (and optionally T2w/PDw) MRI scans
Usage: iterativeN4_multispectral.sh [-e|--exclude <arg>] [-c|--config <arg>] [-l|--logfile <arg>] [-s|--(no-)standalone] [-a|--(no-)autocrop] [--max-iterations <arg>] [--convergence-threshold <arg>] [--classification-prior-weight <arg>] [--(no-)debug] [-v|--verbose] [-h|--help] <input> <output>
        <input>: T1w scan to be corrected
        <output>: Output filename for corrected T1w (also used as basename for other outputs)
        -e, --exclude: Mask file defining regions to exclude from classifcation, region is still corrected (no default)
        -c, --config: Path to an alternative config file defining priors to use, use "auto" to use automatic template selection (no default)
        -l, --logfile: Path to file to log all output (no default)
        -s, --standalone, --no-standalone: Script is run standalone so save all outputs (off by default)
        -a, --autocrop, --no-autocrop: Crop the final output to 10 mm around the head determined by headmask from modelspace (off by default)
        --max-iterations: Maximum number of iterations to run (default: '10')
        --convergence-threshold: Coeffcient of variation limit between two bias field estimates (default: '0.01')
        --classification-prior-weight: How much weight is given to prior classification proabilities during iteration (default: '0.25')
        --debug, --no-debug: Debug mode, increase verbosity further, don't cleanup (off by default)
        -v, --verbose: Set verbose output (can be specified multiple times to increase the effect)
        -h, --help: Prints help
```

**Note: multispectral support is currently not enabled**

## Dependencies

- minc-toolkit-v2 https://bic-mni.github.io/
- BeAST library properly configured https://bic-mni.github.io/#data-files-and-models
- ANTs with MINC support enabled https://github.com/ANTsX/ANTs
- Priors (see below)
- imagemagick for a static QC image
- The webp package from google to get animated QC images: https://developers.google.com/speed/webp/download

## Usage

Currently the functionality of `iterativeN4_multispectral.sh` is mostly automatic.
The most basic operation is:
``iterativeN4_multispectral.sh input.mnc output.mnc``
This will produce a bias-field corrected image.

Since a number of useful output files are generated internally, adding the ``--standalone`` will produce masks, a classification,
and QC images. The ``--autocrop`` option will use an estimated headmask to recrop the output image to remove background
and excess neck data.

If your inputs are from a known population demographic, you can provide customized priors via the ``--config`` option.
Finally, you can use ``--config auto`` for the pipeline to attempt to choose the most similar population template.

# Getting Priors

This pipeline uses the priors available from the MNI at http://nist.mni.mcgill.ca/?page_id=714. The "ants" style priors
are modified versions of https://figshare.com/articles/ANTs_ANTsR_Brain_Templates/915436 to work with this pipeline
and available upon request.
