iterativeN4_multispectral.sh
------------------------------

A preprocessing pipeline which performs iterative N4 bias field correction
and tissue classification until stability of the bias field is achieved.

Utilizes the MNI ICBM priors (brainmask, CSF, GM, WM probabilities) and the
BEaST patch based segmentation tool.

```sh
> iterativeN4_multispectral.sh -h 
iterativeN4_multispectral.sh is script which performs iterative inhomogeneity (bias field) correction and classification on T1w (and optionally T2w/PDw) MRI scans
Usage: ./iterativeN4_multispectral.sh [-e|--exclude <arg>] [--t2 <arg>] [--pd <arg>] [-c|--config <arg>] [-l|--logfile <arg>] [-s|--(no-)standalone] [--max-iterations <arg>] [--convergence-threshold <arg>] [--classification-prior-weight <arg>] [-d|--(no-)debug] [-v|--verbose] [-h|--help] <input> <output>
  <input>: T1w scan to be corrected
  <output>: Output filename for corrected T1w (also used as basename for other outputs)
  -e, --exclude: Mask file defining regions to exclude from classifcation, region is still corrected (no default)
  --t2: T2w scan, will be rigidly registered to T1w and used for classification (no default)
  --pd: PDw scan, will be rigidly registered to T1w and used for classification (no default)
  -c, --config: Path to an alternative config file defining priors to use (no default)
  -l, --logfile: Path to file to log all output (no default)
  -s, --standalone, --no-standalone: Script is run standalone so save all outputs (off by default)
  --max-iterations: Maximum number of iterations to run (default: '10')
  --convergence-threshold: Coeffcient of variation limit between two bias field estimates (default: '0.005')
  --classification-prior-weight: How much weight is given to prior classification proabilities during iteration (default: '0.25')
  -d, --debug, --no-debug: Debug mode, increase verbosity further, don't cleanup (off by default)
  -v, --verbose: Set verbose output (can be specified multiple times to increase the effect)
  -h, --help: Prints help
```

**Note: multispectral support is currently not enabled**

## Configuration

Currently the functionality of `iterativeN4_multispectral.sh` is mostly automatic.

Most the major adjustments of functionality are achieved by providing different priors via the ``--config`` option.

# Getting Priors

This pipeline uses the priors available from the MNI at http://nist.mni.mcgill.ca/?page_id=714
