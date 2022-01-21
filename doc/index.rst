
OXASL_BIDS - Tools for interoperating between BIDS datasets and oxasl/oxford_asl
================================================================================

OXASL_BIDS is a package to enable oxford_asl and oxasl tools to be run more easily
on data in the BIDS format.

The main functions provided are:

 - Generating scripts to run oxford_asl/oxasl with appropriate options on each
   relevant subject/session in a BIDS dataset, by extracting information from the
   BIDS metadata
   
 - Converting the output of an oxford_asl/oxasl run into a BIDS-compatible format
 
 - Converting BIDS-style JSON metadata from an image file into oxford_asl/oxasl
   options independently of its presence in a compliant BIDS format data set (e.g.
   an image/json pair output from DCM2NIIX
   
Current limitations
-------------------

 - Label/control or differenced data only, no support for multiphase, vessel encoded
   or other styles of ASL data
 
 - Output may not be fully BIDS compliant pending better understanding of the standard

 - Not all oxford_asl outputs are included in BIDS transformation

Usage
-----

Generation of an oxford_asl script
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command will generate a script containing a single ```oxford_asl``` command for
each ASL session found in the BIDS dataset::

    oxasl_bids script --bidsdir=<BIDS_DATASET>

Options
^^^^^^^

 - ``--run-fslanat`` : If specified, an additional call to FSL_ANAT will be included for each structural
   image found. The ``oxford_asl`` commands will be directed to the structural data in the output directory
   of this command.
 - ``--output-as-bids`` : If specified, additional commands will be included to convert the oxford_asl output
   back in to BIDS format. This may be merged with the existing BIDS data set as derived data or a new
   BIDS data set can be generated (see ``--merge-output``, ``--bids-output`` in section below)

Converting oxford_asl output to BIDS format
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

BIDS supports two different means of storing derived data:

 1. The derived data is the 'primary' content of the BIDS data set. The source data that it was derived
    from may also be included under the ``sourcedata`` subfolder, but this is not required::

    oxasl_bids bidsout --bidsdir=<SOURCE_DATASET> --bids-output=<DEST_DATASET> --subject=<SUBJECT> --session=<SESSION>

 2. The derived data is merged into the original BIDS data set (or a copy of it). The source data set
    is unchanged apart from the addition of a ``derivatives`` folder in which the derived data is stored::

    oxasl_bids bidsout --bidsdir=<SOURCE_DATASET> --merge-output --subject=<SUBJECT> --session=<SESSION>

Genertion of oxford_asl options from JSON metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This command outputs ``oxford_asl`` options extracted from the metadata of the supplied image
which will be interpreted as the type given::

    oxford_asl args --img=<NIFTI_IMAGE> --img-type=[asl|calib|cblip|struct]

 - ``asl`` - ASL data consisting of label/control pairs or already differenced
 - ``calib`` - M0 image intended for calibration of ASL data
 - ``cblip`` - Phase encoding reversed M0 image intended for distortion correction
 - ``struct`` - Structural image intended for segmentation and transformation of output

The purpose of this is to work with data that is not in BIDS format but does contain BIDS-compatible
metadata sidecar files. An example of this would be data converted using the ``dcm2niix`` tool with the
option ``-b y`` (this is the default in recent versions of ``dcm2niix``). For example freshly converted
ASL session data could be run through ``oxford_asl`` as follows without any manual inspection of the
metadata::

   oxford_asl `oxasl_bids args --img=asldata.nii.gz --img-type=asl` \
            `oxasl_bids args --img=calib.nii.gz --img-type=calib` \
            `oxasl_bids args --img=cblip.nii.gz --img-type=cblip` \
            `oxasl_bids args --img=t1.nii.gz --img-type=struct` \
            -o oxford_asl_out --mc

The only options that need to be specified in this case are those connected with the modelling itself, 
e.g. whether to do motion correction.




