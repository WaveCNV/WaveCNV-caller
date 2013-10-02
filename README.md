WaveCNV-caller
=======

Cancer specific CNV caller for Next Generation sequence


INSTALLATION
Go to the .../WaveCNV/src/ directory and run 'perl Build.PL' to configure the
installation.  Then run './Build installdeps' to install missing perl
dependencies and './Build installexes' to install missing third party
executables.  Finally run './Build install' to install scripts into the
.../WaveCNV/bin/ directory.

To run with the UCSC verion of the human reference genome (HG19) you must run
'./Build hg19' to download it to the data directory.

use the .../WaveCNV/bin/cnv_caller.pl script to annotate existing segmentation
calls.

There are example files to run with in the .../WaveCNV/data/ directory.
