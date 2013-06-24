WaveSeq
=======

Cancer specific CNV caller for Next Generation sequence


INSTALLATION
Go to the .../WaveSeq/src/ directory and run 'perl Build.PL' to configure the
installation.  Then run './Build installdeps' to install missing perl
dependencies and './Build installexes' to install missing third party
executables.  Finally run './Build install' to install scripts into the
.../WaveSeq/bin/ directory.

use the .../WaveSeq/bin/cnv_caller.pl script to annotate existing segmentation
calls.
