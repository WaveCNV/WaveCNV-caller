#!/usr/bin/perl
eval 'exec /usr/bin/env perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;
use FindBin;
use lib "$FindBin::Bin/inc/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/inc/perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../lib"; #libary just needed for configuration
use INSTALL::Build;

chdir($FindBin::Bin); #alway go here for safety

#configure Build
my $build = INSTALL::Build->new(
     module_name => 'CNV_caller',
     dist_author => 'Carson Holt <carson.holt@oicr.on.ca>',
     dist_version_from  => "$FindBin::Bin/version",
     dist_abstract => ('cnv_caller assigns digitized copy numbvers and alleles to CNV segments '.
		       'in cancer genomes using sequencing data and adjusts for xenograft, '.
		       'diploid normal contamination, aneuploidy, and germline events.  Final '.
		       'output is a GVF file with somatic and loss-of-heterozygosity events '.
		       'highlighted.'
		       ),
     dynamic_config => 1,
     license  => 'artistic_2',
     auto_configure_requires => 1,
     requires => {
	 'perl'               => '5.10.0',
	 'Bio::Root::Version' => '1.006901',
	 'Perl::Unsafe::Signals' => '0',
	 'Inline::C'          => '0.48',
	 'forks'              => '0.34',
	 'forks::shared'      => '0.34',
	 'Statistics::Distributions' => 0,
	 'Statistics::KernelEstimation' => 0,
	 'Statistics::Distributions' => 0,
	 'Term::ProgressBar' => 0,
	 'IO::Interactive' => 0,
	 'File::NFSLock' => '0',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     exe_requires => {
	 'samtools' => 'samtools',
	 'vcftools' => 'vcftools',
	 'tabix'    => 'tabix',
     },
     lib_requires => {
	 'samtools bam.h' => 'bam.h',
     },
     build_requires => {
	 'Test::More'         => '0',
	 'ExtUtils::CBuilder' => '0',
	 'File::Which'        => '0',
	 'LWP::Simple'        => '0',
     },
     config_requires => {
	 'File::Which'        => '0',
	 'Term::ReadKey'      => '0',
     },
     recommends => {
	 'Carp'               => '1.24',
	 'GD'                 => '0',
	 #'Other::Module' => '>= 1.2, != 1.5, < 2.0',
     },
     install_base => "$FindBin::Bin/../",
     install_base_relpaths => {
	 arch   => [qw(perl/lib)],
	 lib    => [qw(perl/lib)],
	 bindoc => [qw(perl/man)],
	 libdoc => [qw(perl/man)],
	 data   => [qw(data)],
     },
     PL_files => {'inc/build/Bio-SamTools-mod_of_1.33/Build.PL' => [],
		  'inc/build/Tie-MmapArray-mod_of_0.05/Build.PL' => [],
     },
     create_packlist => 0,
);

#ask for MPI
#my $go = $build->y_n("\ncnv_caller supports distributed parallelization via MPI.\n".
#		     "Would you like to configure MAKER for MPI (This\n".
#		     "requires that you have an MPI client installed)?", 'N');
#$build->feature(mpi_support => $go);
#$build->config_mpi() if($go);

$build->config_exe_loc('samtools');
$build->config_sam_libs();
$build->config_exe_loc('tabix');
$build->config_exe_loc('vcftools');
$build->config_data(PATH => $ENV{PATH});
$build->config_data(build_done => 0);
$build->create_build_script;
$build->log_info("\n\nThe file 'Build' has been created for you to finish installation.\n");
$build->status;
