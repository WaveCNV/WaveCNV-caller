#! /usr/bin/perl
use FindBin;
use lib "$FindBin::Bin/../../lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../../perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../../../lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../../../../perl/lib"; #libary just needed for configuration
use lib "$FindBin::Bin/../../../../lib"; #libary just needed for configuration
use warnings;
use strict;
use Cwd;
use lib "$FindBin::Bin/../../lib/"; #this PL file always run by Build, so cwd should be constant
use Module::Build;
use INSTALL::Build;
use File::Copy;

my $base = INSTALL::Build->resume (
    properties => {
	config_dir => '_build',
    },
   );

chdir($FindBin::Bin);

#=====for Tie::MmapArray Build.PL
my $build = Module::Build->new(
    module_name        => 'Tie::MmapArray',
    dist_version_from  => 'MmapArray.pm',
    dist_author        => 'Andrew Ford <A.Ford@ford-mason.co.uk>',
    dist_abstract      => 'Lets you use mmap to map in a file as a perl array rather than '.
                          ' reading the file into dynamically allocated memory',
    license            => 'artistic',
    extra_linker_flags => [
    ],
    extra_compiler_flags=>[
    ],
    build_requires => {
      'ExtUtils::CBuilder' => 0,
    },
    requires     => {
        'perl'                => '5.10.0',
    },
    pm_files => {
	'MmapArray.pm' => 'lib/Tie/MmapArray.pm'
    },
    xs_files => {
	'MmapArray.xs' => 'lib/Tie/MmapArray.xs'
    },
    install_base => $base->base_dir,
    install_base_relpaths => {
	 bin    => [qw(lib/bin)],
	 script => [qw(lib/bin)],
         arch   => [qw(lib)],
         lib    => [qw(lib)],
         bindoc => [qw(lib/man)],
         libdoc => [qw(lib/man)],
    },
    create_packlist => 0,
);
$build->create_build_script();
$build->dispatch('install');

exit 0;
