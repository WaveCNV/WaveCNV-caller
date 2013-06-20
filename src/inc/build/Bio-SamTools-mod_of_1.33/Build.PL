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
use Module::Build;
use INSTALL::Build;
use File::Copy;

my $base = INSTALL::Build->resume (
    properties => {
	config_dir => '_build',
    },
   );

chdir($FindBin::Bin);

#=====from Bio::DB::Sam Build.PL
my $HeaderFile = "bam.h";
my $LibFile    = "libbam.a";
my $ReadLine;

my ($sam_include,$sam_lib) = find_sam(); # may exit with error here
my $build = Module::Build->new(
    module_name        => 'Bio::SamTools',
    dist_version_from  => 'lib/Bio/DB/Sam.pm',
    dist_author        => 'Lincoln Stein <lincoln.stein@gmail.com>',
    dist_abstract      => 'Perl interface to SamTools library for DNA sequencing',
    license            => 'perl',
    include_dirs       => [$sam_include],
    extra_linker_flags => ["-L$sam_lib",'-lbam','-lz',],

    extra_compiler_flags=>[

        # must match DFLAGS in Samtools Makefile
        '-D_IOLIB=2','-D_FILE_OFFSET_BITS=64',

        # turn off warnings originating in Perl's Newx* calls
        '-Wformat=0',
    ],

    build_requires => {
      'ExtUtils::CBuilder' => 0,
    },
    requires     => {
        'perl'                => '5.008',
        'Bio::Root::Version'  => '1.006009001',
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

sub find_sam {
    my ($sam_include,$sam_lib);

    $sam_include = $base->config_data('sam_include');
    $sam_lib = $base->config_data('sam_lib');

    if (my $samtools = $ENV{SAMTOOLS}) {
	if(!$sam_include){
	    $sam_include = $samtools if(-f "$samtools/$HeaderFile");
	    $sam_include = "$samtools/include" if(-f "$samtools/include/$HeaderFile");
	}
	if(!$sam_lib){
	    $sam_lib     = $samtools if(-f "$samtools/$LibFile");
	    $sam_lib     = "$samtools/lib" if(-f "$samtools/lib/$LibFile");
	}
    }

    my @search_path = qw(/ /usr /usr/share /usr/local);

    unless ($sam_include) {
        for my $p (@search_path) {
            $sam_include ||= "$p/include" if
		-e "$p/include/$HeaderFile";
        }
    }

    unless ($sam_lib) {
        for my $p (@search_path) {
            $sam_lib ||= "$p/lib" if
                -e "$p/lib/$LibFile";
        }
    }


    unless ($sam_include && $sam_lib) {
        print STDOUT "This module requires samtools 0.1.10 or higher (samtools.sourceforge.net).\n";
        my $prompt = "Please enter the location of the bam.h and compiled libbam.a files: ";
        my $found;
        while (!$found) {
            my $path = prompt($prompt);
            print STDOUT "\n";
            last unless $path;
            $sam_include = $path
                if -e "$path/$HeaderFile";
            $sam_include = "$path/include"
                if -e "$path/include/$HeaderFile";
            $sam_lib = $path
                if -e "$path/$LibFile";
            $sam_lib = "$path/lib"
                if -e "$path/lib/$LibFile";
            $found = $sam_include && $sam_lib;
            unless ($found) {
                print STDOUT "That didn't seem to be right.\n";
                $prompt = "Try again, or hit <enter> to cancel: ";
            }
        }
    }
    
    unless ($sam_include && $sam_lib) {
        die "Can\'t find $LibFile and/or $HeaderFile!\n\n".
	    "If you haven\'t done so already, please compile samtools 0.1.10 or\n".
	    "higher from http://samtools.sourceforge.net and set the SAMTOOLS\n".
	    "environment variable to point to a samtools distribution directory\n".
	    "containing the compiled $LibFile and $HeaderFile files.\n";
    }
    
    print STDOUT "Found $sam_include/$HeaderFile and $sam_lib/$LibFile.\n";
    return ($sam_include,$sam_lib);
}


sub prompt {
    my $msg = shift;

    unless (defined $ReadLine) {
        eval "require Term::ReadLine";
        $ReadLine   = Term::ReadLine->can('new') || 0;
        $ReadLine &&= Term::ReadLine->new(\*STDOUT);
        eval {readline::rl_set('TcshCompleteMode','On')};
    }


    unless ($ReadLine) {
        print STDOUT $msg;
        chomp (my $in = <>);
        return $in;
    }

    my $in = $ReadLine->readline($msg);
    chomp $in;
    $in=~ s/\s+//;
    $ReadLine->addhistory($in) if $in =~ /\S/;
    return $in;
}
