#!/usr/bin/perl

#Hack to get around FindBin error where broken Carp is preloaded
BEGIN {
    if(@ARGV == 1 && $ARGV[0] eq 'findbin'){
        eval 'require FindBin';
        print $FindBin::RealBin;
        exit;
    }

    my $Bin = `$0 findbin`;
    eval "use lib '$Bin/../inc/lib'";
    eval "use lib '$Bin/../inc/perl/lib'";
    eval "use lib '$Bin/../../perl/lib'";
    eval "use lib '$Bin/../../lib'";
    eval "use lib '$Bin/../src/inc/lib'";
    eval "use lib '$Bin/../src/inc/perl/lib'";
    eval "use lib '$Bin/../src/lib'";
    eval "use lib '$Bin/../perl/lib'";
    eval "use lib '$Bin/../lib'";
}

use forks;
use forks::shared;

use strict;
use warnings;
use Carp;
use FindBin;
use Getopt::Long;
use File::Which;
use File::Spec;
use File::Basename;
use File::Temp qw(tempfile tempdir);
use POSIX qw(:sys_wait_h ceil);
use Fcntl;
use IPC::Open3;
use Storable;
use URI::Escape;
use Statistics::Distributions;
use Statistics::Regression;
use Statistics::KernelEstimation;
use MLDBM qw(DB_File Storable);
use Term::ProgressBar;
use IO::Interactive qw(is_interactive);
use File::NFSLock;
use File::Copy;
use Bio::DB::Sam;
use Tie::MmapArray;
use Algorithm::KMeans;
use constant PI => 4 * atan2(1, 1); #3.14159265358979323846...

our ($sort_exe, $bgzip_exe, $tabix_exe, $vcftools_exe);
our (%VCF);

BEGIN {
    eval 'require CNV_caller::ConfigData';
    my $edir = "$FindBin::RealBin/../exe/";
    if(!$@){
	$tabix_exe = CNV_caller::ConfigData->config('tabix');
	$vcftools_exe = CNV_caller::ConfigData->config('vcftools');
    }

    $sort_exe = File::Which::which('sort');
    ($tabix_exe) = grep {-f $_} ("$edir/*/tabix", File::Which::which('tabix')) if(!$tabix_exe);
    ($vcftools_exe) = grep {-f $_} ("$edir/*/bin/vcftools", File::Which::which('vcftools')) if(!$vcftools_exe);
    ($bgzip_exe = $tabix_exe) =~ s/tabix$/bgzip/; #comes with tabix

    #put tabix into path or vcf_tools throws a fit
    (my $tab_dir = $tabix_exe) =~ s/[^\/]+$//;
    $ENV{PATH} = "$tab_dir:$ENV{PATH}";

    #fix issue with debugger based restart
    my $last;
    while(1){
	my $stat = waitpid(-1, WNOHANG);
	last if(defined($last) && $stat eq $last);
	$last = $stat;
    }
}

use lib (dirname($vcftools_exe)."/../lib/perl5/site_perl/");
use Vcf;

my ($script) = $0 =~ /([^\/]+)$/;
my $usage = "
Usage:

     $script [options] <vcf_file>

Description:

     This script estimates the amount of mouse contamination in a sample.

Options:

     sid       <STRING>  Sample ID (for mixed VCF and BAM files)

     cpus      <INT>     CPUs to use

     sort                Indicates that VCF needs to be sorted.

";

#--set default values
my %OPT;
$OPT{cpus}  = 1;
$OPT{maf_cov_filt}  = 5;      #min coverage to keep MAF
$OPT{maf_tail_filt} = 0.12;   #threshold for trimming MAF near 0
$OPT{m_aln}         = 0.23;   #estimated fraction alignable by mouse

#--get command line options
GetOptions("cpus=i"      => \$OPT{cpus},
           "sort"        => \$OPT{sort},
           "sid=s"       => \$OPT{sid},
           "help|?"      => sub{print $usage; exit(0);});

$OPT{vcf_file} = shift;

#--print usage
if(! $OPT{vcf_file}){
    print $usage;
    exit(0);
}

#--validate and prepare command line options
my $err = validate_options(\%OPT);
die $err if($err);

#--prepare input files
print STDERR "#Preparing and indexing input files...\n";
$OPT{vcf_file}        = prepare_vcf($OPT{vcf_file});

#--estimate xenograft content
print STDERR "#Estimating mouse DNA content...\n";
($OPT{xcov}, $OPT{xfrac}) = estimate_xeno($OPT{vcf_file}, \%OPT);
print STDERR "#Mean mouse coverage: $OPT{xcov}\n";

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------
#run before forking. clears stored values that may cause core dump
sub prefork {
    undef %VCF;
}

#validates the input options
sub validate_options {
    my $OPT = shift;

    my $vcf_file  = $OPT->{vcf_file};

    my $err;
    $err .= "File $vcf_file does not exist\n" if(! -f $vcf_file);

    return $err;
}

#indexes the vcf files
sub prepare_vcf {
    my @vcf_files = @_;

    foreach my $vcf_file (@vcf_files){
	next if(!$vcf_file);

	my $bgzfile = ($vcf_file =~ /\.gz$/) ? $vcf_file : "$vcf_file.gz";
	if(! -f "$bgzfile"){
	    open(my $OUT, "> $bgzfile.tmp");
	    open(my $IN, "< $vcf_file");
	    my $spid = open3(my $SORT, my $SORTED, '>&STDERR',
			     "$sort_exe -k1,1 -k2,2n");
	    my $bpid = open3(my $BGZ, ">&".fileno($OUT), '>&STDERR', "$bgzip_exe -c");
	    
	    while(my $line = <$IN>){
		if(!$OPT{sort} || $line =~ /^#/){
		    print $BGZ $line;
		}
		else{
		    print $SORT $line;
		}
	    }
	    close($IN);
	    close($SORT);
	    while(my $line = <$SORTED>){
		print $BGZ $line;
	    }
	    close($SORTED);
	    close($BGZ);
	    close($OUT);
	    
	    waitpid($spid, 0);
	    my $sstat = $?;
	    die "ERROR: sort failed with exit status $?\n" if($sstat);
	    
	    waitpid($bpid, 0);
	    my $bstat = $?;
	    die "ERROR: bgzip failed with exit status $bstat\n" if($bstat);
	    
	    #move through hard links
	    link("$bgzfile.tmp", $bgzfile);
	    unlink("$bgzfile.tmp");
	}
	if(! -f "$bgzfile.tbi"){
	    my $tstat = system("$tabix_exe -p vcf $bgzfile");
	    die "ERROR: Tabix failed with exit status $tstat\n" if($tstat);
	}

	$vcf_file = $bgzfile; #replace with new location
    }

    return (wantarray) ? @vcf_files : shift @vcf_files;
}

#estimates mouse xenograft content
sub estimate_xeno {
    my $vcf_file = shift;
    my $OPT      = shift;

    return $OPT->{xcov} if($OPT->{xcov});

    my $sid  = $OPT->{sid};    
    my $cpus = $OPT->{cpus};

    #load mouse markers
    my %marker;
    my $mousef = "$FindBin::RealBin/../data/mouse_markers.txt";
    open(IN, "<$mousef");
    while(my $line = <IN>){
	chomp $line;
	next if(! $line);
	my ($chr, $pos) = split(/\t/, $line);
	$marker{$chr}{$pos}++;
    }
    close(IN);

    my $flag :shared;
    my @work :shared;
    my @threads;
    #--launch threads
    $flag = 1; #thread flag
    prefork(); #prepare for forking
    for(my $i = 1; $i < $cpus; $i++){
	push(@threads, threads->create({context => 'scalar'},
				       \&_estimate_xeno_thread,
				       $vcf_file,
				       $OPT,
				       \@work,
				       \%marker,
				       \$flag));
    }

    #--distribute work to threads
    foreach my $chr (keys %CHRS){
	push(@work, $chr);
    }
    $flag = 0; #allows finished threads to exit
    
    #--jump into the mix
    my $list = _estimate_xeno_thread($vcf_file, $OPT, \@work, \%marker, \$flag);
    foreach my $thr (@threads){
	my $l = $thr->join();
	while(my $key = each %$l){
	    $list->{$key} += $l->{$key};
	}
    }

    #build kernel estimators
    my $ker = Statistics::KernelEstimation->new();
    foreach my $t (keys %$list){
	next unless($t);
	$ker->add_data($t, $list->{$t});
    }

    my $w = ($ker->default_bandwidth() < 2) ? 1 : $ker->default_bandwidth();
    my ($low, $high) = $ker->extended_range();
    $low = 0 if($low < 0);
    my $peak = [0, 0];
    for(my $x = $low; $x <= $high; $x += ($high-$low)/1000) {
	my $y = $ker->pdf($x, $w);
	$peak = [$x, $y] if($y > $peak->[1]);
	last if($x == $high);
    }

    #refine estimate
    $high = $peak->[0] + ($high-$low)/100;
    $low  = $peak->[0] - ($high-$low)/100;
    for(my $x = $low; $x <= $high; $x += 0.001) {
        my $y = $ker->pdf($x, $w);
        $peak = [$x, $y] if($y > $peak->[1]);
        last if($x == $high);
    }

    return sprintf('%.3f', $peak->[0]);
}

sub _estimate_xeno_thread {
    my $vcf_file = shift;
    my $OPT = shift;
    my $work = shift;
    my $marker = shift;
    my $flag = shift;

    my $sid = $OPT->{sid};

    #load VCF files to get expected densities
    my $vcf = Vcf->new(file=>"$vcf_file");
    $vcf->parse_header();
    my ($SID) = grep {/$sid$/} $vcf->get_samples();
    $vcf->set_samples(include=>[$SID]);
    
    my %list;
    while((my $chr = shift @$work) || $$flag){
	next if(!$chr);

	$vcf->open(region=> $chr);
	
	while(my $v = $vcf->next_data_hash()){
	    my $pos = $v->{POS};
	    next unless($marker->{$chr}{$pos});
	    next unless($v->{gtypes}{$SID}{AD});
	    
	    my ($rc, $ac) = split(/,/, $v->{gtypes}{$SID}{AD});
	    next if(!$ac);
	    $list{$ac}++;
	}
    }

    return \%list;
}
