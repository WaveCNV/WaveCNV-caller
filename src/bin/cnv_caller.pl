#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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

our ($sort_exe, $bgzip_exe, $tabix_exe, $vcftools_exe, $hg19, $dummy);
our ($DUM, %VCF, %BAMS, %GEN, %CHRS);

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
    $hg19 = "$FindBin::RealBin/../data/hg19_random.fa";
    $dummy = "$FindBin::RealBin/../data/dummy.bam";

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

     $script [options] <segment_file> <vcf_file>

Description:

     This script produces copy number estimates for specified genomic segments
     using processed sequencing data in VCF and BAM files. It does not produce
     it's own segments. Only a tab-delimited segment and VCF file are required.

     BAM files are optional, but required for calling shorter segments, as the
     SNP density of the human reference genome is only 1 per ~1000bp.

     Supplying a reference sample ID or reference files for tumor samples
     allows for corrections to be made based on variation in the reference
     coverage. Also reference sample data allows CN events to be labeled as
     somatic, germline, or true loss-of-heterozygosity. Use the 'merge' and
     'smooth' options to take full advantage of these corrections.

Options:

     fasta     <PATH>    Reference fasta to use (default is UCSC hg19)

     sid       <STRING>  Sample ID (for mixed VCF and BAM files)

     bam_dir   <PATH>    Directory of BAM files to use

     bam_list  <LIST>    List of BAM files to use (comma separated)

     rid       <STRING>  Reference sample ID

     rvcf_file <PATH>    VCF file to use for reference (if different)

     rbam_dir  <PATH>    Directory of reference BAM files to use (if different)

     rbam_list <LIST>    List of reference BAM files to use (if different)

     xvcf_file <PATH>    VCF file to use for mouse correction

     xbam_dir  <PATH>    Directory of mouse BAM files to use for correction

     xbam_list <LIST>    List of mouse BAM files to use for correction

     lfrag     <INT>     Sequencing fragment length (2x for paired end)

     cell                Estimate copy fraction (akin to cellularity) of sample

     cfrac     <NUM>     Set copy fraction explicitly

     xeno                Estimate xenograft sample contamination

     xcov      <NUM>     Sets xenograft coverage explicitly

     merge               merge neigboring regions of same copy number (recommended)

     smooth              Smooth over small germline events (implies merge)

     seqid     <INT>     Segment seqid column in segment_file (default=1)

     begin     <INT>     Segment start column in segment_file (default=2)

     end       <INT>     Segment end column in segment_file (default=3)

     cpus      <INT>     CPUs to use

     sort                Indicates that VCF needs to be sorted.

     tmp       <PATH>    Alternate temp file directory (default set by system)

";

#--set default values
my %OPT;
$OPT{cfrac}  = 1;
$OPT{rfrac}  = 0;
$OPT{xcov}   = 0;
$OPT{k_bins} = 24;  #bins to use for model distribition or "k"
$OPT{ref}  = {%OPT}; #ref is just a copy of the current sample parameters
$OPT{xeno} = {%OPT}; #xeno is just a copy of the current sample parameters
$OPT{C}     = 1;
$OPT{B}     = 2;
$OPT{E}     = 3;
$OPT{cpus}  = 1;
$OPT{TMP}   = File::Spec->tmpdir();
$OPT{fasta} = $hg19;
$OPT{maf_cov_filt}  = 5;      #min coverage to keep MAF
$OPT{maf_tail_filt} = 0.12;   #threshold for trimming MAF near 0
$OPT{m_aln}         = 0.23;   #estimated fraction alignable by mouse
$OPT{e_rate}        = 0.0005;  #NGS error rate (per bp per sequence)
$OPT{basement}      = 0.05;   #basement expect for model distibutions
$OPT{ref}{cn_max}   = 5;
$OPT{cn_max}        = 5;
$OPT{ref}{models}   = get_models($OPT{cn_max}+1); #i.e. 0:1, 1:2, etc. up to specified CN
$OPT{models}        = get_models($OPT{cn_max}+1); #i.e. 0:1, 1:2, etc. up to specified CN
$OPT{lfrag}         = 1000; #use average SNP density as estimate

#--get command line options
GetOptions("cpus=i"      => \$OPT{cpus},
           "s|seqid=i"   => \$OPT{C},
           "b|begin=i"   => \$OPT{B},
           "e|end=i"     => \$OPT{E},
           "sort"        => \$OPT{sort},
           "sid=s"       => \$OPT{sid},
           "rid=s"       => \$OPT{ref}{sid},
           "xid=s"       => \$OPT{xeno}{sid},
           "xeno"        => \$OPT{use_xeno},
           "xcov=f"      => \$OPT{xcov},
           "cell"        => \$OPT{cell},
           "cfrac=f"     => \$OPT{cfrac},
           "lfrag"       => \$OPT{lfrag},
           "smooth"      => \$OPT{smooth},
           "merge"       => \$OPT{merge},
           "exome"       => \$OPT{exome},
           "tmp=s"       => \$OPT{TMP},
           "fasta=s"     => \$OPT{fasta},
           "fix"         => \$OPT{fix},
           "bam_dir=s"   => \$OPT{bam_dir},
           "bam_list=s"  => \$OPT{bam_list},
           "rvcf_file=s" => \$OPT{ref}{vcf_file},
           "rbam_dir=s"  => \$OPT{ref}{bam_dir},
           "rbam_list=s" => \$OPT{ref}{bam_list},
           "xvcf_file=s" => \$OPT{xeno}{vcf_file},
           "xbam_dir=s"  => \$OPT{xeno}{bam_dir},
           "xbam_list=s" => \$OPT{xeno}{bam_list},
           "t_cov=f"     => \$OPT{t_cov},
           "rt_cov=f"    => \$OPT{ref}{t_cov},
           "base_cov=f"  => \$OPT{base_cov},
           "rbase_cov=f" => \$OPT{ref}{base_cov},
           "ploidy=f"    => \$OPT{ploidy},
           "rploidy=f"   => \$OPT{ref}{ploidy},
           "help|?"      => sub{print $usage; exit(0);});

$OPT{seg_file} = shift;
$OPT{vcf_file} = shift;

#--print usage
if(! $OPT{seg_file} || ! $OPT{vcf_file}){
    print $usage;
    exit(0);
}

#--validate and prepare command line options
print STDERR "#Validating options and preparing input files...\n";
validate_options(\%OPT);

#--process stored settings
$OPT{DB_File} = $OPT{outdir}.'/dbfile';
tie(my %DB, 'MLDBM', $OPT{DB_File}, O_CREAT|O_RDWR, 0640) or die $!;
$DB{ARC}       ||= {ref => {}, xeno => {}}; #archive discoverable values
$DB{disc}      ||= {};
$DB{load}      ||= {};
$DB{merge}     ||= {};
$DB{smooth}    ||= {};

#check to see if the options are changed
my $DBARC = $DB{ARC};
if(my $DBOPT = $DB{OPT}){
    my $disc   = $DB{disc};
    my $load   = $DB{load};
    my $merge  = $DB{merge};
    my $smooth = $DB{smooth};
    
    #rerun of discovery cn and load cn
    my @check1 = qw(base_cov
		    ploidy
		    t_cov);

    #rerun of discovery maf and load maf
    my @check2 = qw(cfrac);
    
    #these force complete rerun
    my @check3 = qw(use_xeno
		    use_ref
		    use_generic
		    vcf_file
		    bam_files
		    fasta
		    sid
		    cell
		    xcov);

    my %effect;
    @effect{@check1} = map {1} @check1;
    @effect{@check2} = map {2} @check2;
    @effect{@check3} = map {3} @check3;
    
    my @comp = ([\%OPT, $DBOPT, $DBARC],
		[ref_param(\%OPT), ref_param($DBOPT), $DBARC->{ref}],
		[xeno_param(\%OPT), xeno_param($DBOPT), $DBARC->{xeno}]);
    
    my $clober = 0;
    foreach my $o (@comp){
	my ($new, $old, $arc) = @$o;
	foreach my $key (@check1, @check2, @check3){
	    next if(!$new->{$key} && !$old->{$key});
	    my $v1 = $new->{$key} || '';
	    my $v2 = $old->{$key} || '';
	    
	    if($key eq 'bam_files'){
		$v1 = join(',', sort map {Cwd::abs_path($_)} @$v1);
		$v2 = join(',', sort map {Cwd::abs_path($_)} @$v2);
	    }
	    next if($v1 eq $v2);
	    if(defined($arc->{$key})){
		if($arc->{$key} eq $v1){
		    next;
		}
		else{
		    $arc->{$key} = $v1;
		}
	    }
	    
	    $clober = $effect{$key} unless($clober > $effect{$key});
	}
    }
    
    if($clober == 1){
	$merge  = {};
	$smooth = {};
	$load->{_hascn} = 0;
	File::Path::rmtree($OPT{outdir}."/merge");
	File::Path::rmtree($OPT{outdir}."/smooth");
    }
    elsif($clober == 2){
	$merge  = {};
	$smooth = {};
	$disc->{_hasmaf} = 0;
	$load->{_premerge} = 0;
	File::Path::rmtree($OPT{outdir}."/merge");
	File::Path::rmtree($OPT{outdir}."/smooth");
    }
    elsif($clober == 3){ #remove everything
	$DBARC  = {};
	$disc   = {};
	$load   = {};
	$merge  = {};
	$smooth = {};
	File::Path::rmtree($OPT{outdir}."/disc");
	File::Path::rmtree($OPT{outdir}."/load");
	File::Path::rmtree($OPT{outdir}."/merge");
	File::Path::rmtree($OPT{outdir}."/smooth");
    }

    $DB{ARC}    = $DBARC;
    $DB{disc}   = $disc;
    $DB{load}   = $load;
    $DB{merge}  = $merge;
    $DB{smooth} = $smooth;
}
$DB{OPT} = \%OPT; #store current options
untie(%DB);

#reuse archived data where possible
$OPT{$_} = $DBARC->{$_} foreach(grep {!/^(ref|xeno)$/} keys %$DBARC);
$OPT{ref}{$_} = $DBARC->{ref}{$_} foreach(keys %{$DBARC->{ref}});
$OPT{xeno}{$_} = $DBARC->{xeno}{$_} foreach(keys %{$DBARC->{xeno}});

#--now read in segments from file
print STDERR "#Reading segment file...\n";
my $segs = read_segment_file($OPT{seg_file}, \%OPT);

#--build discovery segments
print STDERR "#Loading data for initial parameter discovery...\n";
my $disc = discovery_segments(\%OPT);

#--estimate xenograft content
if($OPT{use_xeno}){
    print STDERR "#Estimating mouse DNA content...\n";
    $OPT{xcov} = estimate_xeno($OPT{vcf_file}, \%OPT);

    #use sequenced mouse genome to normalize by region coverage
    if($OPT{xeno}{vcf_file}){
	$OPT{xeno}{xcov} = estimate_xeno($OPT{xeno}{vcf_file}, xeno_param(\%OPT));
	$DBARC->{xeno}{xcov} = $OPT{xeno}{xcov}; #archive the value
	my $obs = $OPT{xeno}{xcov};
	my $exp = $OPT{xcov};
	$OPT{xeno}{ccor_f} = $obs/$exp;
    }
    
    $DBARC->{xcov} = $OPT{xcov}; #archive the value
    tie(%DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
    $DB{ARC} = $DBARC;
    untie(%DB);

    print STDERR "#Mean mouse coverage: $OPT{xcov}\n";
}

#--calculate coverage expects
print STDERR "#Getting expected coverage ...\n";
($OPT{t_cov}, $OPT{ref}{t_cov}) = expected_segment_coverage($disc, \%OPT);

$DBARC->{ref}{t_cov} = $OPT{ref}{t_cov}; #archive the value
$DBARC->{t_cov} = $OPT{t_cov}; #archive the value
tie(%DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
$DB{ARC} = $DBARC;
untie(%DB);

#--estimate tumor contamintation
if($OPT{cell}){
    print STDERR "#Estimating tumor contamination...\n";
    ($OPT{cfrac}, $OPT{cfrac_rsq}, $OPT{hemi_center}) = estimate_copy_fraction($disc, \%OPT);
    $OPT{rfrac}  = 1 - $OPT{cfrac};
    $OPT{k_bins} = ceil($OPT{k_bins}/$OPT{cfrac}); #adust for contamintation
    if($OPT{cfrac} <= 0.30){  #MAF becomes less informative
	$OPT{cn_max} = ($OPT{cfrac} <= 0.15) ? 1 : 4;
	$OPT{models} = get_models($OPT{cn_max}+1);
    }

    $DBARC->{cfrac} = $OPT{cfrac}; #archive the value
    tie(%DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
    $DB{ARC} = $DBARC;
    untie(%DB);

    print STDERR "#Tumor base fraction: $OPT{cfrac}\n";
    print STDERR "#R^2 of fit: $OPT{cfrac_rsq}\n" if(defined $OPT{cfrac_rsq});
}

#--add MAF to discovery segments now I have found contamintation/xenograft etc.
$disc = add_maf_discovery_stats($disc, \%OPT);

#--calculate base coverage
print STDERR "#Getting coverage expects and ploidy...\n";
($OPT{base_cov}, $OPT{ref}{base_cov}) = base_coverage($disc, \%OPT);
$OPT{ref}{chr_expects} = {_ALL => $OPT{ref}{base_cov}};
$OPT{chr_expects} = {_ALL => $OPT{base_cov}};

$DBARC->{ref}{base_cov} = $OPT{ref}{base_cov}; #archive the value
$DBARC->{base_cov} = $OPT{base_cov}; #archive the value
tie(%DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
$DB{ARC} = $DBARC;
untie(%DB);

#--get ploidy
($OPT{ref}{ploidy_ave}, $OPT{ref}{ploidy}) = get_ploidy(ref_seg($disc), ref_param(\%OPT));
($OPT{ploidy_ave}, $OPT{ploidy}) = get_ploidy($disc, \%OPT);
undef $disc; #saves memory

$DBARC->{ref}{ploidy} = $OPT{ref}{ploidy}; #archive the value
$DBARC->{ploidy} = $OPT{ploidy}; #archive the value
tie(%DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
$DB{ARC} = $DBARC;
untie(%DB);

#--report some result statistics
if($OPT{use_ref}){
    print STDERR "\n";
    print STDERR "#Ref expected segment median: $OPT{ref}{t_cov}\n";
    print STDERR "#Ref base coverage (1 copy): $OPT{ref}{base_cov}\n";
    print STDERR "#Reference ploidy mode: $OPT{ref}{ploidy}\n";
    print STDERR "#Reference ploidy mean: $OPT{ref}{ploidy_ave}\n";
}
print STDERR "\n";
print STDERR "#Sample expected segment median: $OPT{t_cov}\n";
print STDERR "#Sample base coverage (1 copy): $OPT{base_cov}\n";
print STDERR "#Sample ploidy mode: $OPT{ploidy}\n";
print STDERR "#Sample ploidy mean: $OPT{ploidy_ave}\n";
print STDERR "\n";

#--fill in the VCF/BAM data for segments
print STDERR "#Loading/processing VCF/BAM data for segments...\n";
my $files = process_segments($segs, \%OPT, 1);
undef $segs; #saves memory

#--smooth out segments and recalculate the copy numbers
if($OPT{merge}){
    print STDERR "#Merging...\n";    
    $files = merge_segments($files, \%OPT);
    $files = smooth_merge_segments($files, \%OPT) if($OPT{smooth});
}

_reassign_cn($files, \%OPT);

#--now write results to file
print STDERR "#Writing GVF...\n";
write_gvf($files, \%OPT);

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------
#run before forking. clears stored values that may cause core dump
sub prefork {
    undef $DUM;
    undef %BAMS;
    undef %VCF;
}

sub _add_maf_stats_thread {
    my $work = shift;
    my $OPT = shift;

    my @files;
    while(my $file = shift @$work){
	my $res = Storable::retrieve($file);
	foreach my $r (@$res){
	    add_maf_stats($r, $OPT) if(!$r->{maf_score});
	}
	Storable::store($res, $file);
        push(@files, $file);
    }

    return @files;
}

sub _add_cn_thread {
    my $work = shift;
    my $OPT = shift;

    my @files;
    while(my $file = shift @$work){
	my $res = Storable::retrieve($file);
	next if(!@$res);

	my $param = ref_param($OPT);
	foreach my $r (@$res){
	    assign_cn($r->{ref}, $param) if($OPT->{use_ref});
	    assign_cn($r, $OPT);
	    label_somatic($r, $OPT) if($OPT->{use_ref});
	}
	Storable::store($res, $file);
        push(@files, $file);
    }

    return @files;
}

sub _premerge_thread {
    my $work = shift;
    my $OPT = shift;

    my @files;
    while(my $file = shift @$work){
	my $res = Storable::retrieve($file);
	next if (!$res);
	delete($_->{maf_score}) foreach(@$res);
	while(1){
	    my $count = @$res;
	    $res = _merge_stdev($res, $OPT, 'ucor', 0);
	    next if($count != @$res);

	    $res = _merge_stdev($res, $OPT, 'ccor', 0);
	    last if($count == @$res);
	}

	Storable::store($res, $file);
        push(@files, $file);
    }

    return @files;
}

#returns a list of models from CN for the run
sub get_models {
    my $cn = shift;

    my @models = '!';
    foreach my $A (0..int($cn/2)){
	foreach my $B ($A..$cn-$A){
	    push(@models, "$A:$B");
	}
    }

    return \@models;
}

#validates the input options
sub validate_options {
    my $OPT = shift;

    my $seg_file  = $OPT->{seg_file};
    my $vcf_file  = $OPT->{vcf_file};
    my $rvcf_file = $OPT->{ref}{vcf_file};
    my $xvcf_file = $OPT->{xeno}{vcf_file};
    my $fasta     = $OPT->{fasta};
    my $sid = $OPT->{sid};
    my $rid = $OPT->{ref}{sid};
    my $xid = $OPT->{xeno}{sid};
    my $bam_dir = $OPT->{bam_dir};
    my $rbam_dir = $OPT->{rbam_dir};
    my $xbam_dir = $OPT->{xbam_dir};
    my $bam_list = $OPT->{bam_list};
    my $rbam_list = $OPT->{rbam_list};
    my $xbam_list = $OPT->{xbam_list};

    #Am I using reference?
    if($OPT{ref}{sid} || $OPT{ref}{vcf_file} || $rbam_dir || $rbam_list){
	$OPT{use_ref} = 1;
	$OPT{use_generic} = 0;
    }
    else{
	$OPT{use_ref} = 0;
	$OPT{use_generic} = 1;
    }
    
    #Am I using reference?
    if($xid || $xvcf_file || $xbam_dir || $xbam_list || $OPT->{xcov}){
	$OPT{use_xeno} = 1;
    }
    else{
	$OPT{use_xeno} = 0;
    }
    
    #check for files
    my $err;
    $err .= "File $seg_file does not exist\n" if(! -f $seg_file);
    $err .= "File $fasta does not exist\n" if(! -f $fasta);
    $err .= "File $vcf_file does not exist\n" if(! -f $vcf_file);
    $err .= "File $rvcf_file does not exist\n" if($rvcf_file && ! -f $rvcf_file);
    $err .= "File $xvcf_file does not exist\n" if($xvcf_file && ! -f $xvcf_file);
    
    #fix columns for chromosome, start, and end
    foreach ($OPT->{C}, $OPT->{B}, $OPT->{E}){
	$_ = int($_-1);
	$err .= "Column selection must be a positive integer > 1\n" if($_ < 0);
    }
    
    #check if BAM files exist
    my @bam_files;
    if($bam_dir){
	if(!-d $bam_dir){
	    $err .= "Directory $bam_dir does not exist\n" if(! -d $bam_dir);
	}
	else{
	    push(@bam_files, <$bam_dir/*.bam>);
	    $err .= "$bam_dir contained no bam file\n" if(! @bam_files);
	}
    }
    if($bam_list){
	my @blist = split(/\,/, $bam_list);
	if(@blist){
	    for(@blist){
		$err .= "File $_ from bam_list does not exist\n" if(!-f $_);
	    }
	    push(@bam_files, @blist);
	}
	else{
	    $err .= "$bam_list contained no bam file\n" if(! @blist);
	}
    }
    
    #check if reference bam files exist
    my @rbam_files;
    if($rbam_dir){
	if(!-d $rbam_dir){
	    $err .= "Directory $rbam_dir does not exist\n" if(! -d $rbam_dir);
	}
	else{
	    push(@rbam_files, <$rbam_dir/*.bam>);
	    $err .= "$rbam_dir contained no bam file\n" if(! @rbam_files);
	}
    }
    if($rbam_list){
	my @blist = split(/\,/, $rbam_list);
	if(@blist){
	    for(@blist){
		$err .= "File $_ from rbam_list does not exist\n" if(!-f $_);
	    }
	    push(@rbam_files, @blist);
	}
	else{
	    $err .= "$rbam_list contained no bam file\n" if(! @blist);
	}
    }
    
    #check if xenograft bam files exist
    my @xbam_files;
    if($xbam_dir){
	if(!-d $xbam_dir){
	    $err .= "Directory $xbam_dir does not exist\n" if(! -d $xbam_dir);
	}
	else{
	    push(@xbam_files, <$xbam_dir/*.bam>);
	    $err .= "$xbam_dir contained no bam file\n" if(! @xbam_files);
	}
    }
    if($xbam_list){
	my @blist = split(/\,/, $xbam_list);
	if(@blist){
	    for(@blist){
		$err .= "File $_ from xbam_list does not exist\n" if(!-f $_);
	    }
	    push(@xbam_files, @blist);
	}
	else{
	    $err .= "$xbam_list contained no bam file\n" if(! @blist);
	}
    }
    
    return $err if($err);
    
    #--prepare FASTA related values
    prepare_fasta($fasta);
    $OPT->{CHRS} = get_chr_lengths($fasta);
    $OPT->{GAPS} = get_chr_gaps($fasta);
    %CHRS = %{$OPT->{CHRS}};
    
    prepare_vcf($vcf_file);
    
    #validate sample ID
    if($vcf_file){
	my $vcf = $VCF{$vcf_file};
	if(!$sid){    
	    die "ERROR: VCF contains multiple samples and you failed to specify a sample ID\n"
		if(@{[$vcf->get_samples()]} > 1);
	    ($sid) = $vcf->get_samples();
	}
	($sid) = grep {/$sid$/} $vcf->get_samples();
    }
    
    #set reference and xenograft vcf if appropriate
    if($rid && !$rvcf_file){
	my $vcf = $VCF{$vcf_file};
	if(grep {/$rid$/} $vcf->get_samples()){
	    $rvcf_file = $vcf_file;
	}
	else{
	    warn "WARNING: You have supplied a sample ID for the reference but no VCF is associated\n";
	}
    }
    if($xid && !$xvcf_file){
	my $vcf = $VCF{$vcf_file};
        if(grep {/$xid$/} $vcf->get_samples()){
	    $xvcf_file = $vcf_file;
        }
	else{
	    warn "WARNING: You have supplied a sample ID for the xenograft but no VCF is associated\n";
	}
    }
    prepare_vcf($rvcf_file);
    prepare_vcf($xvcf_file);
    
    #validate reference ID
    if($rvcf_file){
	my $rvcf = $VCF{$rvcf_file};
	if(!$rid){
	    die "ERROR: Reference VCF contains multiple samples and you failed to specify a sample ID\n"
		if(@{[$rvcf->get_samples()]} > 1);
	    ($rid) = $rvcf->get_samples();
	}
	($rid) = grep {/$rid$/} $rvcf->get_samples();;
    }
    
    #validate xenograft ID
    if($xvcf_file){
	my $xvcf = $VCF{$xvcf_file};
	if(!$xid){
	    die "ERROR: Xenograft VCF contains multiple samples and you failed to specify a sample ID\n"
		if(@{[$xvcf->get_samples()]} > 1);
	    ($xid) = $xvcf->get_samples();
	}
	($xid) = grep {/$xid$/} $xvcf->get_samples();;
    }
    
    #load BAM files
    prepare_bam(@bam_files);
    
    #set reference and xenograft BAM files if appropriate
    if($OPT->{use_ref} && !$rbam_dir && !$rbam_list && @bam_files){
	if(grep {$_ eq $rid} map {keys %{$_->{_samples}}}  @{$BAMS{_files}}{@bam_files}){ #check id
	    @rbam_files = @bam_files;
	}
	else{
	    warn "WARNING: You have supplied bamfiles for you sample but not your reference\n";
	}
    }
    if($OPT->{use_xeno} && !$xbam_dir && !$xbam_list && @bam_files){
	if(grep {$_ eq $xid} map {keys %{$_->{_samples}}}  @{$BAMS{_files}}{@bam_files}){ #check id
            @xbam_files = @bam_files;
	}
	else{
	    warn "WARNING: You have supplied bamfiles for you sample but not your xenograft\n";
	}
    }
    
    prepare_bam(@rbam_files);
    prepare_bam(@xbam_files);
    
    my ($base) = $seg_file =~ /([^\/]+)$/;
    my $outdir = Cwd::cwd()."/$base.cnv.output";
    mkdir($outdir);
    
    #fix dependent values
    $OPT->{fasta} = $fasta;
    $OPT->{sid} = $sid;
    $OPT->{ref}{sid} = $rid;
    $OPT->{xeno}{sid} = $xid;
    $OPT->{vcf_file} = $vcf_file;
    $OPT->{ref}{vcf_file} = $rvcf_file;
    $OPT->{xeno}{vcf_file} = $xvcf_file;
    $OPT->{bam_files} = \@bam_files;
    $OPT->{ref}{bam_files} = \@rbam_files;
    $OPT->{xeno}{bam_files} = \@xbam_files;
    $OPT->{outdir} = $outdir;
    $OPT->{base} = $base;
    $OPT->{merge} = 1 if($OPT->{smooth}); #implies merge
    
    #adust for contamintation
    if($OPT->{cfrac} < 1){
	$OPT->{cell} = 1;
	$OPT->{k_bins} = ceil($OPT->{k_bins}/$OPT->{cfrac});
	$OPT->{cn_max} = ($OPT->{cfrac} <= 0.25) ? 4 : 5;
	$OPT->{models} = get_models($OPT->{cn_max}+1);
    }
    
    return $err;
}

#build fasta index
sub prepare_fasta {
    my @fastas = @_;

    foreach my $fasta (@fastas){
	next if(!$fasta);

	my $bf = $dummy;
	my $bam = Bio::DB::Sam->new(-bam  => $bf,
				    -autoindex => 1,
				    -fasta=> $fasta);
    }

    return (wantarray) ? @fastas : shift @fastas;
}

#indexes the bam files
sub prepare_bam {
    my @bam_files = (@_);

    foreach my $f (@bam_files){
	my $bam = $BAMS{_files}{$f} || Bio::DB::Sam->new(-bam  => $f,
							 -autoindex => 1,
							 -fasta=> $OPT{fasta});
	#get read groups for filtering
	if(!$bam->{_samples}){
	    my @headers = grep {/^\@RG\t/} split(/\n/, $bam->header->text);
	    foreach my $h (@headers){
		my %rg;
		foreach (grep {/^(ID:.*|SM:.*)$/} split(/\t/, $h)){
		    my @F = split(/\:/, $_);
		    $rg{$F[0]} = join(':', @F[1..$#F]);
		}
		next unless(keys(%rg) == 2);
		push(@{$bam->{_samples}{$rg{SM}}}, $rg{ID});
	    }
	}

	if(@bam_files == 1){
	    $BAMS{_files}{$f} = $bam;
	    foreach my $sid (keys %{$bam->{_samples}}){
		$BAMS{$sid}{$_} = $bam foreach(keys %CHRS);
	    }
	}
	else{
	    next unless($f =~ /\.(chr[\dXY]+)\./ && $CHRS{$1});
	    $BAMS{_files}{$f} = $bam;
	    foreach my $sid (keys %{$bam->{_samples}}){
		$BAMS{$sid}{$1} = $bam;
	    }
	}

    }

    return (wantarray) ? @bam_files : shift @bam_files;
}

#indexes the vcf files
sub prepare_vcf {
    my @vcf_files = @_;

    foreach my $vcf_file (@vcf_files){
	next if(!$vcf_file);
	next if($VCF{$vcf_file});

	my $bgzfile = ($vcf_file =~ /\.gz$/) ? $vcf_file : "$vcf_file.gz";
	next if($VCF{$bgzfile});
	if(! -f "$bgzfile"){
	    open(my $OUT, "> $bgzfile.tmp");
	    open(my $IN, "< $vcf_file");
	    my $spid = open3(my $SORT, my $SORTED, '>&STDERR', "$sort_exe -k1,1 -k2,2n");
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

        my $vcf = $VCF{$vcf_file} || $VCF{$bgzfile} || Vcf->new(file=>"$bgzfile");
        $vcf->parse_header();
        $VCF{$vcf_file} = $VCF{$bgzfile} = $vcf;
    }

    return (wantarray) ? @vcf_files : shift @vcf_files;
}

sub get_chr_lengths {
    my $fasta = shift;
    
    my %CHRS;
    if(Cwd::abs_path($fasta) eq Cwd::abs_path($hg19)){
	%CHRS = (chr1 => 249250621,
		 chr2 => 243199373,
		 chr3 => 198022430,
		 chr4 => 191154276,
		 chr5 => 180915260,
		 chr6 => 171115067,
		 chr7 => 159138663,
		 chr8 => 146364022,
		 chr9 => 141213431,
		 chr10 => 135534747,
		 chr11 => 135006516,
		 chr12 => 133851895,
		 chr13 => 115169878,
		 chr14 => 107349540,
		 chr15 => 102531392,
		 chr16 => 90354753,
		 chr17 => 81195210,
		 chr18 => 78077248,
		 chr19 => 59128983,
		 chr20 => 63025520,
		 chr21 => 48129895,
		 chr22 => 51304566,
		 chrX => 155270560,
		 chrY => 59373566);
    }
    else{
	#not implemented yet
    }
    
    return \%CHRS;
}

#used to split chr into smaller chunks (i.e split on centromeres)
sub get_chr_gaps {
    my $fasta = shift;
    
    my %GAPS;
    if(Cwd::abs_path($fasta) eq Cwd::abs_path($hg19)){
	%GAPS = (chr1 => [[121535434, 124535434]],
		 chr2 => [[92326171, 95326171]],
		 chr3 => [[90504854, 93504854]],
		 chr4 => [[49660117, 52660117]],
		 chr5 => [[46405641, 49405641]],
		 chr6 => [[58830166, 61830166]],
		 chr7 => [[58054331, 61054331]],
		 chr8 => [[43838887, 46838887]],
		 chr9 => [[47367679, 50367679]],
		 chrX => [[58632012, 61632012]],
		 chrY => [[10104553, 13104553]],
		 chr10 => [[39254935, 42254935]],
		 chr11 => [[51644205, 54644205]],
		 chr12 => [[34856694, 37856694]],
		 chr13 => [[16000000, 19000000]],
		 chr14 => [[16000000, 19000000]],
		 chr15 => [[17000000, 20000000]],
		 chr16 => [[35335801, 38335801]],
		 chr17 => [[22263006, 25263006]],
		 chr18 => [[15460898, 18460898]],
		 chr19 => [[24681782, 27681782]],
		 chr20 => [[26369569, 29369569]],
		 chr21 => [[11288129, 14288129]],
		 chr22 => [[13000000, 16000000]]);
    }
    else{
	#not implemented yet
    }
    
    return \%GAPS;
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
    prepare_vcf($vcf_file);
    my $vcf = $VCF{$vcf_file};
    $vcf->set_samples(include=>[$sid]);
    
    my %list;
    while((my $chr = shift @$work) || $$flag){
	next if(!$chr);

	$vcf->open(region=> $chr);
	
	while(my $v = $vcf->next_data_hash()){
	    my $pos = $v->{POS};
	    next unless($marker->{$chr}{$pos});
	    next unless($v->{gtypes}{$sid}{AD});
	    
	    my ($rc, $ac) = split(/,/, $v->{gtypes}{$sid}{AD});
	    next if(!$ac);
	    $list{$ac}++;
	}
    }

    return \%list;
}

#read in segments from a file
sub read_segment_file {
    my $seg_file  = shift;
    my $OPT = shift;
    my $C = $OPT->{C};
    my $B = $OPT->{B};
    my $E = $OPT->{E};

    my @segs;
    open(IN, "< $seg_file");
    while(my $line = <IN>){
	chomp $line;
	next if(! $line);
	my @F = map {s/^\s+|\s+$//g; $_} split(/\t/, $line);
	$F[$C] = "chr$F[$C]" if(!$CHRS{$F[$C]} && $CHRS{"chr$F[$C]"});
	next if(!$CHRS{$F[$C]}); #temp
	push(@segs, "$F[$C]\:$F[$B]\-$F[$E]");
    }
    close(IN);

    return \@segs;
}

#create results has from segment coordinates
sub process_segments {
    my $segs     = shift;
    my $OPT      = shift;

    my $flag :shared;
    my $counter :shared;
    my @work :shared;
    my @threads;
    $flag = 1;
    $counter = 0;

    #see what is already finished
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my @files;
    my %seen;
    my $load = $DB{load};
    foreach my $file (keys %$load) {
	next if($file =~ /^_/); #these are parameters

	if(!-f $file){
	    delete $load->{$file};
	    next;
	}
	
	push(@files, $file);
	$seen{$_}++ foreach(@{$load->{$file}});
    }
    $DB{load} = $load;
    untie(%DB);

    #add work stats to status bar and to shared list of what to do
    my $work_count = 0;
    my $finished = 0; #not shared for speed
    my @work_add; #not shared for speed
    foreach my $s (@$segs){
	my ($start, $end) = $s =~ /\:(\d+)\-(\d+)$/;
	$work_count += ($end-$start)+1;
	if($seen{$s}){
	    $finished += ($end-$start)+1;
	    next;
	}
	push(@work_add, $s);
	$load->{_grouped} = 0;
    }
    push(@work, @work_add);
    $flag = 0; #lets threads know everything is loaded
    my $remain = $work_count - $finished;
    
    #progress bar parameters
    my %bparam = (name  => 'Processing',
		  count => $remain,
		  ETA   => 'linear');
    
    #fix terminal size capture in SGE redirect
    if(! is_interactive(*STDERR)){
	$ENV{COLUMNS} = 80;
	$bparam{bar_width}  = 80;
	$bparam{term_width} = 80;
    }
    
    printf STDERR "Total: %ibp\n", $work_count;
    printf STDERR "Processed: %ibp (%.0f%%)\n", $finished, 100*$finished/$work_count;
    printf STDERR "Remaining: %ibp (%.0f%%)\n", $remain, 100*$remain/$work_count;

    #read data into segments and store in files
    if(@work){
	#launch helper threads
	my $cpus = $OPT->{cpus};
	prefork(); #prepare for forking
	for(my $i = 1; $i < $cpus; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_process_segments_thread,
					   \@work,
					   $OPT,
					   \$flag,
					   \$counter,
					   $remain));
	}
	
	#create progress bar object
	my $prog = Term::ProgressBar->new(\%bparam);
	$prog->update($counter);

	#--jump into the mix
	my @ret = _process_segments_thread(\@work,
					   $OPT,
					   \$flag,
					   \$counter,
					   $remain,
					   $prog);    
	push(@files, @ret);

	#--collect thread results
	foreach my $thr (@threads){
	    my @ret = $thr->join;
	    $prog->update($counter) if($prog && $prog->target != $counter);
	    push(@files, @ret);
	}
	$prog->update($work_count) if($prog && $prog->target != $remain);
	$load->{_grouped} = 0;
    }

    #reorder serialized store to represent regions rather than random order
    if(!$load->{_grouped}){
	#split into groups based on gaps
	my $res;
	foreach my $file (@files){
	    my $st = Storable::retrieve($file);
	    push(@$res, @$st);
	}
	$res = _merge_special($res, $OPT, 0);
	my $sets = _split_on_gaps($res, $OPT);
	
	#reserialize based on groups
	my @old = @files;
	undef @files;
	map {delete($load->{$_})} grep {!/^_/} keys %$load;
	foreach my $s (values %$sets){
	    my @ids;
	    push(@ids, @{$_->{contains}}) foreach(@$s);
	    (undef, my $tfile) = tempfile("load_XXXX", DIR => $OPT->{outdir}."/load");
	    Storable::store($s, $tfile);
	    $load->{$tfile} = \@ids;
	    push(@files, $tfile);
	}
        
	#update load status
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	$load->{_grouped} = 1;
	$load->{_premerge} = 0;
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;
	unlink(@old);
    }

    #initial merge of very short segments
    if(!$load->{_premerge}){
	prefork(); #prepare for forking
	my @threads;
	push(@work, @files);
	for(my $i = 1; $i < $OPT->{cpus}; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_premerge_thread,
					   \@work,
					   $OPT));
	}
	_premerge_thread(\@work, $OPT); #jump into the mix
	$_->join foreach(@threads); #gather thread results

	#update load status
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	$load->{_premerge} = 1;
	$load->{_hasmaf} = 0;
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;
    }

    #add MAF to the segments
    if(!$load->{_hasmaf}){
	prefork(); #prepare for forking
	my @threads;
	push(@work, @files);
	for(my $i = 1; $i < $OPT->{cpus}; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_add_maf_stats_thread,
					   \@work,
					   $OPT));
	}
	_add_maf_stats_thread(\@work, $OPT); #jump into the mix
	$_->join foreach (@threads); #gather thread results

	#update load status
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	$load->{_hasmaf} = 1;
	$load->{_hascn} = 0;
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;
    }

    #add CN to each segment
    if(!$load->{_hascn}){
	prefork(); #prepare for forking
	my @threads;
	push(@work, @files);
	for(my $i = 1; $i < $OPT->{cpus}; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_add_cn_thread,
					   \@work,
					   $OPT));
	}
	_add_cn_thread(\@work, $OPT); #jump into the mix
	$_->join foreach (@threads); #gather thread results

	#update load status
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	$load->{_hascn} = 1;
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;
    }

    return \@files;
}

sub _process_segments_thread {
    my $work    = shift;
    my $OPT     = shift;
    my $flag    = shift;
    my $counter = shift;
    my $total   = shift;
    my $prog    = shift;

    my $fasta      = $OPT->{fasta};
    my $sid        = $OPT->{sid};
    my $vcf_file   = $OPT->{vcf_file};
    my $rid        = $OPT->{ref}{sid};
    my $rvcf_file  = $OPT->{ref}{vcf_file};
    my $thr        = $OPT->{maf_tail_filt};
    my $cfrac      = $OPT->{cfrac};
    my $rfrac      = 1 - $cfrac;

    #now process segments
    my @res;
    my @files;
    my $loaded = [];
    my $next_update = 0;
    my $counter_buf = 0; #faster to update than shared $counter
    while((my $seg = shift @$work) || $$flag){
	#update progress bar
	if($counter_buf > $total/1000){
	    $$counter += $counter_buf;
	    $counter_buf = 0;
	    $next_update = $prog->update($$counter) if($prog);
	}

	#nothing to do
	if(! $seg){
	    sleep 1;
	    next;
	}

	#fill in segment hash
	my ($chr, $start, $end) = $seg =~ /^(.*)\:(\d+)\-(\d+)$/;
	my $r = make_seg_hash($chr, $start, $end, $OPT);
	$r = add_vcf_data($r, $OPT);

	#add bam coverage only for short segments
	if($r->{q_length} < 30000 || $r->{maf_count} < 30 || $r->{cov_count}/$r->{q_length} < 0.0005){
	    $r = add_bam_cov($r, $OPT);
	}
	elsif($r->{xeno}){ #always use bam for at xeno match
	    $r->{xeno} = add_bam_cov($r->{xeno}, xeno_param($OPT));
	}

	#add coverage correction
	if(my $ref = $r->{ref}){
	    my $param = ref_param($OPT);
	    my $obs = $ref->{cov_median};
	    my $exp = $param->{t_cov};
	    $r->{ccor_f} = $ref->{ccor_f} = $obs/$exp;
	}

	push(@res, $r); #must have at least some internal info or must be ignored
	push(@$loaded, $seg);
	$counter_buf += ($end-$start)+1; #itterate progress bar

	#keep memory usage lower
	if(@res > 100){
	    mkdir($OPT->{outdir}."/load");
	    (undef, my $tfile) = tempfile("load_XXXX", DIR => $OPT->{outdir}."/load");
	    Storable::store(\@res, $tfile);
	    push(@files, $tfile);

	    my $lock;
	    $lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	    my $load = $DB{load};
	    $load->{$tfile} = $loaded;
	    $DB{load} = $load;
	    untie(%DB);
	    $lock->unlock;

	    undef @res;
	    $loaded = [];
	}
    }

    #keep memory usage lower
    if(@res){
	mkdir($OPT->{outdir}."/load");
	(undef, my $tfile) = tempfile("load_XXXX", DIR => $OPT->{outdir}."/load");
	Storable::store(\@res, $tfile);
	push(@files, $tfile);

	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	my $load = $DB{load};
	$load->{$tfile} = $loaded;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;
	
	undef @res;
	$loaded = [];
    }

    return (wantarray) ? @files : \@files;
    #return (wantarray) ? @res : \@res;
}

sub make_seg_hash {
    my $chr   = shift;
    my $start = shift;
    my $end   = shift;
    my $OPT   = shift;

    my $fasta = $OPT->{fasta};

    #load dummy BAM file (hack to get fasta indexer to work)
    my $bf = $dummy;
    $DUM ||= Bio::DB::Sam->new(-bam  => $bf,
			       -autoindex => 1,
			       -fasta=> $fasta); #static for performance    

    my $segment = $DUM->segment($chr, $start, $end);
    my $dna = $segment->dna;

    my $length = ($end-$start)+1;
    my $mask_count = $dna =~ tr/Nn/Nn/;
    my $q_length = $length - $mask_count; #length adjusted for N's
    my $id = "$chr:$start-$end";

    #put it all in the segment result
    my $res = { chr        => $chr,
		start      => $start,
		end        => $end,
		contains   => [$id],
		length     => $length,
		q_length   => $q_length, #how much queriable
		mask_count => $mask_count,
		cov_count  => 0,
		cov_mean   => 0,
		cov_median => 0,
		cov_bam    => undef,
		cov_vcf    => undef,
		cov_set    => undef,
		maf_count  => 0,
		maf_set    => undef,
		maf_score  => undef};
    
    if($OPT->{use_ref} || $OPT->{use_generic}){
	#need reference
	$res->{ref} = { chr        => $chr,
			start      => $start,
			end        => $end,
			contains   => [$id],
			length     => $length,
			q_length   => $q_length, #how much queriable
			mask_count => $mask_count,
			cov_count  => 0,
			cov_mean   => 0,
			cov_median => 0,
			cov_bam    => undef,
			cov_vcf    => undef,
			cov_set    => undef,
			maf_count  => 0,
			maf_set    => undef,
			maf_score  => undef};
	$res->{ref}{is_generic} = 1 if($OPT->{use_generic});
	$res->{ref}{is_ref} = 1 if($OPT->{use_ref});
    }
    $res->{is_generic} = 1 if($OPT->{is_generic});
    $res->{is_ref} = 1 if($OPT->{is_ref});

    if($OPT->{use_xeno}){
        #need xenograft model
        $res->{xeno} = { chr        => $chr,
			 start      => $start,
			 end        => $end,
			 contains   => [$id],
			 length     => $length,
			 q_length   => $q_length, #how much queriable
			 mask_count => $mask_count,
			 cov_count  => 0,
			 cov_mean   => 0,
			 cov_median => 0,
			 cov_bam    => undef,
			 cov_vcf    => undef,
			 cov_set    => undef,
			 maf_count  => 0,
			 maf_set    => undef,
			 maf_score  => undef};
        $res->{xeno}{is_xeno} = 1;
    }
    $res->{is_xeno} = 1 if($OPT->{is_xeno});
    
    return $res;
}


sub add_vcf_data {
    my $res = shift;
    my $OPT = shift;

    my $sid        = $OPT->{sid};
    my $vcf_file   = $OPT->{vcf_file};
    my $rid        = $OPT->{ref}{sid};
    my $rvcf_file  = $OPT->{ref}{vcf_file};
    my $xid        = $OPT->{xeno}{sid};
    my $xvcf_file  = $OPT->{xeno}{vcf_file};
    my $thr        = $OPT->{maf_tail_filt};

    #load VCF file
    prepare_vcf($vcf_file);
    my $vcf = $VCF{$vcf_file};

    $vcf->set_samples(include=>[$sid]);
   
    #load reference VCF
    my $rvcf;
    if($OPT->{use_ref}){
	if($rvcf_file && $rvcf_file ne $vcf_file){
	    prepare_vcf($rvcf_file);
	    $rvcf = $VCF{$rvcf_file};
	    $rvcf->set_samples(include=>[$rid]);
	}
	else{
	    $vcf->set_samples(include=>[$sid, $rid]);
	}
    }

    #load xeno VCF
    my $xvcf;
    if($OPT->{use_xeno}){
	if($xvcf_file && $xvcf_file ne $vcf_file){
	    prepare_vcf($xvcf_file);
	    $xvcf = $VCF{$xvcf_file};
	    $xvcf->set_samples(include=>[$xid]);
	}
	else{
	    my @ids = ($sid, $xid);
	    push(@ids, $rid) if($OPT->{use_ref} && !$rvcf);
	    $vcf->set_samples(include=>[@ids]);
	}
    }

    #get observed MAF and VCF coverage for segment
    my $cov_count   = 0;
    my $rcov_count  = 0;
    my $xcov_count  = 0;
    my $maf_count   = 0;
    my $rmaf_count  = 0;
    my $xmaf_count  = 0;
    my $rmaf_original = 0;
    my %vcf_set;
    my %rvcf_set;
    my %xvcf_set;
    my %maf_set;
    my %rmaf_set;
    my %xmaf_set;

    #calculate first for separate reference VCF
    my %ref_ok;
    my ($chr, $start, $end) = ($res->{chr}, $res->{start}, $res->{end});
    if($rvcf){
	$rvcf->open(region=> "$chr:$start-$end");
	while(my $v = $rvcf->next_data_hash()){
	    my $pos = $v->{POS};
	    my $rAD = $v->{gtypes}{$rid}{AD};
	    if($rAD && $rAD ne '.'){
		my ($rc, $ac) = split(/,/, $rAD);
		my $cov = $rc+$ac;
		$rvcf_set{$cov}++;
		$rcov_count++;

		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    $rmaf_original++ if($thr <= $maf);
		    if($thr <= $maf && $maf <= 1-$thr){
			$rmaf_set{$maf}{$cov}++;
			$rmaf_count++;
			$ref_ok{$chr}{$pos}++;
		    }
		}
	    }
	}
    }
    
    #calculate from sample VCF
    $vcf->open(region=> "$chr:$start-$end");
    while(my $v = $vcf->next_data_hash()){
	my $pos = $v->{POS};

	#always calculate  reference first if merged in same VCF
	if($OPT->{use_ref} && !$rvcf){
	    my $rAD = $v->{gtypes}{$rid}{AD};
	    if($rAD && $rAD ne '.'){
		my ($rc, $ac) = split(/,/, $rAD);
		my $cov = $rc+$ac;
		$rvcf_set{$cov}++;
		$rcov_count++;

		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($thr <= $maf && $maf <= 1-$thr){
                        $rmaf_set{$maf}{$cov}++;
                        $rmaf_count++;
			$ref_ok{$chr}{$pos}++;
                    }
		}
	    }
	}

	#calculate for sample
	my $sAD = $v->{gtypes}{$sid}{AD};
	if($sAD && $sAD ne '.'){
	    my ($rc, $ac) = split(/,/, $sAD);
	    my $cov = $rc+$ac;
	    $vcf_set{$cov}++;
	    $cov_count++;

	    if($cov >= $OPT->{maf_cov_filt}){
		my $maf = $ac/$cov;
		if($maf >= $thr && ($OPT->{use_generic} || $ref_ok{$chr}{$pos})){
		    $maf_set{$maf}{$cov}++;
		    $maf_count++;
		}
	    }
	}

	#calculate when xeno is merged with sample VCF
	if($res->{use_xeno} && !$xvcf){
	    my $xAD = $v->{gtypes}{$xid}{AD};
	    if($xAD && $xAD ne '.'){
		my ($rc, $ac) = split(/,/, $xAD);
		my $cov = $rc+$ac;
		$xvcf_set{$cov}++;
		$xcov_count++;
		
		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($maf >= $thr){
			$xmaf_set{$maf}{$cov}++;
			$xmaf_count++;
		    }
		}
	    }
	}
    }

    #calculate for separate xeno VCF
    if($xvcf){
	$xvcf->open(region=> "$chr:$start-$end");
	while(my $v = $xvcf->next_data_hash()){
	    my $xAD = $v->{gtypes}{$xid}{AD};
	    if($xAD && $xAD ne '.'){
		my ($rc, $ac) = split(/,/, $xAD);
		my $cov = $rc+$ac;
		$xvcf_set{$cov}++;
		$xcov_count++;
		
		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($maf >= $thr){
			$xmaf_set{$maf}{$cov}++;
			$xmaf_count++;
		    }
		}
	    }
	}
    }

    my ($cov_mean,  $cov_median)  = _cov_stats(\%vcf_set);

    $res->{cov_vcf}    = \%vcf_set;
    $res->{cov_set}    = $res->{vcf_set};
    $res->{cov_count}  = $cov_count;
    $res->{cov_mean}   = $cov_mean;
    $res->{cov_median} = $cov_median;
    $res->{maf_set}    = \%maf_set;
    $res->{maf_count}  = $maf_count;

    if($OPT->{use_ref}){
	my ($rcov_mean, $rcov_median) = _cov_stats(\%rvcf_set);

	$res->{ref}{cov_vcf}    = \%rvcf_set;
	$res->{ref}{cov_set}    = $res->{ref}{vcf_set};
	$res->{ref}{cov_count}  = $rcov_count;
	$res->{ref}{cov_mean}   = $rcov_mean;
	$res->{ref}{cov_median} = $rcov_median;
	$res->{ref}{maf_set}    = \%rmaf_set;
	$res->{ref}{maf_count}  = $rmaf_count;
	$res->{ref}{maf_original} = $rmaf_original;
    }

    if($OPT->{use_xeno} && ($OPT->{xeno}{vcf_file} || $OPT->{xeno}{sid})){
	my ($xcov_mean, $xcov_median) = _cov_stats(\%xvcf_set);

	$res->{xeno}{cov_vcf}    = \%xvcf_set;
	$res->{xeno}{cov_set}    = $res->{xeno}{vcf_set};
	$res->{xeno}{cov_count}  = $xcov_count;
	$res->{xeno}{cov_mean}   = $xcov_mean;
	$res->{xeno}{cov_median} = $xcov_median;
	$res->{xeno}{maf_set}    = \%xmaf_set;
	$res->{xeno}{maf_count}  = $xmaf_count;
    }

    return $res;
}

sub add_maf_stats {
    my $r          = shift;
    my $OPT        = shift;
    my $no_ref     = shift;

    return if($r->{is_xeno}); #never add to xeno

    if($r->{ref} && !$no_ref){
	add_maf_stats($r->{ref}, ref_param($OPT));
    }

    my $q_length   = $r->{q_length};
    my $cov_mean   = $r->{cov_mean};
    my $maf_count  = $r->{maf_count};
    my $maf_set    = $r->{maf_set};
    my $cn_max     = $OPT->{cn_max};
    my $xcov       = $OPT->{xcov};
    my $cfrac      = $OPT->{cfrac};
    my $rfrac      = 1 - $cfrac;

    #get target frequencies (optimization for redundancy)
    my %target_fs;
    foreach my $m (@{$OPT->{models}}){
	next if($m eq '!');
	push(@{$target_fs{get_f_string($m, $cfrac)}}, $m);
    }

    #==rss fit of models with AIC relative liklihoods
    my %maf_score;
    my $obs = maf_bin($maf_set, $OPT->{k_bins}, $OPT); #bin the observed MAF values
    foreach my $t (keys %target_fs){
	next if($t eq '!'); #temp
	my @ms = @{$target_fs{$t}};
	my $m = $ms[0];

	my ($rss, $k, $aic) = (0, 1, 0);
	if($q_length > 30000 && $maf_count){
	    my $expect = maf_expect($maf_set, $m, $cfrac, 0, $OPT); #get expected MAF distribution
	    ($rss, $k) = stat_fit($obs, $expect);
	    my $n = 1 + ($r->{q_length}-1)/$OPT->{lfrag};
	    $n = $maf_count if($maf_count < $n);
	    $n = 100 if($n > 100);
	    $aic = $n*log($rss/$n)+2*$k;
	}

	#set all models with this target freq
	foreach my $m (@ms){
	    $maf_score{$m}{rss} = $rss;
	    $maf_score{$m}{k} = $k;
	    $maf_score{$m}{AIC} = $aic;
	}

	#calculate for xenograft
	next unless($xcov);
	if($q_length > 30000 && $maf_count){
	    my $expect = maf_expect($maf_set, $m, $cfrac, $xcov, $OPT); #get expected MAF distribution
	    ($rss, $k) = stat_fit($obs, $expect);
	    my $n = 1 + ($r->{q_length}-1)/$OPT->{lfrag};
	    $n = $maf_count if($maf_count < $n);
	    $n = 100 if($n > 100);
            $aic = $n*log($rss/$n)+2*$k;
	}

	foreach my $m (@ms){
	    $maf_score{"$m:xeno"}{rss} = $rss;
	    $maf_score{"$m:xeno"}{k} = $k;
	    $maf_score{"$m:xeno"}{AIC} = $aic;
	}
    }

    #determine AIC relative liklihoods
    my @all = map {$maf_score{$_}} keys %maf_score;
    my ($aic_min) = sort {$a <=> $b} map {$_->{AIC}} @all;
    foreach my $s (values %maf_score){
	$s->{L} = exp(($aic_min-$s->{AIC})/2);
    }

    #==mark beast alleles
    _mark_best_alleles(\%maf_score, $OPT);
    $r->{maf_score} = \%maf_score;

    return $r;
}

#Put observed data into hitogram bins
sub maf_bin {
    my $maf_set = shift;
    my $bins    = shift || 100;
    my $OPT     = shift || 0;

    my @obs   = map {0} (0..$bins-1);
    keys %{$maf_set};
    while(defined(my $maf = each %{$maf_set})){
	next if($maf < $OPT{maf_tail_filt}); #filter low maf tail

	keys %{$maf_set->{$maf}};
	while(my ($cov, $num) = each %{$maf_set->{$maf}}){
	    _bin_it($bins,
		    $cov,
		    $num,
		    $maf,
		    \@obs);
	}
    }

    my $sum = 0;
    foreach my $o (@obs){
	$sum  += $o;
    }
    if($sum){
	$_ /= $sum foreach(@obs);
    }

    return \@obs;
}

sub _bin_it {
    my $bins = shift;
    my $cov  = shift;
    my $num  = shift;
    my $maf  = shift;
    my $hist = shift;

    #now get MAF range for coverage
    my $step = ($cov) ? 1/(2*$cov) : 1;
    my $mB = $maf - $step;
    my $mE = $maf + $step;
    $mB = 0 if($mB < 0);
    $mE = 1 if($mE > 1);

    #get bin range for MAF range
    my $bB = $mB * $bins;
    my $bE = $mE * $bins;
    
    #get fraction belonging to each bin
    my $amp = 0;
    foreach my $i (int($bB)..int($bE)){
	next if($i == $bins); #cannot walk past edge
	my $min = $i;
	my $max = $i+1;
	if($min < $bB){
	    $min += ($bB - $min);
	}
	
	if($max > $bE){
	    $max -= ($max - $bE);
	}
	
	#binned by probability density
	$amp = ($max-$min)/(2*$step*$bins) * $num;

	#binned by probability
	#$amp = ($max-$min)/($bE-$bB) * $num;

	$hist->[$i] += $amp;
    }
}

#returns the PDF of the expected model (for the given coverages)
{
my %ARC;
sub maf_expect {
    my $maf_set  = shift;
    my $mod      = shift;
    my $cfrac    = shift;
    my $xcov     = shift;
    my $OPT      = shift;    

    my $rfrac   = 1 - $cfrac;
    my $thr     = $OPT->{maf_tail_filt}; #threshold to trim the edges of the models
    my $bins    = $OPT->{k_bins};
    my $m_aln   = $OPT->{m_aln};
    my $e_rate  = $OPT->{e_rate};

    my @expect   = map {0} (0..$bins-1);
    my $label = "$mod:thr:$thr:e_rate:$e_rate:cfrac:$cfrac";
    $label .= ":is_ref" if($OPT->{is_ref}); #reference includes extremes
    foreach my $mhash (values %{$maf_set}){
	keys %{$mhash};
	while(my ($cov, $num) = each %{$mhash}){
	    #build archive if it does not exist
	    unless($ARC{$label}{$bins}{$cov}){
		my $arc = $ARC{$label}{$bins}{$cov} = [map {0} (0..$bins-1)];
		foreach(my $i = 0; $i <= $cov; $i++){
		    my $maf = $i/$cov;
		    next if($maf < $thr);

		    my $L = _maf_point($i, $cov, $mod, $cfrac, 0, $OPT);

		    _bin_it($bins,
			    $cov,
			    $L,
			    $maf,
			    $arc);
		}
	    }
	    my $arc = $ARC{$label}{$bins}{$cov};

	    #add xenograft to model
	    if($xcov){
		my $mlabel .= "$label:xeno:$xcov:$m_aln";
		unless($ARC{$mlabel}{$bins}{$cov}){
		    my $sf = ($OPT->{use_generic}) ? (1-$m_aln)**2 : 1-$m_aln;
		    my $marc = $ARC{$mlabel}{$bins}{$cov} = [map {$sf * $_} @{$arc}];
		   
		    foreach(my $i = 0; $i <= $cov; $i++){
			my $maf = $i/$cov;
			next if($maf < $thr);

			my $L = _maf_point($i, $cov, $mod, $cfrac, $xcov, $OPT);
			$L *= 1 - $sf; #adjust to mouse percent of total

			_bin_it($bins,
				$cov,
				$L,
				$maf,
				$marc);
		    }
		}
		$arc = $ARC{$mlabel}{$bins}{$cov};
	    }

	    for (my $i = 0; $i < @{$arc}; $i++){
		$expect[$i] += $arc->[$i] * $num;
	    }
	}
    }

    #normalize
    my $count = 0;
    $count += $_ foreach (@expect);
    if($count){
	$_ /= $count foreach (@expect);
    }

    return \@expect;
}
}

#returns the binomial liklihood of a given point of the model
{
my %CUT; #stored cutoff thresholds
sub _maf_point {
    my $k     = shift;
    my $cov   = shift;
    my $mod   = shift;
    my $cfrac = shift;
    my $xcov  = shift;
    my $OPT   = shift;

    my $e_rate   = $OPT->{e_rate};
    my $m_aln    = $OPT->{m_aln};
    my $basement = $OPT->{basement};
    my $rfrac    = 1 - $cfrac;
    my $thr      = $OPT->{maf_tail_filt};

    if($mod eq '!'){
        #I really don't have a model yet
        return 1/($cov+1);
    }

    #get model target frequencies
    my $tot = 0;
    my @p = grep {$_ ne 'xeno'} split(/:/, $mod);

    #create mixed tumor sample targets
    if($rfrac){
        $_ = $rfrac + ($cfrac)*$_ foreach(@p);
    }

    #normalize targets
    $tot += $_ foreach (@p);
    return 1/($cov+1) if(!$tot);
    @p = (0, @p, $tot) if($OPT->{use_generic}); #reference removes non-hets
    @p = map {$_/$tot} @p;

    #create weights
    my %wi;
    if($OPT->{use_generic}){ #expect error derived SNVs
        %wi = ('0' => 999/1000);
        for(my $i = 1; $i < @p; $i++){
            $wi{$p[$i]} += 1/1000 * 1/(@p-1);
        }
    }
    else{ #simpler model (there should only be two peaks)
        $wi{$_} += 1/@p foreach @p;
    }

    #create xenograft target
    if($xcov){
        my $adj = ($cov > $xcov) ? $xcov/$cov : 1;

        #shift these to account for mouse shift
	foreach my $t (@p){
            $t = $t*(1-$adj);
        }

        #adjust weights (not true value, just a normlized estimate)
        undef %wi;
        $wi{$_} += 1/@p foreach @p;
        if($OPT->{use_generic}){ #mouse peak not there when reference filtered
            my $wm = $m_aln/(1 - (1-$m_aln)**2);
            $wi{$_} += 1/@p * (1-$wm) foreach(@p);

            #mouse derived SNV peak
            my $xt = $adj;
            $wi{$xt} += $wm;
            push(@p, $xt);
        }
    }

    #now get likilihoods from binomial
    my %L = (sum => 0);
    for(my $i = 0; $i < @p; $i++){
        next if($L{$p[$i]}++);
        my $t = $p[$i];
        $t += ((1-$p[$i])*$e_rate) - ($p[$i]*$e_rate); #error adjusted target                                            

        #liklihood around observed MAF
        $L{sum} += binomial_pmf($k, $cov, $t)*$wi{$p[$i]};
    }

    #normalization factor for threshold
    my $x = int($cov*$thr);
    $x-- if($x/$cov >= $thr);
    my $pstring = join(':', @p);
    if($thr && !$CUT{$thr}{$pstring}{$cov}){
        my %C = (sum => 0);
        for(my $i = 0; $i < @p; $i++){
            next if($C{$p[$i]}++);
            my $t = $p[$i];
            $t += ((1-$p[$i])*$e_rate) - ($p[$i]*$e_rate); #error adjusted target
            #liklihood around observed MAF
            $C{sum} += binomial_cmf($x, $cov, $t)*$wi{$p[$i]};
        }
        $CUT{$thr}{$pstring}{$cov} = $C{sum};
	
    }
    my $cut = $CUT{$thr}{$pstring}{$cov} || 0;

    #normalize for threshold and basement
    my $L = $L{sum}/(1-$cut);
    $L = $L*(1-$basement);
    $L += $basement/($cov-$x);

    return $L;
}
}

sub _cov_stats {
    my $cov_set = shift;
    my $OPT     = shift;

    my $cov_count = 0;
    my $cov_sum = 0;
    keys %$cov_set;
    while(my ($key, $value) = each %$cov_set){
	$cov_count += $value;
	$cov_sum += $key*$value;
    }

    #calculate coverage mean and median
    my $cov_mean = 0;
    my $cov_median = 0;
    if($cov_count){
	$cov_mean = $cov_sum/$cov_count;

	my $pos = 0;
	my $goal = $cov_count/2;
	foreach my $cov (sort {$a <=> $b} keys %$cov_set){
	    $pos += $cov_set->{$cov};
	    if($pos >= $goal){
		$cov_median = $cov;
		last;
	    }
	}
    }

    return ($cov_mean, $cov_median, $cov_count);
}

sub ref_param{
    my $OPT = shift;

    confess "ERROR: Can't call ref_param without ref" if(!$OPT->{ref});

    return $OPT->{ref} if(!$OPT->{vcf_file});

    my %param = (%$OPT, %{$OPT->{ref}});
    delete($param{ref});
    delete($param{xeno});

    $param{vcf_file}  = $OPT->{vcf_file} if(!$param{vcf_file});
    $param{bam_files} = $OPT->{bam_files} if(!@{$param{bam_files}});
    $param{bam_dir}   = $OPT->{bam_dir} if(!$param{bam_dir});
    $param{bam_list}  = $OPT->{bam_list} if(!$param{bam_list});
    $param{is_generic} = 1 if($OPT->{use_generic});
    $param{is_ref} = 1 if($OPT->{use_ref});
    $param{is_xeno} = 0;
    $param{use_generic} = 0;
    $param{use_ref} = 0;
    $param{use_xeno} = 0;

    return \%param;
}

sub xeno_param{
    my $OPT = shift;

    return $OPT->{xeno} if(!$OPT->{vcf_file});

    my %param = (%$OPT, %{$OPT->{xeno}});
    delete($param{xeno});
    delete($param{ref});

    $param{vcf_file}  = $OPT->{vcf_file} if(!$param{vcf_file});
    $param{bam_files} = $OPT->{bam_files} if(!@{$param{bam_files}});
    $param{bam_dir}   = $OPT->{bam_dir} if(!$param{bam_dir});
    $param{bam_list}  = $OPT->{bam_list} if(!$param{bam_list});
    $param{is_generic} = 0;
    $param{is_ref} = 0;
    $param{is_xeno} = 1;
    $param{use_generic} = 0;
    $param{use_ref} = 0;
    $param{use_xeno} = 0;

    return \%param;
}

sub ref_seg {
    my $res = shift;

    my @ref;
    foreach my $r (@$res){
	push(@ref, $r->{ref}) if($r->{ref});
    }

    return \@ref;
}

sub xeno_seg {
    my $res = shift;

    my @xeno;
    foreach my $r (@$res){
	push(@xeno, $r->{xeno}) if($r->{xeno});
    }

    return \@xeno;
}

sub _mark_best_alleles {
    my $r = shift;
    my $OPT = shift;

    if(ref($r) eq 'HASH' && !$r->{maf_score}){
	$r = {maf_score => $r};
    }

    #already computed so return
    return if($r->{maf_score}{best});

    my $models = $OPT->{models};
    my $thr = $OPT->{maf_tail_filt};
    my $cn_max = $OPT->{cn_max};
    my $xcov = $OPT->{xcov};
    my $cfrac = $OPT->{cfrac};
    my $rfrac = 1-$cfrac;

    #for these everything uses the fine bins
    my $best;
    my %cns;
    foreach my $m (@$models){
	next if($m eq '!');

	my $sc = $r->{maf_score}{$m};
	my $rss = $sc->{rss};
	my $L = $sc->{L};
	my $cn = get_c($m);

	#best in cn category
	$cns{$cn} ||= [$m, $sc];
	my $cmp = ($L <=> $cns{$cn}[1]{L} || $cns{$cn}[1]{rss} <=> $rss);
	$cns{$cn} = [$m, $sc] if($cmp == 1);
	 
	#best overall
	if($cn <= $cn_max && ($m ne '0:0' || $rfrac < $thr)){
	    $best ||= [$m, $sc];
	    my $cmp = ($L <=> $best->[1]{L} || $best->[1]{rss} <=> $rss);
	    $best = [$m, $sc] if($cmp == 1);
	}

	next unless($xcov);

	my $xm = "$m:xeno";
	$sc = $r->{maf_score}{$xm};
	$rss = $sc->{rss};
	$L = $sc->{L};

	#best in cn category
	$cns{$cn} ||= [$xm, $sc];
	$cmp = $L <=> $cns{$cn}[1]{L} || $cns{$cn}[1]{rss} <=> $rss;
	$cns{$cn} = [$xm, $sc] if($cmp == 1);
	 
	#best overall
	if($cn <= $cn_max && ($m ne '0:0' || $rfrac < $thr)){
	    $best ||= [$m, $sc];
	    my $cmp = $L <=> $best->[1]{L} || $best->[1]{rss} <=> $rss;
	    $best = [$xm, $sc] if($cmp == 1);
	}
    }
    $cns{all} = $best;
    $r->{maf_score}{best} = \%cns;

    return;
}

#fills in the coverage of a segment from the bam file
sub add_bam_cov{
    my $r = shift;
    my $OPT = shift;

    my $bam_files = $OPT->{bam_files};
    return $r if(!@$bam_files);

    return $r if($r->{cov_bam}); #already calculated
    return $r if($r->{is_generic}); #not calculatable

    if($r->{ref}){
	add_bam_cov($r->{ref}, ref_param($OPT));
    }

    if($r->{xeno}){
	add_bam_cov($r->{xeno}, xeno_param($OPT));
    }

    #load BAM files
    my $sid = $OPT->{sid};
    my ($chr, $start, $end) = ($r->{chr}, $r->{start}, $r->{end});
    prepare_bam(@$bam_files);
    
    my $bam = $BAMS{$sid}{$chr};
    die "ERROR: There is no bam file for $chr.\n" if(!$bam);
    my $segment = $bam->segment($chr,$start, $end);
    
    #get coverage from BAM
    my $cov_list;
    if(keys %{$bam->{_samples}} == 1){
	my ($obj) = $segment->features(-type => 'coverage');
	$cov_list = $obj->coverage();
    }
    else{
	my ($obj) = $segment->features(-type => 'coverage',
				       -filter => $bam->{_samples}{$sid});
	$cov_list = $obj->coverage();
    }
    $cov_list ||= [];

    #trim for memory optimization
    my %bam_set;
    {#force lexical scope
	my $cov; #initialize outside of foreach loop for optimization
	foreach $cov (@$cov_list){
	    $bam_set{$cov}++;
	}
	undef $cov_list;
    }

    #trim off stored coverages for non-queriable regions
    my $mask_count = $r->{mask_count};
    foreach my $cov (sort {$a <=> $b} keys %bam_set){
	my $count = $bam_set{$cov};
	if($count < $mask_count){
	    $mask_count -= $count;
	    delete($bam_set{$cov});
	}
	else{
	    $count -= $mask_count;
	    $bam_set{$cov} = $count;
	    delete($bam_set{$cov}) if($count == 0);
	    last;
	}
    }

    my ($cov_mean, $cov_median, $cov_count) = _cov_stats(\%bam_set);
    
    $r->{cov_count}  = $cov_count;
    $r->{cov_mean}   = $cov_mean;
    $r->{cov_median} = $cov_median;
    $r->{cov_bam}    = \%bam_set;

    return $r;
}

#fills in the coverage of a segment from the bam file
sub add_generic_cov{
    my $r = shift;
    my $OPT = shift;

    if($r->{ref}){
	add_generic_cov($r->{ref}, ref_param($OPT));
    }

    return $r if($r->{cov_bam}); #already calculated
    return $r if(!$r->{is_generic}); #not calculatable

    #load BAM files
    if(! keys %GEN){
	foreach my $chr (keys %CHRS){
	    my $store_file = "$FindBin::RealBin/generic_background/$chr\_generic.store";
	    my @array;
	    tie @array, 'Tie::MmapArray', $store_file, { template => 'i', nels => $CHRS{$chr}};
	    $GEN{$chr} = \@array;
	}
    }

    my ($chr, $start, $end) = ($r->{chr}, $r->{start}, $r->{end});
    
    #trim for memory optimization
    my %cov_set;
    for (my $i = $start-1; $i < $end; $i++){
	$cov_set{$GEN{$chr}[$i]}++;
    }

    #trim off stored coverages for non-queriable regions
    my $mask_count = $r->{mask_count};
    foreach my $cov (sort {$a <=> $b} keys %cov_set){
	my $count = $cov_set{$cov};
	if($count < $mask_count){
	    $mask_count -= $count;
	    delete($cov_set{$cov});
	}
	else{
	    $count -= $mask_count;
	    $cov_set{$cov} = $count;
	    delete($cov_set{$cov}) if($count == 0);
	    last;
	}
    }

    my ($cov_mean, $cov_median, $cov_count) = _cov_stats(\%cov_set);
    
    $r->{cov_count}  = $cov_count;
    $r->{cov_mean}   = $cov_mean;
    $r->{cov_median} = $cov_median;
    $r->{cov_bam}    = \%cov_set;

    return $r;
}

#returns rss values
sub stat_fit {
    my $obs = shift;
    my $exp = shift;
    my $bins = shift;

    if(!defined($bins) || $bins > @$exp || $bins > @$obs){
	$bins =  (@$exp < @$obs) ? @$exp : @$obs;
    }

    return (0, 0) if($bins < 2);
    
    #if needed find numer of bins to use
    if($bins != @$obs){
	my @obs2;
	my $sum = 0;
	for(my $i = 0; $i < @$obs; $i++){
	    my $j = int(($i * $bins)/@$obs - 1);
	    $j = 0 if($j < 0);
	    $obs2[$j] += $obs->[$i];
	}
	$obs = \@obs2;
    }
    if($bins != @$exp){
	my @exp2;
	for(my $i = 0; $i < @$exp; $i++){
	    my $j = int(($i * $bins)/@$exp - 1);
	    $j = 0 if($j < 0);
	    $exp2[$j] += $exp->[$i];
	}
	$exp = \@exp2;
    }

    my $rss = 0;
    my $n = 0; #bins usable
    for(my $i = 0; $i < @$exp; $i++){
	next if(! $exp->[$i]);
	$n++;
	my $rs = ($obs->[$i] - $exp->[$i])**2;
	$rss += $rs;
    }

    return ($rss, $n);
}

sub find_peaks {
    my $data   = shift;
    my $flank = shift || 1;

    $flank = int(@$data/2) if($flank > @$data);
    $flank = 1 if($flank < 1);

    #for (x,y) type input
    my @index;
    if(ref($data->[0]) eq 'ARRAY'){
	my @data;
	foreach my $d (@$data){
	    push(@index, $d->[0]);
	    push(@data,  $d->[1]);
	}
	$data = \@data;
    }

    #set peaks using a scanning window (absolute max and min within window)
    my %peaks;
    for (my $i = 0; $i < @$data; $i++){
	my $j = $i+$flank;
	$j = $#$data if($j > $#$data);

	my $h = $i-$flank;
	$h = 0 if($h < 0);
	
	my @min;
	my @max;
	for(my $k = $h; $k <= $j; $k++){
	    @min = ($k, $data->[$k]) if(!@min || $min[1] > $data->[$k]);
	    @max = ($k, $data->[$k]) if(!@max || $max[1] < $data->[$k]);
	}
	
	$peaks{$max[0]} =  1 if($max[0] == $i);
	$peaks{$min[0]} = -1 if($min[0] == $i);
    }
    #add edges as peaks
    if(@$data > 1){
	$peaks{0} = ($data->[0] > $data->[1]) ? 1 : -1;
	$peaks{$#$data} = ($data->[$#$data] > $data->[$#$data-1]) ? 1 : -1;
    }

    #now look for missed values between current peaks and valleys
    my @steps = sort {$a <=> $b} (keys %peaks);
    my @last;
    foreach my $j (@steps){
	my $kind = $peaks{$j};
	#search space in between if both are peaks or valleys
	if(@last && $last[1] eq $kind){
	    my $i = $last[0];
	    my @min;
	    my @max;
	    for(my $k = $i; $k <= $j; $k++){
		@min = ($k, $data->[$k]) if(!@min || $min[1] > $data->[$k]);
		@max = ($k, $data->[$k]) if(!@max || $max[1] < $data->[$k]);
	    }
	    
	    $peaks{$max[0]} =  1 if($max[0] != $i && $max[0] != $j);
	    $peaks{$min[0]} = -1 if($min[0] != $i && $min[0] != $j);
	}
	@last = ($j, $kind);
    }

    #convert for (x,y) type input
    if(@index){
	%peaks = map {$index[$_] => $peaks{$_}} keys %peaks;
    }

    return \%peaks;
}

#viualization of MAF curves for debugging
sub view {
    my $dist = shift;
    my $bins = shift;

    my @norm;
    if(ref($dist) eq 'HASH'){
	if($dist->{maf_set}){
	    $dist = maf_bin($dist->{maf_set}, $bins, \%OPT);
	}
    }

    for(my $i = 0; $i < @$dist; $i++){
	$norm[$i] += $dist->[$i]; 
    }

    if($bins && $bins < @norm){
        my @norm2;
        my $sum = 0;
        for(my $i = 0; $i < @norm; $i++){
            my $j = int(($i * $bins)/@norm - 1);
            $j = 0 if($j < 0);
            $norm2[$j] += $norm[$i];
        }
        @norm = @norm2;
    }

    print STDERR "\n";
    for(my $i = 0; $i < @norm; $i++){
	print STDERR "$i";
	print STDERR '*'x(int($norm[$i]*10*@norm));
	print STDERR "\n";
    }
}

#viualization of MAF and coverage in R
sub view2 {
    my $r = shift;
    my $OPT = shift;
    my $flank = shift;

    eval "require Statistics::R"; #load R

    if(! ref($r)){ #if scalar then it's a text based location (from segement)
        my ($chr, $start, $end) = $r =~ /^(.*)\:(\d+)\-(\d+)$/;
        $r = {chr => $chr,
              start => $start,
              end => $end};
    }

    my $sid      = $OPT->{sid};
    my $chr      = $r->{chr};
    my $start    = $r->{start};
    my $end      = $r->{end};
    my $cn       = ($r->{final}) ? $r->{final}{cn} : 'NA';
    my $length   = abs($end-$start)+1;
    my $pad = ($flank) ? $flank : int(0.5*$length);
    my $cs = ($start - $pad > 0) ? $start - $pad : 1;
    my $ce = $end + $pad;
    my ($chrnum) = $chr =~ /(\d+)$/;

    #load VCF file
    my $vcf_file  = $OPT->{vcf_file};
    prepare_vcf($vcf_file);
    my $vcf = $VCF{$vcf_file};
    $vcf->set_samples(include=>[$sid]);

    #load BAM files
    my $bam_files = $OPT->{bam_files};
    prepare_bam(@$bam_files);

    #temp (used to remove xenograft datapoints)
    my $rid = $OPT->{ref}{sid};
    $vcf->set_samples(include=>[$sid, $rid]);

    #get MAF from VCF
    my %vcf_data = (chr => [],
                    start => [],
                    end => [],
                    pos => [],
                    maf=> [],
                    cov => []);
    $vcf->open(region=> "$chr:$cs-$ce");
    while(my $v = $vcf->next_data_hash()){
        my $sAD = $v->{gtypes}{$sid}{AD};

        #temp (used to remove xenograft datapoints)
        my $rAD = $v->{gtypes}{$rid}{AD};
	my $rcov;
        my $rmaf = ($OPT->{use_ref}) ? 0 : 1;
	next unless($rAD && $rAD ne '.');
	my ($rc, $ac) = split(/,/, $rAD);
	$rcov = $rc+$ac;
	next unless($rcov >= $OPT->{maf_cov_filt});
	$rmaf = $ac/$rcov;

        if($sAD && $sAD ne '.'){
            my ($rc, $ac) = split(/,/, $sAD);
            my $cov = $rc+$ac;

            if($cov > $OPT->{maf_cov_filt}){
                my $maf = $ac/$cov;
                #$maf = 0 if(!$rmaf || $rmaf < $OPT->{maf_tail_filt}); #temp

                #create structure for R
                push(@{$vcf_data{chr}},   $chrnum);
                push(@{$vcf_data{start}}, $v->{POS});
                push(@{$vcf_data{end}},   $v->{POS});
                push(@{$vcf_data{pos}},   $v->{POS});
                push(@{$vcf_data{maf}},   $maf);
                push(@{$vcf_data{cov}},   $cov);
                push(@{$vcf_data{rcov}},   $rcov) if($rcov);
            }
        }
    }

    #temp (using reference corrected coverage instead)
    my $ave = 0;
    my $rave = 0;
    for(my $i =0; $i < @{$vcf_data{cov}}; $i++){
        my $cov = $vcf_data{cov}[$i];
        last if(! $vcf_data{rcov});
        my $rcov = $vcf_data{rcov}[$i];

        $ave += $cov/@{$vcf_data{cov}};
        $rave += $rcov/@{$vcf_data{rcov}};
    }
    for(my $i =0; $i < @{$vcf_data{cov}}; $i++){
	my $cov = $vcf_data{cov}[$i];
        last if(! $vcf_data{rcov});
        my $rcov = $vcf_data{rcov}[$i];
        $cov += $rave/$rcov;
        $vcf_data{cov}[$i] = $cov;
    }

    #get coverage from SAM
    my %bam_data = (chr => [],
                    start => [],
                    end => [],
                    pos => [],
                    cov => []);
    if(@$bam_files){
        my $bam = $BAMS{$sid}{$chr};
        my $segment = $bam->segment($chr,$cs, $ce);

        #binned coverage (easier to see on plot)
        my $cov_list;
	my $bins = $length + 2*$pad;
	$bins = 5000 if($bins > 5000);
        if(keys %{$bam->{_samples}} == 1){
            my ($obj) = $segment->features(-type => "coverage:$bins");
            $cov_list = $obj->coverage();
        }
        else{
            my ($obj) = $segment->features(-type => "coverage:$bins",
					   -filter => $bam->{_samples}{$sid});
            $cov_list = $obj->coverage();
        }
        $cov_list ||= [];

        #create structure for R
        my $step = (abs($ce-$cs)+1)/$bins;
        for(my $i = 0; $i < @$cov_list; $i++){
            my $cov = $cov_list->[$i];
            my $B = int($cs + $i*$step);
            my $E = int($cs + ($i+1)*$step) - 1;

            push(@{$bam_data{chr}},   $chrnum);
            push(@{$bam_data{start}}, $B);
            push(@{$bam_data{end}},   $E);
            push(@{$bam_data{pos}},   int(($E+$B)/2));
            push(@{$bam_data{cov}},   $cov);
        }
    }

    #start interpreter
    my $R = Statistics::R->new();

    #inject VCF data
    my ($fh, $tfile) = tempfile();
    print $fh "chr,pos,maf,cov\n";
    for(my $i = 0; $i < @{$vcf_data{pos}}; $i++){
        print $fh "$vcf_data{chr}[$i],$vcf_data{pos}[$i],$vcf_data{maf}[$i],$vcf_data{cov}[$i]\n";
    }
    close($fh);
    $R->run(qq`vcf <- read.csv(file="$tfile",head=TRUE,sep=",")`);
    unlink($tfile);

    #create image plot
    my $pdf = "$sid\_$chr\_$start\-$end\_padding\_$pad.pdf";
    $R->run(qq`pdf("$pdf",width=9,height=7.5,onefile=TRUE,paper='letter',pointsize=12)`);
    if(@$bam_files){
        #inject BAM data
        my ($fh, $tfile) = tempfile();
        print $fh "chr,pos,cov\n";
        for(my $i = 0; $i < @{$bam_data{pos}}; $i++){
            print $fh "$bam_data{chr}[$i],$bam_data{pos}[$i],$bam_data{cov}[$i]\n";
	}
        close($fh);
        $R->run(qq`bam <- read.csv(file="$tfile",head=TRUE,sep=",")`);
        unlink($tfile);

        #plot with VCF and BAM
        $R->run(qq`opar <- par(mfrow=c(2,1));
                   mainstr <- sprintf("chr\%d:\%4.2f-\%4.2f (kb),  length=\%4.2f(kb) ,   cn=\%d",$chrnum,$start/1e3,$end/1e3,($length)/1e3,$cn);
                   plot(bam\$pos,bam\$cov,xlim=c($cs,$ce),pch=20,cex=1,main=mainstr);abline(v=c($start,$end),col=2);
                   plot(vcf\$pos,vcf\$maf,xlim=c($cs,$ce),pch=20,cex=1,ylim=c(0,1),col="black");abline(v=c($start,$end),col=2);
                   par(opar);`);
    }
    else{
        #plot with just VCF                                                                                                                                          
        $R->run(qq`opar <- par(mfrow=c(2,1));
                   mainstr <- sprintf("chr\%d:\%4.2f-\%4.2f (kb),  length=\%4.2f(kb) ,   cn=\%d",$chrnum,$start/1e3,$end/1e3,($length)/1e3,$cn);
                   plot(vcf\$pos,vcf\$cov,xlim=c($cs,$ce),pch=20,cex=1,main=mainstr);abline(v=c($start,$end),col=2);
                   plot(vcf\$pos,vcf\$maf,xlim=c($cs,$ce),pch=20,cex=1,ylim=c(0,1),col="black");abline(v=c($start,$end),col=2);
                   par(opar);`);

    }
    $R->run(q`dev.off()`); #close file and write
    $R->stop(); #call explictly because it is not getting destroyed

    return $pdf;
}

#returns a list of target frequencies for a given model
sub get_f {
    my $mod = shift;
    my $cfrac = shift || 1;
    my $rfrac = 1 - $cfrac;

    return '!' if($mod eq '!');

    my $tot = 0;
    my @p = grep {$_ ne 'xeno'} split(/:/, $mod);
    map {$_ = $rfrac + $cfrac*$_; $tot += $_; $_} @p;

    @p = (0, @p, $tot);
    for (@p){$_ /= $tot if($tot);}
    
    return @p;
}

sub get_f_string {
    return join(':', get_f(@_));
}

#return copy number from allele string
sub get_c {
    my $mod = shift;
    
    return '!' if($mod eq '!');

    my $sum = 0;
    $sum += $_ foreach (grep {$_ ne 'xeno'} split(/:/, $mod));
    
    return $sum;
}

#returns the most common segment median
sub expected_segment_coverage{
    my $res = shift;
    my $OPT = shift;

    my $t_cov  = $OPT->{t_cov} || 0;
    my $rt_cov = $OPT->{ref}{t_cov} || 0;

    #get average coverage for reference
    my $rave = 0 || $rt_cov;
    if(!$rave){
	foreach my $ref (@{ref_seg($res)}){
	    $rave += ($OPT->{exome}) ? $ref->{cov_mean} : $ref->{cov_median};
	}
	$rave /= @$res;
    }

    #make preliminary coverage correction based on average median
    if($rave){
	foreach my $r (@$res){
	    my $ref = $r->{ref};
	    my $obs = ($OPT->{exome}) ? $ref->{cov_mean} : $ref->{cov_median};
	    my $exp = $rave;
	    $r->{ccor_f} = $ref->{ccor_f} = $obs/$exp;
	}
    }

    #get peak for sample using kernel
    my $peak = $t_cov || 0;
    if(!$peak){
	#put data into kernel 
	my $save = 0;
	my $ker = Statistics::KernelEstimation->new();
	foreach my $r (@$res){
	    my $t = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
	    
	    #remove xenograft effect
	    if($r->{xeno} && $OPT->{xeno}{ccor_f}){
		#at very low coverage mean is more informative
		$t -= $r->{xeno}{cov_mean}/$OPT->{xeno}{ccor_f};
	    }
	    
	    $t /= $r->{ref}{ccor_f}; #correction improves resolution
	    
	    $ker->add_data($t);
	    $save += $t;
	}
	$save /= @$res;
    
	#adjust parameters for searching kernel
	my $w = ($ker->default_bandwidth() < 2) ? 1 : $ker->default_bandwidth();
	my ($min, $max) = $ker->extended_range();
	$min = 0 if($min < 0);
	$max = $save*4 if($save*4 < $max);
	my $step = ($max-$min)/100;
	
	#identify kernel peak
	print STDERR "\n";
	my $sum = 0;
	$peak = [0,0]; #expected segment median coverage
	for(my $x = $min; $x <= $max; $x += $step) {
	    my $y = $ker->pdf($x, $w);
	    $sum += $step * $y; #lets me know I've already processed most data
	    print STDERR "$x\t$y\n"; #temp
	    $peak = [$x, $y] if($y >= $peak->[1]);
	    last if($x == $max || $sum >= 0.9);
	}
	$peak = $peak->[0];
    }

    if(!$t_cov || !$rt_cov){
	#select segments matching peak
	my @select;
	my $pad = 8;
	while(!@select && @$res){
	    my ($min, $max) = ($peak-$peak/$pad, $peak+$peak/$pad);
	    foreach my $r (@$res){
		my $t = $r->{cov_median};
		if($r->{xeno} && $OPT->{xeno}{ccor_f}){
		    #at very low coverage mean is more informative
		    $t -= $r->{xeno}{cov_mean}/$OPT->{xeno}{ccor_f};
		}
		$t /= $r->{ref}{ccor_f}; #correction improves resolution
		push(@select, $r) if($min <= $t && $t <= $max);
	    }
	    $pad /= 2;
	}

	#make target for sample and reference be average from same regions
	if(!$rt_cov){	    
	    my $rt_sum = 0;
	    foreach my $r (@select){
		if(my $ref = $r->{ref}){
		    $rt_sum += ($OPT->{exome}) ? $ref->{cov_mean} : $ref->{cov_median};
		}
	    }
	    $rt_cov = $rt_sum/@select;
	}

	if(!$t_cov){	    
	    my $t_sum = 0;
	    foreach my $r (@select){
		my $t = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
		if($r->{xeno} && $OPT->{xeno}{ccor_f}){
		    #at very low coverage mean is more informative
		    $t -= $r->{xeno}{cov_mean}/$OPT->{xeno}{ccor_f};
		}
		$t_sum += $t;
	    }
	    $t_cov = $t_sum/@select;
	}
    }
	
    #add_final coverage correction factors to each segment
    if($rt_cov){
	foreach my $r (@$res){
	    my $ref = $r->{ref};
	    my $obs = $ref->{cov_median};
	    my $exp = $rt_cov;
	    $r->{ccor_f} = $ref->{ccor_f} = $obs/$exp;
	}
    }

    #return results
    return ($t_cov, $rt_cov);
}

#for sorting chromosomes
sub _chrom {
    my $id = shift;

    if($id =~ /^chr(\d+)$/){
	return $1;
    }
    elsif($id eq 'chrX'){
	return 230;
    }
    elsif($id eq 'chrY'){
	return 240;
    }
    elsif($id eq 'chrM'){
	return 500;
    }
    else{
	return 1000;
    }
}

sub add_maf_discovery_stats {
    my $res    = shift;
    my $OPT    = shift;

    mkdir($OPT->{outdir}."/disc");
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my $disc = $DB{disc} || {};
    untie(%DB);

    if(!$disc->{_hasmaf}){
	undef @$res; #don't need these will load from file
	my @work : shared;
	my @files = grep {!/^_/} keys %$disc; 
	push(@work, @files);

	prefork(); #prepare for forking
	my @threads;
	for(my $i = 1; $i < $OPT->{cpus}; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_add_maf_discovery_stats_thread,
					   \@work,
					   $OPT));
	}
	_add_maf_discovery_stats_thread(\@work, $OPT); #jump into the mix
	$_->join foreach (@threads); #gather thread results

	#update load status
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	$disc->{_hasmaf} = 1;
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	$DB{disc} = $disc;
	untie(%DB);
	$lock->unlock;

	#collect merged results from files
	while(my $file = shift @files){
	    my $set = Storable::retrieve($file);
	    push(@$res, @$set);
	}
    }

    return $res;
}

sub _add_maf_discovery_stats_thread {
    my $work = shift;
    my $OPT = shift;

    #add maf_stats for coverage fit
    my @files;
    while(my $file = shift @$work){
        my $res = Storable::retrieve($file);
	foreach my $r (@$res){
	    delete($r->{maf_score}); #just in case
	    add_maf_stats($r, {%$OPT, maf_tail_filt => 0.5});
	    if(my $ref = $r->{ref}){ #add coverage correction if available
		my $param = ref_param($OPT);
		my $obs = $ref->{cov_median};
		my $exp = $param->{t_cov};
		$r->{ccor_f} = $ref->{ccor_f} = $obs/$exp;
	    }
	}
	Storable::store($res, $file);
	push(@files, $file);
    }

    return @files;
}

#get expect by matching the best MAF for segments with that median coverage
sub base_coverage {
    my $res    = shift;
    my $OPT    = shift;

    #get expected coverage
    if(!$OPT->{t_cov} || ! $OPT->{ref}{t_cov}){
	($OPT->{t_cov}, $OPT->{ref}{t_cov}) = expected_segment_coverage($res, $OPT);
    }

    my $param = ref_param($OPT);
    my $rbase = ($param->{base_cov}) ?
	$param->{base_cov} :
	_coverage_fit(ref_seg($res), $param->{t_cov}, $param);
    my $sbase = ($OPT->{base_cov}) ?
	$OPT->{base_cov} :
	_coverage_fit($res, $OPT->{t_cov}, $OPT);

    return ($sbase, $rbase);
}

#returns the average ploidy from all results
sub get_ploidy {
    my $res = shift;
    my $OPT = shift;

    return (2, 2) if($res->[0]->{is_generic});

    my $sum = 0;
    my $cn_sum = 0;
    my @cns;
    foreach my $r (@$res){
	assign_cn($r, $OPT);
        next if($r->{final}{cn} eq '!');
        next if($r->{final}{cn} eq 'N');
        next if($r->{final}{cn} eq '.');
        next if($r->{final}{cn} eq 0 && ($OPT->{is_ref} || $OPT->{is_generic}));

        $sum += $r->{q_length};
        $cn_sum += $r->{q_length} * $r->{final}{cn};
        $cns[$r->{final}{cn}] += $r->{q_length};
    }

    my $p_ave = ($sum) ? $cn_sum/$sum : 0;

    my $p_mode;
    for(my $i = 1; $i < @cns; $i++){
        next if(!$cns[$i]);
        $p_mode = [$i, $cns[$i]] if(!$p_mode || $p_mode->[1] < $cns[$i]);
    }

    return ($p_ave, $OPT->{ploidy}) if($OPT->{ploidy});
    return ($p_ave, $p_mode->[0]);
}

#labels segment results as being somatic or not (requires segment with germline reference)
sub label_somatic {
    my $r = shift;
    my $OPT = shift;

    $r->{is_somatic} = 0; #assume everything is germline
    my $cn = $r->{final}{cn};
    my $rcn = $r->{ref}{cn};
    
    return if($cn eq '!');
    return if($cn eq 'N' || $rcn eq 'N');
    return if($cn eq '.');
    return if($rcn eq 0);
    
    #ploidy is only meaningful if it is calculable
    if($OPT->{base_cov} > 5 && $OPT->{ref}{base_cov}){
	return if($OPT->{ploidy} - $cn == $OPT->{ref}{ploidy} - $rcn); 
	return if($cn == $rcn);
	return if($cn == $r->{ucor}{cn});
	return if($cn == $OPT->{ploidy});
	return if($cn == $OPT->{ref}{ploidy});
	my @v = sort {$a <=> $b} ($r->{ucor}{cn}, $OPT->{ploidy}, $r->{ccor}{cn});
	return if($v[1] == $OPT->{ploidy});
	@v = sort {$a <=> $b} ($r->{ucor}{cn}, $OPT->{ref}{ploidy}, $r->{ccor}{cn});
	return if($v[1] == $OPT->{ref}{ploidy});
    }
    # if just labeling loss and gain
    elsif($OPT->{base_cov} <= 5 && $r->{ref} && sig_somatic($r->{ref}, ref_param($OPT))){
	return if($r->{final}{cn} < $OPT->{ploidy} && $r->{ref}{cn} < $OPT->{ref}{ploidy});
	return if($r->{final}{cn} > $OPT->{ploidy} && $r->{ref}{cn} > $OPT->{ref}{ploidy});
    }
    
    #must be significantly different from ploidy levels
    return if(!sig_somatic($r, $OPT));

    #only label segements that don't have alternate explanations as somatic
    $r->{is_somatic} = 1;
}

#significant difference between expected ploidy levels and assigned copy number
sub sig_somatic {
    my $r   = shift;
    my $OPT = shift;

    my $cfrac = $OPT->{cfrac} || 1;
    my $rfrac = 1 - $cfrac;
    my $bc  = $OPT->{chr_expects}{$r->{chr}} || $OPT->{chr_expects}{_ALL};
    my $rp = ($OPT->{ref} && $OPT->{ref}{ploidy}) ? $OPT->{ref}{ploidy} : 2;
    my $sp = $OPT->{ploidy} || 2;

    my $spc = $sp*$bc+($rp*$bc*$rfrac/$cfrac);
    my $rpc = $rp*$bc+($rp*$bc*$rfrac/$cfrac);

    #calculate if outside length adusted standard deviation of ploidies
    my $med = $r->{cov_median};
    my $n = 1 + ($r->{q_length}-1)/$OPT->{lfrag};
    $n = $r->{cov_count} if($r->{cov_count} < $n);
    $n ||= 1;
    my $sd = sqrt($med/$n) * 3.890592 + 0.5; #e-4 threshold
    return 0 if(abs($med-$rpc) <=  $sd || abs($med-$spc) <=  $sd);

    if($r->{ref} && $r->{ref}{ccor_f}){
	return 0 if($r->{ref}{ccor_f} <= 0.5 || $r->{ref}{ccor_f} >= 2);

	#see if different from reference (ref might not be at ploidy)
	my $r_med = $r->{ref}{cov_median};
	my $r_bc  = $OPT->{ref}{chr_expects}{$r->{chr}} || $OPT->{ref}{chr_expects}{_ALL};
	my $r_n = 1 + ($r->{q_length}-1)/$OPT->{lfrag};
	$r_n = $r->{ref}{cov_count} if($r->{ref}{cov_count} < $r_n);
	$r_n ||= 1;
	my $r_sd = sqrt($r_med/$r_n) * 3.890592 + 0.5;
	return 0 if(abs($med/$bc-$r_med/$r_bc) <=  $sd/$bc ||
		    abs($med/$bc-$r_med/$r_bc) <=  $r_sd/$r_bc);
	return 0 if(abs($med/($sp*$bc)-$r_med/($rp*$r_bc)) <=  $sd/($bc*$sp) ||
		    abs($med/($sp*$bc)-$r_med/($rp*$r_bc)) <=  $r_sd/($r_bc*$rp)); #ploidy adjusted

	#adjust ploidy coverage for coverage correction then check again
	$med /= $r->{ccor_f};
	$sd /= $r->{ccor_f};
	return 0 if(abs($med-$rpc) <=  $sd || abs($med-$spc) <=  $sd);
    }

    return 1;
}

#identify and then merge similar segments in result set
sub merge_segments {
    my $files = shift;
    my $OPT = shift;

    mkdir($OPT->{outdir}."/merge");
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my $load  = $DB{load};
    my $merge = $DB{merge} || {};
    untie(%DB);

    #check files
    foreach my $f (@$files){
	my ($dir, $name) = $f =~ /(.*\/)([^\/]+)$/;
	my $new = $OPT->{outdir}."/merge/$name";
	if(! -f $new){
	    File::Copy::copy($f, $new) or die "Copy failed with: $!";
	    $merge->{$new} = $load->{$f};
	}
	$f = $new;
    }

    #update merge status 
    my $lock;
    $lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
    tie(%DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    $DB{merge} = $merge;
    untie(%DB);
    $lock->unlock;

    my @work : shared;
    push(@work, @$files);
    prefork(); #prepare for forking
    my @threads;
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'list'},
				       \&_merge_thread,
				       \@work,
				       $OPT));
    }
    _merge_thread(\@work, $OPT); #jump into the mix
    $_->join foreach(@threads); #gather thread results

    return $files;
}

#identify and then merge similar segments in result set
sub smooth_merge_segments {
    my $files = shift;
    my $OPT   = shift;

    mkdir($OPT->{outdir}."/smooth");
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my $merge = $DB{merge};
    my $smooth = $DB{smooth} || {};
    untie(%DB);

    #check files
    foreach my $f (@$files){
	my ($dir, $name) = $f =~ /(.*\/)([^\/]+)$/;
	my $new = $OPT->{outdir}."/smooth/$name";
	if(! -f $new){
	    File::Copy::copy($f, $new) or die "Copy failed with: $!";
	    $smooth->{$new} = $merge->{$f};
	}
	$f = $new;
    }

    #update merge status 
    my $lock;
    $lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
    tie(%DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    $DB{smooth} = $smooth;
    untie(%DB);
    $lock->unlock;

    my @work : shared;
    push(@work, @$files);
    prefork(); #prepare for forking
    my @threads;
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'list'},
				       \&_smooth_merge_thread,
				       \@work,
				       $OPT));
    }
    _smooth_merge_thread(\@work, $OPT); #jump into the mix
    $_->join foreach(@threads); #gather thread results

    return $files;
}

sub _split_on_gaps {
    my $res = shift;
    my $OPT = shift;

    my %sets;
    foreach my $r (@$res){
	push(@{$sets{$r->{chr}}}, $r);
    }

    #turn gaps into list of splitting positions
    foreach my $chr (keys %{$OPT->{GAPS}}){
	next unless($sets{$chr});
	my $list = $OPT->{GAPS}{$chr};
	my @mid;
	foreach my $g (@$list){
	    my $s = $g->[0] + int(($g->[1]-$g->[0])/2);
	    push(@mid, $s);
	}
	@mid = sort {$a <=> $b} @mid;
	next if(!@mid);

	#filter the splits for crossing results
	my $c = shift @mid;
	my @keepers;
	@{$sets{$chr}} = sort {$a->{start} <=> $b->{start}} @{$sets{$chr}};
	foreach my $r (@{$sets{$chr}}) {
	    while($c && $c < $r->{start}){
		push(@keepers, $c);
		$c = shift @mid;
	    }
	    last if(!$c);

	    #result crosses split so I can't use it
	    if($r->{start} <= $c && $c <= $r->{end} && $r->{mask_count}/$r->{length} <= 3/4){
		$c = shift @mid;
	    }
	}
	push(@keepers, $c) if($c);
	next if(!@keepers);

	my $cres = $sets{$chr};
	delete($sets{$chr});

	for(my $i = 0; $i < @keepers; $i++){
	    my $s = $keepers[$i];
	    while(my $r = shift @$cres){
		if($r->{end} < $s){
		    push(@{$sets{"$chr:$i"}}, $r);
		}
		else{
		    unshift(@$cres, $r);
		    last;
		}
	    }
	}
	$sets{"$chr:".@keepers} = $cres; #results after last split
    }

    return \%sets;
}

#merge action performed by each thread
#access to results should be copy on write (efficient)
sub _merge_thread {
    my $work = shift; #reference to shared list
    my $OPT   = shift;

    my @files;
    while(my $file = shift @$work){
	my $nostore = ref($file);
	my $res = ($nostore) ? $file : Storable::retrieve($file);
	@$res = sort {$a->{start} <=> $b->{start}} @$res;

	my $total = @$res;
	while(1){
	    my $start = @$res;
	    while(1){
		my $start1 = @$res;

		#merge across stdev
		$res = _merge_stdev($res, $OPT, 'ucor');
		next if($start1 != @$res);
    
		#merge across corrected stdev
		$res = _merge_stdev($res, $OPT, 'ccor');
                next if($start1 != @$res);
		
		#merge based on same cn
		$res = _merge_same_cn($res, $OPT);
                next if($start1 != @$res);
		
		#merge based on very similar neighbor
		#$res = _merge_similar($res, $OPT, 'ucor');
                #next if($start1 != @$res);
		
		#merge based on very similar neighbor
		#$res = _merge_similar($res, $OPT, 'ccor');
                #next if($start1 != @$res);
		
		#merge based on same maf
		#$res = _merge_same_maf($res, $OPT);
		last if($start1 == @$res);
	    }
	    last if($start == @$res);
	}

	if($nostore){
	    push(@files, @$res);
	}
	else{
	    Storable::store($res, $file) unless(@$res == $total);
	    push(@files, $file);
	}
    }

    return @files;
}

#reassigns CN (really just for debugging)
sub _reassign_cn {
    my $work = shift; #reference to shared list
    my $OPT   = shift;
    
    foreach my $file (@$work){
        my $nostore = ref($file);
        my $res = ($nostore) ? $file : Storable::retrieve($file);
        @$res = sort {$a->{start} <=> $b->{start}} @$res;
	assign_cn($_, $OPT) foreach(@$res);
	Storable::store($res, $file);
    }
}

#merge action performed by each thread
#access to results should be copy on write (efficient)
sub _smooth_merge_thread {
    my $work = shift; #reference to shared list
    my $OPT   = shift;

    my @files;
    while(my $file = shift @$work){
	my $nostore = ref($file);
        my $res = ($nostore) ? $file : Storable::retrieve($file);
        @$res = sort {$a->{start} <=> $b->{start}} @$res;
	my $total = @$res;

	while(1){
	    my $start = @$res;

	    #==remove all odd segments
	    my $hold = [];
	    my $to_merge = [];
	    for(my $i = 0; $i < @$res; $i++){
		my $ri  = $res->[$i];
		if($ri->{length} >= 100000){
		    push(@$to_merge, $ri);
		}
		elsif($ri->{final}{cn} eq '.' || $ri->{final}{cn} eq 'N' || $ri->{ref}{ccor_f} <= 0.5 || $ri->{ref}{ccor_f} >= 2){
		    push(@$hold, $ri);
		}
		else{
		    push(@$to_merge, $ri);
		}
	    }
	    undef $res;

	    while(1){
		my $start2 = @$to_merge;
		
		#do standard merge procedure
		@$to_merge = _merge_thread([$to_merge], $OPT);
		next if($start2 != @$to_merge);
		
		#merge systematic rises and falls
		#$to_merge = _merge_systematic_shifts($to_merge, $OPT);
		last if($start2 == @$to_merge);
	    }

	    #now see which removed segments should be added back
	    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
	    @$hold = sort {$a->{start} <=> $b->{start}} @$hold; 
	    my $pos = 0;
	    foreach my $h (@$hold){
		my $bad;
		for(my $i = $pos; $i < @$to_merge; $i++){
		    my $k = $to_merge->[$i];
		    if($k->{end} < $h->{start}){
			$pos = $i;
			next;
		    }
		    
		    if($k->{start} <= $h->{end}){
			#merge in if there is some meaningful data
			if($h->{final}{cn} eq 'N' || ($h->{ref}{cn} ne '!' && $h->{ref}{cn} >= 1)){
			    $to_merge->[$i] = _merge($k, $h, $OPT);
			}
			$bad++;
		    }
		    last if($k->{start} > $h->{end}); #doesn't overlap;
		}
		push(@$res, $h) unless($bad);
	    }
	    @$res = @$to_merge;
	    $to_merge = [];
	    $hold = [];

	    #polish maf data for all those that need it
	    foreach my $r (@$res){
		$r = _polish_merge($r, $OPT) if($r->{_merged});
	    }

	    #==remove segments that appear to be lack of information rather than true CN
	    $hold = [];
	    $to_merge = [];
	    for(my $i = 0; $i < @$res; $i++){
		my $ri  = $res->[$i];
		if($ri->{length} >= 100000){
		    push(@$to_merge, $ri);
		}
		elsif($ri->{final}{cn} eq '.' || $ri->{final}{cn} eq 'N' || $ri->{ref}{ccor_f} <= 0.5){
		    push(@$hold, $ri);
		}
		else{
		    push(@$to_merge, $ri);
		}
	    }
	    undef $res;

	    while(1){
		my $start2 = @$to_merge;
		
		#do standard merge procedure
		@$to_merge = _merge_thread([$to_merge], $OPT);
		next if($start2 != @$to_merge);
		
		#merge systematic rises and falls
		#$to_merge = _merge_systematic_shifts($to_merge, $OPT);
		last if($start2 == @$to_merge);
	    }

	    #now see which removed segments should be added back
	    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
	    @$hold = sort {$a->{start} <=> $b->{start}} @$hold; 
	    $pos = 0;
	    foreach my $h (@$hold){
		my $bad;
		for(my $i = $pos; $i < @$to_merge; $i++){
		    my $k = $to_merge->[$i];
		    if($k->{end} < $h->{start}){
			$pos = $i;
			next;
		    }
		    
		    if($k->{start} <= $h->{end}){
			#merge in if there is some meaningful data
			if($h->{final}{cn} eq 'N' || ($h->{ref}{cn} ne '!' && $h->{ref}{cn} >= 1)){
			    $to_merge->[$i] = _merge($k, $h, $OPT);
			}
			$bad++;
		    }
		    last if($k->{start} > $h->{end}); #doesn't overlap;
		}
		push(@$res, $h) unless($bad);
	    }
	    @$res = @$to_merge;
	    $to_merge = [];
	    $hold = [];

	    #polish maf data for all those that need it
	    foreach my $r (@$res){
		$r = _polish_merge($r, $OPT) if($r->{_merged});
	    }

	    #==remove just the extreme correction bias segments
	    for(my $i = 0; $i < @$res; $i++){
		my $ri  = $res->[$i];
		if($ri->{length} >= 100000){
		    push(@$to_merge, $ri);
		}
		elsif($ri->{final}{cn} eq '.' || $ri->{final}{cn} eq 'N' || $ri->{ref}{ccor_f} <= 0.5){
		    push(@$hold, $ri);
		}
		else{
		    push(@$to_merge, $ri);
		}
	    }
	    undef $res;

	    while(1){
		my $start2 = @$to_merge;
		
		#do standard merge procedure
		@$to_merge = _merge_thread([$to_merge], $OPT);
		next if($start2 != @$to_merge);
		
		#merge systematic rises and falls
		#$to_merge = _merge_systematic_shifts($to_merge, $OPT);
		last if($start2 == @$to_merge);
	    }

	    #now see which removed segments should be added back
	    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
	    @$hold = sort {$a->{start} <=> $b->{start}} @$hold; 
	    $pos = 0;
	    foreach my $h (@$hold){
		my $bad;
		for(my $i = $pos; $i < @$to_merge; $i++){
		    my $k = $to_merge->[$i];
		    if($k->{end} < $h->{start}){
			$pos = $i;
			next;
		    }
		    
		    if($k->{start} <= $h->{end}){
			#merge in if there is some meaningful data
			if($h->{final}{cn} eq 'N' || ($h->{ref}{cn} ne '!' && $h->{ref}{cn} >= 1)){
			    $to_merge->[$i] = _merge($k, $h, $OPT);
			}
			$bad++;
		    }
		    last if($k->{start} > $h->{end}); #doesn't overlap;
		}	
		push(@$res, $h) unless($bad);
	    }
	    @$res = @$to_merge;
	    $to_merge = [];
	    $hold = [];

	    #polish maf data for all those that need it
	    foreach my $r (@$res){
		$r = _polish_merge($r, $OPT) if($r->{_merged});
	    }

	    #===remove just the lack of information segments
	    for(my $i = 0; $i < @$res; $i++){
		my $ri  = $res->[$i];
		if($ri->{length} >= 100000){
		    push(@$to_merge, $ri);
		}
		elsif($ri->{final}{cn} eq '.' || $ri->{final}{cn} eq 'N' || $ri->{ref}{ccor_f} <= 0.5){
		    push(@$hold, $ri);
		}
		else{
		    push(@$to_merge, $ri);
		}
	    }
	    undef $res;

	    while(1){
		my $start2 = @$to_merge;
		
		#do standard merge procedure
		@$to_merge = _merge_thread([$to_merge], $OPT);
		next if($start2 != @$to_merge);
		
		#merge systematic rises and falls
		#$to_merge = _merge_systematic_shifts($to_merge, $OPT);
		last if($start2 == @$to_merge);
	    }

	    #now see which removed segments should be added back
	    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
	    @$hold = sort {$a->{start} <=> $b->{start}} @$hold; 
	    $pos = 0;
	    foreach my $h (@$hold){
		my $bad;
		for(my $i = $pos; $i < @$to_merge; $i++){
		    my $k = $to_merge->[$i];
		    if($k->{end} < $h->{start}){
			$pos = $i;
			next;
		    }
		    
		    if($k->{start} <= $h->{end}){
			#merge in if there is some meaningful data
			if($h->{final}{cn} eq 'N' || ($h->{ref}{cn} ne '!' && $h->{ref}{cn} >= 1)){
			    $to_merge->[$i] = _merge($k, $h, $OPT);
			}
			$bad++;
		    }
		    last if($k->{start} > $h->{end}); #doesn't overlap;
		}	
		push(@$res, $h) unless($bad);
	    }
	    @$res = @$to_merge;
	    $to_merge = [];
	    $hold = [];

	    #polish maf data for all those that need it
	    foreach my $r (@$res){
		$r = _polish_merge($r, $OPT) if($r->{_merged});
	    }

	    #now do standard merge
	    while(1){
		my $start2 = @$res;
		
		#do standard merge procedure
		@$res = _merge_thread([$res], $OPT);
		last if($start2 == @$res);
	    }
	    last if($start == @$res)
	}

        if($nostore){
            push(@files, @$res);
        }
        else{
	    Storable::store($res, $file) unless(@$res == $total);
	    push(@files, $file);
        }
    }

    return @files;
}

#only merges special types used for splitting
sub _smooth_merge_special {
    my $res = shift;
    my $OPT = shift;
    my $polish = shift;

    $polish = 1 if(! defined($polish));
    @$res = sort {$a->{start} <=> $b->{start}} @$res;
    
    #remove segments that appear to be lack of information rather than true CN
    my $hold = [];
    my $to_merge = [];
    for(my $i = 0; $i < @$res; $i++){
	my $ri  = $res->[$i];
	if($ri->{length} >= 100000){
	    push(@$to_merge, $ri);
	}
	elsif($ri->{final}{cn} eq '.' || $ri->{final}{cn} eq 'N' || $ri->{ref}{ccor_f} <= 0.5){
	    push(@$hold, $ri);
	}
	else{
	    push(@$to_merge, $ri);
	}
    }
    undef @$res;

    #now run the special merge (no polish for now)
    $to_merge = _merge_special($to_merge, $OPT, 0);

    #now see which removed segments should be added back
    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
    @$hold = sort {$a->{start} <=> $b->{start}} @$hold; 
    my $pos = 0;
    foreach my $h (@$hold){
	my $bad;
	for(my $i = $pos; $i < @$to_merge; $i++){
	    my $k = $to_merge->[$i];
	    if($k->{end} < $h->{start}){
		$pos = $i;
		next;
	    }
	    
	    if($k->{start} <= $h->{end}){
		#merge in if there is some meaningful data
		if($h->{final}{cn} eq 'N' || ($h->{ref}{cn} ne '!' && $h->{ref}{cn} >= 1)){
		    $to_merge->[$i] = _merge($k, $h, $OPT);
		}
		$bad++;
	    }
	    last if($k->{start} > $h->{end}); #doesn't overlap;
	}	
	push(@$res, $h) unless($bad);
    }
    push(@$res, @$to_merge);

    #polish maf data for all those that need it
    foreach my $r (@$res){
	$r = _polish_merge($r, $OPT) if($r->{_merged} && $polish);
    }

    return $res;
}

#only merges special types used for splitting
sub _merge_special {
    my $res = shift;
    my $OPT = shift;
    my $polish = shift;

    $polish = 1 if(! defined($polish));

    #quick inital sort
    @$res = sort {$a->{chr} cmp $b->{chr} || $a->{start} <=> $b->{start}} @$res;
    
    #merge special types only first ('N' and '!' - used for splitting)
    my $act;
    my @keepers;
    for (my $i = 0; $i < @$res; $i++){
	my $r0 = $res->[$i];
	if(!$act){
            $act = $r0;
            next;
        }

	if($r0->{chr} ne $act->{chr}){
            push(@keepers, $act);
            $act = $r0;
	}
	elsif(abs($r0->{start} - $act->{end}) + 1 >= 100000){
	    push(@keepers, $act);
            $act = $r0;
	}
	elsif(!$r0->{final} || !$act->{final}){
	    if($r0->{mask_count}/$r0->{length} > 3/4 && #these are all N's
		$act->{mask_count}/$act->{length} > 3/4
	      ){
		$act = _merge($act, $r0, $OPT);
		$act->{final} = $r0->{final};
	    }
	    else{
		push(@keepers, $act);
		$act = $r0;
	    }
	}
	elsif($r0->{final}{cn} eq 'N' && $act->{final}{cn} eq 'N'){
	    $act = _merge($act, $r0, $OPT);
	    $act->{final} = $r0->{final};
	}
	elsif($r0->{final}{cn} eq '!' && $act->{final}{cn} eq '!'){
	    $act = _merge($act, $r0, $OPT);
	    $act->{final} = $r0->{final};
	}
	else{
            push(@keepers, $act);
            $act = $r0;
	}
    }
    push(@keepers, $act) if($act);

    #polish maf data for all those that need it
    foreach my $r (@keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged} && $polish);
    }

    return \@keepers;
}

#merges if cn is same and allele is same or empty
sub _merge_same_cn {
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{start} <=> $b->{start}} @$s;

    my $cfrac = $OPT->{cfrac};

    my $keepers = [];
    my $act; #active working segment
    for(my $i = 0; $i < @$s; $i++){
	my $rm = $s->[$i-1] if($i-1 >= 0); #merger minus 1 (ancestral)
	my $r0  = $s->[$i]; #this is merger

	if(!$act){
	    $act = $r0;
	    next;
	}

        die "Logic error\n" if($act->{chr} ne $r0->{chr});

	#make sure copy numbers match
	my $match = 0;
	if($act->{final}{cn} eq '.' || $r0->{final}{cn} eq '.'){ #matching '.' may be because of reference
	    $match = 0;
	}
	elsif(abs($r0->{start} - $act->{end}) + 1 >= 100000){
	    $match = 0;
	}
	elsif($act->{final}{cn} eq $r0->{final}{cn} && _alleles_match($act, $r0, $OPT)){
	    $match = 1;
	}
	elsif($act->{ucor}{cn} eq $r0->{ucor}{cn} && _alleles_match($act, $r0, $OPT)){
	    $match = 1;
	}

	if($match){
	    my $new = _merge($act, $r0, $OPT);

	    #decide when to recompute (makes merging faster)
	    my ($big, $small) =  ($r0->{q_length} >= $act->{q_length}) ? ($r0, $act) : ($act, $r0);
	    my $size = $big->{_merged} || $big->{q_length};
	    my $diff = $big->{q_length} - $size;
	    if($new->{q_length} > 30000 && $big->{q_length} <= 30000){ #maf threshold recompute
		$new = _polish_merge($new, $OPT);
	    }
	    elsif($size/10 > $diff + $small->{q_length}){ #if much bigger don't recompute yet
		$new->{_merged} = $size; #keeps track of size used to compute maf so far
		$new->{final} = $big->{final};
		$new->{ucor} = $big->{ucor};
		$new->{maf_score} = $big->{maf_score}; #approximation
		if($new->{ref}){
		    $new->{ref}{final} = $big->{ref}{final};
		    $new->{ref}{ucor} = $big->{ref}{ucor};
		    $new->{ref}{maf_score} = $big->{ref}{maf_score};
		}
	    }
	    else{
		$new = _polish_merge($new, $OPT);
	    }

	    $act = $new;
	    next;
	}
	else{
	    #no more merging so add stats and recall CN
	    push(@$keepers, $act);
	    $act = $r0;
	    next;
	}
    }
    push(@$keepers, $act) if($act);

    #polish maf data for all those that need it
    foreach my $r (@$keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return $keepers;
}

#merges if cn is same and allele is same or empty
sub _merge_same_maf {
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{start} <=> $b->{start}} @$s;

    my $cfrac = $OPT->{cfrac};
    my $cn_max = $OPT->{cn_max};

    my $keepers = [];
    my $act; #active working segment
    for(my $i = 0; $i < @$s; $i++){
	my $r0  = $s->[$i]; #this is merger

	if(!$act){
	    $act = $r0;
	    next;
	}

        die "Logic error\n" if($act->{chr} ne $r0->{chr});

	my $match;
	if(abs($r0->{start} - $act->{end}) + 1 >= 100000){
	    $match = 0
	}
	elsif($act->{q_length} < 1000000 || $r0->{q_length} < 1000000){
	    $match = 0;
	}
	elsif($act->{final}{cn} eq 'N' || $act->{final}{cn} eq '!' || $act->{final}{cn} eq '.' ||
	      $r0->{final}{cn} eq 'N' || $r0->{final}{cn} eq '!' || $r0->{final}{cn} eq '.'
	){
	    $match = 0;
	}
        elsif($act->{ref}{loh} || $r0->{ref}{loh}){
	    $match = 0;
	}
	elsif($act->{final}{model} eq '.' || $r0->{final}{model} eq '.'){ #no MAF to compare
	    $match = 0;
	}
	else{
	    my $obs_a = maf_bin($act->{maf_set}, $OPT->{k_bins}, $OPT); #bin the observed MAF values
	    my $obs_r = maf_bin($r0->{maf_set},  $OPT->{k_bins}, $OPT); #bin the observed MAF values
	    my ($rss, $k) = stat_fit($obs_a, $obs_r);
	    my $rss_a = $act->{maf_score}{$act->{final}{model}}{rss};
	    my $rss_r = $r0->{maf_score}{$r0->{final}{model}}{rss};

	    $match = 1 if($rss < $rss_a || $rss < $rss_r);
	}

	if($match){
	    my $new = _merge($act, $r0, $OPT);
	    $new = _polish_merge($new, $OPT);   
	    $act = $new;
	    next;
	}
	else{
	    #no more merging so add stats and recall CN
	    push(@$keepers, $act);
	    $act = $r0;
	    next;
	}
    }
    push(@$keepers, $act) if($act);

    #polish maf data for all those that need it
    foreach my $r (@$keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return $keepers;
}

#merges if more like neighbor than expected CN
sub _merge_similar {
    my $s = shift;
    my $OPT = shift;
    my $tag  = shift || 'ucor';

    @$s = sort {$a->{start} <=> $b->{start}} @$s;

    my $cfrac = $OPT->{cfrac};
    my $cn_max = $OPT->{cn_max};

    my $keepers = [];
    my $act; #active working segment
    for(my $i = 0; $i < @$s; $i++){
	my $r0  = $s->[$i]; #this is merger

	if(!$act){
	    $act = $r0;
	    next;
	}

        die "Logic error\n" if($act->{chr} ne $r0->{chr});

	my $match;
	if(abs($r0->{start} - $act->{end}) + 1 >= 100000){
	    $match = 0
	}
	elsif($act->{final}{cn} eq 'N' && $r0->{final}{cn} eq 'N'){
	    $match = 1;
	}
	elsif($act->{final}{cn} eq 'N' || $r0->{final}{cn} eq 'N'){
	    $match = 0;
	}
	elsif($act->{ref}{ccor_f} <= 0.5 && $r0->{ref}{ccor_f} <= 0.5 ){
	    $match = 1;
	}
	elsif($act->{ref}{ccor_f} <= 0.5 || $r0->{ref}{ccor_f} <= 0.5 ){
	    $match = 0;
	}
	elsif($act->{ref}{ccor_f} >= 2 && $r0->{ref}{ccor_f} >= 2 ){
	    $match = 1;
	}
	elsif($act->{ref}{ccor_f} >= 2 || $r0->{ref}{ccor_f} >= 2 ){
	    $match = 0;
	}
	elsif(abs($r0->{cov_mean} - $act->{cov_mean}) >= $OPT->{base_cov}){
	    $match = 0
	}
	elsif(abs($r0->{cov_median} - $act->{cov_median}) >= $OPT->{base_cov}){
	    $match = 0
	}
        elsif($act->{ref}{loh} && $r0->{ref}{loh}){
	    $match = 1;
	}
        elsif($act->{ref}{loh} || $r0->{ref}{loh}){
	    $match = 0;
	}
	elsif(_matches_neighbor($act, $r0, $OPT, $tag)){
	    $match = 1;
	}
	else{
	    $match = 0;
	}

	if($match){
	    my $new = _merge($act, $r0, $OPT);
	    $new = _polish_merge($new, $OPT);   
	    $act = $new;
	    next;
	}
	else{
	    #no more merging so add stats and recall CN
	    push(@$keepers, $act);
	    $act = $r0;
	    next;
	}
    }
    push(@$keepers, $act) if($act);

    #polish maf data for all those that need it
    foreach my $r (@$keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return $keepers;
}

sub _polish_merge {
    my $r = shift;
    my $OPT = shift;

    $r = add_maf_stats($r, $OPT);

    if(my $ref = $r->{ref}){
	my $param = ref_param($OPT);
	assign_cn($ref, $param);
	delete($ref->{_merged});
    }

    assign_cn($r, $OPT);
    label_somatic($r, $OPT) if($OPT->{use_ref});
    delete($r->{_merged});

    return $r;
}

#merges short segments within 1 SD od the ciurrent level
sub _merge_stdev {
    my $s = shift;
    my $OPT = shift;
    my $tag  = shift || 'ucor';
    my $polish = shift;

    $polish = 1 if(! defined($polish));
    return $s if(@$s <= 1);

    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};

    my @all = sort {$a->{chr} cmp $b->{chr} || $a->{start} <=> $b->{start}} @$s;
    return @all if(!@all);

    my $sfactor = (3.890592/($OPT->{base_cov} - 0.5))**2; #for getting min seg length
    my $rfactor = (3.890592/($OPT->{ref}{base_cov} - 0.5))**2; #for getting min seg length

    while(1){
	my @order = sort {$all[$a]{q_length} <=> $all[$b]{q_length}} (0..$#all);
	my @seen;
	
	my $min = $OPT->{CHRS}{$all[0]{chr}};
	my $start = @all; #how many I started with
	while(defined(my $i = shift @order)){ #works from smallest to largest
	    my $r0 = $all[$i]; #this is merger
	    $seen[$i] = 1;
	    next if(!$all[$i]);
	    next if($r0->{mask_count}/$r0->{length} > 3/4); #equivilent of cn=N
	    next if($tag eq 'ccor' && $r0->{ref}{ccor_f} <= 0.5 || $r0->{ref}{ccor_f} >= 2); #equivilent of cn=!

	    my $n_min     = $r0->{cov_median} * $sfactor;
	    my $n_min_ref = ($r0->{ref}) ? $r0->{ref}{cov_median} * $rfactor : $n_min;

	    my $h = $i-1; #merger minus 1
	    $h-- while($h > 0 && !$all[$h]); #keep moving upstream
	    my $rm = $all[$h] if($h >= 0);

	    my $j = $i+1; #merger plus 1
	    $j++ while($j < $#all && !$all[$j]); #keep moving downstream
	    my $rp = $all[$j] if($j <= $#all);

	    die "Logic error\n" if($rm && $rm->{chr} ne $r0->{chr});
	    die "Logic error\n" if($rp && $rp->{chr} ne $r0->{chr});


	    next if($r0->{final} && $r0->{final}{cn} eq '!');
	    next if($r0->{final} && $r0->{final}{cn} eq 'N');
	    next if($r0->{ref} && ($r0->{ref}{ccor_f} <= 0.5 || $r0->{ref}{ccor_f} >= 2));
	    undef $rm if($rm && $rm->{final} && $rm->{final}{cn} eq '!');
	    undef $rm if($rm && $rm->{final} && $rm->{final}{cn} eq 'N');
	    undef $rm if($rm && $rm->{ref} && ($rm->{ref}{ccor_f} <= 0.5 || $rm->{ref}{ccor_f} >= 2));
	    undef $rm if($rm && abs($r0->{start} - $rm->{end}) >= 100000);
	    undef $rp if($rp && $rp->{final} && $rp->{final}{cn} eq '!');
	    undef $rp if($rp && $rp->{final} && $rp->{final}{cn} eq 'N');
	    undef $rp if($rp && $rp->{ref} && ($rp->{ref}{ccor_f} <= 0.5 || $rp->{ref}{ccor_f} >= 2));
	    undef $rp if($rp && abs($rp->{start} - $r0->{end}) >= 100000);
	    next if(!$rm && ! $rp);

	    #calculate number of independent datapoints
	    my ($n, $n_ref);
	    $n = $n_ref = 1 + ($r0->{q_length}-1)/$OPT->{lfrag};
	    $n = $r0->{cov_count} if($r0->{cov_count} < $n);
	    $n_ref = $r0->{ref}{cov_count} if($r0->{ref} && $r0->{ref}{cov_count} < $n_ref);
	    $n ||= 1;
	    $n_ref ||= 1;

	    next if(!$polish && $n > $n_min && $n_ref > $n_min_ref);

	    #calculate length adusted standard deviation
	    my $med = $r0->{cov_median};
	    my $sd = 3.890592*sqrt($med/$n) + 0.5; #e-4 threshold
	    if($tag eq 'ccor'){
		my $rx_adj = ($r0->{xeno}) ? $r0->{xeno}{cov_median}/$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		my $med2 = $med - $rx_adj;
		$med2 = 0 if($med <= 0);
		$sd += 3.890592*(sqrt($med)-sqrt($med2));
		$med2 /= $r0->{ref}{ccor_f};
		$sd /= $r0->{ref}{ccor_f};
		$med = $med2;
	    }
	    
	    my $dist_m;
	    if($rm){
		my $m_med = $rm->{cov_median};
		if($tag eq 'ccor'){
		    my $mx_adj = ($rm->{xeno}) ? $rm->{xeno}{cov_median}/$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		    $m_med = $m_med - $mx_adj;
		    $m_med = 0 if($m_med <= 0);
		    $m_med /= $rm->{ref}{ccor_f};
		}
		
		$dist_m = abs($med - $m_med);
		undef $dist_m if($dist_m > $sd);
	    }
	    
	    my $dist_p;
	    if($rp){
		my $p_med = $rp->{cov_median};
		if($tag eq 'ccor'){
		    my $px_adj = ($rp->{xeno}) ? $rp->{xeno}{cov_median}/$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		    $p_med = $p_med - $px_adj;
		    $p_med = 0 if($p_med <= 0);
		    $p_med /= $rp->{ref}{ccor_f};
		}
		
		$dist_p = abs($med - $p_med);
		undef $dist_p if($dist_p > $sd);
	    }

	    if($polish){
		undef $dist_m if(defined($dist_m) && ! _alleles_match($rm, $r0, $OPT));
		undef $dist_p if(defined($dist_p) && ! _alleles_match($rp, $r0, $OPT));
	    }

	    if(defined($dist_m) && (!defined($dist_p) || $dist_m <= $dist_p)){
		my $new = _merge($r0, $rm, $OPT);

		#decide when to recompute (makes merging faster)
		my ($big, $small) =  ($r0->{q_length} >= $rm->{q_length}) ? ($r0, $rm) : ($rm, $r0);
		my $size = $big->{_merged} || $big->{q_length};
		my $diff = $big->{q_length} - $size;
		if(!$polish){
		    $new->{_merged} = $size;
		    $new->{final} = $big->{final};
		    $new->{ucor} = $big->{ucor};
		    $new->{maf_score} = $big->{maf_score};
		    if($new->{ref}){
			$new->{ref}{final} = $big->{ref}{final};
			$new->{ref}{ucor} = $big->{ref}{ucor};
			$new->{ref}{maf_score} = $big->{ref}{maf_score};
		    }
		}
		elsif($new->{q_length} > 30000 && $big->{q_length} <= 30000){ #maf threshold recompute
		    $new = _polish_merge($new, $OPT);
		}
		elsif($size/10 > $diff + $small->{q_length}){ #if much bigger don't recompute yet
		    $new->{_merged} = $size; #keeps track of size used to compute maf so far
		    $new->{final} = $big->{final};
		    $new->{ucor} = $big->{ucor};
		    $new->{maf_score} = $big->{maf_score}; #approximation
		    if($new->{ref}){
			$new->{ref}{final} = $big->{ref}{final};
			$new->{ref}{ucor} = $big->{ref}{ucor};
			$new->{ref}{maf_score} = $big->{ref}{maf_score};
		    }
		}
		else{
		    $new = _polish_merge($new, $OPT);
		}

		$all[$i] = $new;
		$all[$h] = undef;
		push(@order, $i); #end of list
		$seen[$i] = 0;
		$min = $new->{q_length} if($new->{q_length} < $min);

		$h-- while($h > 0 && !$all[$h]); #find upstream
		if($all[$h] && $seen[$h]){
		    push(@order, $h);
		    $seen[$h] = 0;
		    $min = $all[$h]->{q_length} if($all[$h]->{q_length} < $min);
		}

		$j++ while($j < $#all && !$all[$j]); #find downstream
		if($all[$j] && $seen[$j]){
		    push(@order, $j);
		    $seen[$j] = 0;
		    $min = $all[$j]->{q_length} if($all[$j]->{q_length} < $min);
		}

		#see if I need to resort the @order list
		my $prox = $order[0];
		while(!$all[$prox]){
		    shift @order; #throw away
		    last if(!@order);
		    $prox = $order[0];
		}
		if($all[$prox]->{q_length} >= $min){
		    @order = sort {$all[$a]{q_length} <=> $all[$b]{q_length}} grep {$all[$_]} @order;
		    $min = $OPT->{CHRS}{$new->{chr}}; #reset
		}
	    }
	    elsif(defined($dist_p)){
		my $new = _merge($r0, $rp, $OPT);

		#decide when to recompute (makes merging faster)
		my ($big, $small) =  ($r0->{q_length} >= $rp->{q_length}) ? ($r0, $rp) : ($rp, $r0);
		my $size = $big->{_merged} || $big->{q_length};
		my $diff = $big->{q_length} - $size;
                if(!$polish){
                    $new->{_merged} = $size;
                    $new->{final} = $big->{final};
		    $new->{ucor} = $big->{ucor};
                    $new->{maf_score} = $big->{maf_score};
		    if($new->{ref}){
			$new->{ref}{final} = $big->{ref}{final};
			$new->{ref}{ucor} = $big->{ref}{ucor};
			$new->{ref}{maf_score} = $big->{ref}{maf_score};
		    }
                }
		elsif($new->{q_length} > 30000 && $big->{q_length} <= 30000){ #maf threshold recompute
		    $new = _polish_merge($new, $OPT);
		}
		elsif($size/10 > $diff + $small->{q_length}){ #if much bigger don't recompute yet
		    $new->{_merged} = $size; #keeps track of size used to compute maf so far
		    $new->{final} = $big->{final};
		    $new->{ucor} = $big->{ucor};
		    $new->{maf_score} = $big->{maf_score}; #approximation
		    if($new->{ref}){
			$new->{ref}{final} = $big->{ref}{final};
			$new->{ref}{ucor} = $big->{ref}{ucor};
			$new->{ref}{maf_score} = $big->{ref}{maf_score};
		    }
		}
		else{
		    $new = _polish_merge($new, $OPT);
		}

		$all[$i] = $new;
		$all[$j] = undef;
		push(@order, $i); #end of list
		$seen[$i] = 0;
		$min = $new->{q_length} if($new->{q_length} < $min);

		$h-- while($h > 0 && !$all[$h]); #find upstream
		if($all[$h] && $seen[$h]){
		    push(@order, $h);
		    $seen[$h] = 0;
		    $min = $all[$h]->{q_length} if($all[$h]->{q_length} < $min);
		}

		$j++ while($j < $#all && !$all[$j]); #find downstream
		if($all[$j] && $seen[$j]){
		    push(@order, $j);
		    $seen[$j] = 0;
		    $min = $all[$j]->{q_length} if($all[$j]->{q_length} < $min);
		}

		#see if I need to resort the @order list
		$min = $new->{q_length} if($new->{q_length} < $min);
		my $prox = $order[0];
		while(!$all[$prox]){
		    shift @order; #throw away
		    last if(!@order);
		    $prox = $order[0];
		}
		if($all[$prox]->{q_length} >= $min){
		    @order = sort {$all[$a]{q_length} <=> $all[$b]{q_length}} grep {$all[$_]} @order;
		    $min = $OPT->{CHRS}{$new->{chr}}; #reset
		}
	    }
	}

	@all = grep {$_} @all; #filter out empty entries that were merged
	last if(@all == $start);
    }

    #polish anything that needs it
    foreach my $r (@all){
	next unless($r->{_merged});

	if($polish){
	    $r = _polish_merge($r, $OPT);
	}
	else{
	    delete($r->{_merged});
	}
    }

    return \@all;
}

sub _alleles_match {
    my $r1 = shift;
    my $r2 = shift;
    my $OPT = shift;

    return 1 if($r1->{is_generic});

    my $cn_max = $OPT->{cn_max};
    my $cfrac  = $OPT->{cfrac};
    my $rfrac  = 1 - $cfrac;

    my $a1 = $r1->{final}{model};
    my $a2 = $r2->{final}{model};
    _mark_best_alleles($r1, $OPT);
    _mark_best_alleles($r2, $OPT);

    my $match;
    if($r1->{final}{model} eq $r2->{final}{model}){
	$match = 1;
    }
    elsif($r1->{final}{model} eq 'N'  || $r2->{final}{model} eq 'N'){
	$match = 0;
    }
    elsif($r1->{final}{model} eq '!'  || $r2->{final}{model} eq '!'){
	$match = 0;
    }
    elsif($r1->{final}{model} eq '.'  || $r2->{final}{model} eq '.'){
	$match = 1;
    }
    elsif($r1->{q_length} < 30000 || $r2->{q_length} < 30000){
	$match = 1;
    }
    elsif($r1->{ref}{final}{loh} && $r2->{ref}{final}{loh}){ #if ref is LOH that's all that's possible
	$match = 1;
    }
    elsif($r1->{maf_count} < 15 && !$r1->{ref}{final}{loh}){ #not enough info
	$match = 1;
    }
    elsif($r2->{maf_count} < 15 && !$r2->{ref}{final}{loh}){ #not enough info
	$match = 1;
    }
    elsif($r1->{ref}{final}{loh} && $r1->{q_length} >= 1000000){ #big reference LOH
	$match = 0;
    }
    elsif($r2->{ref}{final}{loh} && $r2->{q_length} >= 1000000){ #big reference LOH
	$match = 0;
    }
    elsif($r1->{maf_score}{best}{all}[0] eq $r2->{maf_score}{best}{all}[0]){ #sure why not
	$match = 1;
    }
    elsif(!$rfrac && $r1->{final}{loh} && $r2->{final}{loh}){
	$match = 1;
    }
    elsif(!$rfrac && $r1->{final}{cn} <= 1 && $r2->{final}{cn} <= 1){
	$match = 1;
    }
    elsif($r1->{final}{cn} > $cn_max  || $r2->{final}{cn} > $cn_max){ #can't say the alleles don't match
	$match = 1;
    }
    elsif($r2->{ref}{final}{loh} || $r1->{ref}{final}{loh}){ #lack of info going back to ref (these are short)
	$match = 1;
    }
    else{ #now compare relative liklihood of alternate allele
	my $alt_a1 = $a2;
	my $alt_a2 = $a1;
	my $aic1 = $r1->{maf_score}{$a1}{AIC};
	my $aic2 = $r2->{maf_score}{$a2}{AIC};
	my $alt_aic1 = $r1->{maf_score}{$alt_a1}{AIC};
	my $alt_aic2 = $r2->{maf_score}{$alt_a2}{AIC};

	#scaled relative liklihood
	my ($aic_min1) = ($aic1 < $alt_aic1) ? $aic1 : $alt_aic1;
	my ($aic_min2) = ($aic2 < $alt_aic2) ? $aic2 : $alt_aic2;
	my $alt_rL1 = exp(($aic_min1-$alt_aic1)/2);
	my $alt_rL2 = exp(($aic_min2-$alt_aic2)/2);

	#compare relative liklihoods
	if($alt_rL1 > 0.05 || $alt_rL2 > 0.05){
	    $match = 1;
	}
	else{ #see if directly match each other
            my $obs_1 = maf_bin($r1->{maf_set}, $OPT->{k_bins}, $OPT); #bin the observed MAF values
            my $obs_2 = maf_bin($r2->{maf_set}, $OPT->{k_bins}, $OPT); #bin the observed MAF values
            my ($rss, $k) = stat_fit($obs_1, $obs_2);

            my $n1 = 1 + ($r1->{q_length}-1)/$OPT->{lfrag};
            $n1 = $r1->{maf_count} if($r1->{maf_count} < $n1);
            $n1 = 100 if($n1 > 100);
            $alt_aic1 = $n1*log($rss/$n1)+2*$k;

            my $n2 = 1 + ($r2->{q_length}-1)/$OPT->{lfrag};
            $n2 = $r2->{maf_count} if($r2->{maf_count} < $n2);
            $n2 = 100 if($n2 > 100);
            $alt_aic2 = $n2*log($rss/$n2)+2*$k;

	    ($aic_min1) = ($aic1 < $alt_aic1) ? $aic1 : $alt_aic1;
	    ($aic_min2) = ($aic2 < $alt_aic2) ? $aic2 : $alt_aic2;
	    $alt_rL1 = exp(($aic_min1-$alt_aic1)/2);
	    $alt_rL2 = exp(($aic_min2-$alt_aic2)/2);

	    $match = ($alt_rL1 > 0.05 || $alt_rL2 > 0.05) ? 1 : 0;
	}
    }

    return $match;
}

sub _merge {
    my $r1 = shift;
    my $r2 = shift;    
    my $OPT = shift;

    ($r2, $r1) = ($r1, $r2) if($r1->{start} > $r2->{start});

    #split off reference segments and merge separately
    my $ref;
    if($r1->{ref} ){ 
	my $rr1 = $r1->{ref};
	my $rr2 = $r2->{ref};
	my $param = ref_param($OPT);
	$ref = _merge($rr1, $rr2, $param);
    }

    #split off xenograft segments and merge separately
    my $xeno;
    if($r1->{xeno} ){ 
	my $rr1 = $r1->{xeno};
	my $rr2 = $r2->{xeno};
	my $param = xeno_param($OPT);
	$xeno = _merge($rr1, $rr2, $param);
    }

    my $chr = $r1->{chr};
    my $start = ($r1->{start} < $r2->{start}) ? $r1->{start} : $r2->{start};
    my $end = ($r1->{end} > $r2->{end}) ? $r1->{end} : $r2->{end};

    #get middle if mergers are separated
    #my $mid;
    #if(($r2->{start} - $r1->{end}) - 1 >= 100){ #ksseg max separation
    #	my %mparam = (%$OPT, ref => {}); #force to skip ref
    #	$mid = _get_mid($r1, $r2, \%mparam);
    #}

    #copy coverage and MAF data from mergers
    my $q_length = 0;
    my %vcf_set;
    my $cov_count = 0;
    my $cov_sum   = 0;
    my %maf_set;
    my $maf_count = 0;    
    foreach my $rr ($r1, $r2){
	next if(!$rr);

	$q_length += $rr->{q_length};

	#merge VCF coverage info
	keys %{$rr->{cov_vcf}};
	while(my ($cov, $count) = each %{$rr->{cov_vcf}}){
	    $vcf_set{$cov} += $count;
	    $cov_count += $count;
	    $cov_sum += $cov*$count;
	}
	
	#merge VCF MAF info
	keys %{$rr->{maf_set}};
	while(my ($maf, $chash) = each %{$rr->{maf_set}}){
	    keys %$chash;
	    while(my ($cov, $count) = each %$chash){
		$maf_count += $count;
		$maf_set{$maf}{$cov} += $count;
	    }
	}
    }

    #get coverage and MAF statistics
    my ($cov_mean, $cov_median) = _cov_stats(\%vcf_set);

    #make merged result
    my $res = make_seg_hash($chr, $start, $end, $OPT);
    $res->{contains} = [];
    $res->{gene} = $r1->{gene} if($r1->{gene});

    #just in case
    delete($res->{ref}) if($r1->{is_ref} || $r1->{is_generic});
    $res->{is_ref} = 1 if($r1->{is_ref});
    $res->{is_generic} = 1 if($r1->{is_generic});
    delete($res->{xeno}) if($r1->{is_xeno});
    $res->{is_xeno} = 1 if($r1->{is_xeno});

    $res->{q_length}   = $q_length;
    $res->{mask_count} = $res->{length} - $res->{q_length}; #will include ignored bases 
    $res->{cov_count}  = $cov_count;
    $res->{cov_mean}   = $cov_mean;
    $res->{cov_median} = $cov_median;
    $res->{cov_vcf}    = \%vcf_set;
    $res->{cov_set}    = $res->{cov_vcf}; #reused address
    $res->{maf_count}  = $maf_count;
    $res->{maf_set}    = \%maf_set;
    $res->{_merged} = 1;

    #decide if BAM coverage is necessary
    my $need_bam = 0;
    if($res->{q_length} < 30000 || $res->{maf_count} < 30 ||
       $res->{cov_count}/$res->{q_length} < 0.0005 || $res->{is_xeno}
        ){
	$need_bam = 1;
    }
    
    #add BAM coverage
    if($need_bam){
	#make use of existing info where possible
	add_bam_cov($r1, $OPT);
	add_bam_cov($r2, $OPT);

	if(@{$OPT->{bam_files}}){ #must actually have bam to calculate stats
	    my %bam_set;
	    $cov_count = 0;
	    $cov_sum = 0;
	    foreach my $rr ($r1, $r2){
		keys %{$rr->{cov_bam}}; #reset
		while(my ($cov, $count) = each %{$rr->{cov_bam}}){
		    $bam_set{$cov} += $count;
		    $cov_count += $count;
		    $cov_sum += $cov*$count;
		}
	    }
	    
	    ($cov_mean, $cov_median) = _cov_stats(\%bam_set);
	    $res->{cov_bam}    = \%bam_set;
	    $res->{cov_set}    = $res->{cov_bam}; #reused address
	    $res->{cov_count}  = $cov_count;
	    $res->{cov_mean}   = $cov_mean;
	    $res->{cov_median} = $cov_median;
	}
    }

    #add merged reference fragment back to segment
    $res->{ref} = $ref if($ref);
    $res->{xeno} = $xeno if($xeno);
    
    #keep track of where I came from
    foreach my $id (@{$r1->{contains}}, @{$r2->{contains}}){
	push(@{$res->{contains}}, $id)
    }

    #get new correction factor
    if($ref){
	my $param = ref_param($OPT);
	my $obs = $ref->{cov_median};
	my $exp = $param->{t_cov};
	$res->{ccor_f} = $ref->{ccor_f} = $obs/$exp;
    }

    return $res;
}


#looks for systematic rise and fall patterns in the sample and refereence
sub _merge_systematic_shifts {
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{start} <=> $b->{start}} @$s;

    my $keepers = [];
    my $act; #active working segment
    foreach my $r (@$s){
	if(!$act){
	    $act = $r;
	    next;
	}

        if($act->{final}{cn} =~ /^[\!N\.]$/){
            push(@$keepers, $act);
            $act = $r;
	    next;
        }

        if($r->{final}{cn} =~ /^[\!N\.]$/){
            push(@$keepers, $act);
            $act = $r;
	    next;
        }
	
	#make sure copy numbers match
        my $match;
        $match = 1 if($act->{final}{cn} - $r->{final}{cn} == $act->{ref}{cn} - $r->{ref}{cn});
	$match = 0 if(abs($r->{start} - $act->{end}) + 1 >= 100000);
	$match = _alleles_match($act, $r, $OPT) if($match);

	#make sure alleles match for long segments
	if($match){
	    $act = _merge($act, $r, $OPT);
	    $act = _polish_merge($act, $OPT);
	}
	else{
	    push(@$keepers, $act);
	    $act = $r;
	    next;
	}
    }
    push(@$keepers, $act) if($act);

    return $keepers;
}

sub write_gvf {
    my $files  = shift;
    my $OPT = shift;

    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my $data  = $DB{load};
    $data     = $DB{merge}  if($OPT->{merge});
    $data     = $DB{smooth} if($OPT->{smooth});
    untie(%DB);

    #make index to easily order files
    my %index;
    foreach my $f (@$files){
	my $seg = $data->{$f}[0];
	my ($chr, $start, $end) = $seg =~ /^(.*)\:(\d+)\-(\d+)$/;
	my %h = (chr   => $chr,
		 start => $start,
		 end   => $end);
	$index{$f} = \%h;
    }

    #sort by position
    @$files = sort {_chrom($index{$a}{chr}) <=> _chrom($index{$b}{chr}) ||
			$index{$a}{start} <=> $index{$b}{start}} @$files;

    my $OUT;
    if(my $outfile = $OPT->{outfile}){
	open($OUT, ">$outfile")
    }
    else{
	open($OUT, ">&STDOUT")
    }

    print $OUT "##gvf-version 1.06\n";
    print $OUT "#xenograft_coverage=".$OPT->{xcov}."\n";
    print $OUT "#copy fraction (sample/total)=".$OPT->{cfrac}."\n";
    print $OUT "#target_expect=".$OPT->{t_cov}."\n";
    print $OUT "#base_coverage=".$OPT->{base_cov}."\n";
    print $OUT "#ploidy=".$OPT->{ploidy}."\n";

    foreach my $f (@$files){
	my $res = Storable::retrieve($f);
	@$res = sort {$a->{start} <=> $b->{start}} @$res;
	print $OUT _gvf_line($_) foreach(@$res);
    }
    close($OUT);
}


#formats one line of GVF from a result
{
my $ID_COUNT = 1;
sub _gvf_line {
    my $r = shift;

    my @data;
    push(@data, $r->{chr});
    push(@data, '.');
    push(@data, 'copy_number_variation');
    push(@data, $r->{start});
    push(@data, $r->{end});
    push(@data, $r->{cov_median});
    push(@data, '.');
    push(@data, '.');

    my $att = 'ID='.$ID_COUNT++;
    #my $att = 'Name=CNA:'.$r->{chr}.':'.$r->{start}.'-'.$r->{end};

    #add ref cn
    my $add = '';
    if($r->{ref}){
	my $cn     = $r->{ref}{cn};
	my $allele = $r->{ref}{model};
	my $L      = $r->{ref}{rL};
	my $aL     = $r->{ref}{rL_maf};
	my $ccor_f = $r->{ref}{ccor_f};
	my $ratio = $r->{ref}{maf_count}/$r->{q_length} * 10000;
	$add .= ";ref_cn=$cn";
	$add .= ";ref_cn_L=".sprintf('%.2f', $L);
	$add .= ";ref_cn_allele=$allele";
	$add .= ";ref_cn_allele_L=$aL";
	$add .= ";ref_ccor_f=".sprintf('%.2f', $ccor_f);
	$add .= ";ref_maf_ratio=".sprintf('%.2f', $ratio);
    }
    if($r->{ccor}){
	my $cn     = $r->{ccor}{cn};
	my $allele = $r->{ccor}{model};
	my $L      = $r->{ccor}{rL};
	my $aL     = $r->{ccor}{rL_maf};
	$add .= ";ccor_cn=$cn";
	$add .= ";ccor_cn_L=".sprintf('%.2f', $L);
	$add .= ";ccor_cn_allele=$allele";
	$add .= ";ccor_cn_allele_L=$aL";
    }
    if($r->{ucor}){
	my $cn     = $r->{ucor}{cn};
	my $allele = $r->{ucor}{model};
	my $L      = $r->{ucor}{rL};
	my $aL     = $r->{ucor}{rL_maf};
	$add .= ";ucor_cn=$cn";
	$add .= ";ucor_cn_L=".sprintf('%.2f', $L);
	$add .= ";ucor_cn_allele=$allele";
	$add .= ";ucor_cn_allele_L=$aL";
    }

    #add last part to attributes
    if($r->{final}{cn} eq '.'){
        $data[2] = 'region_lacks_data';
        $att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} eq '!'){
        $data[2] = 'uncallable_region';
        $att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} eq 'N'){
        $data[2] = 'masked_region';
        $att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} == $OPT{ploidy} && $r->{final}{loh}){
	$data[2] = 'LOH_region';
	$att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} == $OPT{ploidy} && $r->{final}{loh}){
	$data[2] = 'LOH_region';
	$att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} eq $OPT{ploidy}){
        $data[2] = 'region';
        $att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} > $OPT{ploidy}){
	$data[2] = 'copy_number_gain';
	$att .= ';Reference_seq=~;Variant_seq=~';
    }
    elsif($r->{final}{cn} < $OPT{ploidy}){
	$data[2] = 'copy_number_loss';
	$att .= ';Reference_seq=~;Variant_seq=-';
    }
    else{
	die "ERROR: Logic error";
    }

    #add final cn
    my $f_cn     = $r->{final}{cn};
    my $f_allele = $r->{final}{model};
    my $f_L      = $r->{final}{rL};
    my $f_aL     = $r->{final}{rL_maf};
    $att .= ";Variant_copy_number=$f_cn";
    $att .= ";liklihood=".sprintf('%.2f', $f_L);
    $att .= ";allele=$f_allele";
    $att .= ";allele_liklihood=$f_aL";
    $att .= ';is_somatic='.$r->{is_somatic} if($r->{is_somatic});
    $att .= ';is_loh='.$r->{final}{loh} if($r->{final}{loh});
    $att .= $add;
    push(@data, $att);  
    return join("\t", @data)."\n"
}
}

sub p_cov {
    my $r = shift;
    my $m = shift;
    my $points = shift;

    my $x = $r;
    my $spm = ($m < 0) ? abs($m/2): $m; #special context for expect of 0
    my $rpm = 0;
    if(ref($r)){
        if(!$points){
            $points = 1 + ($r->{q_length}-1)/$OPT{lfrag};;
            $points = $r->{cov_count} if($r->{cov_count} < $points);
	    $points ||= 1;
        }
        $x = $r->{cov_median};

        $rpm = $r->{ref}{cov_median} if($r->{ref}); #aproximates reference poisson mean
    }
    $points /= 4; #temp
    $points = 1 if(!$points || $points < 1);

    $x += abs($m/2) if($m < 0); #special context for expect of 0 since most of these are really 1
    $x   *= $points; #aproximates observed within poisson
    $spm *= $points; #aproximates sample poisson mean
    $rpm *= $points;
    if($rpm){
        $x /= $rpm;
        return c_transform($x, $spm, $rpm);
    }
    else{
        return poisson_pmf($x, $spm);
    }

    #this essentially gets skipped
    my $s = sqrt(abs($m)/$points);
    $m = 0 if($m < 0); #special context for expect of 0

    #approximates with normal (poor fit but ~4x SD works empirically)
    return gaus_pdf($x, $m, 4*$s);
}

# Normal distribution PDF function
sub gaus_pdf {
    my ($x, $m, $s) = @_; #position, mean, standard deviation

    if($s == 0){
	return ($x == $m) ? 1 : 0;
    }

    my $z = ($x - $m)/$s;
    return exp(-0.5*$z*$z)/($s*sqrt(2.0*PI));
}

# Normal distribution CDF function
# Abramowitz & Stegun method
sub gaus_cdf {
    my ($x, $m, $s) = @_; #position, mean, standard deviation

    if($s == 0){
	return ($x >= $m) ? 1 : 0;
    }

    my $z = abs($x - $m)/$s;
    my $t = 1.0/(1.0 + 0.2316419*$z);
    my $y = $t*( 0.319381530
                 + $t*( -0.356563782
                        + $t*( 1.781477937
                               + $t*( -1.821255978 + $t*1.330274429 ))));
    my $pdf = exp(-0.5*$z*$z)/($s*sqrt(2.0*PI));

    if( $x >= $m ) {
        return 1.0 - $pdf*$y*$s;
    } else {
        return $pdf*$y*$s;
    }
}

sub assign_cn {
    my $r   = shift;
    my $OPT = shift;

    my $cn_max = $OPT->{cn_max};
    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};
    my $cfrac  = $OPT->{cfrac};
    my $rfrac  = 1 - $cfrac;

    my $bc = $OPT->{chr_expects}{$r->{chr}} || $OPT->{chr_expects}{_ALL};

    my $ucn = best_cn($r, 0, $bc, $OPT);
    my $cn  = ($r->{ccor_f}) ? best_cn($r, 0, $bc*$r->{ccor_f}, $OPT) : $ucn;
    my ($fin) = sort {$b->{rL} <=> $a->{rL} || $b->{L} <=> $a->{L}} ($cn, $ucn);
    $fin = $ucn if(!$r->{ref}); #temp
    if($r->{ref}){
	my $loi = {L => 0,
		   L_cov => 0,
		   L_maf => 0,
		   cn => '.',
		   e_cov => 0,
		   model => '.',
		   rL => 0,
		   rL_cov => 0,
		   rL_maf => 0,
		   score => undef};
	my $bad = {L => 0,
		   L_cov => 0,
		   L_maf => 0,
		   cn => '!',
		   e_cov => 0,
		   model => '!',
		   rL => 0,
		   rL_cov => 0,
		   rL_maf => 0,
		   score => undef};

	if($r->{ref}{ccor_f} >= 2){ #reference spikes make theses calls unreliable
	    $fin = $bad;
	}
	elsif($r->{ref}{ccor_f} <= 0.5){ #there is not enough information
	    $fin = $bad;
	}
	elsif($r->{q_length}/$r->{length} < 0.5){ #most of the sequence is unqueriable
	    $fin = $loi;
	}
	elsif($rfrac && $cn->{cn} == 0){ #check internal ref of primary for lack of coverage
	    my $rbc = $bc * $rfrac/$cfrac;
	    $rbc *= $r->{ccor_f} if($r->{ccor_f});
	    my $rcn = best_cn($r, 0, $rbc, {%$OPT, cfrac => 1, rfrac => 0});

	    #if ref internal to CNV is not diploid then everything is wrong
	    $fin = $loi if($rcn->{cn} <= 1);
	}
	else{ #is seg long enough to distinguish the ref copy number as different than 0
	    #relationship must always be true --> n' > (c_o * a^2)/(d - 0.5)^2
	    my $median = $r->{ref}{cov_median};
	    my $min_len = $OPT->{lfrag}*$median*3.890592**2/($median-0.5)**2;
	    $fin = $loi if($r->{ref}{q_length} < $min_len);
	}
    }

    #it's a masked segment
    if($r->{mask_count}/$r->{length} > 3/4){
	$fin = {L => 0,
		L_cov => 0,
		L_maf => 0,
		cn => 'N',
		e_cov => 0,
		model => 'N',
		rL => 0,
		rL_cov => 0,
		rL_maf => 0,
		score => undef};
    }

    keys  %$cn;
    while(my ($key, $value) = each %$cn){
	$r->{ccor}{$key} = $value;
    }

    keys  %$ucn;
    while(my ($key, $value) = each %$ucn){
	$r->{ucor}{$key} = $value;
    }

    keys  %$fin;
    while(my ($key, $value) = each %$fin){
	$r->{final}{$key} = $value;
    }

    #mark loss of heterozygosity
    if($OPT->{use_ref} && $r->{q_length} > 30000 && ($r->{maf_count} > 15 || $r->{q_length} > 1000000)){
	if($r->{ref}{final}{model} =~  /^0/){
	    if($r->{final}{model} =~ /^[1-9]/){ #something odd
		if($r->{ref}{ccor}{model} !~ /^0/){
		    $r->{ref}{model} = '.';
		    $r->{ref}{final}{model} = '.';
		    $r->{ref}{final}{loh} = 0;
		    $r->{ref}{loh} = 0;
		    $r->{final}{loh} = 0;
		}
		else{
		    $r->{final}{model} = '.';
		    $r->{ref}{final}{loh} = 1;
		    $r->{ref}{loh} = 1;
		    $r->{final}{loh} = 0;
		}
	    }
	    else{
		$r->{ref}{final}{loh} = 1;
		$r->{ref}{loh} = 1;
		$r->{final}{loh} = 0;
	    }
	}
	else{
            $r->{ref}{final}{loh} = 0;
            $r->{ref}{loh} = 0;
	    $r->{final}{loh} = ($r->{final}{model} =~ /^0/) ? 1 : 0;
	}
    }
    else{
	$r->{ref}{final}{loh} = 0 if($r->{ref});
	$r->{ref}{loh}        = 0 if($r->{ref});
	$r->{final}{loh}      = 0;
    }

    #shortcuts for quick access
    $r->{cn} = $r->{final}{cn};
    $r->{rL} = $r->{final}{rL};
    $r->{model} = $r->{final}{model};
    $r->{rL_maf} = $r->{final}{rL_maf};
    $r->{loh} = $r->{final}{loh};

    return;
}

#returns log liklihood ratio between two coverage models
sub logLr_cov {
    my ($k, $m1, $m2) = @_;

    if(ref($k)){
	$k = $k->{cov_median};
    }

    return "INF"+1 if($m2 == 0); #perl's infinity value
    return -1*("INF"+1) if($m1 == 0); #perl's infinity value

    #equation before log --> ($m1/$m2)**int($k) * exp($m2-$m1);
    return ($m2-$m1) + int($k)*log($m1/$m2);
}

sub estimate_copy_fraction {
    my $res = shift;
    my $OPT = shift;
    
    die unless($OPT->{use_ref});

    #--find global MAF peaks
    my @select;
    my $rker = Statistics::KernelEstimation->new();
    my $sker = Statistics::KernelEstimation->new();
    my $fker = Statistics::KernelEstimation->new();
    my $cker = Statistics::KernelEstimation->new(); #for coverage
    foreach my $r (@$res){
        #ignore if reference is not as expected
        next if($r->{ref}{maf_target} < 0.46 || $r->{ref}{maf_target} > 0.54);
        #sample is closer than ref
        next if(abs(0.5-$r->{maf_target}) < abs(0.5-$r->{ref}{maf_target}));

    	$rker->add_data($r->{ref}{maf_target});
    	$sker->add_data($r->{maf_target});

	next if($r->{maf_target} >= $r->{ref}{maf_target});
	
    	#only add different from reference
    	$fker->add_data($r->{maf_target});
	my $cov = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
    	$cker->add_data($cov/$r->{ref}{ccor_f});
	push(@select, $r);
    }
    
    #get global peaks of MAF
    my @rdis;
    my @sdis;
    my @fdis;
    my $w = 0.005; #for gausian
    my ($min, $max) = (0, 500); #0 to 0.5 when divided by 1000
    for(my $x = $min; $x <= $max; $x++) {
    	my $y = $rker->pdf($x/1000, $w);
    	push(@rdis, $y);
    	$y = $sker->pdf($x/1000, $w);
    	push(@sdis, $y);
    	$y = $fker->pdf($x/1000, $w);
    	push(@fdis, $y);
    }
    my $rpeaks = find_peaks(\@rdis, 10);
    my $speaks = find_peaks(\@sdis, 10);
    my $fpeaks = find_peaks(\@fdis, 10);

    #filter all MAF peaks
    my @rpeaks = sort {$a <=> $b} grep {$rpeaks->{$_} == 1} keys %$rpeaks;
    my @speaks = sort {$a <=> $b} grep {$speaks->{$_} == 1} keys %$speaks;
    my @fpeaks = sort {$a <=> $b} grep {$fpeaks->{$_} == 1} keys %$fpeaks;
    my ($rtop) = sort {$rdis[$b] <=> $rdis[$a]} @rpeaks;
    my ($stop) = sort {$sdis[$b] <=> $sdis[$a]} @speaks;
    @rpeaks = grep {$rdis[$_] >= 0.1*$rdis[$rtop]} @rpeaks; #ignore very small peaks
    @speaks = grep {$sdis[$_] >= 0.1*$sdis[$stop]} @fpeaks; #filtered fpeaks using sdis

    #get global peaks of coverage
    my @cdis;
    ($min, $max) = (0, $cker->range()+10);
    for(my $x = $min; $x <= $max; $x++) {
        push(@cdis, $cker->pdf($x, 1));
    }
    my $cpeaks = find_peaks(\@cdis, 10);
    
    #filter coverage peaks
    my @cpeaks = sort {$a <=> $b} grep {$cpeaks->{$_} == 1} keys %$cpeaks;
    my ($ctop) = sort {$cdis[$b] <=> $cdis[$a]} @cpeaks;
    @cpeaks = grep {$cdis[$_] >= 0.1*$cdis[$ctop]} @cpeaks; #ignore very small peaks    

    if(!@speaks){
    	@speaks = grep {$sdis[$_] >= 0.05*$sdis[$stop]} @fpeaks;
    	return (0, 0) if(!@speaks);
    }
    return (1, 0) if($speaks[0] <= 1); #peak at standard MAF extreme

    #precluster MAF to normalize or use as final clusters
    my @divs = map {$_/1000} sort {$a <=> $b} grep {$fpeaks->{$_} == -1} keys %$fpeaks;
    shift @divs; #all results are always past first cutoff extreme
    my $pos = 0;
    my @clusters;
    @select = sort {$a->{maf_target} <=> $b->{maf_target}} @select;
    foreach my $r (@select){
	$pos++ if($r->{maf_target} > $divs[$pos]);
	push(@{$clusters[$pos]{data}}, $r);
    }

    #cluster by MAF and coverage
    if($OPT->{exome} || (@speaks <= 2 && @cpeaks == 1)){ #base primarilly on MAF if exome (too much coverage variance)
	#build convenient cluster hash (center of cluster is peak's mean MAF and coverage)
	foreach my $c (@clusters){
	    my $maf_ave;
	    my $cov_ave;
	    foreach my $r (@{$c->{data}}){
		$maf_ave += $r->{maf_target};
		my $cov = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
		$cov_ave += $cov/$r->{ref}{ccor_f};
	    }
	    $maf_ave /= @{$c->{data}};
	    $cov_ave /= @{$c->{data}};

	    $c->{maf} = $maf_ave;
	    $c->{cov} = $cov_ave;
	}
    }
    else{ #use k-means clustering
	#get MAF variance
	my $maf_var;
	foreach my $c (@clusters){
	    my $mean;
	    foreach my $r (@{$c->{data}}){
		$mean += $r->{maf_target};
	    }
	    $mean /= @{$c->{data}};
	    
	    foreach my $r (@{$c->{data}}){
		$maf_var += ($r->{maf_target} - $mean)**2;
	    }
	}
	$maf_var /= @select;

	#precluster on coverage to normalize in that axis
	@divs = sort {$a <=> $b} grep {$cpeaks->{$_} == -1} keys %$cpeaks;
	shift @divs; #all results are always past first cutoff extreme
	$pos = 0;
	undef @clusters;
	if($OPT->{exome}){
	    @select = sort {$a->{cov_mean}/$a->{ref}{ccor_f} <=> $b->{cov_mean}/$b->{ref}{ccor_f}} @select;
	}
	else{
	    @select = sort {$a->{cov_median}/$a->{ref}{ccor_f} <=> $b->{cov_median}/$b->{ref}{ccor_f}} @select;
	}
	foreach my $r (@select){
	    if(!defined($divs[$pos])){ #just add to last cluster
		push(@{$clusters[$pos-1]{data}}, $r);
		next;
	    }

	    my $cov = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
	    $cov /= $r->{ref}{ccor_f};
	    
	    $pos++ if($cov > $divs[$pos]);
	    push(@{$clusters[$pos]{data}}, $r);
	}
    
	#get coverage variance
	my $cov_var;
	foreach my $c (@clusters){
	    my $mean;
	    foreach my $r (@{$c->{data}}){
		my $cov = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
		$cov /= $r->{ref}{ccor_f};
		$mean += $cov;
	    }
	    $mean /= @{$c->{data}};
	    
	    foreach my $r (@{$c->{data}}){
		my $cov = ($OPT->{exome}) ? $r->{cov_mean} : $r->{cov_median};
		$cov /= $r->{ref}{ccor_f};
		$cov_var += ($cov - $mean)**2;
	    }
	}
	$cov_var /= @select;

	#make normalized datafile 
	mkdir($OPT->{outdir}."/tmp");
	my ($tfh, $tfile) = tempfile(DIR => $OPT->{outdir}."/tmp", CLEANUP => 1);
	my %index;
	my $maf_sd = sqrt($maf_var);
	my $cov_sd = sqrt($cov_var);
	foreach my $r (@select){
	    my $id = $r->{chr}.':'.$r->{start}.'-'.$r->{end};
	    my $cov = $r->{cov_median}/($r->{ref}{ccor_f}*$cov_sd);
	    my $maf = $r->{maf_target}/$maf_sd;
	    print $tfh "$id\t$maf\t$cov\n";
	    $index{$id} = $r;
	}
	close($tfh);

	#use Algorithm::ExpectationMaximization;
	#my $em = Algorithm::ExpectationMaximization->new( datafile            => $tfile,
	#						  mask                => 'N11',
	#						  K                   => scalar(@fpeaks),
	#						  max_em_iterations   => 300,
	#						  seeding             => 'kmeans',
	#						  terminal_output     => 1,
	#						  debug               => 0,);
	#$em->read_data_from_file();
	#my ($ids, $centers) = $em->EM();

	#--do k-means clustering to identify peaks
	my $Kmin = (scalar(@fpeaks) > 1) ? scalar(@fpeaks) : scalar(@cpeaks); #at least as many peaks as the kernel
	my $kmeans = Algorithm::KMeans->new(datafile        => $tfile,
					    mask            => 'N11', #ignore coverage for exome
					    K               => 0, #discoverable
					    Kmin            => $Kmin,
					    cluster_seeding => 'smart',
					    terminal_output => 0,
					    debug           => 0);
	$kmeans->read_data_from_file();
	my ($ids, $centers) = $kmeans->kmeans();

	eval "require Graphics::GnuplotIF"; #temp
	$kmeans->visualize_clusters('11') if(!$@); #temp

	#make convenient clusters hash (correct normalization in centers)
	undef @clusters;
	for (my $i = 0; $i < @$ids; $i++){
	    my %c = (data => [@index{@{$ids->[$i]}}],
		     maf => $centers->[$i][0]*$maf_sd,
		     cov => $centers->[$i][1]*$cov_sd);
	    push(@clusters, \%c);
	}
    }
	
    #filter out non-LOH clusters
    my @loh;
    @clusters = sort {$a->{cov} <=> $b->{cov} || $a->{maf} <=> $b->{maf}} @clusters;
    while(my $l = shift @clusters){
	@clusters = grep {$_->{cov} >= $l->{cov} && $_->{maf} <= $l->{maf}} @clusters;
	push(@loh, $l);
    }
    return (0,0) if(! @loh);

    #do linear regression of outer two peaks to fit CN and celularity
    #format --> (MAF, CN)
    my @data;
    my ($low) = sort {$a->{maf} <=> $b->{maf}} grep {$_->{cov} < 2*$loh[0]{cov}} @loh;
    my @loh1 = ([$low, 1]); #CN 1
    push(@data, \@loh1);
    
    if(@loh >= 2){
	my @loh2 = ([$loh[0], 1], [$loh[1], 2]); #CN 1-2
	push(@data, \@loh2);
    }
    
    if(@loh == 2 && @fpeaks >= 5){
	my @loh2b = ([$loh[0], 2], [$loh[1], 3]);  #CN 2-3
	push(@data, \@loh2b);
    }
    
    if(@loh >= 3){
	my @loh3 = ([$loh[0], 1], [$loh[1], 2], [$loh[2], 3]); #CN 1-3
	push(@data, \@loh3);
    }
    
    foreach my $d (@data){
	#now force interecept for final slope
	my $reg = Statistics::Regression->new( "Cellularity Regression", ["Intercept", "Slope"]);
	$reg->include(1/$_->[0]{maf} - 2, [0, $_->[1]]) foreach(@$d); #-2 to force intercept of 0
	$reg->include(1/0.5 - 2, [0, 0]); #add point at 0 and 0.5
	
	my $m = $reg->theta()->[1]; #get slope for celularity
	my $cell = $m/(1+$m);
	my $rsq = $reg->rsq();
	$d = [$cell, $rsq, $d->[0][0], $d->[0][1]];
    }

    my ($cell, $rsq, $loh, $cn);
    my $one = shift(@data);
    @data = sort {$b->[1] <=> $a->[1]} grep {$_->[1] > 0.95} @data;
    if(@data){
	($cell, $rsq, $loh, $cn) = @{$data[0]};
    }
    else{
	($cell, $rsq, $loh, $cn) = @{$one};
    }

    #initial base coverage estimate based on one copy loh
    delete($loh->{data});
    $loh->{cn} = $cn;

    return wantarray ? ($cell, $rsq, $loh) : $cell;    

}

sub discovery_segments {
    my $OPT    = shift;

    die unless($OPT->{use_ref}); #temp

    #see what is already finished
    mkdir($OPT->{outdir}."/disc");
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my @files;
    my %seen;
    my $disc = $DB{disc};
    foreach my $file (keys %$disc) {
	next if($file =~ /^_/); #these are parameters

	if(!-f $file){
	    delete $disc->{$file};
	    next;
	}
	
	push(@files, $file);
	$seen{$_}++ foreach(@{$disc->{$file}});
    }
    $DB{disc} = $disc;
    untie(%DB);

    #make segments to see what still needs to be finished
    my @work :shared;
    foreach my $chr (keys %CHRS){
	next if($chr =~ /[XYM]$/);
	my $end = $OPT->{CHRS}{$chr};
	my $gaps = $OPT->{GAPS}{$chr};

	my $pos = 1;
	foreach my $g (@$gaps){
	    my ($B, $E) = @$g;
	    if($B <= $pos){
		$pos = $E+1;
		next
	    }

	    my $id = "$chr:$pos-".($B-1);
	    push(@work, $id) if(!$seen{$id});
	    $pos = $E+1;
	}
	if($pos < $end){
	    my $id = "$chr:$pos-$end";
	    push(@work, $id) if(!$seen{$id});
	}
    }

    #process segments into files
    if(@work){
	#--launch threads
	my @threads;
	prefork(); #prepare for forking
	for(my $i = 1; $i < $OPT->{cpus}; $i++){
	    push(@threads, threads->create({context => 'list'},
					   \&_discovery_segments_thread,
					   \@work,
					   $OPT));
	}
	
	#--jump into the mix
	my @ret = _discovery_segments_thread(\@work, $OPT);
	push(@files, @ret);

	#--collect thread results
	foreach my $thr (@threads){
	    my @ret = $thr->join;
	    push(@files, @ret);
	}
	$disc->{_hasmaf} = 0;
    }

    my @res;
    while(my $file = shift @files){
	my $set = Storable::retrieve($file);
        push(@res, @$set);
    }

    return \@res;
}

sub _discovery_segments_thread {
    my $work = shift;
    my $OPT    = shift;

    my $sid        = $OPT->{sid};
    my $vcf_file   = $OPT->{vcf_file};
    my $rid        = $OPT->{ref}{sid};
    my $rvcf_file  = $OPT->{ref}{vcf_file};
    my $thr        = $OPT->{maf_tail_filt};

    #load VCF file
    prepare_vcf($vcf_file);
    my $vcf = $VCF{$vcf_file};
    $vcf->set_samples(include=>[$sid]);
    
    #load reference VCF
    my $rvcf;
    if($OPT->{use_ref}){
	if($rvcf_file && $rvcf_file ne $vcf_file){
	    prepare_vcf($rvcf_file);
	    $rvcf = $VCF{$rvcf_file};
	    $rvcf->set_samples(include=>[$rid]);
	}
	else{
	    $vcf->set_samples(include=>[$sid, $rid]);
	}
    }

    #read in each region of work
    my @files;
    while (my $seg = shift @$work){
	my @sres;
	
	#calculate first for separate reference VCF
	my %ref_ok;
	if($rvcf){
	    $rvcf->open(region=> $seg);
	    while(my $v = $rvcf->next_data_hash()){
		my $pos = $v->{POS};
		my $rAD = $v->{gtypes}{$rid}{AD};
		if($rAD && $rAD ne '.'){
		    my ($rc, $ac) = split(/,/, $rAD);
		    my $cov = $rc+$ac;
		    
		    if($cov >= $OPT->{maf_cov_filt}){
			my $maf = $ac/$cov;
			if($thr <= $maf && $maf <= 1-$thr){
			    $ref_ok{$seg}{$pos} = [$maf, $cov];
			}
		    }
		}
	    }
	}
	
	#calculate from sample VCF
	$vcf->open(region => $seg);
	my @bin;
	my @rbin;
	my $start;
	while(my $v = $vcf->next_data_hash()){
	    my $chr = $v->{CHROM};
	    my $pos = $v->{POS};
	    $start = $pos if(! $start);
	    
	    #always calculate  reference first if merged in same VCF
	    if($OPT->{use_ref} && !$rvcf){
		my $rAD = $v->{gtypes}{$rid}{AD};
		if($rAD && $rAD ne '.'){
		    my ($rc, $ac) = split(/,/, $rAD);
		    my $cov = $rc+$ac;
		    
		    if($cov >= $OPT->{maf_cov_filt}){
			my $maf = $ac/$cov;
			if($thr <= $maf && $maf <= 1-$thr){
			    $ref_ok{$seg}{$pos} = [$maf, $cov];
			}
		    }
		}
	    }
	    
	    #calculate for sample
	    my $sAD = $v->{gtypes}{$sid}{AD};
	    if($sAD && $sAD ne '.'){
		my ($rc, $ac) = split(/,/, $sAD);
		my $cov = $rc+$ac;
		
		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($maf >= $thr && $ref_ok{$seg}{$pos}){
			push(@bin, [$maf, $cov]);
			push(@rbin, $ref_ok{$seg}{$pos});
		    }
		}
	    }
	    
	    #in groups of 1000 SNVs
	    my $num = ($OPT->{exome}) ? 200 : 1000;
	    if(@bin && @bin % $num == 0){
		next if(abs($pos - $start)+1 < 2000000); #minimum 2 Mb segment
		
		#get MAF subset of reference
		my %rmaf_set;
		my $rcov_mean;
		my $rmaf_count;
		foreach my $s (@rbin){
		    my ($maf, $cov) = @$s;
		    $rcov_mean += $cov;
		    
		    next if($maf < 0.5);
		    $rmaf_set{$maf}{$cov}++;
		    $rmaf_count++;
		}
		next if($rmaf_count < $num); #too few datapoints
		
		$rcov_mean /= @rbin;
		@rbin = sort {$a->[1] <=> $b->[1]} @rbin;
		my $rcov_median = $rbin[int(@rbin/2)][1];
		
		if(abs($rcov_mean - $rcov_median) >= 2 && !$OPT->{exome}){ #not poisson
		    undef @rbin;
		    undef @bin;
		    undef $start;
		    next;
		}
		
		#fit reference
		my $robs = maf_bin(\%rmaf_set, 100);
		my $rbest;
		for(my $i = 500; $i >= 0; $i-=5){
		    my $maf = $i/1000; #make frequency
		    my $exp = maf_expect(\%rmaf_set,
					 "$maf:".(1-$maf),
					 1,
					 0,
					 {%OPT,
					  k_bins => 100,
					  maf_tail_filt => 0.5});
		    my ($rss, $k) = stat_fit($robs, $exp);
		    if(!$rbest || $rbest->[1] > $rss){
			$rbest = [$maf, $rss];
		    }
		    else{
			$rbest->[2]++;
			last if($rbest->[2] >= 10);
		    }
		}
		my $r_m = $rbest->[0];
		
		#now get MAF subset of sample
		my %smaf_set;
		my $scov_mean;
		my $smaf_count;
		foreach my $s (@bin){
		    my ($maf, $cov) = @$s;
		    $scov_mean += $cov;
		    
		    next if($maf < 0.5);
		    $smaf_set{$maf}{$cov}++;
		    $smaf_count++;
		}
		next if($smaf_count < $num); #too few datapoints

		$scov_mean /= @bin;
		@bin = sort {$a->[1] <=> $b->[1]} @bin;
		my $scov_median = $bin[int(@bin/2)][1];
		
		if(abs($scov_mean - $scov_median) >= 2 && !$OPT->{exome}){ #not poisson
		    undef @rbin;
		    undef @bin;
		    undef $start;
		    next;
		}
		
		#fit the sample
		my $sobs = maf_bin(\%smaf_set, 100);
		my $sbest;
		for(my $i = 500; $i >= 0; $i-=5){
		    my $maf = $i/1000; #make frequency
		    my $exp = maf_expect(\%smaf_set,
					 "$maf:".(1-$maf),
					 1,
					 0,
					 {%OPT,
					  k_bins => 100,
					  maf_tail_filt => 0.5});
		    my ($rss, $k) = stat_fit($sobs, $exp);
		    if(!$sbest || $sbest->[1] > $rss){
			$sbest = [$maf, $rss];
		    }
		    else{
			$sbest->[2]++;
			last if($sbest->[2] >= 10);
		    }
		}		
		my $s_m = $sbest->[0];
		
		#make a standard segment object
		my $r = make_seg_hash($chr, $start, $pos, $OPT);
		$r->{maf_target}      = $s_m;
		$r->{maf_count}       = $smaf_count;
		$r->{maf_set}         = \%smaf_set;
		$r->{cov_count}       = @bin;
		$r->{cov_mean}        = $scov_mean;
		$r->{cov_median}      = $scov_median;
		$r->{ref}{maf_target} = $r_m;
		$r->{ref}{maf_count}  = $rmaf_count;
		$r->{ref}{maf_set}    = \%rmaf_set;
		$r->{ref}{cov_count}  = @rbin;
		$r->{ref}{cov_mean}   = $rcov_mean;
		$r->{ref}{cov_median} = $rcov_median;
		
		#must add xenograft info if available
		if($r->{xeno}){
		    $r->{xeno} = add_vcf_data($r->{xeno}, xeno_param($OPT));
		    $r->{xeno} = add_bam_cov($r->{xeno}, xeno_param($OPT));
		}
		
		push(@sres, $r);
		
		undef @rbin;
		undef @bin;
		undef $start;
	    }
	}

	#save into file
	mkdir($OPT->{outdir}."/disc");
	(undef, my $tfile) = tempfile("disc_XXXX", DIR => $OPT->{outdir}."/disc");
	Storable::store(\@sres, $tfile);
	push(@files, $tfile);
	
	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	my $disc = $DB{disc};
	$disc->{$tfile} = [$seg];
	$DB{disc} = $disc;
	untie(%DB);
	$lock->unlock;
    }

    return @files;
}

#fits the coverage expect to all segments
sub _coverage_fit {
    my $res    = shift;
    my $t_cov  = shift;
    my $OPT    = shift;

    my $cn_max = $OPT->{cn_max};
    my $cfrac  = $OPT->{cfrac};
    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};
    my $thr    = $OPT->{maf_tail_filt};
    my $rfrac  = 1 - $cfrac;

    #estimated target base coverages
    my @bcs = (undef, map {$t_cov * $cfrac/(2*$rfrac + $_*$cfrac)} (1..4));
    my $n = $bcs[$#bcs]/(abs($bcs[$#bcs-1]-$bcs[$#bcs])/8)**2; #optimal points for smoother curve

    #get coverage that generages min RSS and max L
    #my $max = (($t_cov+sqrt($t_cov)) * $cfrac/(2*$rfrac + 1*$cfrac))); #if cn 1 - SD
    #my $min = (($t_cov-sqrt($t_cov)) * $cfrac/(2*$rfrac + 4*$cfrac))); #if cn 4 + SD

    #make list to test with (makes curve)
    my @test = 1..$#bcs;
    for(my $i = 0.6; $i <= $#bcs+0.2; $i += ($#bcs+0.2-0.6)/(51-@bcs)){
	push(@test, $i);
    }

    #map the coverage and MAF liklihood distributions
    my @curve;
    my @stats;
    my $L_min;
    my $L_max;
    my $maf_L_min;
    my $maf_L_max;
    my $cov_L_min;
    my $cov_L_max;
    foreach my $i (sort {$a <=> $b} @test){
	my $t1 = $t_cov * $cfrac/(2*$rfrac + $i*$cfrac); #temp
	my $cov_L_sum = 0;
	my $maf_L_sum = 0;
	my $L_sum     = 0;
	my $count = 0;
	foreach my $r (@$res){
	    my $bc = $t1;
	    if($r->{ref}){
		$bc *= $r->{ref}{ccor_f};
	    }
	    elsif($r->{ccor_f}){
		$bc *= $r->{ccor_f};
	    }

	    my $best = best_cn($r, $n, $bc, $OPT, 0);
	    next if($best->{cn} > $cn_max);

	    $cov_L_sum += $best->{L_cov}*$best->{rL_cov}*$r->{maf_count}; #length weighted
	    $L_sum     += $best->{L}*$best->{rL}*$r->{maf_count}; #length weighted
	    $maf_L_sum += $best->{L_maf}*$best->{rL_maf}*$r->{maf_count}; #length weighted
	    $count += $r->{maf_count};
	}

	#average out sum to smooth out low coverage lack of MAF
	next if(!$count);
	$L_sum     /= $count;
	$cov_L_sum /= $count;
	$maf_L_sum /= $count;

	#find mix/max for normaization
	$cov_L_min = $cov_L_sum if(!defined($cov_L_min) || $cov_L_sum < $cov_L_min);
	$cov_L_max = $cov_L_sum if(!defined($cov_L_max) || $cov_L_sum > $cov_L_max);
	$maf_L_min = $maf_L_sum if(!defined($maf_L_min) || $maf_L_sum < $maf_L_min);
	$maf_L_max = $maf_L_sum if(!defined($maf_L_max) || $maf_L_sum > $maf_L_max);
	$L_min = $L_sum if(!defined($L_min) || $L_sum < $L_min);
	$L_max = $L_sum if(!defined($L_max) || $L_sum > $L_max);

	#make curve
	push(@curve, [$i, $t1, $cov_L_sum, $maf_L_sum, $L_sum]);

	next unless($i == int($i) && $bcs[$i] && !$stats[$i]);
	$stats[$i] = [$i, $t1, $cov_L_sum, $maf_L_sum, $L_sum];
    }
    @stats = grep {$_} @stats; #always missing 0

    print STDERR "\n";
    foreach my $s (@curve){
    	my ($c, $t1, $cov_L_sum, $maf_L_sum, $L_sum) = @$s;
    
    	#scale the liklihood
    	$cov_L_sum = (1 * ($cov_L_sum-$cov_L_min)/($cov_L_max-$cov_L_min));
    	$maf_L_sum = (1 * ($maf_L_sum-$maf_L_min)/($maf_L_max-$maf_L_min));
    	$L_sum = (1 * ($L_sum-$L_min)/($L_max-$L_min));
    
    	print STDERR "$t1\t$cov_L_sum\t$maf_L_sum\t$L_sum\n"; #temp
    }

    #sort and filter candidates
    shift @stats if($stats[0][0] == 1 && $cfrac < 0.15); #remove this one if resolution is too low
    pop @stats if($stats[-1][0] == 4 && $cfrac < 0.15); #remove this one if resolution is too low
    @stats = sort {$b->[4] <=> $a->[4]} @stats; #most likly first
    @stats = grep {$_->[4]/$stats[0][4] > 0.1} @stats; #remove low liklihood
    @stats = grep {$_->[4]/$stats[0][4] > 0.5} @stats; #slightly more strict threshold

    my $base;
    if($stats[0][0] == 1 || $stats[0][0] == 3){ #these are unambiguous
	$base = $stats[0];
    }
    elsif(@stats > 1){
	shift @stats if($stats[0][0] == 4 && $stats[1][0] == 2);
	shift @stats if(@stats > 1 && $stats[0][0] == 2 && $stats[1][0] == 1 && $cfrac > 0.3);
	$base = $stats[0];
    }
    else{
	$base = $stats[0];
    }

    #now adjust to maximum
    my $step = 0.01;
    for(my $t1 = $base->[1]+$step; 1; $t1 += $step){
	my $cov_L_sum = 0;
	my $maf_L_sum = 0;
	my $L_sum     = 0;
	foreach my $r (@$res){
	    my $bc = $t1;
	    if($r->{ref}){
		$bc *= $r->{ref}{ccor_f};
	    }
	    elsif($r->{ccor_f}){
		$bc *= $r->{ccor_f};
	    }

	    my $best = best_cn($r, $n, $bc, $OPT);
	    $cov_L_sum += $best->{L_cov} * $r->{maf_count}; #length weighted
	    $maf_L_sum += $best->{L_maf} * $r->{maf_count}; #length weighted
	    $L_sum     += $best->{L} * $r->{maf_count}; #length weighted
	}
	
	#find closer match
	if($L_sum > $base->[4]){
	    $base = [$base->[0], $t1, $cov_L_sum, $maf_L_sum, $L_sum];
	}
	elsif($step > 0){ #search the opposite way
	    $step = -$step;
	    $t1 = $base->[1]+$step;
	}
	else{
	    last;
	}
    }

    return $base->[1];
}

sub _matches_neighbor {
    my $r1 = shift;
    my $r2 = shift;
    my $OPT = shift;
    my $tag  = shift || 'ucor';

    my $match = 0;

    my $o_cov1 = $r1->{cov_mean};
    my $o_cov2 = $r2->{cov_mean};

    #I should use Skellam to estimate difference probability
    #but I'll aproximate as normal for now
    my $n1 = 1 + ($r1->{q_length}-1)/$OPT->{lfrag};
    $n1 = $r1->{cov_count} if($r1->{cov_count} < $n1);

    my $n2 = 1 + ($r2->{q_length}-1)/$OPT->{lfrag};
    $n2 = $r2->{cov_count} if($r2->{cov_count} < $n2);

    my $n = ($n1 < $n2) ? $n1 : $n2;
    $n = 1 if ($n < 1);
    $n = 100 if($n > 100);

    my $m = $n*$OPT->{base_cov}; #expected difference
    my $s = sqrt($n*($o_cov1+$o_cov2));
    my $x = abs($n*($o_cov1-$o_cov2));
    my $p_diff = gaus_pdf($x, $m, $s) + gaus_pdf(-$x, $m, $s); #get point probability

    $m = 0;
    $s = sqrt($n*($o_cov1+$o_cov2));
    $x = abs($n*($o_cov1-$o_cov2));
    my $p_same = gaus_pdf($x, $m, $s) + gaus_pdf(-$x, $m, $s); #aproximate skellam 

    $match = 1 if($p_same/$p_diff > 10000);
    $match = _alleles_match($r1, $r2, $OPT) if($match);

    return $match;
}

sub best_cn {
    my $r    = shift;
    my $n    = shift;
    my $bc   = shift;
    my $OPT  = shift;
    my $use_maf = shift;

    $use_maf = 1 if(! defined($use_maf));

    if(! $n){
	$n = 1 + ($r->{q_length}-1)/$OPT->{lfrag};
	$n = $r->{cov_count} if($r->{cov_count} < $n);
	$n ||= 1;
    }

    my $cn_max = $OPT->{cn_max};
    my $cfrac  = $OPT->{cfrac};
    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};
    my $thr    = $OPT->{maf_tail_filt};
    my $rfrac  = 1 - $cfrac;
    my $x_adj = ($r->{xeno}) ? $r->{xeno}{cov_median}/$OPT->{xeno}{ccor_f} : 0;
    
    my $cne = ($r->{cov_median} - 2*$bc * $rfrac/$cfrac - $x_adj)/$bc;
    $cne = 0 if($cne < 0);
    
    my $cn0 = int($cne);
    my $cn1 = int($cne)+1;
    
    my $e_cov0 = $cn0*$bc + 2*$bc * $rfrac/$cfrac + $x_adj;
    my $e_cov1 = $cn1*$bc + 2*$bc * $rfrac/$cfrac + $x_adj;
    
    my $pc0 = ($e_cov0) ? p_cov($r, $e_cov0, $n) : p_cov($r, -$e_cov1, $n); #expect of 0 causes issues
    my $pc1 = p_cov($r, $e_cov1, $n);
    
    #scaled relative coverage liklihood
    if($pc0 == $pc1){
    	my $rel = logLr_cov($r, $e_cov0, $e_cov1);
	if($rel <= -700){
    	    ($pc0, $pc1) = (0, 1);
    	}
    	elsif($rel >= 700){
    	    ($pc0, $pc1) = (1, 0);
    	}
    	elsif($rel >= 0){
    	    ($pc0, $pc1) = (1, 1/exp($rel));
    	}
    	elsif($rel < 0){
    	    ($pc0, $pc1) = (exp($rel), 1);
    	}
    }

    #scale relative coverge liklihood
    my $rc0 = $pc0/($pc0+$pc1);
    my $rc1 = $pc1/($pc0+$pc1);

    #use MAF information if available
    my $ratio = ($r->{ref}) ? $r->{ref}{maf_count}/$r->{q_length} : $r->{maf_count}/$r->{q_length};
    my ($m0, $pm0, $rm0, $sc0) = ('.', 1, 0.5, undef);
    my ($m1, $pm1, $rm1, $sc1) = ('.', 1, 0.5, undef);
    if($cn0 <= $cn_max && $r->{maf_count} > 15 && $ratio > 0.0001){
	($m0, $sc0) = @{$r->{maf_score}{best}{"$cn0"}};
	($m1, $sc1) = @{$r->{maf_score}{best}{"$cn1"}};

	#MAF liklihoods
	$pm0 = $sc0->{L};
	$pm1 = $sc1->{L};
	
	#scaled relative liklihood
	my ($aic_min) = ($sc0->{AIC} < $sc1->{AIC}) ? $sc0->{AIC} : $sc1->{AIC};
	$rm0 = exp(($aic_min-$sc0->{AIC})/2);
	$rm1 = exp(($aic_min-$sc1->{AIC})/2);

	my $norm = ($rm0+$rm1);
	$rm0 /= $norm;
	$rm1 /= $norm;
    }

    #combined liklihood
    my $L0 = $pc0*$pm0;
    my $L1 = $pc1*$pm1;

    #combined relative liklihood
    my $rL0 = $rc0*$rm0;
    my $rL1 = $rc1*$rm1;

    my $norm = ($rL0+$rL1);
    $rL0 /= $norm;
    $rL1 /= $norm;
    
    my $res0 = {rL     => $rL0,
		rL_cov => $rc0,
		rL_maf => $rm0,
		L      => $L0,
		L_cov  => $pc0,
		L_maf  => $pm0,
		e_cov  => $e_cov0,
		cn     => $cn0,
		model  => $m0,
		score => $sc0};
    my $res1 = {rL     => $rL1,
		rL_cov => $rc1,
		rL_maf => $rm1,
		L      => $L1,
		L_cov  => $pc1,
		L_maf  => $pm1,
		e_cov  => $e_cov1,
		cn     => $cn1,
		model  => $m1,
		score => $sc1};

    #fix as LOH if the information is very thin (reference derived LOH)
    if($ratio <= 0.0001 && $r->{q_length} > 100000){
	$res0->{model} = '0:'.$res0->{cn};
	$res1->{model} = '0:'.$res1->{cn};
    }

    my $best;

    if($use_maf){
	($best) = sort { $b->{rL} <=> $a->{rL} ||
			 $b->{rL_cov} <=> $a->{rL_cov} ||
			 $b->{rL_maf} <=> $a->{rL_maf} } ($res0, $res1);
    }
    else{
	($best) = sort { $b->{rL_cov} <=> $a->{rL_cov} ||
			 $b->{rL_maf} <=> $a->{rL_maf} ||
			 $b->{rL} <=> $a->{r} } ($res0, $res1);
    }

    return $best;
}

sub poisson_pmf {
    my ($x, $m) = @_;

    if($m == 0){
        return ($x == 0) ? 1 : 0;
    }

    if($m > 1000){#normal aproximation                                          
        return gaus_pdf(int($x), $m, sqrt($m));
    }

    my $pmf = 1;
    my $e = exp(1);
    my $r = $m;
    for(my $i = 1; $i <= $x; $i++){
        $pmf *= $m/$i;

        while($pmf > 1 && $r > 0){
            $pmf /= $e;
            $r--;
        }
    }

    return $pmf/exp($r);
}

sub poisson_cmf {
    my ($x, $m) = @_;

    if($m == 0){
        return ($x >= $m) ? 1 : 0;
    }

    #aproximate with gausian for big numbers
    if($m > 1000){
        return gaus_cdf(int($x)+0.5, $m, sqrt($m));
    }

    my $cmf = 0;
    for(my $i = 0; $i <= $x; $i++){
        $cmf += poisson_pmf($i, $m);
    }

    $cmf = 1 if($cmf > 1);
    return $cmf;
}

#binomial distribution pmf function
sub binomial_pmf {
    my $k = shift;
    my $n = shift;
    my $p = shift;

    if($n > 100 && $n*$p > 5 && $n*(1-$p) > 5){#normal aproximation
        my $m = $n*$p;
        return gaus_pdf($k, $m, sqrt($m*(1-$p)));
    }

    my $d = $n-$k;
    return binomial_coef($n, $k)*($p**$k)*((1-$p)**$d);
}

#binomial distribution CMF function
sub binomial_cmf {
    my $k = shift;
    my $n = shift;
    my $p = shift;

    if($n > 100){#normal aproximation
        my $m = $n*$p;
        return gaus_cdf(int($k)+0.5, $m, sqrt($m*(1-$p)));
    }

    my $cmf = 0;
    foreach (my $i = 0; $i <= $k; $i++){
        $cmf += binomial_pmf($i, $n, $p);
    }

    return $cmf;
}

{
my @BUF; #holds previously calculated coefficients
sub binomial_coef {
    my $n = shift;
    my $k = shift;

    return 0 if($k < 0 || $k > $n);

    $k = $n - $k if($k > $n - $k); # take advantage of symmetry

    return $BUF[$n][$k] if($BUF[$n][$k]);

    my $c = 1;
    for(my $i = 1; $i <= $k; $i++){
        $c *= ($n-($k-$i));
	$c /= $i;
    }

    $BUF[$n][$k] = $c;
    return $c;
}
}

#Statistics::Distributions has counterintuitive variable order. use
#this instead --> chisqr(a, k); where a is the left tail probability
sub chisqr {
    return Statistics::Distributions::chisqrdistr(int($_[1]), 1-$_[0]);
}

#variation of the Geary–Hinkley transformation to get distribution
#ratio of poisson over poisson. Basically scaled diference of mean
#divided by scaled standard deviation to transform into normal
sub c_transform{
    my $x  = shift;
    my $n0 = shift; #numerator distribution mean
    my $d0 = shift; #denominator distribution mean

    #ratio of distribution means
    my $r = $n0/$d0;

    #confidence interval of 1 SD of mean
    #tranformation assume 1 SD
    my $crit = (1-.682689492);

    #poisson confidence innterval eqution for numerator
    #normal one SD equivilent for either side
    my $nm = int(2*$n0) ? 0.5 * chisqr($crit/2, 2*$n0) : $n0;
    my $np = 0.5 * chisqr(1-$crit/2, 2*$n0+2);

    #poisson confidence innterval eqution for denominator
    #normal one SD equivilent for either side
    my $dm = int(2*$d0) ? 0.5 * chisqr($crit/2, 2*$d0) : $d0;
    my $dp = 0.5 * chisqr(1-$crit/2, 2*$d0+2);
    
    #perform transform
    my $nsd = ($x <= $r) ? $nm : $np; #nummerator SD equivilent
    my $dsd = ($x <= $r) ? $dp : $dm; #denominator SD equivilent
    my $t = ($d0*$x-$n0)/sqrt($dsd* $x**2+$nsd);

    return gaus_pdf($t, 0, 1);
}
