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
use Bio::DB::Sam;
use Tie::MmapArray;
use constant PI => 4 * atan2(1, 1); #3.14159265358979323846...

our ($sort_exe, $bgzip_exe, $tabix_exe, $vcftools_exe, $hg19, $dummy);
our ($DUM, %VCF, %BAMS, %RG, %GEN, %CHRS);

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

     cell                Estimate cellularity of mixed tumor sample

     cfrac     <NUM>     Set cellularity explicitly

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
$OPT{cn_max}        = 5;
$OPT{models}        = get_models($OPT{cn_max}); #i.e. 0:1, 1:2, etc. up to specified CN
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
           "t_cov=f"     => \$OPT{user_t_cov},
           "rt_cov=f"    => \$OPT{ref}{user_t_cov},
           "base_cov=f"  => \$OPT{user_base_cov},
           "rbase_cov=f" => \$OPT{ref}{user_base_cov},
           "ploidy=f"    => \$OPT{user_ploidy},
           "rploidy=f"   => \$OPT{ref}{user_ploidy},
           "help|?"      => sub{print $usage; exit(0);});

$OPT{seg_file} = shift;
$OPT{vcf_file} = shift;

#--print usage
if(! $OPT{seg_file} || ! $OPT{vcf_file}){
    print $usage;
    exit(0);
}

#--validate and prepare command line options
my $err = validate_options(\%OPT);
die $err if($err);

#--prepare input files
print STDERR "#Preparing and indexing input files...\n";
$OPT{vcf_file}        = prepare_vcf($OPT{vcf_file});
$OPT{ref}{vcf_file}   = prepare_vcf($OPT{ref}{vcf_file});
$OPT{xeno}{vcf_file}  = prepare_vcf($OPT{xeno}{vcf_file});
$OPT{bam_files}       = [prepare_bam(@{$OPT{bam_files}})];
$OPT{ref}{bam_files}  = [prepare_bam(@{$OPT{ref}{bam_files}})];
$OPT{xeno}{bam_files} = [prepare_bam(@{$OPT{xeno}{bam_files}})];
$OPT{fasta}           = prepare_fasta($OPT{fasta});
$OPT{CHRS}            = get_chr_lengths($OPT{fasta});
$OPT{GAPS}            = get_chr_gaps($OPT{fasta});

%CHRS = %{$OPT{CHRS}};

#--process stored settings
$OPT{DB_File} = $OPT{outdir}.'/dbfile';
tie(my %DB, 'MLDBM', $OPT{DB_File}, O_CREAT, 0640) or die $!;

#check to see if the options are changed
if(my $DBOPT = $DB{OPT}){
    my @check = qw(use_xeno
		   use_ref
		   use_generic
		   vcf_file
		   bam_files
		   fasta
		   sid
		   xcov
		   cell
		   cfrac);

    my @comp = ([\%OPT,$DBOPT],
		[ref_param(\%OPT), ref_param($DBOPT)],
		[xeno_param(\%OPT), xeno_param($DBOPT)]);

    my $clober;
    foreach my $o (@comp){
	my ($new, $old) = @$o;
	foreach my $key (@check){
	    next if(!$new->{$key} && !$old->{$key});
	    my ($v1, $v2) = ($new->{$key}, $old->{$key});
	    if($key eq 'bam_files'){
		$v1 = join(',', sort map {Cwd::abs_path($_)} @$v1);
		$v2 = join(',', sort map {Cwd::abs_path($_)} @$v2);
	    }
	    next if($v1 eq $v2);

	    $clober = 1;
	    last;
	}
	if($clober){
	    untie(%DB);
	    File::Path::rmtree($OPT{outdir});
	    mkdir($OPT{outdir});
	    tie(my %DB, 'MLDBM', $OPT{DB_File}, O_CREAT, 0640) or die $!;
	    last;
	}
    }
}
$DB{OPT} = \%OPT; #store current options
$DB{load} ||= {};
$DB{ARC}  ||= {};
my $DBARC = $DB{ARC}; #archived parameters for hard to calculate values
untie(%DB);

#--estimate xenograft content
if($OPT{use_xeno}){
    print STDERR "#Estimating mouse DNA content...\n";
    $OPT{xcov} = estimate_xeno($OPT{vcf_file}, \%OPT) if(!$OPT{xcov});

    #use sequenced mouse genome to normalize by region coverage
    if($OPT{xeno}{vcf_file}){
	$OPT{xeno}{xcov} = estimate_xeno($OPT{xeno}{vcf_file}, xeno_param(\%OPT));
	$OPT{xeno}{ccor_f} = $OPT{xcov}/$OPT{xeno}{xcov};
    }

    print STDERR "#Mean mouse coverage: $OPT{xcov}\n";
}

#--estimate tumor cellularity
if($OPT{cell}){
    print STDERR "#Estimating tumor cellularity...\n";
    ($OPT{cfrac}, $OPT{cfrac_rsq}) = estimate_cellularity(\%OPT);
    $OPT{rfrac}  = 1 - $OPT{cfrac};
    $OPT{k_bins} = ceil($OPT{k_bins}/$OPT{cfrac}); #adust for cellularity
    $OPT{cn_max} = ($OPT{cfrac} <= 0.25) ? 3 : 5;
    $OPT{models} = get_models($OPT{cn_max});

    #archive the value
    $DBARC->{cfrac} = $OPT{cfrac};
    tie(my %DB, 'MLDBM', $OPT{DB_File}, O_RDWR, 0640) or die $!;
    $DB{ARC} = $DBARC;
    untie(%DB);

    print STDERR "#Tumor cellularity: $OPT{cfrac}\n";
    print STDERR "#R^2 of fit: $OPT{cfrac_rsq}\n" if(defined $OPT{cfrac_rsq});
}

#--now read in segments from file
print STDERR "#Reading segment file...\n";
my $segs = read_segment_file($OPT{seg_file}, \%OPT);

#--fill in the VCF/BAM data for segments
print STDERR "#Loading VCF/BAM data into segments...\n";
my $res;
if($OPT{smooth} && -f "$OPT{outdir}/smooth_merge.store"){
    $res = Storable::retrieve("$OPT{outdir}/smooth_merge.store");
}
elsif($OPT{merge} && -f "$OPT{outdir}/merge.store"){
    $res = Storable::retrieve("$OPT{outdir}/merge.store");
}
else{
    $res = process_segments($segs, \%OPT, 1);
}

#--don't know if I need these but they may come in useful
$OPT{n50} = n50($res);
$OPT{seg_mean}= seg_mean($res);

#--get coverage expects and add copy number info to reference
my @refs;
foreach my $r (@$res){
    push(@refs, $r->{ref});
}

#--calculate coverage expects
my $param = ref_param(\%OPT);
print STDERR "#Getting expected coverage for reference...\n" if($OPT{use_ref});
($OPT{ref}{t_cov}, $OPT{ref}{base_cov}) = fit_expected_coverage(\@refs, $param);
print STDERR "#Ref expected segment median: $OPT{ref}{t_cov}\n" if($OPT{use_ref});
print STDERR "#Ref base coverage (1 copy): $OPT{ref}{base_cov}\n" if($OPT{use_ref});

#--now assign reference copy numbers
$param = ref_param(\%OPT); #reset
print STDERR "#Assigning copy numbers to reference...\n";
$OPT{ref}{chr_expects} = {_ALL => $OPT{ref}{base_cov}};
assign_copy_number(\@refs, $OPT{ref}{chr_expects}, $param);

#--get reference ploidy
$param = ref_param(\%OPT); #reset
($OPT{ref}{ploidy_ave}, $OPT{ref}{ploidy}) = get_ploidy(\@refs, $param);
print STDERR "#Reference ploidy mode: $OPT{ref}{ploidy}\n" if($OPT{use_ref});
print STDERR "#Reference ploidy mean: $OPT{ref}{ploidy_ave}\n" if($OPT{use_ref});

#--calculate coverage correction based on referece
print STDERR "#Calculating coverage corrections from reference...\n";
foreach my $r (@refs){
    my $exp = $OPT{ref}{chr_expects}{$r->{chr}};
    my $obs = $r->{cov_median}/$OPT{ref}{ploidy};
    $r->{ccor_f} = ($r->{ucor}{cn}) ? $obs/$exp : 1;
}

#--now maximize liklihood to identify expected coverage
print STDERR "#Getting expected coverage for sample...\n";
($OPT{t_cov}, $OPT{base_cov}) = fit_expected_coverage($res, \%OPT);
print STDERR "#Expected segment median: $OPT{t_cov}\n";
print STDERR "#Base coverage (1 copy): $OPT{base_cov}\n";

#--now assign copy numbers to each allele
print STDERR "#Assigning copy numbers...\n";
$OPT{chr_expects} = {_ALL => $OPT{base_cov}};
assign_copy_number($res, $OPT{chr_expects}, \%OPT);

#--report some result statistics
($OPT{ploidy_ave}, $OPT{ploidy}) = get_ploidy($res, \%OPT);
print STDERR "#Sample ploidy mode: $OPT{ploidy}\n";
print STDERR "#Sample ploidy mean: $OPT{ploidy_ave}\n";
label_somatic($res, \%OPT);

#--smooth out segments and recalculate the copy numbers
my $count = @$res;
while($OPT{merge}){
    print STDERR "#Merging...\n";
    my $count2 = @$res;
    $res = merge_segments($res, {%OPT, smooth => 0}); #force no smoothing at first
    if($count2 != @$res){
	#assign all final values
	print STDERR "#Reassigning copy numbers...\n";
	recall_cn($res, \%OPT);
	Storable::store($res, "$OPT{outdir}/merge.store");
	next;
    }
    last;
}

while($OPT{smooth}){
    print STDERR "#Smoothing and remerging...\n";
    my $count2 = @$res;
    $res = merge_segments($res, \%OPT);
    if($count2 != @$res){
	#assign all final values
	print STDERR "#Reassigning copy numbers...\n";
	recall_cn($res, \%OPT);
	Storable::store($res, "$OPT{outdir}/smooth_merge.store");
	next;
    }
    last;
}

#print final statistics
print STDERR "\n#Final statistics...\n";
if($OPT{use_ref}){
    print STDERR "#Reference segment median: $OPT{ref}{t_cov}\n";
    print STDERR "#Reference base coverage (1 copy): $OPT{ref}{base_cov}\n";
    print STDERR "#Reference ploidy mode: $OPT{ref}{ploidy}\n";
    print STDERR "#Reference ploidy mean: $OPT{ref}{ploidy_ave}\n";
    my $exp = $OPT{ref}{chr_expects};
    foreach my $chr (sort {_chrom($a) <=> _chrom($b)} keys %$exp){
	my $bc = $exp->{$chr};
	print STDERR "#Reference adjust $chr base coverage: $bc\n";
    }
}
print STDERR "#Sample segment median: $OPT{t_cov}\n";
print STDERR "#Sample base coverage (1 copy): $OPT{base_cov}\n";
print STDERR "#Sample ploidy mode: $OPT{ploidy}\n";
print STDERR "#Sample ploidy mean: $OPT{ploidy_ave}\n\n";
my $exp = $OPT{chr_expects};
foreach my $chr (sort {_chrom($a) <=> _chrom($b)} keys %$exp){
    my $bc = $exp->{$chr};
    print STDERR "#Sample Adjusted $chr base coverage: $bc\n";
}

#--now write results to file
print STDERR "#Writing GVF...\n";
write_gvf($res);

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------
#run before forking. clears stored values that may cause core dump
sub prefork {
    undef $DUM;
    undef %BAMS;
    undef %VCF;
}

sub fit_maf_kernel {
    my $maf_set = shift;

    my $ker = Statistics::KernelEstimation->new();
    foreach my $maf (keys %$maf_set){
	foreach my $cov (keys %{$maf_set->{$maf}}){
	    my $count = $maf_set->{$maf}{$cov};
	    $ker->add_data($maf, $count, 1/($cov+1));
	}
    }

    my $w = $ker->default_bandwidth;
    my $best = [0.5, 0];
    for(my $x = 0; $x <= 100; $x++) { #0 to 1 when dividing by 100
	#my $y = $ker->pdf_width_from_data($x/100);
	my $y = $ker->pdf($x/100, $w);
	$best = [$x/100, $y] if($best->[1] < $y);
    }

    #refine
    my $i = int($best->[0]*1000+0.5); #round just incase
    for(my $x = $i-10; $x <= $i+10; $x++) {
	#my $y = $ker->pdf_width_from_data($x/1000);
	my $y = $ker->pdf($x/1000, $w);
	$best = [$x/1000, $y] if($best->[1] < $y);
    }

    return $best->[0];
}

#runs after estimating cellularity to recalculate MAF stats
sub recalc_maf {
    my $OPT = shift;

    my @work :shared;
    tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
    my $load = $DB{load};
    untie(%DB);

    my @files;
    my @work_add; #not shared for speed
    foreach my $file (keys %$load){
	if($file =~ /\.cell$/){
	    push(@files, $file);
	}
	else{
	    push(@work_add, $file);
	}
    }
    push(@work, @work_add);
    
    #split off threads
    prefork(); #prepare for forking
    my @threads;
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'list'},
				       \&_recalc_maf_thread,
				       \@work,
				       $OPT));
    }

    #jump into the mix
    my @ret = _recalc_maf_thread(\@work, $OPT);
    push(@files, @ret);

    #gather thread results
    foreach my $thr (@threads){
	my @ret = $thr->join;
	push(@files, @ret);
    }

    #retrieve serialized results
    my @res;
    foreach my $file (@files){
	my $st = Storable::retrieve($file);
	push(@res, @$st);
    }

    return \@res;
}

sub _recalc_maf_thread {
    my $work = shift;
    my $OPT = shift;

    my @files;
    while(my $file = shift @$work){
	my $res = Storable::retrieve($file);
	add_maf_stats($_, $OPT, 1) foreach(@$res); #all MAF scores must be fixed
	Storable::store($res, "$file.cell");

	my $lock;
	$lock = File::NFSLock->new($OPT->{DB_File}, 'EX', 10, 30) while(!$lock);
	tie(my %DB, 'MLDBM', $OPT->{DB_File}, O_RDWR, 0640) or die $!;
	my $load = $DB{load};
	my $list = $load->{$file};
	delete($load->{$file});
	$load->{"$file.cell"} = $list;
	$DB{load} = $load;
	untie(%DB);
	$lock->unlock;

	unlink($file);
        push(@files, "$file.cell");
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
    my $fasta     = $OPT->{fasta};

    my $err;
    $err .= "File $seg_file does not exist\n" if(! -f $seg_file);
    $err .= "File $fasta does not exist\n" if(! -f $fasta);
    $err .= "File $vcf_file does not exist\n" if(! -f $vcf_file);
    $err .= "File $rvcf_file does not exist\n" if($rvcf_file && ! -f $rvcf_file);
    
    foreach ($OPT->{C}, $OPT->{B}, $OPT->{E}){
	$_ = int($_-1);
	$err .= "Column selection must be a positive integer > 1\n" if($_ < 0);
    }
    
    my @bam_files;
    if(my $bam_dir = $OPT->{bam_dir}){
	if(!-d $bam_dir){
	    $err .= "Directory $bam_dir does not exist\n" if(! -d $bam_dir);
	}
	else{
	    push(@bam_files, <$bam_dir/*.bam>);
	    $err .= "$bam_dir contained no bam file\n" if(! @bam_files);
	}
    }
    if(my $bam_list = $OPT->{bam_list}){
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
    my @rbam_files;
    if(my $rbam_dir = $OPT->{ref}{bam_dir}){
	if(!-d $rbam_dir){
	    $err .= "Directory $rbam_dir does not exist\n" if(! -d $rbam_dir);
	}
	else{
	    push(@rbam_files, <$rbam_dir/*.bam>);
	    $err .= "$rbam_dir contained no bam file\n" if(! @rbam_files);
	}
    }
    if(my $rbam_list = $OPT->{ref}{bam_list}){
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
    my @xbam_files;
    if(my $xbam_dir = $OPT->{xeno}{bam_dir}){
	if(!-d $xbam_dir){
	    $err .= "Directory $xbam_dir does not exist\n" if(! -d $xbam_dir);
	}
	else{
	    push(@xbam_files, <$xbam_dir/*.bam>);
	    $err .= "$xbam_dir contained no bam file\n" if(! @xbam_files);
	}
    }
    if(my $xbam_list = $OPT->{xeno}{bam_list}){
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

    my ($base) = $seg_file =~ /([^\/]+)$/;
    my $outdir = Cwd::cwd()."/$base.cnv.output";
    mkdir($outdir);

    #fix dependent values
    $OPT->{bam_files}  = \@bam_files;
    $OPT->{ref}{bam_files} = \@rbam_files;
    $OPT->{xeno}{bam_files} = \@xbam_files;
    $OPT->{use_xeno}   = 1 if($OPT->{xcov});
    $OPT->{merge}      = 1 if($OPT->{smooth}); #implies merge
    $OPT->{outdir} = $outdir;
    $OPT->{base}   = $base;

    #adust for cellularity
    if($OPT->{cfrac} < 1){
	$OPT->{cell}   = 1;
	$OPT->{k_bins} = ceil($OPT->{k_bins}/$OPT->{cfrac});
	$OPT->{cn_max}   = ($OPT->{cfrac} <= 0.25) ? 3 : 5;
	$OPT->{models}   = get_models($OPT->{cn_max});
    }


    if($OPT{ref}{sid} || $OPT{ref}{vcf_file}){
	$OPT{use_ref} = 1;
    }
    else{
	$OPT{use_generic} = 1;
    }
    if($OPT{xeno}{sid} || $OPT{xeno}{vcf_file}){
	$OPT{use_xeno} = 1;
    }

    return $err;
}

#build fasta index
sub prepare_fasta {
    my @fastas = @_;

    foreach my $fasta (@fastas){
	next if(!$fasta);

	my $bf = $dummy;
	my $sam = Bio::DB::Sam->new(-bam  => $bf,
				    -autoindex => 1,
				    -fasta=> $fasta);
    }

    return (wantarray) ? @fastas : shift @fastas;
}

#indexes the bam files
sub prepare_bam {
    my @bam_files = @_;

    foreach my $f (@bam_files){
	next if(!$f);

	Bio::DB::Sam->new(-bam  => $f,
			  -autoindex => 1);
    }

    return (wantarray) ? @bam_files : shift @bam_files;
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
				       \&_estimate_xeno,
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
    my $list = _estimate_xeno($vcf_file, $OPT, \@work, \%marker, \$flag);
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
    $w = 1 if($w < 1);
    my ($low, $high) = $ker->extended_range();
    $low = 0 if($low < 0);
    my $peak = [0, 0];
    for(my $x = $low; $x <= $high; $x += ($high-$low)/1000) {
	my $y = $ker->pdf($x, $w);
	#print "$x\t$y\n"; #temp
	$peak = [$x, $y] if($y > $peak->[1]);
	last if($x == $high);
    }

    $high = $peak->[0] + ($high-$low)/100;
    $low  = $peak->[0] - ($high-$low)/100;
    for(my $x = $low; $x <= $high; $x += 0.001) {
        my $y = $ker->pdf($x, $w);
        $peak = [$x, $y] if($y > $peak->[1]);
        last if($x == $high);
    }

    return sprintf('%.2f', $peak->[0]);
}

sub _estimate_xeno {
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

sub n50 {
    my $segs = shift;

    my $sum = 0;
    my @lengths;
    foreach my $s (@$segs){
	my $chr;
	my $start;
	my $end;

	if(ref($s)){
	    ($chr, $start, $end) = ($s->{chr}, $s->{start}, $s->{end});
	}
	else{
	    ($chr, $start, $end) = $s =~ /^(.*)\:(\d+)\-(\d+)$/;
	}

	push(@lengths, ($end - $start)+1);
	$sum += $lengths[-1];
    }

    my $half = 0;
    foreach my $l (sort {$a <=> $b} @lengths){
	$half += $l;
	if($half/$sum >= 0.5){
	    return $l;
	}
    }
}

sub seg_mean {
    my $segs = shift;

    my $sum = 0;
    my $count = 0;
    foreach my $s (@$segs){
	my $chr;
	my $start;
	my $end;

	if(ref($s)){
	    ($chr, $start, $end) = ($s->{chr}, $s->{start}, $s->{end});
	}
	else{
	    ($chr, $start, $end) = $s =~ /^(.*)\:(\d+)\-(\d+)$/;
	}

	$sum += ($end - $start)+1;
	$count++;
    }
    
    return ($count) ? $sum/$count : 0;
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
	if(!-f $file){
	    delete $load->{$file};
	    next;
	}
	
	push(@files, $file);
	$seen{$_}++ foreach(@{$load->{$file}});
    }
    $DB{load} = $load;
    untie(%DB);

    #distribute the work
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
    }
    push(@work, @work_add);
    $flag = 0; #lets threads know everything is loaded
    my $remain = $work_count - $finished;

    printf STDERR "Total: %ibp\n", $work_count;
    printf STDERR "Processed: %ibp (%.0f%%)\n", $finished, 100*$finished/$work_count;
    printf STDERR "Remaining: %ibp (%.0f%%)\n", $remain, 100*$remain/$work_count;

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

	#--progress bar configuration
	#create progress bar parameters
	my %bparam = (name  => 'Processing',
		      count => $remain,
		      ETA   => 'linear');
	
	#fix terminal size capture in SGE redirect
	if(! is_interactive(*STDERR)){
	    $ENV{COLUMNS} = 80;
	    $bparam{bar_width}  = 80;
	    $bparam{term_width} = 80;
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
    }

    #retrieve serialized results
    my @res;
    foreach my $file (@files){
	my $st = Storable::retrieve($file);
	push(@res, @$st);
    }

    return \@res;
}

sub _process_segments_thread {
    my $work    = shift;
    my $OPT     = shift;
    my $flag    = shift;
    my $counter = shift;
    my $total   = shift;
    my $prog    = shift;

    my $fasta      = $OPT->{fasta};
    my $models     = $OPT->{models};
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
	$r = add_maf_stats($r, $OPT);
	$r = add_generic_cov($r, $OPT) if($OPT->{use_generic});

	#add bam coverage only for short segments
	if($r->{q_length} < 30000 || $r->{maf_count} < 30 ||
	   $r->{q_length}/$r->{cov_count} > 2000
	){
	    $r = add_bam_cov($r, $OPT);
	}
	#always use bam for at xeno match
	elsif($r->{xeno}){
	    $r->{xeno} = add_bam_cov($r->{xeno}, xeno_param($OPT));
	}
	
	push(@res, $r);
	push(@$loaded, $seg);
	$counter_buf += ($end-$start)+1; #itterate progress bar

	#keep memory usage lower
	if(@res > 500){
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

    #put it all in the segment result
    my $res = { chr        => $chr,
		start      => $start,
		end        => $end,
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

    if($OPT->{use_xeno} && ($OPT->{xeno}{vcf_file} || $OPT->{xeno}{sid})){
        #need xenograft model
        $res->{xeno} = { chr        => $chr,
			 start      => $start,
			 end        => $end,
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
    my $vcf;
    if($VCF{$vcf_file}){
	$vcf = $VCF{$vcf_file};
    }
    else{
	$vcf = Vcf->new(file=>"$vcf_file");
	$vcf->parse_header();
	$VCF{$vcf_file} = $vcf;
    }

    my ($SID) = grep {/$sid$/} $vcf->get_samples() if($sid);
    $vcf->set_samples(include=>[$SID]) if($sid);
   
    #load reference VCF
    my $rvcf;
    my $RID;
    if($rvcf_file && $rvcf_file ne $vcf_file){
	if($VCF{$rvcf_file}){
	    $rvcf = $VCF{$rvcf_file};
	}
	else{
	    $rvcf = Vcf->new(file=>"$rvcf_file");
	    $rvcf->parse_header();
	    $VCF{$rvcf_file} = $rvcf;
	}
	($RID) = grep {/$rid$/} $rvcf->get_samples() if($rid);
	$rvcf->set_samples(include=>[$RID]) if($rid);
    }
    elsif($rid){
	($RID) = grep {/$rid$/} $vcf->get_samples();
	$vcf->set_samples(include=>[$SID, $RID]);
    }

    #load xeno VCF
    my $xvcf;
    my $XID;
    if($xvcf_file && $xvcf_file ne $vcf_file){
	if($VCF{$xvcf_file}){
	    $xvcf = $VCF{$xvcf_file};
	}
	else{
	    $xvcf = Vcf->new(file=>"$xvcf_file");
	    $xvcf->parse_header();
	    $VCF{$xvcf_file} = $xvcf;
	}
	($XID) = grep {/$xid$/} $xvcf->get_samples() if($xid);
	$xvcf->set_samples(include=>[$XID]) if($xid);
    }
    elsif($xid){
	($XID) = grep {/$xid$/} $vcf->get_samples();
	my @ids = grep {$_} ($SID, $RID, $XID);
	$vcf->set_samples(include=>[@ids]);
    }

    #get observed MAF and VCF coverage for segment
    my $cov_count   = 0;
    my $rcov_count  = 0;
    my $xcov_count  = 0;
    my $maf_count   = 0;
    my $rmaf_count  = 0;
    my $xmaf_count  = 0;
    my %vcf_set;
    my %rvcf_set;
    my %xvcf_set;
    my %maf_set;
    my %rmaf_set;
    my %xmaf_set;

    #calculate first for separate reference VCF
    my %ref_ok;
    my ($chr, $start, $end) = ($res->{chr}, $res->{start}, $res->{end});
    if($rvcf && $OPT->{use_ref}){
	$rvcf->open(region=> "$chr:$start-$end");
	while(my $v = $rvcf->next_data_hash()){
	    my $pos = $v->{POS};
	    my $rAD = ($RID) ? $v->{gtypes}{$RID}{AD} : $v->{gtypes}{AD};
	    if($rAD && $rAD ne '.'){
		my ($rc, $ac) = split(/,/, $rAD);
		my $cov = $rc+$ac;
		$rvcf_set{$cov}++;
		$rcov_count++;

		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($thr <= $maf){
			$rmaf_set{$maf}{$cov}++;
			$rmaf_count++;
			if($maf <= 1-$thr){
			    $ref_ok{$chr}{$pos}++;
			}
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
	if($OPT->{use_ref} && $RID && !$rvcf){
	    my $rAD = $v->{gtypes}{$RID}{AD};
	    if($rAD && $rAD ne '.'){
		my ($rc, $ac) = split(/,/, $rAD);
		my $cov = $rc+$ac;
		$rvcf_set{$cov}++;
		$rcov_count++;

		if($cov >= $OPT->{maf_cov_filt}){
		    my $maf = $ac/$cov;
		    if($thr <= $maf){
                        $rmaf_set{$maf}{$cov}++;
                        $rmaf_count++;
			if($maf <= 1-$thr){
			    $ref_ok{$chr}{$pos}++;
			}
                    }
		}
	    }
	}

	#calculate for sample
	my $sAD = ($SID) ? $v->{gtypes}{$SID}{AD} : $v->{gtypes}{AD};
	if($sAD && $sAD ne '.'){
	    my ($rc, $ac) = split(/,/, $sAD);
	    my $cov = $rc+$ac;
	    $vcf_set{$cov}++;
	    $cov_count++;

	    if($cov >= $OPT->{maf_cov_filt}){
		my $maf = $ac/$cov;
		if($maf >= $thr && (!$OPT->{use_ref} || $ref_ok{$chr}{$pos})){
		    $maf_set{$maf}{$cov}++;
		    $maf_count++;
		}
	    }
	}

	#calculate when xeno is merged with sample VCF
	if($res->{use_xeno} && $XID && !$xvcf){
	    my $xAD = $v->{gtypes}{$XID}{AD};
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
    if($xvcf && $OPT->{use_xeno}){
	$xvcf->open(region=> "$chr:$start-$end");
	while(my $v = $xvcf->next_data_hash()){
	    my $xAD = ($XID) ? $v->{gtypes}{$XID}{AD} : $v->{gtypes}{AD};
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
    }

    if($OPT->{use_xeno}){
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
    my $xcov       = $OPT->{xcov};
    my $cfrac      = $OPT->{cfrac};
    my $rfrac      = 1 - $cfrac;

    #get target frequencies (optimization for redundancy)
    my %target_fs;
    foreach my $m (@{$OPT->{models}}){
	next if($m eq '!');
	push(@{$target_fs{get_f_string($m, $cfrac)}}, $m);
    }

    #==for long segments sample around the mean to build expect
    my $exp_set;
    if($maf_count > 100 && defined($cov_mean)){
	my $max = (3*$cov_mean)+1;
	$max = 100 if($max < 100);
	my %all;
	foreach my $rhash (values %$maf_set){
            foreach my $key (keys %{$rhash}){
		next unless($rhash->{$key});
		$all{$key} += $rhash->{$key};
            }
        }
	my @keys = keys %all;
	if(@keys > 20){
	    my @greater = sort {$a <=> $b} grep {$_ > int($cov_mean) && $_ < $max} @keys;
	    my @less = sort {$b <=> $a} grep {$_ <= int($cov_mean)} @keys;
	    undef @keys;
	    for (my $i = 0; $i < @greater && $i < 20; $i++){
		push(@keys, $greater[$i]) if($i % 2);
	    }

	    for (my $i = 0; $i < @less && $i < 20; $i++){
		push(@keys, $less[$i]) if(!$i % 2);
	    } 
	}
	my %new = (map {$_ => $all{$_}} @keys);
	$exp_set = {0.5 => \%new}; #new smaller maf_set

	#revise maf_count
	$maf_count = 0;
	while(my ($cov, $num) = each %new){
	    $maf_count += $num;
	}
    }
    else{
	$exp_set = $maf_set; #all data
    }

    #==score with normalized composite liklihood
    #my %t_scores;
    #foreach my $t (keys %target_fs){ #initialize score for each target
    #	next if($t eq '!'); #temp 
    #	$t_scores{$t} = 1;
    #	$t_scores{"$t:xeno"} = 1 if($xcov);
    #}
    #
    #combine score for all mafs
    #keys %{$maf_set};
    #while(defined(my $maf = each %{$maf_set})){
    #	keys %{$maf_set->{$maf}};
    #   while(my ($cov, $num) = each %{$maf_set->{$maf}}){
    #	    #calculate individual liklihoods
    #	    my $max = 0;
    #	    my %Ls;
    #	    foreach my $t (keys %t_scores){
    #		next if($t =~ /xeno$/); #would be redundant
    #		my $i = int($maf*$cov+0.5); #round just in case
    #		my $mod = $target_fs{$t}[0];
    #		$Ls{$t} = _maf_point($i, $cov, $mod, $cfrac, 0, $OPT);
    #		$max = $Ls{$t} if($Ls{$t} > $max);
    #
    #		next unless($xcov);
    #		$Ls{"$t:xeno"} = _maf_point($i, $cov, $mod, $cfrac, $xcov, $OPT);
    #		$max = $Ls{"$t:xeno"} if($Ls{"$t:xeno"} > $max);
    #	    }
    #	    next if($max == 0);
    # 
    #       #scale individual liklihoods
    #       foreach my $t (keys %t_scores){
    #           $t_scores{$t} += ($Ls{$t}/$max)*$num; #composite liklihood
    #       }
    #   }
    #}
    #
    #scale combined liklihoods
    #my $maxt = 0;
    #foreach my $t (keys %t_scores){
    #    $maxt = $t_scores{$t} if($t_scores{$t} > $maxt);
    #}
    #next if($maxt == 0);
    #normalizes and scales
    #scaling still underestimates differences in liklihood
    #$_ = ($_/$maxt)**$r->{maf_count} foreach (values %t_scores);
    #
    #add to final output hash
    #my %maf_score;
    #foreach my $t (keys %t_scores){
    #	if($t =~ s/\:xeno$//){
    #	    my @ms = @{$target_fs{$t}};
    #	    $maf_score{"$_:xeno"}{L} = $t_scores{"$t:xeno"} foreach(@ms);
    #	}
    #	else{
    #	    my @ms = @{$target_fs{$t}};
    #	    $maf_score{$_}{L} = $t_scores{$t} foreach(@ms);
    #	}
    #}
    
    #==rss fit of models with AIC relative liklihoods
    my %maf_score;
    my $obs = maf_bin($maf_set, $OPT->{k_bins}, $OPT); #bin the observed MAF values
    foreach my $t (keys %target_fs){
	next if($t eq '!'); #temp
	my @ms = @{$target_fs{$t}};
	my $m = $ms[0];

	my ($rss, $k, $aic) = (0, 1, 0);
	if($q_length > 30000 && $maf_count){
	    my $expect = maf_expect($exp_set, $m, $cfrac, 0, $OPT); #get expected MAF distribution
	    ($rss, $k) = stat_fit($obs, $expect);
	    $aic = $maf_count*log($rss/$maf_count)+2*$k;
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
	    my $expect = maf_expect($exp_set, $m, $cfrac, $xcov, $OPT); #get expected MAF distribution
	    ($rss, $k) = stat_fit($obs, $expect);
	    $aic = $maf_count*log($rss/$maf_count)+2*$k;
	}

	foreach my $m (@ms){
	    $maf_score{"$m:xeno"}{rss} = $rss;
	    $maf_score{"$m:xeno"}{k} = $k;
	    $maf_score{"$m:xeno"}{AIC} = $aic;
	}
    }

    #determine AIC relative liklihoods
    my ($aic_min) = sort {$a <=> $b} (map {$_->{AIC}} values %maf_score);
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
    my $dup     = shift || 0;

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

	    #forced symetry
	    if($dup){
		$count += $num;
		_bin_it($bins,
			$cov,
			$num,
			1-$maf,
			\@obs);
	    }
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
    my $dup      = shift || 0;

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

		    #force symetry
		    if($dup){
			_bin_it($bins,
				$cov,
				$L,
				1-$maf,
				$arc);
		    }
		}
	    }
	    my $arc = $ARC{$label}{$bins}{$cov};

	    #add xenograft to model
	    if($xcov){
		my $mlabel .= "$label:xeno:$xcov:$m_aln";
		unless($ARC{$mlabel}{$bins}{$cov}){
		    my $sf = ($OPT->{use_ref}) ? 1-$m_aln : (1-$m_aln)**2;
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

			#force symetry
			if($dup){
			    _bin_it($bins,
				    $cov,
				    $L,
				    1-$maf,
				    $marc);
			}
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
    @p = (0, @p, $tot) if(!$OPT->{use_ref}); #reference removes non-hets
    @p = map {$_/$tot} @p;

    #create weights
    my %wi;
    if(!$OPT->{use_ref}){ #expect error derived SNVs
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
        if(!$OPT->{use_ref}){ #mouse peak not there when reference filtered
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
    my $cut = $CUT{$thr}{$pstring}{$cov};

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
	if($cns{$cn}){
	    my $cmp = $L <=> $cns{$cn}[1]{L} || $cns{$cn}[1]{rss} <=> $rss;
	    $cns{$cn} = [$m, $sc] if($cmp == 1);
	}
	else{
	    $cns{$cn} = [$m, $sc];
	}
	 
	#best overall
	if($best){
            my $cmp = $L <=> $best->[1]{L} || $best->[1]{rss} <=> $rss;
	    $cmp = 0 if($m eq '0:0' && $rfrac < $thr); #0 has no maf to match
            $best = [$m, $sc] if($cmp == 1);
        }
        else{
            $best = [$m, $sc] unless($m eq '0:0' && $rfrac < $thr); #0 has no maf to match;
        }

	next unless($xcov);

	my $xm = "$m:xeno";
	$sc = $r->{maf_score}{$xm};
	$rss = $sc->{rss};
	$L = $sc->{L};

	#best in cn category
	if($cns{$cn}){
	    my $cmp = $L <=> $cns{$cn}[1]{L} || $cns{$cn}[1]{rss} <=> $rss;
	    $cns{$cn} = [$xm, $sc] if($cmp == 1);
	}
	else{
	    $cns{$cn} = [$xm, $sc];
	}
	 
	#best overall
	if($best){
            my $cmp = $L <=> $best->[1]{L} || $best->[1]{rss} <=> $rss;
	    $cmp = 0 if($m eq '0:0' && $rfrac < $thr); #0 has no maf to match
            $best = [$xm, $sc];
        }
        else{
            $best = [$xm, $sc];
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
    my $tag = $OPT->{bam_files}[0];
    if(! keys %{$BAMS{$tag}}){
	if(@$bam_files == 1){
	    foreach my $chr (keys %CHRS){
		$BAMS{$tag}{$chr} = $bam_files->[0];
	    }
	}
	else{
	    foreach my $f (@$bam_files){
		next unless($f =~ /\.(chr[\dXY]+)\./ && $CHRS{$1});
		$BAMS{$tag}{$1} = $f;
	    }
	}
	foreach my $chr (keys %CHRS){
	    next if(!$BAMS{$tag}{$chr});
	    $BAMS{$tag}{$chr} = Bio::DB::Sam->new(-bam  => $BAMS{$tag}{$chr},
						  -autoindex => 1,
						  -fasta=> $OPT->{fasta});
	}
    }

    my ($chr, $start, $end) = ($r->{chr}, $r->{start}, $r->{end});
    
    my $bam = $BAMS{$tag}{$chr};
    die "ERROR: There is no bam file for $chr.\n" if(!$bam);
    my $segment = $bam->segment($chr,$start, $end);

    #process BAM files
    my $sid = $OPT->{sid};
    if($sid && !$RG{$tag}{$chr}){ #get read groups for filtering
	my @headers = grep {/^\@RG\t/} split(/\n/, $bam->header->text);
	foreach my $h (@headers){
	    my %rgs;
	    foreach (grep {/^(ID:.*|SM:.*)$/} split(/\t/, $h)){
		my @F = split(/\:/, $_);
		$rgs{$F[0]} = $F[1];
	    }
	    next unless( keys(%rgs) == 2);
	    push(@{$RG{$tag}{$chr}{$rgs{SM}}}, $rgs{ID});
	}
	$RG{$tag}{$chr} ||= {};
    }
    
    #get coverage from BAM
    my $cov_list;
    if(!$sid || !keys %{$RG{$tag}{$chr}}){
	my ($obj) = $segment->features(-type => 'coverage');
	$cov_list = $obj->coverage();
    }
    else{
	my ($SID) = grep {/$sid$/} keys %{$RG{$tag}{$chr}};
	my ($obj) = $segment->features(-type => 'coverage',
				       -filter => $RG{$tag}{$chr}{$SID});
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
    my $vcf;
    my $vcf_file  = $OPT->{vcf_file};
    if($VCF{$vcf_file}){
        $vcf = $VCF{$vcf_file};
    }
    else{
        $vcf = Vcf->new(file=>"$vcf_file");
        $vcf->parse_header();
        $VCF{$vcf_file} = $vcf;
    }
    my ($SID) = grep {/$sid$/} $vcf->get_samples() if($sid);
    $vcf->set_samples(include=>[$SID]) if($sid);

    #load BAM files
    my $bam_files = $OPT->{bam_files};
    my $btag = $OPT->{bam_files}[0];
    if(! keys %BAMS){
	if(@$bam_files == 1){
            foreach my $chr (keys %CHRS){
                $BAMS{$btag}{$chr} = $bam_files->[0];
            }
        }
        elsif(@$bam_files){
            foreach my $f (@$bam_files){
                next unless($f =~ /\.(chr[\dXY]+)\./ && $CHRS{$1});
                $BAMS{$btag}{$1} = $f;
            }
        }
        else{
            foreach my $chr (keys %CHRS){
                $BAMS{$btag}{$chr} = $dummy;
            }
        }
        foreach my $chr (keys %CHRS){
            next if(!$BAMS{$btag}{$chr});
            $BAMS{$btag}{$chr} = Bio::DB::Sam->new(-bam  => $BAMS{$btag}{$chr},
                                            -autoindex => 1,
                                            -fasta=> $OPT->{fasta});
        }
    }

    #temp (used to remove xenograft datapoints)
    my $rid = $OPT->{ref}{sid};
    my ($RID) = grep {/$rid$/} $vcf->get_samples() if($rid);
    $vcf->set_samples(include=>[$SID, $RID]) if($sid && $rid);

    #get MAF from VCF
    my %vcf_data = (chr => [],
                    start => [],
                    end => [],
                    pos => [],
                    maf=> [],
                    cov => []);
    $vcf->open(region=> "$chr:$cs-$ce");
    while(my $v = $vcf->next_data_hash()){
        my $sAD = ($SID) ? $v->{gtypes}{$SID}{AD} : $v->{gtypes}{AD};

        #temp (used to remove xenograft datapoints)
        my $rAD = $v->{gtypes}{$RID}{AD} if($RID);
	my $rcov;
        my $rmaf;
        if($RID){
            next unless($rAD && $rAD ne '.');
            my ($rc, $ac) = split(/,/, $rAD);
            $rcov = $rc+$ac;
            next unless($rcov >= $OPT->{maf_cov_filt});
            $rmaf = $ac/$rcov;
        }

        if($sAD && $sAD ne '.'){
            my ($rc, $ac) = split(/,/, $sAD);
            my $cov = $rc+$ac;

            if($cov > $OPT->{maf_cov_filt}){
                my $maf = $ac/$cov;
                $maf = 0 if(!$rmaf || $rmaf < $OPT->{maf_tail_filt}); #temp

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
        my $sam = $BAMS{$btag}{$chr};
        my $segment = $sam->segment($chr,$cs, $ce);

        #get BAM sample ID andf read groups
        if($sid && !$RG{$btag}{$chr}){ #get read groups for filtering
            my @headers = grep {/^\@RG\t/} split(/\n/, $sam->header->text);
            foreach my $h (@headers){
                my %tags;
                foreach (grep {/^(ID:.*|SM:.*)$/} split(/\t/, $h)){
                    my @F = split(/\:/, $_);
                    $tags{$F[0]} = $F[1];
                }
                next unless( keys(%tags) == 2);
                push(@{$RG{$btag}{$chr}{$tags{SM}}}, $tags{ID});
            }
            $RG{$btag}{$chr} ||= {};
        }

        #binned coverage (easier to see on plot)
        my $cov_list;
        my $bins = int($length/$OPT->{lfrag}) || 1;
        if(!$sid || !keys %{$RG{$btag}{$chr}}){
            my ($obj) = $segment->features(-type => "coverage:$bins");
            $cov_list = $obj->coverage();
        }
        else{
            my ($SID) = grep {/$sid$/} keys %{$RG{$btag}{$chr}};
            my ($obj) = $segment->features(-type => "coverage:$bins",
					   -filter => $RG{$btag}{$chr}{$SID});
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
    $R->run(qq`pdf("test.pdf",width=9,height=7.5,onefile=TRUE,paper='letter',pointsize=12)`);
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
                   plot(vcf\$pos,vcf\$maf,xlim=c($cs,$ce),pch=20,cex=1,ylim=c(0,1),col="green");abline(v=c($start,$end),col=2);
                   par(opar);`);
    }
    else{
        #plot with just VCF                                                                                                                                          
        $R->run(qq`opar <- par(mfrow=c(2,1));
                   mainstr <- sprintf("chr\%d:\%4.2f-\%4.2f (kb),  length=\%4.2f(kb) ,   cn=\%d",$chrnum,$start/1e3,$end/1e3,($length)/1e3,$cn);
                   plot(vcf\$pos,vcf\$cov,xlim=c($cs,$ce),pch=20,cex=1,main=mainstr);abline(v=c($start,$end),col=2);
                   plot(vcf\$pos,vcf\$maf,xlim=c($cs,$ce),pch=20,cex=1,ylim=c(0,1),col="green");abline(v=c($start,$end),col=2);
                   par(opar);`);

    }
    $R->run(q`dev.off()`); #close file and write
    $R->stop(); #call explictly because it is not getting destroyed
}

sub target_v {
    my $t_cov = shift;
    my $bc    = shift;
    my $res   = shift;
    my $OPT   = shift;

    my $cfrac  = $OPT->{cfrac};
    my $xcov   = $OPT->{xcov};
    my $models = $OPT->{models};
    my $m_aln  = $OPT->{m_aln};
    my $rfrac  = 1 - $cfrac;

    my $V = 0;
    my $count = 0;
    foreach my $r (@$res){
        my $cn = $r->{final}{cn};
        next if($cn eq '!');
        next if($cn eq '.');
        next if($cn eq 'N');
        next if($rfrac && $cn == 0);

        my $e_cov = $cn*$bc+$xcov*$m_aln+(2*$bc*$rfrac)/($cfrac);

        $V += ($e_cov-$r->{cov_median})**2 * $r->{q_length};
        $count += $r->{q_length};
    }
    $V /= $count if($count);

    return $V;
}

#returns a list of target frequencies for a given model
sub get_f {
    my $mod = shift;
    my $cfrac = shift || 1;
    my $rfrac = 1 - $cfrac;

    return '!' if($mod eq '!');

    my $tot = 0;
    my @p = grep {$_ ne 'xeno'} split(/:/, $mod);
    map {$_ = $rfrac + ($cfrac)*$_; $tot += $_; $_} @p;

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
sub expected_segment_median{
    my $res = shift;
    my $OPT = shift;

    return 100 if($res->[0]->{is_generic});

    my %medians;
    foreach my $r (@$res){
	my $t = $r->{cov_median};

	#remove xenograft effect
	if($r->{xeno} && $OPT->{xeno}{ccor_f}){
	    $t -= $r->{xeno}{cov_median}*$OPT->{xeno}{ccor_f};
	}

	#use corrected coverage when available
	if($r->{ref} && $r->{ref}{ccor_f}){
	    $t = $t/$r->{ref}{ccor_f}; #round to nearest
	}

	$t = int($t+0.5); #always round to nearest
	next unless($t > 0);

        my $points = $r->{q_length}/$OPT->{lfrag};
        $points = $r->{cov_count} if($r->{cov_count} < $points);
	$medians{$t} += $points;
    }

    my $sum = 0;
    my $count = 0;
    my $ker = Statistics::KernelEstimation->new();
    foreach my $t (keys %medians){
	next unless($t);
	$ker->add_data($t, $medians{$t});
	$sum += $t*$medians{$t};
	$count += $medians{$t};   
    }
    
    my $w = ($ker->default_bandwidth() < 2) ? 1 : $ker->default_bandwidth();
    $w = 1 if($w < 1);
    my $ave = $sum/$count;
    my ($min, $max) = $ker->extended_range();
    $min = 0 if($min < 0);
    $max = $ave*4 if($ave*4 < $max);
    my $step = ($max-$min)/100;
    if($step > 0.1){
	($min, $max) = (int($min), int($max));
	$step = 0.1;
    }
    my $t_cov = [0,0]; #expected segment mean coverage
    $sum = 0;
    $count = 0;
    for(my $x = $min; $x <= $max; $x += $step) {
	my $y = $ker->pdf($x, $w);
	$sum += 0.1 * $y; #lets me know I've already processed most data
	print STDERR "$x\t$y\n" if($count++ % 50 == 0);
	$t_cov = [$x, $y] if($y >= $t_cov->[1]);
	last if($x == $max || $sum >= 0.9);
    }
    $t_cov = $t_cov->[0];
    
    return $t_cov;
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

#get expect by matching the best MAF for segments with that median coverage
sub fit_expected_coverage {
    my $res    = shift;
    my $OPT    = shift;

    my $models = $OPT->{models};

    if($res->[0]->{is_generic}){
	return (100, 50);
    }

    #get median of most frequent segment
    my $t_cov;
    if($OPT->{user_t_cov}){
        $t_cov = $OPT->{user_t_cov};
    }
    else{
        $t_cov = expected_segment_median($res, $OPT);
    }

    if($OPT->{user_base_cov}){
	return ($t_cov, $OPT->{user_base_cov});
    }

    #test between coverages for copy numbers 1-4
    my @fits = _coverage_fit($res, $t_cov, $models, $OPT, 50);
    my $fit2 = _fit_targets($res, $t_cov, $models, $OPT); #match t_cov
    return ($t_cov, $fits[0]) if(@fits == 1);
    return ($t_cov, $fit2) if(!@fits && $fit2);
    return ($t_cov, $t_cov/2) if(!@fits && !$fit2);

    #if the top fit is really close, just return it
    my @close = map {[abs($_-$fit2), $_]} @fits;
    my $top = shift @close;
    @close = sort {$a->[1] <=> $b->[1]} @close;
    return ($t_cov, $fits[0]) if($top->[0] <= $close[0][0]);

    #check if supported by other fitting data
    my $cn2 = cn_for_cov($t_cov, $fit2, $OPT);
    if($cn2 != 1 || $fit2 > 5){
	foreach my $fit1 (@fits){	    
	    my $cn1 = cn_for_cov($t_cov, $fit1, $OPT);
	    
	    return ($t_cov, $fit1) if($cn1 == $cn2);
	    return ($t_cov, $fit1) if($cn1 == $cn2/2 && $cn1 != 1);
	}
    }

    #throw out cn=4 (only allowed if fit1 and fit2 don't match)
    if(@fits > 1){
        @fits = grep {cn_for_cov($t_cov, $_, $OPT) < 4} @fits;
    }

    if(@fits && $fit2 && ($cn2 != 1 || $fit2 > 5)){
	my ($closest)= sort {abs($fit2-$a) <=> abs($fit2-$b)} @fits;
	return ($t_cov, $closest);
    }
    elsif(@fits){
	return ($t_cov, $fits[0]);
    }
    else{
	return ($t_cov, $fit2);
    }
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
	next if(!$r->{final});
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

    return ($p_ave, $OPT->{user_ploidy}) if($OPT->{user_ploidy});
    return ($p_ave, $p_mode->[0]);
}

#labels segment results as being somatic or not (requires segment with germline reference)
sub label_somatic {
    my $res = shift;
    my $OPT = shift;

    foreach my $r (@$res){
	$r->{is_somatic} = 0; #assume everything is germline
	my $cn = $r->{final}{cn};
        my $rcn = $r->{ref}{cn};

	next if($cn eq '!');
	next if($cn eq 'N');
	next if($rcn eq 'N'); #why?
        next if($cn eq '.');
        next if($rcn eq 0); #skip if matches reference after correction

        #can only be LOH if ref is not
	if($r->{q_length} > 100000 && $r->{final}{loh} && !$r->{ref}{loh}){
            $r->{is_somatic} = 1;
            next;
        }

	#ploidy is only meaningful if it is calculable
	if($OPT->{base_cov} > 5){
	    next if($OPT->{ref}{ploidy} - $cn == $OPT->{ref}{ploidy} - $rcn);
	    next if($OPT->{ploidy} - $cn == $OPT->{ref}{ploidy} - $rcn); 
	    next if($cn == $rcn);
	    next if($cn == $OPT->{ploidy});
	    next if($cn == $OPT->{ref}{ploidy});
	}

	#cases where correction leaves wide uncertainty around ploidy
	my ($cn0, $cn1) = sort {$a <=> $b} ($r->{ucor}{cn}, $r->{ccor}{cn});
	next if($cn0 <= $OPT->{ploidy} && $OPT->{ploidy} <= $cn1);
	next if($cn0 <= $OPT->{ref}{ploidy} && $OPT->{ref}{ploidy} <= $cn1);
	
	#must be significantly different from ploidy levels
	next if(!sig_somatic($r, $OPT));

	# if just labeling loss and gain
	if($OPT->{base_cov} <= 5 && $r->{ref} && sig_somatic($r->{ref}, ref_param($OPT))){
	    next if($r->{final}{cn} < $OPT->{ploidy} && $r->{ref}{cn} < $OPT->{ref}{ploidy});
	    next if($r->{final}{cn} > $OPT->{ploidy} && $r->{ref}{cn} > $OPT->{ref}{ploidy});
	}

        #be more strict on what to accept for short segments
        #if($r->{q_length} < 10000){
        #    next if(abs(($OPT->{ploidy} - $cn) - ($OPT->{ref}{ploidy} - $rcn)) <= 1);
        #    next if(abs(($OPT->{ref}{ploidy} - $cn) - ($OPT->{ref}{ploidy} - $rcn)) <= 1);
        #}

        #only label segements that don't have alternate explanations as somatic
        $r->{is_somatic} = 1;
    }
}

#significant difference between expected ploidy levels and assigned copy number
sub sig_somatic {
    my $r   = shift;
    my $OPT = shift;

    my $cfrac = $OPT->{cfrac} || 1;
    my $rfrac = 1 - $cfrac;
    my $bc  = $OPT->{chr_expects}{$r->{chr}} || $OPT->{chr_expects}{_ALL};
    my $rp = ($OPT->{ref} && $OPT->{ref}{ploidy}) ? $OPT->{ref}{ploidy} : 2;

    my $spc = $OPT->{ploidy}*$bc+($rp*$bc*$rfrac)/($cfrac);
    my $rpc = ($OPT->{chr_expects}{_ALL} > 5) ?
	$rp*$bc+($rp*$bc*$rfrac)/($cfrac) : $spc; #ignored if sample ploidy not calculable

    #calculate length adusted standard deviation
    my $r_med = $r->{cov_median};
    my $r_n = $r->{q_length}/$OPT->{lfrag};
    $r_n = $r->{cov_count} if($r->{cov_count} < $r_n);
    $r_n = 1 if($r_n < 1);
    my $r_sd = sqrt($r_med/$r_n) * 3.890592 + 0.5; #e-4 threshold

    #adjust for uncertainty in target frerquency
    my $t_sd = sqrt($spc/(0.001*$OPT->{n50}));
    my $sd = sqrt($r_sd**2+$t_sd**2);

    if(abs($r_med-$rpc) <=  $sd || abs($r_med-$spc) <=  $sd){
	return 0;
    }

    if($r->{ref} && $r->{ref}{ccor_f}){
	$r_med /= $r->{ref}{ccor_f};
	$r_sd /= $r->{ref}{ccor_f};
	$sd = sqrt($r_sd**2+$t_sd**2);
	if(abs($r_med-$rpc) <=  $sd || abs($r_med-$spc) <=  $sd){
	    return 0;
	}
    }

    return 1;
}

#identify and then merge similar segments in result set
sub merge_segments {
    my $res = shift;
    my $OPT = shift;

    $res = _merge_special($res, $OPT); #merge special CN types first

    my $count = @$res;
    my ($sets, $hold) = _split_on_gaps($res, $OPT); #split into reonable groups

    #make stack for threads to draw from
    my @stack :shared;
    mkdir($OPT->{outdir}."/tmp");
    foreach my $tag (keys %$sets) {
	push(@stack, $tag);
	(undef, my $tfile) = tempfile(DIR => $OPT->{outdir}."/tmp", CLEANUP => 1);
	Storable::store($sets->{$tag}, $tfile);
	$sets->{$tag} = $tfile;
    }
    undef @$res; #for memory optimization

    #split off threads
    my @threads;
    prefork(); #prepare for forking
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'scalar'},
                                       \&_merge_thread,
                                       \@stack, #shared
				       $sets, 
				       $OPT));
    }

    #collect results
    my $merged = _merge_thread(\@stack,$sets, $OPT);
    foreach my $t (@threads){
	my $new = $t->join();
	push(@$merged, @$new)
    }
    push(@$merged, @$hold); #add back uncallable segments

    #merge special classes across the lack of data intervals
    $merged = _merge_around_special($merged, $OPT) if($OPT->{smooth});
    $merged = _merge_uncallable($merged, $OPT) if($OPT->{smooth});

    #recall the copy numbers
    if($count != @$merged){
	recall_cn($merged, $OPT);
    }
    
    return $merged;
}

sub _split_on_gaps {
    my $res = shift;
    my $OPT = shift;

    my %sets;
    my @hold;
    foreach my $r (@$res){
	#add big mask sections to the gap list
	if($r->{final}{cn} eq 'N' && $r->{mask_count} > 100000){
	    my $gaps = $OPT->{GAPS}{$r->{chr}};
	    my $rep = grep {$_->[0] == $r->{start} && $_->[0] == $r->{end}} @{$gaps};
	    push(@{$gaps}, [$r->{start}, $r->{end}]) unless($rep);
	    push(@hold, $r);
	}
	elsif($r->{final}{cn} eq '!'){ #add uncallable sections as gaps as well
	    my $gaps = $OPT->{GAPS}{$r->{chr}};
	    my $rep = grep {$_->[0] == $r->{start} && $_->[0] == $r->{end}} @{$gaps};
	    push(@{$gaps}, [$r->{start}, $r->{end}]) unless($rep);
	    push(@hold, $r);
	}
	else{ #otherwise add to chromosome set
	    push(@{$sets{$r->{chr}}}, $r);
	}
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
	    if($r->{start} <= $c && $c <= $r->{end} && $r->{final}{cn} !~ /N|\!/){
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

    return (\%sets, \@hold);
}

#merge action performed by each thread
#access to results should be copy on write (efficient)
sub _merge_thread {
    my $stack = shift; #reference to shared list
    my $set   = shift; #reference to cromosome divided sets
    my $OPT   = shift;

    my $keepers = [];
    while(my $tag = shift @$stack){
	my $file = $set->{$tag};
	my $res = (ref($file)) ? $file : Storable::retrieve($file);
	@$res = sort {$a->{start} <=> $b->{start}} @$res;

	while(1){
	    my $start = @$res;
	    my $use_ref = 1 if($OPT->{use_ref});
	    my $smooth = 1 if($OPT->{smooth});
	    while(1){
		my $start1 = @$res;

		#merge across stdev
		$res = _merge_stdev($res, $OPT, 'ucor');
		next if($start1 != @$res);
    
		#merge across corrected stdev
		$res = _merge_stdev($res, $OPT, 'ccor');
                next if($start1 != @$res);
		
		#merge based on same copy number
		$res = _merge_same_cn($res, $OPT, 'ucor');
		next if($start1 != @$res);

		#merge based on ref corrected cn
		$res = _merge_same_cn($res, $OPT, 'ccor');
                next if($start1 != @$res);

		#merge based on ref corrected cn
		$res = _merge_same_cn($res, $OPT, 'final');
                next if($start1 != @$res);

		#merge across regions where data was insufficient to make a call
		$res = _merge_around_missing($res, $OPT) if($smooth);
                next if($start1 != @$res);

		#merge systematic rises and falls
		$res = _merge_systematic_shifts($res, $OPT) if($smooth);
                next if($start1 != @$res);

		last;
	    }

	    next if($start != @$res);

	    last;
	}

	push(@$keepers, @$res);
    }

    return $keepers;
}

#only merges special types used for splitting
sub _merge_special {
    my $res = shift;
    my $OPT = shift;

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
	
	if($r0->{final}{cn} ne $act->{final}{cn}){
	    push(@keepers, $act);
	    $act = $r0;
	}
	elsif($r0->{chr} ne $act->{chr}){
            push(@keepers, $act);
            $act = $r0;
	}
	elsif($r0->{final}{cn} !~ /N|\!/){
            push(@keepers, $act);
            $act = $r0;
	}
	else{
	    $act = _merge($act, $r0, $OPT);
	    $act = _polish_merge($act, $OPT);
	}
    }
    push(@keepers, $act) if($act);

    return \@keepers;
}

sub _merge_around_special {
    my $res = shift;
    my $OPT = shift;

    #split by chromosome
    my %sets;
    foreach my $r (@$res){
	push(@{$sets{$r->{chr}}}, $r)
    }

    #now process
    my @keepers;
    foreach my $s (values %sets){
	@$s = sort {$a->{start} <=> $b->{start}} @$s;
	#remove segments that appear to be lack of information rather than true CN
	my $to_merge = [];
	my @hold;
	for(my $i = 0; $i < @$s; $i++){
	    my $ri  = $s->[$i];
	    if($ri->{final}{cn} eq '.' && $ri->{length} < 100000){
		push(@hold, $ri);
	    }
	    elsif($ri->{final}{cn} eq 'N' && $ri->{length} < 100000){
		push(@hold, $ri);
	    }
	    else{
		push(@$to_merge, $ri);
	    }
	}

	if(@$s == @$to_merge){ #nothing changed
	    push(@keepers, @$s);
	    next;
	}

	#sanity check
	die "Logic error for _merge_somatic\n" if(@$s != @hold + @$to_merge);
	
	#run standard merge with newly introduced gaps
	my $chr = $s->[0]{chr};
	$to_merge = _merge_special($to_merge, $OPT);

	if(@$s == @$to_merge + @hold){ #nothing changed
	    push(@keepers, @$s);
	    next;
	}
	
	#now see which removed segments should be added back
	@$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
	@hold = sort {$a->{start} <=> $b->{start}} @hold; 
	push(@keepers, @$to_merge);
	my $pos = 0;
	foreach my $h (@hold){
	    my $bad;
	    for(my $i = $pos; $i < @$to_merge; $i++){
		my $k = $to_merge->[$i];
		if($k->{end} < $h->{start}){
		    $pos = $i;
		    next;
		}
		
		if($k->{start} <= $h->{end}){
		    $bad++;
		    
		    #merge in if there is meaningful data
		    $k = _merge($k, $h, $OPT) if($h->{ucor}{cn} > 1);
		}
		last if($k->{start} > $h->{end}); #doesn't overlap;
	    }	
	    push(@keepers, $h) unless($bad);
	}
    }

    foreach my $r (@keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return \@keepers;
}

sub _merge_exome{
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{gene} cmp $b->{gene} || $a->{start} <=> $b->{start}} @$s;

    my $keepers = [];
    my $act;
    for(my $i = 0; $i < @$s; $i++){
        my $r0  = $s->[$i]; #merger

        if(!$act){
            $act = $r0;
            next;
        }

        if($act->{gene} ne $r0->{gene}){
            push(@$keepers, $act);
            $act = $r0;
            next;
        }

        elsif ($act->{gene} eq $r0->{gene}){
            $act = _merge($act, $r0, $OPT);
            next;
        }
    }
    push(@$keepers, $act) if($act);

    foreach my $r (@$keepers){
	$r = add_maf_stats($r, $OPT) 
    }

    return $keepers;
}

#merges if cn is same and allele is same or empty
sub _merge_same_cn {
    my $s = shift;
    my $OPT = shift;
    my $tag = shift || '';

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
	my $match;
	$match = 1 if($act->{$tag}{cn} eq $r0->{$tag}{cn});
	$match = 0 if($act->{final}{cn} eq '!' && $r0->{final}{cn} ne '!');
	$match = 0 if($r0->{final}{cn} eq '!' && $act->{final}{cn} ne '!');
	$match = 0 if($act->{final}{cn} eq 'N' && $r0->{final}{cn} ne 'N');
	$match = 0 if($r0->{final}{cn} eq 'N' && $act->{final}{cn} ne 'N');
	$match = 0 if($tag eq 'final' && $r0->{final}{cn} eq '.'); #matching '.' may be reference

	#make sure alleles match
	$match = _alleles_match($act, $r0, $OPT) if($match);

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
		$new->{ccor} = $big->{ccor};
		$new->{ucor} = $big->{ucor};
		$new->{final} = $big->{final};
		$new->{maf_score} = $big->{maf_score}; #approximation
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
sub _merge_uncallable {
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{chr} cmp $b->{chr} || $a->{start} <=> $b->{start}} @$s;

    my $keepers = [];
    my $act; #active working segment
    for(my $i = 0; $i < @$s; $i++){
	my $rm = $s->[$i-1] if($i-1 >= 0); #merger minus 1 (ancestral)
	my $r0  = $s->[$i]; #this is merger

	if(!$act){
	    $act = $r0;
	    next;
	}

	if($act->{chr} ne $r0->{chr}){
	    push(@$keepers, $act);
	    $act = $r0;
	    next;
	}

	if($act->{final}{cn} ne '.' || $act->{ucor}{cn} > 1){
            push(@$keepers, $act);
            $act = $r0;
            next;
        }

	if($r0->{final}{cn} ne '.' || $r0->{ucor}{cn} > 1){
            push(@$keepers, $act);
            $act = $r0;
            next;
        }

	#make new but recompute CN later
	my $new = _merge($act, $r0, $OPT);
	$new->{ccor} = $act->{ccor};
	$new->{ucor} = $act->{ucor};
	$new->{final} = $act->{final};
	$act = $new;
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

    if($r->{ref}){
	my $param = ref_param($OPT);
	my $rc = $param->{chr_expects}{$r->{chr}};
	_assign_cn($r->{ref}, $rc, $param);
	delete($r->{ref}{_merged});
	my $obs = $r->{ref}{cov_median}/$param->{ploidy};
	$r->{ref}{ccor_f} = ($r->{ref}{ucor}{cn}) ? $obs/$rc : 1;
    }

    my $rc = $OPT->{chr_expects}{$r->{chr}};
    _assign_cn($r, $rc, $OPT);
    delete($r->{_merged});

    label_somatic([$r], $OPT);

    return $r;
}

#merges short segments within 1 SD od the ciurrent level
sub _merge_stdev {
    my $s = shift;
    my $OPT = shift;
    my $tag = shift || '';

    return $s if(@$s <= 1);

    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};

    my @all = sort {$a->{chr} cmp $b->{chr} || $a->{start} <=> $b->{start}} @$s;
    while(1){
	my @order = sort {$all[$a]{q_length} <=> $all[$b]{q_length}} (0..$#all);

	my $min = $OPT->{CHRS}{$all[0]{chr}};
	my $start = @all; #how many I started with
	while(defined(my $i = shift @order)){ #works from smallest to largest
	    my $r0 = $all[$i]; #this is merger
	    next if(!$all[$i]);
	    next if($r0->{final}{cn} eq '!' || $r0->{final}{cn} eq 'N');
	    
	    my $h = $i-1;
	    $h-- while($h > 0 && !$all[$h]);
	    my $rm = $all[$h] if($h >= 0); #merger minus 1

	    my $j = $i+1;
	    $j++ while($j < $#all && !$all[$j]);
	    my $rp = $all[$j] if($j <= $#all); #merger plus 1

	    if($rm){
		die "Logic error\n" if($rm->{chr} ne $r0->{chr});
		undef $rm if($rm->{final}{cn} eq '!' ||
			     $rm->{final}{cn} eq 'N' ||
			     abs($r0->{start} - $rm->{end}) >= 100000);

	    }
	    if($rp){
		die "Logic error\n" if($rp->{chr} ne $r0->{chr});
		undef $rp if($rp->{final}{cn} eq '!' ||
			     $rp->{final}{cn} eq 'N' ||
			     abs($rp->{start} - $r0->{end}) >= 100000);
	    }
	    next if(!$rm && ! $rp);
	    
	    #calculate length adusted standard deviation
	    my $r_med = $r0->{cov_median};
	    my $r_n = $r0->{q_length}/$OPT->{lfrag};
	    $r_n = $r0->{cov_count} if($r0->{cov_count} < $r_n);
	    $r_n = 1 if($r_n < 1);
	    my $r_sd = sqrt($r_med/$r_n) * 3.890592 + 0.5; #e-4 threshold
	    if($tag eq 'ccor'){
		my $rx_adj = ($r0->{xeno}) ? $r0->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		$r_med = $r_med - $rx_adj;
		$r_med = 0 if($r_med <= 0);
		$r_med /= $r0->{ref}{ccor_f};
	    }
	    
	    my $dist_m;
	    if($rm){
		my $m_med = $rm->{cov_median};
		my $m_n = $rm->{q_length}/$OPT->{lfrag};
		$m_n = $rm->{cov_count} if($rm->{cov_count} < $m_n);
		$m_n = 1 if($m_n < 1);
		if($tag eq 'ccor'){
		    my $mx_adj = ($rm->{xeno}) ? $rm->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		    $m_med = $m_med - $mx_adj;
		    $m_med = 0 if($m_med <= 0);
		    $m_med /= $rm->{ref}{ccor_f};
		}
		
		$dist_m = abs($r_med - $m_med);
		undef $dist_m if($dist_m > $r_sd);
	    }
	    
	    my $dist_p;
	    if($rp){
		my $p_med = $rp->{cov_median};
		my $p_n = $rp->{q_length}/$OPT->{lfrag};
		$p_n = $rp->{cov_count} if($rp->{cov_count} < $p_n);
		$p_n = 1 if($p_n < 1);
		if($tag eq 'ccor'){
		    my $px_adj = ($rp->{xeno}) ? $rp->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		    $p_med = $p_med - $px_adj;
		    $p_med = 0 if($p_med <= 0);
		    $p_med /= $rp->{ref}{ccor_f};
		}
		
		$dist_p = abs($r_med - $p_med);
		undef $dist_p if($dist_p > $r_sd);
	    }

	    undef $dist_m if(defined($dist_m) && ! _alleles_match($rm, $r0, $OPT));
	    undef $dist_p if(defined($dist_p) && ! _alleles_match($rp, $r0, $OPT));

	    if(defined($dist_m) && (!defined($dist_p) || $dist_m <= $dist_p)){
		my $new = _merge($r0, $rm, $OPT);

		#decide when to recompute (makes merging faster)
		my ($big, $small) =  ($r0->{q_length} >= $rm->{q_length}) ? ($r0, $rm) : ($rm, $r0);
		my $size = $big->{_merged} || $big->{q_length};
		my $diff = $big->{q_length} - $size;
		if($new->{q_length} > 30000 && $big->{q_length} <= 30000){ #maf threshold recompute
		    $new = _polish_merge($new, $OPT);
		}
		elsif($size/10 > $diff + $small->{q_length}){ #if much bigger don't recompute yet
		    $new->{_merged} = $size; #keeps track of size used to compute maf so far
		    $new->{ccor} = $big->{ccor};
		    $new->{ucor} = $big->{ucor};
		    $new->{final} = $big->{final};
		    $new->{maf_score} = $big->{maf_score}; #approximation
		}
		else{
		    $new = _polish_merge($new, $OPT);
		}

		$all[$i] = $new;
		$all[$h] = undef;
		push(@order, $i); #end of list

		#see if I need to resort the @order list
		$min = $new->{q_length} if($new->{q_length} < $min);
		my $prox = $order[0];
		while(!$all[$prox]){
		    shift @order; #throw away
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
		if($new->{q_length} > 30000 && $big->{q_length} <= 30000){ #maf threshold recompute
		    $new = _polish_merge($new, $OPT);
		}
		elsif($size/10 > $diff + $small->{q_length}){ #if much bigger don't recompute yet
		    $new->{_merged} = $size; #keeps track of size used to compute maf so far
		    $new->{ccor} = $big->{ccor};
		    $new->{ucor} = $big->{ucor};
		    $new->{final} = $big->{final};
		    $new->{maf_score} = $big->{maf_score}; #approximation
		}
		else{
		    $new = _polish_merge($new, $OPT);
		}

		$all[$i] = $new;
		$all[$j] = undef;
		push(@order, $i); #end of list

		#see if I need to resort the @order list
		$min = $new->{q_length} if($new->{q_length} < $min);
		my $prox = $order[0];
		while(!$all[$prox]){
		    shift @order; #throw away
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
        $r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return \@all;

}

sub _alleles_match {
    my $r1 = shift;
    my $r2 = shift;
    my $OPT = shift;

    return 1 if($r1->{is_generic});

    my $a1 = $r1->{final}{allele};
    my $a2 = $r2->{final}{allele};
    _mark_best_alleles($r1, $OPT);
    _mark_best_alleles($r2, $OPT);

    my $match;
    if($r1->{final}{cn} eq 'N' && $r2->{final}{cn} eq 'N'){
	$match = 1;
    }
    elsif($r1->{final}{cn} eq 'N'  || $r2->{final}{cn} eq 'N'){
	$match = 0;
    }
    elsif($r1->{final}{cn} eq '!' && $r2->{final}{cn} eq '!'){
	$match = 1;
    }
    elsif($r1->{final}{cn} eq '!'  || $r2->{final}{cn} eq '!'){
	$match = 0;
    }
    elsif($r1->{maf_score}{best}{all}[0] eq $r1->{maf_score}{best}{all}[0]){ #sure why not
	$match = 1;
    }
    elsif($r1->{q_length} < 30000 || $r2->{q_length} < 30000){
	$match = 1;
    }
    elsif($r1->{maf_count} < 15 || $r2->{maf_count} < 15){
	$match = 1;
    }
    elsif($r1->{final}{allele} eq '.' || $r2->{final}{allele} eq '.'){ #no allele, so I can't say they don't match
	$match = 1;
    }
    elsif($r1->{final}{cn} <= 1 && $r2->{final}{cn} <= 1){
	$match = 1;
    }
    elsif($r1->{final}{cn} >= 7  || $r2->{final}{cn} >= 7){ #can't say the alleles don't match
	$match = 1;
    }
    elsif($r2->{ref}{loh} || $r1->{ref}{loh}){ #lack of info going back to ref
	$match = 1;
    }
    else{ #now compare the alleles to one another if cn jumping
	$match = 1;
	my $alt1 = $r2->{final}{allele};
	my $alt2 = $r1->{final}{allele};
	my $L1 = $r1->{final}{allele_L};
	my $L2 = $r2->{final}{allele_L};
	my $altL1 = $r1->{maf_score}{$alt1}{L};
	my $altL2 = $r2->{maf_score}{$alt2}{L};
	
	#compare relative liklihoods
	if($altL2/$L2 > 0.1){
	    $match = 1;
	}
	elsif($altL1/$L1 > 0.1){
	    $match = 1;
	}
	else{
	    $match = 0;
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

	#get new correction factor
	if(defined($OPT->{ref}{ploidy})){
	    my $exp = $OPT->{ref}{chr_expects}{$ref->{chr}};
	    my $obs = $ref->{cov_median}/$OPT->{ref}{ploidy};
	    $ref->{ccor_f} = ($ref->{ucor}{cn}) ? $obs/$exp : 1;
	}
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
       $res->{q_length}/$res->{cov_count} > 1500 || $res->{is_xeno}
        ){
	$need_bam = 1;
    }
    
    #add BAM coverage
    if(@{$OPT->{bam_files}}){
	#make use of existing info where possible
	if($res->{is_generic}){
	    add_generic_cov($r1, $OPT);
	    add_generic_cov($r2, $OPT);
	}
	else{
	    add_bam_cov($r1, $OPT);
	    add_bam_cov($r2, $OPT);
	}
	
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

    #add merged reference fragment back to segment
    $res->{ref} = $ref if($ref);
    $res->{xeno} = $xeno if($xeno);

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

        if($act->{final}{cn} eq '!' || $act->{final}{cn} eq 'N' || $act->{final}{cn} eq '.'){
            push(@$keepers, $act);
            $act = $r;
	    next;
        }

        if($r->{final}{cn} eq '!' || $r->{final}{cn} eq 'N' || $r->{final}{cn} eq '.'){
            push(@$keepers, $act);
            $act = $r;
	    next;
        }
	
	#make sure copy numbers match
        my $match;
        $match = 1 if($act->{final}{cn} - $r->{final}{cn} == $act->{ref}{cn} - $r->{ref}{cn});
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

sub _merge_around_missing {
    my $s = shift;
    my $OPT = shift;

    @$s = sort {$a->{start} <=> $b->{start}} @$s;

    #remove segments that appear to be lack of information rather than true CN
    my $to_merge = [];
    my @hold;
    for(my $i = 0; $i < @$s; $i++){
	my $ri  = $s->[$i];
	if($ri->{final}{cn} eq '.' && $ri->{length} < 100000){
	    push(@hold, $ri);
	}
	elsif($ri->{final}{cn} eq 'N' && $ri->{length} < 100000){
	    push(@hold, $ri);
	}
	else{
	    push(@$to_merge, $ri);
	}
    }

    return $s if(@$s == @$to_merge); #nothing changed

    #sanity check
    die "Logic error for _merge_somatic\n" if(@$s != @hold + @$to_merge);

    #run standard merge with newly introduced gaps
    my $chr = $s->[0]{chr};
    $to_merge = _merge_thread([$chr], {$chr => $to_merge}, $OPT);

    return $s if(@$s == @$to_merge + @hold); #nothing changed

    #now see which removed segments should be added back
    @$to_merge = sort {$a->{start} <=> $b->{start}} @$to_merge;
    @hold = sort {$a->{start} <=> $b->{start}} @hold; 
    my @keepers;
    push(@keepers, @$to_merge);
    my $pos = 0;
    foreach my $h (@hold){
	my $bad;
	for(my $i = $pos; $i < @$to_merge; $i++){
	    my $k = $to_merge->[$i];
	    if($k->{end} < $h->{start}){
		$pos = $i;
		next;
	    }

	    if($k->{start} <= $h->{end}){
		$bad++;
		
		#merge in if there is meaningful data
		$k = _merge($k, $h, $OPT) if($h->{ucor}{cn} > 1);
	    }
	    last if($k->{start} > $h->{end}); #doesn't overlap;
	}	
	push(@keepers, $h) unless($bad);
    }

    foreach my $r (@keepers){
	$r = _polish_merge($r, $OPT) if($r->{_merged});
    }

    return \@keepers;
}

sub maybe_loh{
    my $r = shift;
    my $OPT = shift;

    return 0 if($r->{is_generic});

    _mark_best_alleles($r, $OPT);
    foreach my $key (keys %{$r->{maf_score}{best}}){
	next if($key eq '0' || $key eq '1');
	return 1 if($r->{maf_score}{best}{$key}[0] =~ /^0/);
    }

    return 0;
}

sub write_gvf {
    my $res = shift;
    my $file = shift;

    #divide by chromosome
    my %sets;
    foreach my $r (@$res){
	push(@{$sets{$r->{chr}}}, $r);
    }    

    #sort by position
    foreach my $s (values %sets){
	@$s = sort {$a->{start} <=> $b->{start}} @$s;
    }

    my $OUT;
    if($file){
	open($OUT, ">$file")
    }
    else{
	open($OUT, ">&STDOUT")
    }

    print $OUT "#xenograft_coverage=".$OPT{xcov}."\n";
    print $OUT "#celularity=".$OPT{cfrac}."\n";
    print $OUT "#target_expect=".$OPT{t_cov}."\n";
    print $OUT "#base_coverage=".$OPT{base_cov}."\n";
    print $OUT "#ploidy=".$OPT{ploidy}."\n";

    foreach my $k (sort {_chrom($a) <=> _chrom($b)} keys %sets){
	my $set = $sets{$k};
	foreach my $r (@$set){
	    print $OUT _gvf_line($r);
	}
    }
    close($OUT);
}


#formats one line of GVF from a result
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

    my $att = 'ID=CNA:'.$r->{chr}.':'.$r->{start}.'-'.$r->{end};
    my $add = '';

    #add uncorrected cn
    my $u_cn     = $r->{ucor}{cn};
    my $u_allele = $r->{ucor}{allele};
    my $u_L      = $r->{ucor}{L};
    my $u_aL     = $r->{ucor}{allele_L};
    my $u_sc     = $r->{maf_score}{best}{$u_cn}[1];

    $add .= ";ucor_cn=$u_cn";
    $add .= ";ucor_cn_L=".sprintf('%.2f', $u_L);
    if($u_allele ne '.' && $u_aL >= 0.9){
	$add .= ";ucor_cn_allele=$u_allele";
	$add .= ";ucor_cn_allele_L=$u_aL";
    }
    
    #add corrected cn
    my $c_cn     = $r->{ccor}{cn};
    my $c_allele = $r->{ccor}{allele};
    my $c_L      = $r->{ccor}{L};
    my $c_aL     = $r->{ccor}{allele_L};
    my $c_sc     = $r->{maf_score}{best}{$c_cn}[1];
    
    $add .= ";ccor_cn=$c_cn";
    $add .= ";ccor_cn_L=".sprintf('%.2f', $c_L);
    if($c_allele ne '.' && $c_aL >= 0.9){
	$add .= ";ccor_cn_allele=$c_allele";
	$add .= ";ccor_cn_allele_L=$c_aL";
    }
    
    #add ref cn
    my $r_cn     = $r->{ref}{cn};
    my $r_allele = $r->{ref}{allele};
    my $r_L      = $r->{ref}{L};
    my $r_aL     = $r->{ref}{allele_L};
    my $r_sc     = $r->{ref}{maf_score}{best}{$r_cn}[1];
    
    $add .= ";ref_cn=$r_cn";
    $add .= ";ref_cn_L=".sprintf('%.2f', $r_L);
    if($r_allele ne '.' && $r_aL >= 0.9){
	$add .= ";ref_cn_allele=$r_allele";
	$add .= ";ref_cn_allele_L=$r_aL";
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
    elsif($r->{final}{cn} eq $OPT{ploidy} && !$r->{is_somatic}){
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
    elsif($r->{final}{cn} == $OPT{ploidy} && $r->{final}{loh} && $r->{is_somatic}){
	$data[2] = 'LOH_region';
	$att .= ';Reference_seq=~;Variant_seq=~';
    }
    else{
	$data[2] = 'region';
	$att .= ';Reference_seq=~;Variant_seq=~';
    }

    #add genes affected
    if($r->{gene}){
	$att .= ";gene=$r->{gene}"
    }

    #add final cn
    my $f_cn     = $r->{final}{cn};
    my $f_allele = $r->{final}{allele};
    my $f_L      = $r->{final}{L};
    my $f_aL     = $r->{final}{allele_L};
    my $f_sc     = $r->{maf_score}{best}{$f_cn}[1];

    $att .= ";Variant_copy_number=$f_cn";
    $att .= ";liklihood=".sprintf('%.2f', $f_L);
    if($f_allele ne '.' && $f_aL >= 0.9){
	$att .= ";allele=$f_allele";
	$att .= ";allele_liklihood=$f_aL";
    }
    $att .= ';is_somatic='.$r->{is_somatic} if($r->{is_somatic});
    $att .= ';is_loh='.$r->{final}{loh} if($r->{final}{loh});


    $att .= $add;
    push(@data, $att);  

    return join("\t", @data)."\n"
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
            $points = $r->{q_length}/$OPT{lfrag};
            $points = $r->{cov_count} if($r->{cov_count} < $points);
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


#assigns digital copy numbers to the results
sub assign_copy_number {
    my $res   = shift;
    my $exp   = shift;
    my $OPT  = shift;

    #assign preliminary copy numbers
    foreach my $r (@$res){
	my $cn = _best_cn($r, $exp->{_ALL}, $OPT);
	keys %$cn;
	while(my ($key, $value) = each %$cn){
	    $r->{final}{$key} = $value;
	}
    }

    #divide by chromosome
    my %sets;
    foreach my $r (@$res){
	push(@{$sets{$r->{chr}}}, $r);
    }
    
    #refine coverage by chromosome for final cn
    foreach my $chr (keys %sets){
	my $set = $sets{$chr};
	my $rc = refine_coverage($set, $exp->{_ALL}, $OPT->{models}, $OPT);
	$exp->{$chr} = $rc;
	foreach my $r (@$set){
	    _assign_cn($r, $rc, $OPT);
	}
    }
}

sub _assign_cn {
    my $r   = shift;
    my $bc  = shift;
    my $OPT = shift;

    my $cn     = _best_cn($r, $bc, $OPT);
    my $cn_max = $OPT->{cn_max};
    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};
    my $cfrac  = $OPT->{cfrac};
    my $rfrac  = 1 - $cfrac;

    my $fin = $cn;
    if($r->{ref}){
	my $rc = $r->{ref}{ccor_f}*$bc || $bc;
	my $ccn = _best_cn($r, $rc, $OPT);

	my $points = $r->{q_length}/$OPT->{lfrag};
        $points = $r->{maf_count} if($r->{maf_count} < $points);
        $points = 1 if($points < 1);

	if($cn->{cn} <= $cn_max && $ccn->{cn} <= $cn_max && $r->{maf_count} > 15){
	    if($ccn->{allele_L} > $cn->{allele_L} && $cn->{allele_L}/$ccn->{allele_L} < 0.1){
		$fin = $ccn;
	    }
	    elsif($cn->{allele_L} > $ccn->{allele_L} && $ccn->{allele_L}/$cn->{allele_L} < 0.1){
		$fin = $cn;
	    }
	    else{
		$fin = $ccn;
	    }
	}
	else{
	    $fin = $ccn; #gotta go with the corrected value
	}

	#reference spikes make the calls unreliable
	my $d_cn = abs($cn->{cn} - $ccn->{cn});
	if($r->{ref}{ccor_f} >= 2.5){
	    $fin = {cn       => '!',
		    L        => 0,
		    allele   => '!',
		    allele_L => 0,
		    allele_rss => 1000};
	}
	elsif($r->{ref}{cn} eq 0){
	    $fin = {cn       => '.',
                    L        => 0,
                    allele   => '.',
                    allele_L => 0,
		    allele_rss => 1000};
	}
	elsif($r->{q_length}/$r->{length} < 0.5){ #almost entirely unqueriable
	    $fin = {cn       => '.',
                    L        => 0,
                    allele   => '.',
                    allele_L => 0,
		    allele_rss => 1000};
	}
	elsif($d_cn >= 2 && $d_cn/(0.5*($cn->{cn}+$ccn->{cn})) >= 0.3){ #too much variation
	    $fin = {cn       => '.',
                    L        => 0,
                    allele   => '.',
                    allele_L => 0,
		    allele_rss => 1000};
	}
	elsif($cfrac < 1 && $cn->{cn} == 0){ #check internal ref of primary for lack of coverage
	    my $rbc = $bc * $rfrac/$cfrac;
	    my $x_adj = ($r->{xeno}) ? $r->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
	    my $cne = ($r->{cov_median}-$x_adj)/$bc;
	    $cne = 0 if($cne < 0);

	    my $e_cov0 = int($cne)*$bc+$x_adj+(2*$bc*$rfrac)/($cfrac);
	    my $e_cov1 = (int($cne)+1)*$bc+$x_adj+(2*$bc*$rfrac)/($cfrac);

	    my $p0 = (int($cne) == 0) ? p_cov($r, -$e_cov1) : p_cov($r, $e_cov0); #expect of 0 causes issues
	    my $p1 = p_cov($r, $e_cov1);

	    my ($L0, $cn0, $L1, $cn1) = ($p0, int($cne), $p1, int($cne)+1);
	    if($L0 > $L1){
		$cne = $cn0;
	    }
	    elsif($L1 > $L0){
		$cne = $cn1;
	    }
	    elsif($L0 == $L1){
		if(logLr_cov($r, $e_cov0, $e_cov1) >= 0){
		    $cne = $cn0;
		}
		else{
		    $cne = $cn1;
		}
	    }

	    #if ref internal to CNV is not diploid then everything is wrong
	    if($cne != 2){
		$fin = {cn       => '.',
			L        => 0,
			allele   => '.',
			allele_L => 0,
			allele_rss => 1000};
	    }
	}
	else{ #is seg long enough to distinguish the ref copy number as different than 0
	    my $e_cov0 = $bc;

	    #relationship must alwauys be true --> n' > (c_o * a^2)/(d - 0.5)^2
	    my $median = $r->{ref}{cov_median};
	    my $minL = $OPT->{lfrag}*$median*3.890592**2/($median-0.5)**2;
	    if($r->{ref}{q_length} < $minL){
		$fin = {cn       => '.',
			L        => 0,
			allele   => '.',
			allele_L => 0,
			allele_rss => 1000};
	    }
	}

	keys %$ccn;
	while(my ($key, $value) = each %$ccn){
	    $r->{ccor}{$key} = $value;
	}
    }

    #it's a masked segment
    if($r->{mask_count}/$r->{length} > 1/2){
	$fin = {cn       => 'N',
		L        => 0,
		allele   => 'N',
		allele_L => 0,
		allele_rss => 1000};
    }

    keys  %$cn;
    while(my ($key, $value) = each %$cn){
	$r->{ucor}{$key} = $value;
    }

    keys  %$fin;
    while(my ($key, $value) = each %$fin){
	$r->{final}{$key} = $value;
    }

    #mark loss of heterozygosity
    if($OPT->{use_ref} && $r->{maf_count} > 15 && $r->{q_length} > 30000){
	if($r->{ref}{allele} =~  /^0/){
	    $r->{ref}{final}{loh} = 1;
	    $r->{ref}{loh} = 1;
	    $r->{ucor}{loh}  = 0;
	    $r->{ccor}{loh}  = 0;
	    $r->{final}{loh} = 0;
	}
	else{
            $r->{ref}{final}{loh} = 0;
            $r->{ref}{loh} = 0;
	    $r->{ucor}{loh}  = ($r->{ucor}{allele} =~ /^0/) ? 1 : 0;
	    $r->{ccor}{loh}  = ($r->{ccor}{allele} =~ /^0/) ? 1 : 0;
	    $r->{final}{loh} = ($r->{final}{allele} =~ /^0/) ? 1 : 0;
	}	
    }
    else{
	$r->{ref}{final}{loh} = 0 if($r->{ref});
	$r->{ref}{loh}        = 0 if($r->{ref});
	$r->{ucor}{loh}       = 0;
	$r->{ccor}{loh}       = 0;
	$r->{final}{loh}      = 0;
    }

    #shortcuts for quick access
    $r->{cn} = $r->{final}{cn};
    $r->{L} = $r->{final}{L};
    $r->{allele} = $r->{final}{allele};
    $r->{allele_L} = $r->{final}{allele_L};
    $r->{allele_rss} = $r->{final}{allele_rss};

    return;
}

sub _best_cn {
    my $r = shift;
    my $bc = shift;
    my $OPT = shift;

    my $cfrac  = $OPT->{cfrac};
    my $xcov   = $OPT->{xcov};
    my $models = $OPT->{models};
    my $m_aln  = $OPT->{m_aln};
    my $t_cov  = $OPT->{t_cov};
    my $cn_max = $OPT->{cn_max};
    my $thr    = $OPT->{maf_tail_filt};
    my $rfrac  = 1 - $cfrac;

    my $x_adj = ($r->{xeno}) ? $r->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;

    my $cne = ($r->{cov_median}-$x_adj-(2*$bc*$rfrac/($cfrac)))/$bc;
    $cne = 0 if($cne < 0);
    
    my $ret = {};
    my $maf_ignore;
    $maf_ignore = 1 if($r->{q_length} < 30000 || $r->{maf_count} < 15);
    $maf_ignore = 1 if($r->{cov_median} <= 1);
    $maf_ignore = 1 if($rfrac && $r->{ref} && $r->{ref}{allele} =~ /^0/); #messes up MAF expect
    $maf_ignore = 1 if($r->{cov_median} < $t_cov/8); #if t_cov is 4 copy then 1/8 is closer to 0
    if($cne >= $cn_max || $maf_ignore){
	my $e_cov0 = int($cne)*$bc+$x_adj+(2*$bc*$rfrac)/($cfrac);
	my $e_cov1 = (int($cne)+1)*$bc+$x_adj+(2*$bc*$rfrac)/($cfrac);
	my $p0 = (int($cne) == 0) ? p_cov($r, -$e_cov1) : p_cov($r, $e_cov0); #expect of 0 causes issues
	my $p1 = p_cov($r, $e_cov1);
	my ($L0, $cn0, $L1, $cn1) = ($p0, int($cne), $p1, int($cne)+1);
	if($L0 > $L1){
	    $cne = $cn0;
	}
	elsif($L1 > $L0){
	    ($L0, $L1) = ($L1, $L0);
	    ($cn0, $cn1) = ($cn1, $cn0);
	}
	elsif($L0 == $L1){
	    if(logLr_cov($r, $e_cov0, $e_cov1) >= 0){
		($L0, $L1) = ($p0, $p1);
	    }
	    else{
		($L0, $L1) = ($p1, $p0);
		($cn0, $cn1) = ($cn1, $cn0);
	    }
	}
	
	if($cn0 > $cn_max || $maf_ignore){
	    $ret->{cn}         = $cn0;
	    $ret->{L}          = $L0;
	    $ret->{allele}     = '.';
	    $ret->{allele_L}   = 0;
	    $ret->{allele_rss} = 1000;

	    return $ret;
	}
    }

    #separate alleles by copy number
    my @cnsc;
    _mark_best_alleles($r, $OPT);
    keys %{$r->{maf_score}{best}};
    while(my ($key, $value) = each %{$r->{maf_score}{best}}){
	next if($key eq 'all');
	next if($key eq '!'); #temp
	
	my $cn = $key;
	my ($m, $sc) = @{$value};
	my $e_cov = $cn*$bc+(2*$bc*$rfrac)/($cfrac);
	my $a_cov = ($cn+1)*$bc+(2*$bc*$rfrac)/($cfrac) if($e_cov == 0); #expect of 0 causes issues
	$e_cov += $x_adj if($m =~ /xeno/);
	my $pt = ($a_cov) ? p_cov($r, -$a_cov) : p_cov($r, $e_cov);
	next if($pt <= 0 && abs($cne - $cn) >= 1.5);
	my $L = $pt*$sc->{L};
	$cnsc[$cn] = [$L, $sc, $m, $cn];
    }
    my ($best0, $best1) = sort {$b->[0] <=> $a->[0] ||
				    $a->[1]{rss} <=> $b->[1]{rss} ||
				    abs($cne-$a->[3]) <=> abs($cne-$b->[3])} grep {$_} @cnsc;

    #can't assign these based on MAF alone
    #if($best0->[0] == 0 && $best0->[3] <= 1 && int($cne) <= 1 && $rfrac < $thr){
    #	my $e_0 = 0*$bc+(2*$bc*$rfrac)/($cfrac) + $x_adj; #max cov for 0
    #	my $e_1 = 1*$bc+(2*$bc*$rfrac)/($cfrac); #min cov for 1
    #	my $d0 = abs($e_0 - $r->{cov_median});
    #	my $d1 = abs($e_1 - $r->{cov_median});
    #	$best0 = ($d1 <= $d0) ? $cnsc[1] : $cnsc[0];
    #}
    #elsif($r->{maf_score}{best}{all}[0] !~ /^0/ && $best1->[2] !~ /^0/){
    #	my $m0 = $best0->[2];
    #	my $m1 = $best1->[2];
    #
    #	$best0 = $best1 if($r->{maf_score}{$m1}{rss} < $r->{maf_score}{$m0}{rss});
    #}

    #best
    my $L0  = $best0->[0];
    my $pc0 = $best0->[1]{L};
    my $rs0 = $best0->[1]{rss};
    my $m0  = $best0->[2];
    my $cn0 = $best0->[3];

    #fill in return value
    $ret->{cn}         = $cn0;
    $ret->{L}          = $L0;
    $ret->{allele}     = $m0;
    $ret->{allele_L}   = $pc0;
    $ret->{allele_rss} = $rs0;
    
    return $ret;
}

#given a coverage and expect returns the best CN for that coverage
sub cn_for_cov {
    my $o_cov = shift;
    my $bc = shift;
    my $OPT = shift;

    my $cfrac = $OPT->{cfrac} || 1;
    my $xcov  = $OPT->{xcov}  || 0;
    my $m_aln = $OPT->{m_aln} || 0;
    my $rfrac = 1 - $cfrac;

    my $cne = int(($o_cov -($xcov*$m_aln) -(2*$bc*$rfrac/($cfrac)))/$bc);
    $cne = 0 if($cne < 0);
    my $e_cov0 = $cne*$bc+$xcov*$m_aln+(2*$bc*$rfrac)/($cfrac);
    my $e_cov1 = ($cne+1)*$bc+$xcov*$m_aln+(2*$bc*$rfrac)/($cfrac);

    #around 0 SD is closer to 1 copy so simulate with normal
    #I don't need to be exact, just need to know which is greater
    if($cne == 0){
	my $p0 = p_cov($o_cov, -$e_cov1, 1);
	my $p1 = p_cov($o_cov, $e_cov1, 1);
	return ($p0 > $p1) ? $cne : $cne+1;
    }

    return (logLr_cov($o_cov, $e_cov0, $e_cov1) >= 0) ? $cne : $cne+1;
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

sub _fit_targets{
    my $res = shift;
    my $t_cov  = shift;
    my $models = shift;
    my $OPT    = shift;

    return 50 if($res->[0]->{is_generic});

    my $cfrac = $OPT->{cfrac} || 1;
    my $xcov  = $OPT->{xcov}  || 0;
    my $m_aln = $OPT->{m_aln} || 0;
    my $thr   = $OPT->{maf_tail_filt};
    my $rfrac = 1 - $cfrac;

    #test only those where best MAF fit is not LOH
    my @test;
    my $sd = sqrt($t_cov);
    my @mod = grep {!/^(0|\!)/} @$models;
    @mod = grep {!/^(0\:0|\!)/} @$models if($rfrac > $thr); #include loh in mixed models
    foreach my $r (@$res){
	next if($r->{maf_count} < 15 || $r->{q_length} < 30000 || $r->{cov_median} <= 1);
	next if(abs($t_cov - $r->{cov_median}) > $sd); #too far away from target
	#next if($rfrac < $thr && maybe_loh($r, $OPT));
	_mark_best_alleles($r, $OPT);
	my ($best) = $r->{maf_score}{best}{all}[0];
	next if(get_c($best) >= 5);

	my $sc = $r->{maf_score}{$best};
	#next if($sc->{L} <= 0.9); #don't even try if it's a poor fit
	next if($sc->{rss}/$sc->{k} > 0.001); #don't even try if it's a poor fit
	push(@test, [$r, $sc, $best]);
    }

    while(@test > 200){
	#filter by distance from target frequency
	@test = sort {abs($t_cov-$a->[0]{cov_median}) <=> abs($t_cov-$b->[0]{cov_median})} @test;
	@test = @test[0..int(@test/2)];
	last if(@test < 200);

	#filter by longest
	@test = sort {$b->[0]{q_length} <=> $a->[0]{q_length}} @test;
	@test = @test[0..int(@test/2)];
	last if(@test < 200);

	#filter for closeness
	@test = sort {$a->[1]{rss} <=> $b->[1]{rss}} @test;
	@test = @test[0..int(@test/2)];
	last if(@test < 200);
    }

    #separate by target frequency (max 100 per bin)
    my %sets;
    foreach my $s (@test){
	my $tf = get_f_string($s->[-1], $OPT->{cfrac});
	$sets{$tf} ||= [];
	push(@{$sets{$tf}}, $s);
    }
    return $t_cov/2 if(!@test); #assume cn of 2
    
    #fit closest for each target frequency
    my %peaks;
    foreach my $k (keys %sets){
	my $size = 0;
	$size += $_->[0]{q_length} foreach (@{$sets{$k}});

	my @v = (@{$sets{$k}} > 100) ? @{$sets{$k}}[0..99] : @{$sets{$k}};
	@v = @v[0..int(@v/10)];

	my $rss_sum = 0;
	my $k_sum = 0;
	my %medians;
	foreach my $d (@v){
	    my ($r, $sc, $mod) = @$d;
	    my $t = $r->{cov_median};
	    next unless($t);
	    $rss_sum += $sc->{rss};
	    $k_sum += $sc->{k};
	    $medians{$t} += $r->{maf_count};
	}
	my $rss_ave = $rss_sum/$k_sum;

	my $ker = Statistics::KernelEstimation->new();
	foreach my $t (keys %medians){
	    next unless($t);
	    $ker->add_data($t, $medians{$t});
	}

	my $w = ($ker->default_bandwidth() < 2) ? 1 : $ker->default_bandwidth();
	$w = 1 if($w < 1);
	my ($low, $high) = $ker->extended_range();
	$low = 0 if($low < 0);
	my $o_cov = [0,0]; #expected segment mean coverage
	for(my $x = $low; $x <= $high; $x += ($high-$low)/100) {
	    my $y = $ker->pdf($x, $w);
	    $o_cov = [$x, $y] if($y >= $o_cov->[1]);
	    last if($x == $high);
	}
	$o_cov = $o_cov->[0];

	my $pt = p_cov($o_cov, $t_cov, $OPT->{n50}/$OPT->{lfrag});
	$peaks{$k} = [$o_cov, $pt, $rss_ave, $size];
    }

    #get cns for target frequencies
    my %target_fs;
    foreach my $m (@mod){
	push(@{$target_fs{get_f_string($m, $cfrac)}}, get_c($m));
    }

    #select the best matching target frequency
    my @select = keys %peaks;
    my $big;
    foreach my $s (@select){
	$big = $s if(!$big || $peaks{$s}[3] > $peaks{$big}[3]);
    }
    return _fit_targets($res, 2*$t_cov, $models, $OPT) if($big =~ /^0\:0\:/); #target is LOH
    @select = grep {!/^0\:0\:/} @select; #filter LOH
    @select = grep {$peaks{$_}[3]/$peaks{$big}[3] >= 0.5} @select; #get most abundant subset first
    my ($tf) = sort {$peaks{$a}[2] <=> $peaks{$b}[2]} @select; #closest MAF fit next

    #now get target copy number and base coverage
    my $tcn = ($tf) ? $target_fs{$tf}[0] : 2; #assume 2 copy if unclear
    $t_cov *= (($cfrac)*$tcn)/($rfrac*2+($cfrac)*$tcn);
    my $bc = $t_cov/$tcn;

    return $bc;
}

sub estimate_cellularity {
    my $OPT    = shift;

    die unless($OPT->{use_ref});

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

	    push(@work, "$chr:$pos-".($B-1));
	    $pos = $E+1;
	}
	if($pos < $end){
	    push(@work, "$chr:$pos-$end");
	}
    }  

    #--launch threads
    my @threads;
    prefork(); #prepare for forking
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'list'},
				       \&_estimate_cellularity_thread,
				       \@work,
				       $OPT));
    }
    
    #--jump into the mix
    my @res = _estimate_cellularity_thread(\@work, $OPT);

    #--collect results from threads
    foreach my $thr (@threads){
	my @ret = $thr->join();
	push(@res, @ret);
    }

    #--find global MAF peaks
    my $tcov;
    foreach my $s (@res){
	$tcov += $s->[2];
    }
    $tcov /= @res;

    my $rker = Statistics::KernelEstimation->new();
    my $sker = Statistics::KernelEstimation->new();
    my $fker = Statistics::KernelEstimation->new();
    foreach my $s (@res){
	$rker->add_data($s->[0]);
	$sker->add_data($s->[1]);
	#only add different from reference
	$fker->add_data($s->[1]) if($s->[0] != $s->[1]);
	$s->[2] /= $tcov;
	$s->[3] *= $s->[2];
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

    #filter all peaks
    my @rpeaks = sort {$a <=> $b} grep {$rpeaks->{$_} == 1} keys %$rpeaks;
    my @speaks = sort {$a <=> $b} grep {$speaks->{$_} == 1} keys %$speaks;
    my @fpeaks = sort {$a <=> $b} grep {$fpeaks->{$_} == 1} keys %$fpeaks;
    my ($rtop) = sort {$rdis[$b] <=> $rdis[$a]} @rpeaks;
    my ($stop) = sort {$sdis[$b] <=> $sdis[$a]} @speaks;
    @rpeaks = grep {$rdis[$_] >= 0.1*$rdis[$rtop]} @rpeaks; #ignore very small peaks
    @speaks = grep {$sdis[$_] >= 0.1*$sdis[$stop]} @fpeaks; #filtered fpeaks using sdis
    my ($filt) = grep {$_ <= 470} (@rpeaks, 470); #0.47 is ~10% cellularity

    return (1, 0) if($speaks[0] <= 1);
    if(!@speaks){
	@speaks = grep {$sdis[$_] >= 0.05*$sdis[$stop]} @fpeaks;
	return (0, 0) if(!@speaks);
    }

    #remove those near reference
    my @keep = grep {$_ < $filt} @speaks;
    if(!@keep && @speaks == 1){
	@speaks = grep {$sdis[$_] >= 0.05*$sdis[$stop]} @fpeaks;
	@keep = grep {$_ < $filt} @speaks;
    }
    @speaks = @keep if(@keep);

    #too difficult to calculate rear the edges
    $_ /= 1000 foreach(@speaks); #convert to decimal
    $_ /= 1000 foreach(@rpeaks); #convert to decimal
    $rtop /= 1000;
    $stop /= 1000;

    #do linear regression of outer two peaks to fit CN and celularity
    my $cell;
    my $rsq;
    if(($speaks[0] < 0.47 && @speaks > 2) || ($speaks[0] <= 0.35 && @speaks == 2)){
	#format --> (MAF, CN)
	my @data0 = ([$speaks[0], 2],
		     [$speaks[1], 1]);
	my @data1 = ([$speaks[0], 3],
		     [$speaks[1], 2]);

	#format --> (regression object, base CN key)
	my @reg = ([Statistics::Regression->new( "Base CN 1", ["Intercept", "Slope"]), \@data0],
		   [Statistics::Regression->new( "Base CN 2", ["Intercept", "Slope"]), \@data1]);

	#minus 2 in Y for expected intercept at 0
	$reg[0][0]->include(1/$_->[0] - 2, [1, $_->[1]]) foreach(@{$reg[0][1]});
	$reg[1][0]->include(1/$_->[0] - 2, [1, $_->[1]]) foreach(@{$reg[1][1]});

	if(@speaks > 3){
	    my @data2 = ([$speaks[2], 3],
			 [$speaks[1], 4],
			 [$speaks[0], 2]);

	    #peak in between will be CN 4 (occurs between LOH 1 and 2)
	    my @data3 = ([$speaks[0], 2],
			 [$speaks[2], 1]);

	    my @add = ([Statistics::Regression->new( "Base CN 3", ["Intercept", "Slope"]), \@data2],
		       [Statistics::Regression->new( "Base CN 2", ["Intercept", "Slope"]), \@data3]);

	    $add[0][0]->include(1/$_->[0] - 2, [1, $_->[1]]) foreach(@{$add[0][1]});
	    $add[1][0]->include(1/$_->[0] - 2, [1, $_->[1]]) foreach(@{$add[1][1]});
	    push(@reg, @add);
	}

	#the object and key closest to expected intercept
	my ($reg) = sort {abs($a->[0]->theta()->[0]) <=>
			  abs($b->[0]->theta()->[0])} @reg;

	#now force interecept for final slope
	my $data = $reg->[1]; #get key for best fitting model
	$reg = Statistics::Regression->new( "Title", ["Intercept", "Slope"]);
	$reg->include(1/$_->[0] - 2, [0, $_->[1]]) foreach(@$data);
	$reg->include(2 - 2, [0, 0]); #add point at 0 and 0.5 (inverse is 2)

	my $m = $reg->theta()->[1]; #get slope for celularity
	$cell = $m/(1+$m);
	$rsq = $reg->rsq();
    }

    #assume LOH is CN 1 and just uses outer peak
    my $min_rsq = (@speaks == 2) ? 0.98 : 0.95;
    if(!$rsq || $rsq < $min_rsq){
	my $reg = Statistics::Regression->new( "Title", ["Intercept", "Slope"]);
	if(@speaks <= 2){ #assume CN is 1 for primary peak
	    $reg->include(1/$speaks[0] - 2, [0, 1]); #0 for contant and -2 at Y forces the intercept
	}
	else{ #otherwise assume peak is CN 2
	    $reg->include(1/$speaks[0] - 2, [0, 2]);
	}
	$reg->include(2 - 2, [0, 0]); #add point at 0 and 0.5 (inverse is 2)
	my $m = $reg->theta()->[1]; #get slope for celularity
	$cell =$m/(1+$m);
	$rsq = $reg->rsq();
    }

    return wantarray ? ($cell, $rsq) : $cell;    

}

sub _estimate_cellularity_thread {
    my $work = shift;
    my $OPT    = shift;

    my $sid        = $OPT->{sid};
    my $vcf_file   = $OPT->{vcf_file};
    my $rid        = $OPT->{ref}{sid};
    my $rvcf_file  = $OPT->{ref}{vcf_file};
    my $thr        = $OPT->{maf_tail_filt};

    #load VCF file
    my $vcf;
    if($VCF{$vcf_file}){
	$vcf = $VCF{$vcf_file};
    }
    else{
	$vcf = Vcf->new(file=>"$vcf_file");
	$vcf->parse_header();
	$VCF{$vcf_file} = $vcf;
    }

    my ($SID) = grep {/$sid$/} $vcf->get_samples() if($sid);
    $vcf->set_samples(include=>[$SID]) if($sid);
   
    #load reference VCF
    my $rvcf;
    my $RID;
    if($rvcf_file && $rvcf_file ne $vcf_file){
	if($VCF{$rvcf_file}){
	    $rvcf = $VCF{$rvcf_file};
	}
	else{
	    $rvcf = Vcf->new(file=>"$rvcf_file");
	    $rvcf->parse_header();
	    $VCF{$rvcf_file} = $rvcf;
	}
	($RID) = grep {/$rid$/} $rvcf->get_samples() if($rid);
	$rvcf->set_samples(include=>[$RID]) if($rid);
    }
    elsif($rid){
	($RID) = grep {/$rid$/} $vcf->get_samples();
	$vcf->set_samples(include=>[$SID, $RID]);
    }
    
    #calculate first for separate reference VCF
    my @res;
    mkdir($OPT->{outdir}."/cell");
    while (my $seg = shift @$work){
	my $safe = uri_escape($seg);
	my $file = $OPT->{outdir}."/cell/$safe";

	my $sres = [];
	if(-f $file){
	    $sres = Storable::retrieve($file);
	    push(@res, @$sres);
	    next;
	}

	my %ref_ok;
	if($rvcf && $OPT->{use_ref}){
	    $rvcf->open(region=> $seg);
	    while(my $v = $rvcf->next_data_hash()){
		my $pos = $v->{POS};
		my $rAD = ($RID) ? $v->{gtypes}{$RID}{AD} : $v->{gtypes}{AD};
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
	while(my $v = $vcf->next_data_hash()){
	    my $pos = $v->{POS};
	    
	    #always calculate  reference first if merged in same VCF
	    if($OPT->{use_ref} && $RID && !$rvcf){
		my $rAD = $v->{gtypes}{$RID}{AD};
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
	    my $sAD = ($SID) ? $v->{gtypes}{$SID}{AD} : $v->{gtypes}{AD};
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
	    
	    #for 1000 SNVs
	    if(@bin >= 5000){ #1-2 Mb segment
		my %rmaf_set;
		my $rcov_mean;
		foreach my $s (@rbin){
		    my ($maf, $cov) = @$s;
		    next if($maf < 0.5);
		    $rmaf_set{$maf}{$cov}++;
		    $rcov_mean += $cov;
		}
		$rcov_mean /= @rbin;

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
			last if($maf < 0.48); #threshold I use to filter
		    }
		    else{
			last;
		    }
		}
		my $r_m = $rbest->[0];

		#ignore if reference is not as expected
		if($r_m < 0.48 || $r_m > 0.52){
		    undef @bin;
		    undef @rbin;
		    next;
		}

		#now fit sample
		my %smaf_set;
		my $scov_mean;
		foreach my $s (@bin){
		    my ($maf, $cov) = @$s;
		    next if($maf < 0.5);
		    $smaf_set{$maf}{$cov}++;
		    $scov_mean += $cov;
		}
		$scov_mean /= @bin;

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
			last;
		    }
		}		
		my $s_m = $sbest->[0];

		undef @bin;
		undef @rbin;
                next if(abs(0.5-$s_m) < abs(0.5-$r_m)); #sample is closer than ref

		push(@$sres, [$r_m, $s_m, $rcov_mean, $scov_mean]);
	    }
	}

	#save copy for recovery on failure or restart
	Storable::store($sres, $file);
	push(@res, @$sres);
    }

    return @res;
}

sub _add_maf_peak {
    my $r = shift;

    #kernel for MAF of segment
    my $ker = Statistics::KernelEstimation->new();
    my $maf_set = $r->{maf_set};
    foreach my $maf (keys %$maf_set){
	foreach my $cov (keys %{$maf_set->{$maf}}){
	    my $count = $maf_set->{$maf}{$cov};
	    $ker->add_data($maf, $count, 1/($cov+1));
	    #$ker->add_data(1-$maf, $count, 1/($cov+1)); #for symetry
	}
    }

    my @dis;
    my ($min, $max) = (0, 0.5);
    for(my $x = $min; $x <= $max; $x += 0.01) {
	my $y = $ker->pdf_width_from_data($x);
	push(@dis, $y);
    }
    my $peaks = find_peaks(\@dis, 3);

    #filter peaks
    my @peaks = grep {$peaks->{$_} == 1} keys %$peaks;
    @peaks = sort {$dis[$b] <=> $dis[$a]} @peaks;
    @peaks = grep {$dis[$_] >= 0.1*$dis[0]} @peaks; #ignore very small peaks
    
    if($r->{ref} && $r->{ref}{is_ref}){
	$r->{maf_peak} = $peaks[0];
    }
    else{
	$r->{maf_peak} = ($peaks[0] == 0 && $peaks[1]) ? $peaks[1] : $peaks[0];
    }
}

#fits the coverage expect to all segments
sub _coverage_fit {
    my $res = shift;
    my $t_cov = shift;
    my $models = shift;
    my $OPT    = shift;
    my $points = shift;

    return 50 if($res->[0]->{is_generic});

    my $cfrac = $OPT->{cfrac};
    my $xcov  = $OPT->{xcov};
    my $m_aln = $OPT->{m_aln};
    my $thr   = $OPT->{maf_tail_filt};
    my $rfrac = 1 - $cfrac;

    #test only those where best MAF fit is not LOH 
    my @test;
    my @mod = grep {!/^(0|\!)/} @$models;
    @mod = grep {!/^(0\:0|\!)/} @$models if($rfrac > $thr); #include loh in mixed models
    foreach my $r (@$res){
	next if($r->{maf_count} < 15 || $r->{q_length} < 30000 || $r->{cov_median} <= 1);
	next if($rfrac < $thr && maybe_loh($r, $OPT));
	_mark_best_alleles($r, $OPT);
	my ($best) = $r->{maf_score}{best}{all}[0];
        next if(get_c($best) >= 5);
	next if($best =~ /^0/ && $rfrac <= $thr); #don't use loh models

	my $sc = $r->{maf_score}{$best};
	#next if($sc->{L} <= 0.9); #don't even try if it's a poor fit
	next if($sc->{rss}/$sc->{k} > 0.001); #don't even try if it's a poor fit
	push(@test, [$r, $sc, $best]);
    }

    #now select max 200 segments for testing
    if(@test > 200){
	#filter by distance from target frequency
	@test = sort {abs($t_cov-$a->[0]{cov_median}) <=> abs($t_cov-$b->[0]{cov_median})} @test;
	@test = @test[0..int(@test*0.75)];
    }
    while(@test > 200){
	#filter by longest
	@test = sort {$b->[0]{q_length} <=> $a->[0]{q_length}} @test;
	@test = @test[0..int(@test*0.75)];
	last if(@test <= 200);

	#filter by closest matches with RSS
	@test = sort { $a->[1]{rss} <=> $b->[1]{rss} } @test;
	@test = @test[0..int(@test*0.75)];
    }
    $_ = $_->[0] foreach(@test);

    return if(!@test);

    #get coverage that generages min RSS and max L (approximate)
    my @stats;
    my $L_max;
    my $L_min;
    my $rss_min;
    my $rss_max;
    my $aicL_min;
    my $aicL_max;
    my $max = ($t_cov+sqrt($t_cov)*0.75) * ($cfrac)/($rfrac*2 + ($cfrac));
    my $min = ($t_cov * (4*$cfrac)/(2*$rfrac + (4*$cfrac)))/4; #if cn 4
    my $step = ($max/$min - $max/$max)/$points; #freq increase with each step
    for(my $f = $max/$max; $f <= $max/$min; $f+=$step){
	my $t1 = $max/$f;
	my $L_sum = 0;
	my $rss_sum = 0;
	my $aicL_sum = 0;
	foreach my $r (@test){
	    my @best;
	    foreach my $m (@mod){
		my $cn = get_c($m);
		next if($cn eq '0');
		
		my $e_cov = $cn*$t1+(2*$t1*$rfrac)/($cfrac);
		my $a_cov = ($cn+1)*$t1+(2*$t1*$rfrac)/($cfrac) if($e_cov == 0); #expect of 0 causes issues
		my $pt = ($a_cov) ? p_cov($r, -$a_cov) : p_cov($r, $e_cov);
		my $sc = $r->{maf_score}{$m};
		my $cL = $pt;
		if(!@best || $cL > $best[0] || ($cL == $best[0] && $sc->{rss} < $best[1]->{rss})){
		    @best = ($cL, $sc, $m);
		}
		
		next unless($xcov);
		my $x_adj =($r->{xeno}) ? $r->{xeno}{cov_median}*$OPT->{xeno}{ccor_f} : $xcov*$m_aln;
		$e_cov += $x_adj;
		$pt = p_cov($r, $e_cov);
		$sc = $r->{maf_score}{"$m:xeno"};
		$cL = $pt;
		if(!@best || $cL > $best[0] || ($cL == $best[0] && $sc->{rss} < $best[1]->{rss})){
		    @best = ($cL, $sc, $m);
		}
	    }
	    $L_sum   += $best[0] * $r->{maf_count}; #length weighted liklihood
	    $rss_sum += $best[1]->{rss}/$best[1]->{k} * $r->{maf_count}; #length weighted
	    $aicL_sum += $best[1]->{L} * $r->{maf_count}; #length weighted
	}
	$L_min = $L_sum if(!defined($L_min) || $L_sum < $L_min);
	$L_max = $L_sum if(!defined($L_max) || $L_sum > $L_max);
	$rss_min = $rss_sum if(!defined($rss_min) || $rss_sum < $rss_min);
	$rss_max = $rss_sum if(!defined($rss_max) || $rss_sum > $rss_max);
	$aicL_min = $aicL_sum if(!defined($aicL_min) || $aicL_sum < $aicL_min);
	$aicL_max = $aicL_sum if(!defined($aicL_max) || $aicL_sum > $aicL_max);
	push(@stats, [$t1, $L_sum, $rss_sum, $aicL_sum]);
    }

    my @bcs;
    my @bcs2;
    my @dist;
    my @dist2;
    my %index;
    my $last = 0;
    my $pass = 0;
    my $last2 = 0;
    my $pass2 = 0;
    my $max_sep;
    foreach my $s (@stats){
	my ($t1, $L_sum, $rss_sum, $aicL_sum) = @$s;

	#scale the liklihood
	$L_sum = (1 * ($L_sum-$L_min)/($L_max-$L_min));
	$aicL_sum = (1 * ($aicL_sum-$aicL_min)/($aicL_max-$aicL_min));

	#scale the rss
	$rss_sum = (1 * ($rss_sum-$rss_min)/($rss_max-$rss_min));

	my $sep = $L_sum - $rss_sum;
	$max_sep = [$t1, $sep] if(!$max_sep || $sep > $max_sep->[1]);

	print STDERR "$t1\t$L_sum\t$rss_sum\t$sep\n"; #temp

	#check occilation of sep peaks to get list of best
	my $cross = ($sep <=> $rss_sum);
	if($last && $last != $cross){
	    $pass++;
	}
	$last = $cross;

	#check occilation of lilihood peaks to get list
	my $cross2 = ($L_sum <=> $rss_sum);
	if($last2 && $last2 != $cross2){
	    $pass2++;
	}
	$last2 = $cross2;

	push(@dist, [$t1, $L_sum]);
	push(@dist2, [$t1, $rss_sum]); #second distribution for missing peaks
	$index{$t1} = [$L_sum, $rss_sum, $aicL_sum, $sep];
	if($cross > 0){
	    if(!$bcs[$pass] || $bcs[$pass][1] < $sep){
		$bcs[$pass] = [$t1, $sep];
	    }
	}
	if($cross2 > 0){
	    if(!$bcs2[$pass2] || $bcs2[$pass2][1] < $sep){
                $bcs2[$pass2] = [$t1, $sep];
            }
	}
    }
    @bcs = grep {$_} (@bcs, @bcs2);
    @bcs = $max_sep if(!@bcs);

    #add missing peaks
    my $peaks2 = find_peaks(\@dist2, int($points/2));
    my @alt = map {[$_, $index{$_}[3]]} grep {$peaks2->{$_} == -1} keys %$peaks2;
    if(@alt > @bcs){
	foreach my $o (@alt){
	    next if(grep {$_->[0] == $o->[0]} @bcs);
	    push(@bcs, $o);
	}
    }

    #now adjust to the max L peak near value
    my $peaks = find_peaks(\@dist, int($points/20));
    my @order = sort {$a <=> $b} keys %$peaks;
    @bcs = map {$_->[0]} sort {$b->[1] <=> $a->[1]} @bcs;
    foreach my $bc (@bcs){
	for(my $i = 0; $i < @order - 1; $i++){
	    my $j = $i + 1;
	    my $pi = $order[$i];
	    my $pj = $order[$j];
	    next unless($pi <= $bc && $bc <= $pj);
	    
	    #adjust base coverage to be at max liklihood
	    $bc = ($peaks->{$pi} == 1) ? $pi : $pj;
	}
    }

    #don't match to extremes (testing error)
    ($min, $max) = (sort {$a <=> $b} keys %$peaks)[0,-1];
    @bcs = grep {$_-$min >= 1e-13 && ($max-$_) >= 1e-13} @bcs; #1e-13 handles floating point precision

    #uniq
    my %uniq;
    @bcs = grep{! $uniq{$_}++} @bcs;

    #sort based on combined cov/MAF liklihood
    @bcs = sort {$index{$b}[0]*$index{$b}[2] <=> $index{$a}[0]*$index{$a}[2]} @bcs;

    #values lower in both coverage and liklihood are multiples of the base coverage
    @bcs = grep {$_ >= $bcs[0]} @bcs;

    return (wantarray) ? @bcs : $bcs[0];
}

sub refine_coverage {
    my $res = shift;
    my $base = shift;
    my $models = shift;
    my $OPT    = shift;

    return 50 if($res->[0]->{is_generic});

    my $cfrac  = $OPT->{cfrac};
    my $xcov   = $OPT->{xcov};
    my $m_aln  = $OPT->{m_aln};
    my $cn_max = $OPT->{cn_max};
    my $thr    = $OPT->{maf_tail_filt};
    my $rfrac  = 1 - $cfrac;

    #test only those where best MAF fit is not LOH
    my $before = 0;
    my $after = 0;
    my @cns = map {[]} (0..$cn_max);
    my @loh = map {[]} (0..$cn_max);
    foreach my $r (@$res){
	next if($r->{final}{cn} eq '.');
	next if($r->{final}{cn} eq '!');
	next if($r->{final}{cn} eq 'N');
	next if($r->{final}{cn} eq '0');
	$before += $r->{q_length};
	next if($r->{final}{cn} > 5);
	next if($r->{maf_count} < 15 || $r->{q_length} < 30000);
	if($rfrac > $thr && $r->{final}{allele} =~ /^0:/){
	    push(@{$loh[$r->{final}{cn}]}, $r);
	}
	else{
	    push(@{$cns[$r->{final}{cn}]}, $r);
	}
	$after += $r->{q_length};
    }
    return $base if(!$before || $after/$before < 0.5); #is not a representative sample

    #add set with good MAF values up to CN 3 (higher CN leads to uncertainty)
    my @test;
    my $content = 0;
    for(my $i = 2; $i < @cns && $i <= 3; $i++){
	push(@test, @{$cns[$i]});
	$content += $_->{q_length} foreach(@{$cns[$i]});
	last if($content/$before >= 0.5);
    }

    #add CN 1
    if($content/$before < 0.5 &&  @{$loh[1]}){
	push(@test, @{$loh[1]});
	$content += $_->{q_length} foreach(@{$loh[1]});
    }

    #add CN 4 with MAF
    if($content/$before < 0.5 && @{$cns[4]}){
	push(@test, @{$cns[4]});
	$content += $_->{q_length} foreach(@{$cns[4]});
    }

    #add other LOH regions if necessary
    if($content/$before < 0.5){
	for(my $i = 2; $i < @loh && $i <= 4; $i++){
	    push(@test, @{$loh[$i]});
	    $content += $_->{q_length} foreach(@{$loh[$i]});
	    last if($content/$before >= 0.5);
	}
    }

    #last resort add CN 5
    if($content/$before < 0.3 && @{$cns[5]}){
	push(@test, @{$cns[5]});
        $content += $_->{q_length} foreach(@{$cns[5]});
    }

    #CN 5 LOH (should I even do this?)
    if($content/$before < 0.3 && @{$loh[5]}){
	push(@test, @{$loh[5]});
        $content += $_->{q_length} foreach(@{$loh[5]});
    }

    my %medians;
    foreach my $r (@test){
	my $cn = $r->{final}{cn};
	my $rploidy = $OPT->{ref}{ploidy} || 2;
	my $t = ($r->{cov_median} * ($cn*$cfrac)/($rploidy*$rfrac + ($cn*$cfrac)))/$cn;
        next unless($t);

	#adjust for copy correction factor if available
	if($r->{ref} && $r->{ref}{ccor_f}){
	    $t /= $r->{ref}{ccor_f};
	}

        my $points = $r->{q_length}/$OPT->{lfrag};
        $points = $r->{cov_count} if($r->{cov_count} < $points);
        $points = 1 if($points < 1);
        $medians{$t} += $points;
    }

    my $sum = 0;
    my $count = 0;
    my $ker = Statistics::KernelEstimation->new();
    foreach my $t (keys %medians){
        next unless($t);
        $ker->add_data($t, $medians{$t});
        $sum += $t*$medians{$t};
        $count += $medians{$t};    
    }

    my $w = ($ker->default_bandwidth() < 2) ? 1 : $ker->default_bandwidth();
    $w = 1 if($w < 1);
    my $ave = $sum/$count;
    my ($min, $max) = $ker->extended_range();
    $min = 0 if($min < 0);
    $max = $ave*4 if($ave*4 < $max);
    my $step = ($max-$min)/100;
    if($step > 0.1){
	($min, $max) = (int($min), int($max));
	$step = 0.1;
    }

    my $t_cov = [0,0]; #expected segment mean coverage
    for(my $x = $min; $x <= $max; $x += $step) {
        my $y = $ker->pdf($x, $w);
        $t_cov = [$x, $y] if($y >= $t_cov->[1]);
        last if($x == $max);
    }
    $t_cov = $t_cov->[0];

    return ($t_cov) ? $t_cov : $base;
}

#redo the steps for estimating target coverage and copy number
sub recall_cn {
    my $res = shift; 
    my $OPT = shift;

    $OPT->{n50}= n50($res);

    if($res->[0]{ref}){
	my @refs;
	push(@refs, $_->{ref}) foreach(@$res);
	my $param = ref_param($OPT);
	($OPT->{ref}{t_cov}, $OPT->{ref}{base_cov}) = fit_expected_coverage(\@refs, $param);
	$param = ref_param($OPT); #reset
	$OPT->{ref}{chr_expects} = {_ALL => $OPT->{ref}{base_cov}};
	assign_copy_number(\@refs, $OPT->{ref}{chr_expects}, $param);
	$param = ref_param($OPT); #reset
	($OPT->{ref}{ploidy_ave}, $OPT->{ref}{ploidy}) = get_ploidy(\@refs, $param);
	foreach my $ref (@refs){
	    my $exp = $OPT->{ref}{chr_expects}{$ref->{chr}};
	    my $obs = $ref->{cov_median}/$OPT->{ref}{ploidy};
	    $ref->{ccor_f} = ($ref->{ucor}{cn}) ? $obs/$exp : 1;
	}
    }
    
    ($OPT->{t_cov}, $OPT->{base_cov}) = fit_expected_coverage($res, $OPT);
    $OPT->{chr_expects} = {_ALL => $OPT->{base_cov}};
    assign_copy_number($res, $OPT->{chr_expects}, $OPT);
    ($OPT->{ploidy_ave}, $OPT->{ploidy}) = get_ploidy($res, $OPT);
    label_somatic($res, $OPT);
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
    my $nm = 0.5 * chisqr($crit/2, 2*$n0);
    my $np = 0.5 * chisqr(1-$crit/2, 2*$n0+2);

    #poisson confidence innterval eqution for denominator
    #normal one SD equivilent for either side
    my $dm = 0.5 * chisqr($crit/2, 2*$d0);
    my $dp = 0.5 * chisqr(1-$crit/2, 2*$d0+2);
    
    #perform transform
    my $nsd = ($x <= $r) ? $nm : $np; #nummerator SD equivilent
    my $dsd = ($x <= $r) ? $dp : $dm; #denominator SD equivilent
    my $t = ($d0*$x-$n0)/sqrt($dsd* $x**2+$nsd);

    return gaus_pdf($t, 0, 1);
}
