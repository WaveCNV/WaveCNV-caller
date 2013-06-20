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
use Term::ProgressBar;
use IO::Interactive qw(is_interactive);
use File::NFSLock;
use Bio::DB::Sam;
use constant PI => 4 * atan2(1, 1); #3.14159265358979323846...

our ($sort_exe, $bgzip_exe, $tabix_exe, $vcftools_exe, $ksseg_exe, $hg19, $dummy);
our ($DUM, %VCF, %BAMS, %RG, %GEN);

BEGIN {
    eval 'require CNV_caller::ConfigData';
    my $edir = "$FindBin::RealBin/../exe/";
    if(!$@){
	$tabix_exe = CNV_caller::ConfigData->config('tabix');
	$vcftools_exe = CNV_caller::ConfigData->config('vcftools');
	$ksseg_exe = CNV_caller::ConfigData->config('ksseg');
    }

    $sort_exe = File::Which::which('sort');
    $ksseg_exe = "/.mounts/labs/lakshmilab/private/scripts/src/KSseg.code/KSseg" if(!$ksseg_exe);
    ($ksseg_exe) = grep {-f $_} ("$edir/*/ksseg", File::Which::which('KSseg')) if(! -f $ksseg_exe);
    ($tabix_exe) = grep {-f $_} ("$edir/*/tabix", File::Which::which('tabix')) if(!$tabix_exe);
    ($vcftools_exe) = grep {-f $_} ("$edir/*/bin/vcftools", File::Which::which('vcftools')) if(!$vcftools_exe);
    ($bgzip_exe = $tabix_exe) =~ s/tabix$/bgzip/; #comes with tabix

    $hg19 = "$FindBin::RealBin/../data/hg19_random.fa";
    $dummy = "$FindBin::RealBin/../data/dummy.bam";

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

     $script [options]

Description:

     Script does blah

Options:

     sid       <STRING>  Sample ID (for mixed VCF and BAM files)

     bam_dir   <PATH>    Directory of BAM files to use

     bam_list  <LIST>    List of BAM files to use (comma separated)

     vcf_file  <PATH>    VCF file to use (instead of BAM)

     sort                Indicates that the VCF needs to be sorted.

     fasta     <PATH>    Reference fasta to use (default is UCSC hg19)

     cpus      <INT>     CPUs to use

     tmp       <PATH>    Alternate temp file directory (default set by system)

";

#--set default values
my %OPT;
$OPT{cpus}  = 1;
$OPT{TMP}   = File::Spec->tmpdir();
$OPT{fasta} = $hg19;

#--get command line options
GetOptions("sid=s"       => \$OPT{sid},
	   "fasta=s"     => \$OPT{fasta},
	   "vcf_file=s"  => \$OPT{vcf_file},
           "bam_dir=s"   => \$OPT{bam_dir},
	   "bam_list=s"  => \$OPT{bam_list},
	   "cpus=i"      => \$OPT{cpus},
           "tmp=s"       => \$OPT{TMP},
           "sort"        => \$OPT{sort},
           "help|?"      => sub{print $usage; exit(0);});

#--validate and prepare command line options
my $err = validate_options(\%OPT);
die $err if($err);

#--print usage
if(!@{$OPT{bam_files}} && !$OPT{vcf_file}){
    print $usage;
    exit(0);
}

$OPT{vcf_file}  = prepare_vcf($OPT{vcf_file});
$OPT{bam_files} = [prepare_bam(@{$OPT{bam_files}})];
$OPT{fasta}     = prepare_fasta($OPT{fasta});
$OPT{CHRS}      = get_chr_lengths($OPT{fasta});
$OPT{GAPS}      = get_chr_gaps($OPT{fasta});    

#read in data to convert from BAM or VCF
my $segs;
if(@{$OPT{bam_files}}){
    $segs = bam2ksseg(\%OPT);
}
elsif($OPT{vcf_file}){
    $segs = vcf2ksseg(\%OPT);
}

#------------------------------------------------------------------------------
#-------------------------------- SUBROUTINES ---------------------------------
#------------------------------------------------------------------------------
sub ksseg {
    my $param = shift;

    return [] if(!$param->{lines});

    my $infile  = $param->{file};
    my $mapfile = $param->{map};
    my $lines   = $param->{lines};
    my $thr     = $param->{thr};
    my (undef, $outfile) = tempfile(CLEANUP => 1, TMPDIR => 1);

    my $cmd = "$ksseg_exe $lines $infile $outfile $thr";
    my $pid = open3(my $IN, '>&STDERR', '>&STDERR', $cmd);
    waitpid $pid, 0;
    my $tstat = $?;
    die "ERROR: KSSEG failed with exit status $tstat\n" if($tstat);

    open(IN, "< $outfile");
    open(MAP, "< $mapfile");
    <IN>; #drop header line
    my @last;
    my $seg;
    my @res;
    while(my $line = <IN>){
	my $mline = <MAP>;
	chomp $line;
	chomp $mline;
	
	my ($chr, $ratio, $mean, $std, $median) = split("\t", $line);
	my ($start, $end) = split("\t", $mline);
	
	if(!$seg || $seg->[0] ne $chr || $seg->[3] != $mean ||
	   $seg->[4] != $std || $seg->[5] != $median
	  ){
	    if (@last){
		$seg->[2] = $last[2];
		push(@res, $seg);
	    }
	    $seg = [$chr, $start, $end, $mean, $std, $median];
	}
	@last = ($chr, $start, $end, $ratio, $mean, $std, $median);
    }
    close(IN);
    close(MAP);
    if($seg){
	$seg->[2] = $last[2];
	push(@res, $seg);
    }
    unlink($outfile);

    return \@res;
}

#run before forking. clears stored values that may cause core dump
sub prefork {
    undef $DUM;
    undef %BAMS;
    undef %VCF;
}

#validates the input options
sub validate_options {
    my $OPT = shift;

    my $err;
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
    $OPT->{bam_files}  = \@bam_files;

    my $vcf_file  = $OPT->{vcf_file};
    my $fasta     = $OPT->{fasta};
    $err .= "File $fasta does not exist\n" if(! -f $fasta);
    $err .= "File $vcf_file does not exist\n" if($vcf_file && ! -f $vcf_file);
    $err .= "You cannot specify both BAM and VCF files\n" if($vcf_file && @bam_files);

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

sub vcf2ksseg{
    my $OPT = shift;

    my $vcf_file = $OPT->{vcf_file};
    return [] if(!$vcf_file);

    my @work :shared;
    my $CHRS = $OPT->{CHRS};
    my $GAPS = $OPT->{GAPS};
    foreach my $chr (keys %$CHRS){
	next if($chr =~ /[XYM]$/);
	my $end = $CHRS->{$chr};
	my $gaps = $GAPS->{$chr};
	my $pos = 1;
	foreach my $g (@$gaps){
	    my ($B, $E) = @$g;
	    push(@work, "$chr:$pos-".($B-1)) if($B > $pos);
	    $pos = $E+1;
	}
	push(@work, "$chr:$pos-$end") if($pos < $end);
    }

    #--launch threads
    my @threads;
    prefork(); #prepare for forking
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'scalar'},
				       \&_vcf2ksseg_thread,
				       \@work,
				       $OPT));
    }
    
    #--jump into the mix
    my $segs = _vcf2ksseg_thread(\@work, $OPT);

    #--collect results from threads
    foreach my $thr (@threads){
	my $seg = $thr->join();
	push(@$segs, @$seg);
    }

    #sort the results
    @$segs = sort {$a->[0] <=> $b->[0] ||
		       $a->[1] <=> $b->[1]} @$segs;
    print "CHR\tSTART\tEND\tMEAN\tSTDEV\tMEDIAN\n";
    foreach my $s (@$segs){
	print join("\t", @$s)."\n";
    }
}

sub _vcf2ksseg_thread {
    my $work = shift;
    my $OPT = shift;

    my $sid        = $OPT->{sid};
    my $vcf_file   = $OPT->{vcf_file};

    #load VCF file into object
    my $vcf;
    if($VCF{$vcf_file}){
	$vcf = $VCF{$vcf_file};
    }
    else{
	$vcf = Vcf->new(file=>"$vcf_file");
	$vcf->parse_header();
	$VCF{$vcf_file} = $vcf;
    }

    #parse only sample of interest
    my ($SID) = grep {/$sid$/} $vcf->get_samples() if($sid);
    $vcf->set_samples(include=>[$SID]) if($sid);
   
    #get coverage for each SNV
    my @segs;
    my @map;
    while(my $region = shift @$work){
	my ($kfh, $kfilename) = tempfile(CLEANUP => 1, TMPDIR => 1);
	my ($mfh, $mfilename) = tempfile(CLEANUP => 1, TMPDIR => 1);

	my ($chr, $start, $end) = $region =~ /^(.*)\:(\d+)\-(\d+)$/;
	$vcf->open(region=> $region);
	my $count = 0;
	my @kdata;
	my %cov;
	while(my $v = $vcf->next_data_hash()){
	    my $pos = $v->{POS};
	    my $sAD = ($SID) ? $v->{gtypes}{$SID}{AD} : $v->{gtypes}{AD};
	    next if(!$sAD || $sAD eq '.');
	    my ($rc, $ac) = split(/,/, $sAD);
	    next if(!defined($rc) || !defined($ac));
	    my $c = $rc+$ac;
	    push(@kdata, [_chrom($chr), $c]);
	    $cov{$c}++;
	    print $mfh "$pos\t$pos\n";
	    $count++;
	}

	my $median;
	my $place = 0;
	foreach my $c (sort {$a <=> $b} keys %cov){
	    $place += $cov{$c};
	    if($place > ($count/2)){
		$median = $c;
		last;
	    }
	}
	foreach my $d (@kdata){
	    $d->[1] = log($d->[1]/$median);	    
	    print $kfh $d->[0]."\t".$d->[1]."\n";
	}
	close($kfh);
	close($mfh);

	my $seg = ksseg({chr   => _chrom($chr),
			 start => $start,
			 end   => $end,
			 lines => $count,
			 thr   => 1e-10,
			 file  => $kfilename,
			 map   => $mfilename});
	push(@segs, @$seg);
	unlink($kfilename, $mfilename); #cleanup
    }

    return \@segs;
}

#fills in the coverage of a segment from the bam file
sub bam2ksseg{
    my $OPT = shift;

    my $bam_files = $OPT->{bam_files};
    return [] if(!@$bam_files);

    my @work :shared;
    my $CHRS = $OPT->{CHRS};
    my $GAPS = $OPT->{GAPS};
    foreach my $chr (keys %$CHRS){
	next if($chr =~ /[XYM]$/);
	my $end = $CHRS->{$chr};
	my $gaps = $GAPS->{$chr};
	my $pos = 1;
	foreach my $g (@$gaps){
	    my ($B, $E) = @$g;
	    push(@work, "$chr:$pos-".($B-1)) if($B > $pos);
	    $pos = $E+1;
	}
	push(@work, "$chr:$pos-$end") if($pos < $end);
    }

    #--launch threads
    my @threads;
    prefork(); #prepare for forking
    for(my $i = 1; $i < $OPT->{cpus}; $i++){
	push(@threads, threads->create({context => 'scalar'},
				       \&_bam2ksseg_thread,
				       \@work,
				       $OPT));
    }
    
    #--jump into the mix
    my $segs = _bam2ksseg_thread(\@work, $OPT);

    #--collect results from threads
    foreach my $thr (@threads){
	my $seg = $thr->join();
	push(@$segs, @$seg);
    }

    #sort the results
    @$segs = sort {$a->[0] <=> $b->[0] ||
		       $a->[1] <=> $b->[1]} @$segs;
    print "CHR\tSTART\tEND\tMEAN\tSTDEV\tMEDIAN\n";
    foreach my $s (@$segs){
        print join("\t", @$s)."\n";
    }
}

sub _bam2ksseg_thread {
    my $work = shift;
    my $OPT = shift;

    my $bam_files = $OPT->{bam_files};
    my $sid = $OPT->{sid};    
    my $CHRS = $OPT->{CHRS};
    my $GAPS = $OPT->{GAPS};

    #load BAM files into objects
    my $tag = $OPT->{bam_files}[0];
    if(! keys %{$BAMS{$tag}}){
	if(@$bam_files == 1){
	    foreach my $chr (keys %$CHRS){
		$BAMS{$tag}{$chr} = $bam_files->[0];
	    }
	}
	else{
	    foreach my $f (@$bam_files){
		next unless($f =~ /\.(chr[\dXY]+)\./ && $CHRS->{$1});
		$BAMS{$tag}{$1} = $f;
	    }
	}
	foreach my $chr (keys %$CHRS){
	    next if(!$BAMS{$tag}{$chr});
	    $BAMS{$tag}{$chr} = Bio::DB::Sam->new(-bam  => $BAMS{$tag}{$chr},
						  -autoindex => 1,
						  -fasta=> $OPT->{fasta});
	}
    }

    my @segs;
    my @map;
    foreach my $region (@$work){
        my ($kfh, $kfilename) = tempfile(CLEANUP => 1, TMPDIR => 1);
        my ($mfh, $mfilename) = tempfile(CLEANUP => 1, TMPDIR => 1);

	my ($chr, $start, $end) = $region =~ /^(.*)\:(\d+)\-(\d+)$/;
	my $bam = $BAMS{$tag}{$chr};
	die "ERROR: There is no bam file for $chr.\n" if(!$bam);

	#process BAM files
	my $count = 0;
        my @kdata;
        my %cov;
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

	#refined region edge
	my $length = ($end-$start)+1;
	my $int = int($length/101);
	my $target = 101*$int;

	#get segment
	my $segment = $bam->segment($chr,$start, $start+$target+1);

	#get overhange of region
	my $over;
	if($target < $length){
	    $over = $bam->segment($chr, $start+$target+2, $end);
	}
	
	#get coverage from BAM
	my $c;
	if(!$sid || !keys %{$RG{$tag}{$chr}}){
	    ($c) = $segment->features(-type => "coverage:$int");
	    ($over) = $over->features(-type => "coverage:1") if($over);
	}
	else{
	    my ($SID) = grep {/$sid$/} keys %{$RG{$tag}{$chr}};
	    ($c) = $segment->features(-type => "coverage:$int",
					-filter => $RG{$tag}{$chr}{$SID});
	    ($over) = $over->features(-type => "coverage:$int",
				      -filter => $RG{$tag}{$chr}{$SID}) if($over);
	}

	$c = $c->coverage() || [];
	for(my $i = 0; $i <@$c; $i++){
	    my $c = $c->[$i];
	    my $B = $start+$i*101+1;
	    my $E = $start+($i+1)*101;
            push(@kdata, [_chrom($chr), $c]);
            $cov{int($c)}++;
	    print $mfh "$start\t$E\n";
	    $count++;
	}

	#get coverage for overhang
	if($over){
	    my $B = $over->start;
	    my $E = $over->end;
	    my $c = ($over->coverage)[0];
            push(@kdata, [_chrom($chr), $c]);
            $cov{int($c)}++;
	    print $mfh "$start\t$E\n";
	    $count++;
	}
	close($kfh);
	close($mfh);

        my $median;
        my $place = 0;
        foreach my $c (sort {$a <=> $b} keys %cov){
            $place += $cov{$c};
            if($place > ($count/2)){
                $median = $c;
                last;
            }
        }
        foreach my $d (@kdata){
	    $d->[1] = log($d->[1]/$median);
             print $kfh $d->[0]."\t".$d->[1]."\n";
        }

	my $seg = ksseg({chr   => _chrom($chr),
                         start => $start,
                         end   => $end,
                         lines => $count,
                         thr   => 1e-10,
                         file  => $kfilename,
                         map   => $mfilename});
        push(@segs, @$seg);
        unlink($kfilename, $mfilename); #cleanup
    }

    return \@segs;
}


#integer for chromosome
sub _chrom {
    my $id = shift;

    if($id =~ /^chr(\d+)$/){
	return $1;
    }
    elsif($id eq 'chrX'){
	return 23;
    }
    elsif($id eq 'chrY'){
	return 24;
    }
    elsif($id eq 'chrXY'){
	return 25;
    }
    elsif($id eq 'chrM'){
	return 26;
    }
    else{
	return 1000;
    }
}

sub parse {
    my $OPT    = shift;

    die unless($OPT->{use_ref});


}
