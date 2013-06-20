#!/usr/bin/perl
use threads;
use threads::shared;
use strict;
use warnings;
use FindBin;
use Storable;
use Bio::DB::Sam;
use feature "state";
use Tie::MmapArray;
use File::Copy;

my $bam_base = '/.mounts/labs/spbprojects/PCSI/data/analysis/somaticSNVtoDate6/';
$bam_base = "/Users/cholt/hn1/lab_share/carson/pcsi_data/data2/" if(! -d $bam_base);
my @samp = qw(PCSI0002 PCSI0005 PCSI0006 PCSI0024);
my %CHRS = (chr1 => 249250621,
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

my %files;
foreach my $s (@samp){
    my $dir = "$bam_base/$s/wgs/bam_files";
    my @bam_files = <$dir/*.bam>;
    $files{$s} = \@bam_files;
}

foreach my $chr (keys %CHRS){
    my $file = "$chr\_generic.store";
    open(my $MM, ">$file"); #create empty file
    close($MM);
}

my @list :shared;
foreach my $s (@samp){
    foreach my $chr (keys %CHRS){
	push(@list, "$chr\t$s");
    }
}

for (my $c = 1; $c < 10; $c++){
    #threads->create(\&for_thread);
}
for_thread();
$_->join() foreach(threads->list());

sub for_thread{
    while (my $job = shift @list){
	my ($chr, $s) = split(/\t/, $job);
	my $max = $CHRS{$chr};
	my $bam = $files{$s};
	my $sid = $s.'X';
	my $store_file = "$sid\_$chr\_generic.store";
	if(!-f $store_file){
	    open(MM, ">$store_file") && close(MM); #create file
	}

	my @array;
	tie @array, 'Tie::MmapArray', $store_file, { template => 'i', nels => $max};
	for(my $i = 1; $i < $max; $i += 1000000){
	    my $start = $i;
	    my $end = ($i + 1000000) - 1;
	    $end = $max if($end > $max);    
	    next if($array[$start-1] || $array[$end-1]);
    
	    my $cov = get_bam_coverage($chr, $start, $end, $bam, $sid);
	    #push(@array, @$cov); #push doesn't work on this tie yet
	    my $offset = $start-1;
	    $array[$offset++] = $_ foreach(@$cov);
	}
	
	my $mean = 0;
	$mean += $_ foreach(@array);
	$mean /= @array;

	my $store_file2 = "$chr\_generic.store";
	my @array2;
	tie @array2, 'Tie::MmapArray', $store_file2, { template => 'i', nels => $max};
	for(my $i = 0; $i < @array; $i++){
	    $array2[$i] += $array[$i]/$mean * 1/@samp;
	}
    }
}

sub get_bam_coverage{
    my $chr = shift;
    my $start = shift;
    my $end = shift;
    my $bam_files = shift;
    my $sid = shift;

    state %BAMS;
    state %RG;

    if(!$BAMS{$sid}){
	if(@$bam_files == 1){
	    foreach my $chr (keys %CHRS){
		$BAMS{$sid}{$chr} = $bam_files->[0];
	    }
	}
	else{
	    foreach my $f (@$bam_files){
		next unless($f =~ /\.(chr[\dXY]+)\./ && $CHRS{$1});
		$BAMS{$sid}{$1} = $f;
	    }
	}
	foreach my $chr (keys %CHRS){
	    next if(!$BAMS{$sid}{$chr});
	    my $fasta = "$FindBin::Bin/../overlap_package/lib/Feature/hg19_random.fa";
	    $BAMS{$sid}{$chr} = Bio::DB::Sam->new(-bam  => $BAMS{$sid}{$chr},
						  -autoindex => 1,
						  -fasta=> $fasta);
	}
    }

    my $sam = $BAMS{$sid}{$chr};
    die "ERROR: There is no bam file for $chr.\n" if(!$sam);
    my $segment = $sam->segment($chr,$start, $end);

    if($sid && !$RG{$sid}{$chr}){ #get read groups for filtering
        my @headers = grep {/^\@RG\t/} split(/\n/, $sam->header->text);
        foreach my $h (@headers){
            my %tags;
            foreach (grep {/^(ID:.*|SM:.*)$/} split(/\t/, $h)){
                my @F = split(/\:/, $_);
                $tags{$F[0]} = $F[1];
            }
            next unless( keys(%tags) == 2);
            push(@{$RG{$sid}{$chr}{$tags{SM}}}, $tags{ID});
        }
        $RG{$sid}{$chr} ||= {};
    }

    my $cov_list;
    if(!$sid || !$RG{$sid}{$chr}){
        my ($obj) = $segment->features(-type => 'coverage');
        $cov_list = $obj->coverage();
    }
    else{
        my ($SID) = grep {/$sid$/} keys %{$RG{$sid}{$chr}};
        my ($obj) = $segment->features(-type => 'coverage',
                                       -filter => $RG{$sid}{$chr}{$SID});
        $cov_list = $obj->coverage();
    }
    $cov_list ||= [];

    return $cov_list;
}
