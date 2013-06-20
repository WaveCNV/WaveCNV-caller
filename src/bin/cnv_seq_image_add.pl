#!/usr/bin/perl
use strict;
use warnings;
use GD;
use FindBin;
use Bio::DB::Sam;

my $usage = "
USAGE:
     $0 <vcf_file> <sample_type_id> <regions>

";

my $file = shift;
my ($name) = $file =~ /([^\/]+)$/;

my %CHR = (chr1 => 249250621,
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
	   chrY => 59373566,);

my ($h, $w) = (303, 3030);
my $s = shift || 'X';
my $seg_file = shift;

if(! $file){
    print $usage;
    exit();
}

my %segs;
if($seg_file){
    open(IN, "< $seg_file");
    while(my $line = <IN>){
	next if($line =~ /^\#/);
	chomp($line);
	next if(! $line);
	my @F = split(/\t/, $line);
	next if(!@F);

	my $x1 = int((($F[3]-1)/($CHR{$F[0]}-1)) * ($w-1));
	my $x2 = int((($F[4]-1)/($CHR{$F[0]}-1)) * ($w-1));
	my @att = split(/\;/, $F[8]);
	my ($cn1) = grep {/^Variant_copy_number=/} @att;
	$cn1 =~ s/Variant_copy_number=//;
	$cn1 .= '*' if(grep {$_ eq "is_somatic=1"} @att);
	my ($al1) = grep {/^allele=/} @att;
	if($al1){
	    $al1 =~ s/allele=//;
	    $al1 =~ s/\:xeno//;
	}
	push(@{$segs{$F[0]}}, [$x1, $x2, $cn1, $al1]);
    }
    close(IN);
}

while(my $chr = each %segs){
    my $png = "$name.$s.$chr.png";
    next if(! -f $png); #skip non-existant

    #create image pallet for chromosome
    my $x = $w;
    my $y = ($h * 2) + 1; #adjust to fit both images
    my $im1 =  GD::Image->newFromPng($png);
    my $red = $im1->colorClosest(255,0,0);
    
    #write lines for segments
    if($segs{$chr}){
	foreach my $d (@{$segs{$chr}}){
	    my ($x1, $x2, $cn1, $al1) = @$d;
	    if(abs($x2-$x1)+1 > 10){
		my $x = int((abs($x2-$x1)+1)/2)-2+$x1;
		my $y = int(($h-1)/10);
		$y = 5 if($y < 5);
		$im1->string(gdSmallFont,$x,$y,"$cn1",$red);
		$y += 10;
		$im1->string(gdSmallFont,$x,$y,"$al1",$red) if($al1);
	    }
	    
	    foreach my $x ($x1, $x2){
		$im1->line($x, 0, $x, $h-1, $red);
	    }
	}
    }
	
    open(OUT, "> $name.$s.$chr.section.png");
    binmode OUT;
    print OUT $im1->png;
    close(OUT);
}
