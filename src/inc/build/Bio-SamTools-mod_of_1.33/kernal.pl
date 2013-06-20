#!/usr/bin/perl

use strict;
use warnings;
use Statistics::KernelEstimation;

my $s = Statistics::KernelEstimation->new();

my $sum = 0;
my %data;
while(<STDIN>){
    chomp($_);
    my ($pos, $len, $cov) = split(/\t/, $_);
    $cov = int($cov * 10)/10;
    $data{$cov} += $len;

}

foreach my $cov (keys %data){    
    next unless($data{$cov});
    $s->add_data($cov, $data{$cov});
}

my $bw;
#$bw = $s->optimal_bandwidth();
$bw = 1 if(!$bw || $bw < 0);

my ($min, $max) = $s->extended_range();
for( my $x=$min; $x<=$max; $x+=($max-$min)/100 ) {
    print $x, "\t", $s->pdf( $x, $bw ), "\n";
}
