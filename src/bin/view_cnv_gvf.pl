#!/usr/bin/perl 

use strict;
use warnings;

my $USAGE = "
USAGE:
     $0 <gvf_file>
     cat <gvf_file> | $0

";

my $IN;
if(@ARGV){
    my $file = shift @ARGV;
    open($IN, "<$file");
}
elsif(! -t){
    open($IN, "<&STDIN");
}
else{
    print $USAGE;
    exit;
}

while(my $line = <$IN>){
    next if($line =~ /^\#/);
    chomp $line;
    next if(! $line);
    my @F = split("\t", $line);
    my %at = %{parse_att($F[8])};
    my $l = ($F[4]-$F[3]) +1;
    print "$F[0]\t$F[3]\t$F[4]\t$F[5]\t$l".
	"\tcn=$at{Variant_copy_number}".
	"\tu_cn=$at{ucor_cn}".
	"\tc_cn=$at{ccor_cn}".
	"\tref_cn=$at{ref_cn}".
	"\tsom=$at{is_somatic}".
	"\tallele=$at{allele}".
	#"\tacn_L=$at{alt_cn_liklihood}".
	#"\tccor_L=$at{ccor_cn_liklihood}".
	#"\taccor_L=$at{ccor_alt_cn_liklihood}".
	#"\taL=$at{allele_liklihood}".
	#"\ta_aL=$at{alt_allele_liklihood}".
	#"\tc_aL=$at{ccor_allele_liklihood}".
	#"\tca_aL=$at{ccor_alt_allele_liklihood}".
	"\n";
}

sub parse_att {
    my $nine = shift;
    my %at = map {my @s = split(/\=/, $_); @s[0..1]} grep {$_} split(/\;/, $nine);
    #$_ = [split(/\,/, $_)] foreach (values %att);
    $at{Variant_copy_number} = '' if(!defined $at{Variant_copy_number});
    $at{ucor_cn}  = '' if(!defined $at{ucor_cn});
    $at{ccor_cn} = '' if(!defined $at{ccor_cn});
    $at{ref_cn}  = '' if(!defined $at{ref_cn});
    $at{liklihood} = '' if(!defined $at{liklihood});
    $at{is_somatic} = '' if(!defined $at{is_somatic});
    $at{allele} = '' if(!defined $at{allele});
    $at{ccor_cn_liklihood} = '' if(!defined $at{ccor_cn_liklihood});
    $at{ucor_cn_liklihood} = '' if(!defined $at{ucor_cn_liklihood});
    

    return \%at;
}
