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

use strict;
use warnings;
use Statistics::Distributions;
use constant PI => 4 * atan2(1, 1);

my $n0 = 8000; #numerator
my $d0 = 500; #denominator
my $r = $n0/$d0; #ratio
my $crit = (1-.682689492); #confidence interval of 1 SD of mean

#confidence innterval around numerator
my $nm = 0.5 * chisqr($crit/2, 2*$n0);
my $np = 0.5 * chisqr(1-$crit/2, 2*$n0+2);

#confidence innterval around denominator
my $dm = 0.5 * chisqr($crit/2, 2*$d0);
my $dp = 0.5 * chisqr(1-$crit/2, 2*$d0+2);

#perform transform
for(my $i=0; $i < 3*$r; $i += 3*$r/100){
    my $nsd = ($i <= $r) ? $nm : $np;
    my $dsd = ($i <= $r) ? $dp : $dm;
    my $t = ($d0*$i-$n0)/sqrt($dsd* $i**2+$nsd);

    print "$i\t$t\n";
}

#Statistics::Distributions has counterintuitive chisquared function
sub chisqr {
    return Statistics::Distributions::chisqrdistr($_[1], 1-$_[0]);
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
