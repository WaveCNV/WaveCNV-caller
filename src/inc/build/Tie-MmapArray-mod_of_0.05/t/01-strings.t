#!/usr/bin/perl

use Test;
use strict;
use lib qw(blib/lib blib/arch ../blib/lib ../blib/arch);

BEGIN { plan tests => 2 };

use Tie::MmapArray;

my $file = "testfile";

my $string = pack("iii", 42, 43, 23);
my @array;
my $failed;


open(FILE, ">$file") or die "cannot create testfile\n";
binmode FILE;
print FILE $string;
close FILE;

tie @array, 'Tie::MmapArray', $file, { template => 'i', nels => 3 };

print "hello world";
