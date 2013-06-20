#!/usr/bin/perl
use strict;
use warnings;
use GD;
use FindBin;
use Bio::DB::Sam;

my $usage = "
USAGE:
     cnv_seq_imager.pl <vcf_file> <sample_type_id> <optional_segments>

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
my $cov = 8;
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

my @sids;
my %bins1;
my %bins2;
my $lc = 0;
my $last;
open(IN, "<$file");
while(my $line = <IN>){
    next if($line =~ /^\#\#/);
    chomp($line);
    next if(! $line);
    my @F = split(/\t/, $line);
    if($line =~ /^\#CHROM/){
	$F[0] =~ s/^\#//;
	@sids = @F[9..$#F] if($#F > 8);
	($s) = grep {/$s/} @sids;
	next;
    }

    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @E) = @F;
    last if($last && $chrom ne $last);
    $last = $chrom;

    my %samp;
    if(@sids){
	@samp{@sids} = @E;
    }
    elsif(@F == 3){
	@samp{qw(P R X)} = @E;
    }
    elsif(@F == 4){
	@samp{qw(C P R X)} = @E;
    }

    foreach my $v (values %samp){
	my %hash;
	@hash{split(/\:/, $format)} = map {[split(/\,/, $_)]} split(/\:/, $v);
	$v = \%hash;
    }

    next unless($samp{$s}{AD} && @{$samp{$s}{AD}});

    my $cr = $samp{$s}{AD}[0];
    my $ca = $samp{$s}{AD}[1];
    next unless(defined($cr) && defined($ca));
    my $tot = $cr + $ca;

    next unless($tot);

    my $maf = $ca/$tot;
    next if(!$CHR{$chrom});
    my $x = int((($pos-1)/($CHR{$chrom}-1)) * ($w-1));
    my $y = int($maf * ($h - 1));

    next if($tot < $cov);
    push(@{$bins2{$x}}, $tot);

    next if($maf <= 0.05 || $maf >= 0.95);
    $bins1{$x}{$y}++;

    $lc++;
    last if($lc > 100000);
}
close(IN);

my @counts1;
foreach my $x (keys %bins1){
    push(@counts1, values %{$bins1{$x}});
}

#make max intensity be at median 
@counts1 = sort {$a <=> $b} @counts1;
my $fac1 = int(255/($counts1[int($#counts1/2)]));
$fac1 = 1 if($fac1 < 1);

#get expected cov
my @counts2;
foreach my $x (keys %bins2){
    push(@counts2, @{$bins2{$x}});
}
@counts2 = sort {$a <=> $b} @counts2;
my $e = $counts2[int($#counts2/2)];

#now normalize to image size
my $lrmax = 2; #cn 8
my $lrmin = int(log(1/$e)); #around cn 0

print "min=$lrmin\tmax=$lrmax\n";

undef @counts2;
foreach my $x (keys %bins2){
    my %norm;
    foreach my $y (@{$bins2{$x}}){
	my $y = (log($y/$e) - $lrmin)/($lrmax-$lrmin);
	$y = 1 if($y > 1);
	$y = 0 if($y < 0);
	$y = int($y*$h);
	$norm{$y}++;
    }
    $bins2{$x} = \%norm;
    push(@counts2, values %{$bins2{$x}});
}
@counts2 = sort {$a <=> $b} @counts2;
my $fac2 = int(255/($counts2[int($#counts2 * 3/4)]));
$fac2 = 1 if($fac2 < 1);

my %images;
my @smooth;
my $exit = 0;
undef $last;
open(IN, "<$file");
while(my $line = <IN>){
    last if($exit);
    next if($line =~ /^\#\#/);
    chomp($line);
    next if(! $line);

    #get sample ids from the headers
    my @F = split(/\t/, $line);
    if($line =~ /^\#/){
	$F[0] =~ s/^\#//;
	@sids = @F[9..$#F] if($#F > 8);
	($s) = grep {/$s/} @sids;
	next;
    }

    #get standard column values
    my ($chrom, $pos, $id, $ref, $alt, $qual, $filter, $info, $format, @E) = @F;
    next if(!$CHR{$chrom}); #skip unmapped chromosomes
    next if(-f "$name.$s.$chrom.png"); #skip finished chromosomes

    #split meta data by sample
    my %samp;
    if(@sids){
	@samp{@sids} = @E;
    }
    elsif(@F == 3){
	@samp{qw(P R X)} = @E;
    }
    elsif(@F == 4){
	@samp{qw(C P R X)} = @E;
    }
    foreach my $v (values %samp){
	my %hash;
	@hash{split(/\:/, $format)} = map {[split(/\,/, $_)]} split(/\:/, $v);
	$v = \%hash;
    }

    #get MAF and coverage
    next unless($samp{$s}{AD} && @{$samp{$s}{AD}});
    my $cr = $samp{$s}{AD}[0];
    my $ca = $samp{$s}{AD}[1];
    next unless(defined($cr) && defined($ca));
    my $tot = $cr + $ca;
    next unless($tot);

    my $maf = $ca/$tot;
    #my $cov = $tot;

    #create image pallet for chromosome
    if(!$images{$chrom}){
	my $x = $w;
	my $y = ($h * 2) + 1; #adjust to fit both images
	my $im1 = new GD::Image($x,$y);
	my $white = $im1->colorAllocate(255,255,255);
	my $red = $im1->colorAllocate(255,0,0);
	my $black = $im1->colorAllocate(0,0,0);
	$images{$chrom}[0] = $im1;
	
	$im1->line(0, $h, $w, $h, $black); #line to separate images

	#write file as I advance to next chromosome
	if($last){
	    #prepare smoothed line for coverage
	    foreach my $p (@smooth) {
		$p = ($p) ? $p->[0]/$p->[1] : 0;
	    }
	    my $sum = 0;
	    my $count = 0;
	    
	    #init window
	    for(my $x = 0; $x < @smooth; $x++){
		$sum += $smooth[$x];
		$count++;
		last if($count == 5);
	    }
	    
	    #slide window across for running average
	    my $last_xy;
	    my $red = $images{$last}[0]->colorExact(255,0,0);
	    for(my $x = 0; $x < @smooth; $x++){
		my $i = $x - 6; #to drop
		my $j = $x + 5; #to keep
		
		if($i >= 0 && $count > 10){
		    $sum -= $smooth[$i];
		    $count--;
		}
		
		if($j < @smooth){
		    $sum += $smooth[$j];
		    $count++;
		}
		
		#paint smoothed line
		my $y = int($sum/$count);
		$y = $h-1 if($y >= $h);
		$y = $h - $y - 1; #flip for image as it's upsidedown
		
		if($last_xy){
		    $images{$last}[0]->line(@$last_xy, $x, $y, $red);
		    $last_xy = [$x, $y];
		}
		else{
		    $images{$last}[0]->setPixel($x, $y, $red);
		    $last_xy = [$x, $y];
		}
	    }
	    undef @smooth;

	    #write lines for segments
	    if($segs{$last}){
		foreach my $d (@{$segs{$last}}){
		    my ($x1, $x2, $cn1, $al1) = @$d;
		    if(abs($x2-$x1)+1 > 10){
			my $x = int((abs($x2-$x1)+1)/2)-2+$x1;
			my $y = int(($h-1)/10);
			$y = 5 if($y < 5);
			$images{$last}[0]->string(gdSmallFont,$x,$y,"$cn1",$red);
			$y += 10;
			$images{$last}[0]->string(gdSmallFont,$x,$y,"$al1",$red);
		    }
		
		    foreach my $x ($x1, $x2){
			$images{$last}[0]->line($x, 0, $x, $h-1, $red);
		    }
		}
	    }

	    open(OUT, "> $name.$s.$last.png");
	    binmode OUT;
	    print OUT $images{$last}[0]->png;
	    close(OUT);
	}
    }
    $last = $chrom;

    #get image pallete for this feature
    my $im1 = $images{$chrom}[0];

    #the x and y position for each image
    my $x = int((($pos-1)/($CHR{$chrom}-1)) * ($w-1));
    my $y1 = int($maf * ($h - 1));
    my $y2 = $tot/$e;

    #ranges for widening pixels
    my $rx = int(250000 * $w/$CHR{$chrom} + 1);
    my $ry1 = int(($h/$tot - 1)/2 + 1);
    $ry1 = int($ry1*.5);

    #range to widen pixels for image 2 (is log based so differnt above vs below)
    my $ry2min = $y2 - 1/(2*$e);
    my $ry2max = $y2 + 1/(2*$e);
    my $ry2 = 1/$e;
    foreach($y2, $ry2min, $ry2max){$_ = int((log($_)-$lrmin)/($lrmax-$lrmin) * ($h-1))};
    $ry2min++ if($ry2min < $y2);

    #get all pixels to fill for the ranges
    my @xs1 = grep {$_ >= 0} ($x-$rx)..($x+$rx);
    my @ys1 = map {$h - $_ - 1} grep {0 <= $_ && $_ < $h} ($y1-$ry1)..($y1+$ry1);
    my @ys2 = map {$h - $_ - 1} grep {0 <= $_ && $_ < $h} ($ry2min..$ry2max);
    @ys2 = $h-1 if(!@ys2 && $y2 >= $h);
    @ys2 = 0 if(!@ys2 && $y2 < 0);

    #determine smoothed average line for coverage in image 2
    $smooth[$x] ||= [];
    $smooth[$x][0] += $y2; #sum
    $smooth[$x][1]++; #count

    #intensity ajustment for image 1 pixels
    my $o1 = int($fac1/((4*2+1)*($ry1*2+1)));
    $o1 = 1 if($o1 < 1);

    #intensity ajustment for image 2 pixels
    my $o2 = int($fac2/((4*2+1)*($ry2+1+1)));
    $o2 = 1 if($o1 < 1);

    #now paint
    foreach my $x1 (@xs1){
	#paint image 1
	if($tot >= $cov){
	    foreach my $y (@ys1){
		my $y1 = $y + $h; #adjust to take up bottom frame
		my ($red,$green,$blue) = $im1->rgb($im1->getPixel($x1, $y1));
		if($red == 0 && $green == 0 && $blue == 0){
		    next;
		}
		else{
		    $red -= $o1;
		    $red = 0 if($red < 0);
		    $green = $blue = $red;		
		    my $grey = $im1->colorResolve($red,$green,$blue);
		    $grey = $im1->colorClosest($red,$green,$blue) if($grey == -1);
		    $im1->setPixel($x1, $y1, $grey);
		}
	    }
	}

	#paint image 2
	foreach my $y2 (@ys2){
	    my ($red,$green,$blue) = $im1->rgb($im1->getPixel($x1, $y2));
	    if($red == 0 && $green == 0 && $blue == 0){
		next;
	    }
	    else{
		$red -= $o2;
		$red = 0 if($red < 0);
		$green = $blue = $red;
		my $grey = $im1->colorResolve($red,$green,$blue);
		$grey = $im1->colorClosest($red,$green,$blue) if($grey == -1);
		$im1->setPixel($x1, $y2, $grey);
	    }
	}
    }
}
close(IN);

#handle last image in buffer
if($last){
    #prepare smoothed line for coverage
    foreach my $p (@smooth) {
	$p = ($p) ? $p->[0]/$p->[1] : 0;
    }
    my $sum = 0;
    my $count = 0;
    
    #init window
    for(my $x = 0; $x < @smooth; $x++){
	$sum += $smooth[$x];
	$count++;
	last if($count == 5);
    }
    
    #slide window across for running average
    my $last_xy;
    my $red = $images{$last}[0]->colorExact(255,0,0);
    for(my $x = 0; $x < @smooth; $x++){
	my $i = $x - 6; #to drop
	my $j = $x + 5; #to keep
	
	if($i >= 0 && $count > 10){
	    $sum -= $smooth[$i];
	    $count--;
	}
	
	if($j < @smooth){
	    $sum += $smooth[$j];
	    $count++;
	}
	
	#paint smoothed line
	my $y = int($sum/$count);
	$y = $h-1 if($y >= $h);
	$y = $h - $y - 1; #flip for image as it's upsidedown

	if($last_xy){
	    $images{$last}[0]->line(@$last_xy, $x, $y, $red);
            $last_xy = [$x, $y];
	}
	else{
	    $images{$last}[0]->setPixel($x, $y, $red);
	    $last_xy = [$x, $y];
	}
    }
    undef @smooth;

    #write lines for segments
    if($segs{$last}){
	foreach my $d (@{$segs{$last}}){
	    my ($x1, $x2, $cn1, $al1) = @$d;
	    if(abs($x2-$x1)+1 > 10){
		my $x = int((abs($x2-$x1)+1)/2)-2+$x1;
		my $y = int(($h-1)/10);
		$y = 5 if($y < 5);
		$images{$last}[0]->string(gdSmallFont,$x,$y,"$cn1",$red);
		$y += 10;
		$images{$last}[0]->string(gdSmallFont,$x,$y,"$al1",$red);
	    }
	    
	    foreach my $x ($x1, $x2){
		$images{$last}[0]->line($x, 0, $x, $h-1, $red);
	    }
	}
    }

    open(OUT, "> $name.$s.$last.png");
    binmode OUT;
    print OUT $images{$last}[0]->png;
    close(OUT);    
}

