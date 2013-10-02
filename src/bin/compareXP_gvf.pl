#!/usr/bin/perl 
use strict;

my $xfile = shift;
my $pfile = shift;
my $tag = shift || 'Variant_copy_number';

my $xploidy = 2;
my $pploidy = 2;

open(my $xIN, "< $xfile");
while(my @d = next_line($xIN, 'Variant_copy_number')){
    next if($d[3] eq '!' || $d[3] eq 'N' || $d[3] eq '.');
    if($d[4] eq 'copy_number_increase' || $d[4] eq 'copy_number_gain'){	
	$xploidy-- while($d[3] <= $xploidy);
    }
    elsif($d[4] eq 'copy_number_decrease' || $d[4] eq 'copy_number_loss'){
	$xploidy++ while($d[3] >= $xploidy);
    }
    elsif($d[4] eq 'region' || $d[4] eq 'LOH_region'){
	$xploidy = $d[3];
    }
}
close($xIN);

open(my $pIN, "< $pfile");
while(my @d = next_line($pIN, 'Variant_copy_number')){
    next if($d[3] eq '!' || $d[3] eq 'N' || $d[3] eq '.');
    if($d[4] eq 'copy_number_increase' ||  $d[4] eq 'copy_number_gain'){	
	$pploidy-- while($d[3] <= $pploidy);
    }
    elsif($d[4] eq 'copy_number_decrease' || $d[4] eq 'copy_number_loss'){
	$pploidy++ while($d[3] >= $pploidy);
    }
    elsif($d[4] eq 'region' || $d[4] eq 'LOH_region'){
	$pploidy = $d[3];
    }
}
close($pIN);

my %xdata;
my %pdata;
open(my $xIN, "< $xfile");
while(my @d = next_line($xIN, $tag)){
    next if($d[3] eq '!' || $d[3] eq 'N' || $d[3] eq '.' || $d[3] eq '?');
    push(@{$xdata{$d[0]}}, \@d);
}
close($xIN);

open(my $pIN, "< $pfile");
while(my @d = next_line($pIN, $tag)){
    next if($d[3] eq '!' || $d[3] eq 'N' || $d[3] eq '.'|| $d[3] eq '?');
    push(@{$pdata{$d[0]}}, \@d);
}
close($pIN);

foreach my $chr (keys %xdata,grep {!$xdata{$_}} keys %pdata){
    $xdata{$chr} ||= [];
    $pdata{$chr} ||= [];
}

my $ptotal = 0; #total bp in primary events
my $xtotal = 0; #total bp in xenograft events
my $pseg   = 0; #total segments in priamry events
my $xseg   = 0; #total segments in xenograft events
my $xbpok   = 0; #xenograft breakpoints confirmed
my $pbpok   = 0; #primary breakpoints confirmed

my $plseg = 0; #primary loss segment count
my $pgseg = 0; #primary gain segment count
my $plsegm = 0; #primary loss segment match count
my $pgsegm = 0; #primary gain segment match count
my $psegc = 0; #primary segment crossmatch count
my $xlseg = 0; #xenograft loss segment count
my $xgseg = 0; #xenograft gain segment count
my $xlsegm = 0; #xenograft loss segment match count
my $xgsegm = 0; #xenograft gain segment match count
my $xsegc = 0; #xenograft segment crossmatch count

my $pbmatch = 0; #primary boundary match
my $xbmatch = 0; #xenograft boundary match

my $xgsum = 0; #xenograft gain bp count
my $pgsum = 0; #primary gain bp count
my $xlsum = 0; #xenograft loss bp count
my $plsum = 0; #primary loss bp count

my $bgsum = 0; #total match gain bp count
my $blsum = 0; #total match loss bp count
my $bcsum = 0; #total crossmatch bp count

my $xcount = 0; #per segment xenograft bp count for matching loss/gain
my $pcount = 0; #per segment primary bp count for matching loss/gain

my $zxcount = 0; #per segment xenograft bp count for LOH overlap
my $zpcount = 0; #per segment primary bp count for LOH overlap
foreach my $chr (keys %xdata){
    my $xd;
    my $pd;
    my ($xC, $xB, $xE, $xcn, $xt, $xz, $xl, $xm, $xlb, $xrb);
    my ($pC, $pB, $pE, $pcn, $pt, $pz, $pl, $pm, $plb, $prb);

    while (1){
	if(!$xd || !@{$pdata{$chr}}){
	    $xd = shift @{$xdata{$chr}};
	    if($xl){
		my $f = $xcount/$xl if($xl);
		my $f2 = $zxcount/$xl if($xl);
		print "$xfile\t$xC\t$xB\t$xE\t$xl\t$xcn\t$xz\t$f\t$f2\n" if($xC);
		undef $xl;
	    }
	    $xcount = 0;
	    $zxcount = 0;

	    if($xd){
		($xC, $xB, $xE, $xcn, $xt, $xz) = @$xd;
		$xl = abs($xE-$xB)+1;
		$xm = 0; #reset
		$xlb = 0;
		$xrb = 0;
		if($xcn > $xploidy){
		    $xt = 1;
		    $xgsum += $xl;
		    $xtotal += $xl;
		    $xgseg++;
		    $xseg++;
		}
		elsif($xcn < $xploidy){
		    $xt = -1;
		    $xlsum += $xl;
		    $xtotal += $xl;
		    $xlseg++;
		    $xseg++;
		}
		else{
		    #undef $xd;
		    #next;
		}
	    }
	}

	if(!$pd || !$xd){
	    $pd = shift @{$pdata{$chr}};
	    if($pl){
		my $f = $pcount/$pl if($pl);
		my $f2 = $zpcount/$pl if($pl);
		print "$pfile\t$pC\t$pB\t$pE\t$pl\t$pcn\t$pz\t$f\t$f2\n" if($pC);
		undef $pl;
	    }
	    $pcount = 0;
	    $zpcount = 0;

	    if($pd){
		($pC, $pB, $pE, $pcn, $pt, $pz) = @$pd;
		$pl= abs($pE-$pB)+1;
		$pm = 0; #reset
		$plb = 0;
		$prb = 0;
		if($pcn > $pploidy){
		    $pt = 1;
		    $pgsum += $pl;
		    $ptotal += $pl;
		    $pgseg++;
		    $pseg++;
		}
		elsif($pcn < $pploidy){
		    $pt = -1;
		    $plsum += $pl;
		    $ptotal += $pl;
		    $plseg++;
		    $pseg++;
		}
		else{
		    #undef $pd;
		    #next;
		}
	    }
	}

	last if(!$xd && !$pd);
	next if(!$xd || !$pd);

	#check breakpoints
	my $dis = 1000;
	if(abs($xB-$pB) <= $dis){
	    $xbpok++ if(!$xlb);
	    $pbpok++ if(!$plb);
	    $xlb = $plb = 1;
	}
	if(abs($xE-$pE) <= $dis){
	    $xbpok++ if(!$xrb);
	    $pbpok++ if(!$prb);
	    $xrb = $prb = 1;
	}
	if(abs($xE-$pB) <= $dis){
	    $xbpok++ if(!$xrb);
	    $pbpok++ if(!$plb);
	    $xrb = $plb = 1;
	}
	if(abs($xB-$pE) <= $dis){
	    $xbpok++ if(!$xlb);
	    $pbpok++ if(!$prb);
	    $xlb = $prb = 1;
	}

	#get overlap
	my @o = sort {$a <=> $b} ($xB, $xE, $pB, $pE);
	if($o[1] == $xB || $o[1] == $pB){
	    my $l = abs($o[2]-$o[1])+1;
	    if($xt > 0 && $pt > 0){ #match gain
		$bgsum += $l;
		$xcount += $l;
		$pcount += $l;
		if(! $xm && $l/$xl >= 0.5){ #at least 50% overlap
		    $xgsegm++;
		    $xm = 1;
		}
		if(!$pm && $l/$pl >= 0.5){ #at least 50% overlap
		    $pgsegm++;
		    $pm = 1;
		}
	    }
	    elsif($xt < 0 && $pt < 0){ #match loss
		$blsum += $l;
		$xcount += $l;
		$pcount += $l;
		if(! $xm && $l/$xl >= 0.5){ #at least 50% overlap
		    $xlsegm++;
		    $xm = 1;
		}
		if(!$pm && $l/$pl >= 0.5){ #at least 50% overlap
		    $plsegm++;
		    $pm = 1;
		}
	    }
	    elsif(($xt > 0 && $pt < 0) || ($xt < 0 && $pt > 0)){ #crossmatch
		$bcsum += $l;
		if(! $xm && $l/$xl >= 0.5){ #at least 50% overlap
		    $xsegc++;
		    $xm = 1;
		}
		if(!$pm && $l/$pl >= 0.5){ #at least 50% overlap
		    $psegc++;
		    $pm = 1;
		}
	    }

	    #consistency of zygosity
	    if($xz == $pz){
		$zxcount += $l;
                $zpcount += $l;
	    }
	}

	if($xE <= $pE){
	    undef $xd;
	}
	elsif($pE < $xE){
	    undef $pd;
	}
	else{
	    die "Logic error";
	}
    }
}

my $mtotal = ($bgsum+$blsum);
my $xmsegtotal = $xgsegm+$xlsegm;
my $pmsegtotal = $pgsegm+$plsegm;
my $xmp = $mtotal/$xtotal;
my $pmp = $mtotal/$ptotal;
my $xgp = ($xgsum) ? $bgsum/$xgsum : 0;
my $pgp = ($pgsum) ? $bgsum/$pgsum : 0;
my $xlp = ($xlsum) ? $blsum/$xlsum : 0;
my $plp = ($plsum) ? $blsum/$plsum : 0;
my $xcp = $bcsum/$xtotal;
my $pcp = $bcsum/$ptotal;

print "SEG_ALL\tSEG_AMP\tSEG_DEL\tBASE_ALL\tBASE_AMP\tBASE_DEL\tMATCH_ALL\tMATCH_AMP\tMATCH_DEL\tCROSS\tP_MATCH\tP_MATCH_AMP\tP_MATCH_DEL\tP_CROSS\tSEG_MATCH_ALL\tSEG_MATCH_AMP\tSEG_MATCH_DEL\tSEG_CROSSMATCH\tBP_MATCH\n";
print "$xseg\t$xgseg\t$xlseg\t$xtotal\t$xgsum\t$xlsum\t$mtotal\t$bgsum\t$blsum\t$bcsum\t$xmp\t$xgp\t$xlp\t$xcp\t$xmsegtotal\t$xgsegm\t$xlsegm\t$xsegc\t$xbpok\n";
print "$pseg\t$pgseg\t$plseg\t$ptotal\t$pgsum\t$plsum\t$mtotal\t$bgsum\t$blsum\t$bcsum\t$pmp\t$pgp\t$plp\t$pcp\t$pmsegtotal\t$pgsegm\t$plsegm\t$psegc\t$pbpok\n";

sub next_line {
    my $IN = shift;
    my $tag = shift || 'Variant_copy_number';

    while(my $line = <$IN>){
	next if($line =~ /^#/);
	chomp($line);
	next if(!$line);
	my @F = split(/\t/, $line);
	my $loh = ($F[8] =~ /is_loh=1/)? 1 : 0;
	if($F[8] =~ /[\;\t]$tag=([^\;\n]+)/){
	    my @d = ($F[0], $F[3], $F[4], $1, $F[2], $loh);
	    if($d[3] eq '.' || $d[3] eq '!' || $d[3] eq 'N'){#temp
		if($d[3] eq 'N'){
		    #$d[3] = 0;
		}
		elsif($d[3] eq '.' && $F[8] =~ /[\;\t]ucor_cn=([^\;\n]+)/){
		    #$d[3] = $1;
		}
	    	elsif($F[8] =~ /[\;\t]ucor_cn=([^\;\n]+)/){
	    	    #$d[3] = $1;
	    	}
	    }
	    return @d;
	}
	else{
	    next;
	}
    }

    return;
}
