use strict;
use warnings;

my $genesize=1000;
my $spreadstep=1000;
my $failrate=0.3;
my $it=100;
my $g26no=20;

my %dis; my %motif;
for my $j (1..$it) {
	print $j,"\n";
	for my $i (1..1000) {
		&assign_26G ($j,\%dis,\%motif);
	}
}

foreach my $p (sort {$b <=> $a} keys %dis) {
	foreach my $i (1..$it) {
		$dis{$p}{$i}+=0;
	}
}

open (hand1,">ERGO1_26G_distribution_simulation_expression_$failrate\_$g26no.txt");
foreach my $p (sort {$a <=> $b} keys %dis) {
	my $pp=$p*(-1);
	my @val=values %{$dis{$p}};
	my $count=$#val+1; my $tot=0;
	foreach my $val (@val) {
		$tot+=$val;
	}
	my $ave=$tot/$count;
	my $var=0;
	foreach my $val (@val) {
		$var+=($val-$ave)**2;
	}
	my $sterror=($var/($count-1))**0.5/$count**0.5;
	$ave=int($ave+0.5); $sterror=int($sterror*100+0.5)/100;
	print hand1 "$pp\t$ave\t$sterror\n";
}
close hand1;

open (hand1,">ERGO1_26G_distribution_simulation_motif_$failrate\_$g26no.txt");
print hand1 "pos\tA\tC\tG\tT\n";
foreach my $p (sort {$b <=> $a} keys %motif) {
	foreach my $nt (qw (T G C A)) {
		foreach my $i (1..$it) {
			$motif{$p}{$nt}{$i}+=0;
		}
	}
}

foreach my $p (sort {$b <=> $a} keys %motif) {
	my $pp=$p*(-1);
	print hand1 $pp;
	foreach my $nt (qw (T G C A)) {
		my @val=values %{$motif{$p}{$nt}};
		my $count=$#val+1; my $tot=0;
		foreach my $val (@val) {
			$tot+=$val;
		}
		my $ave=$tot/$count;
		my $var=0;
		foreach my $val (@val) {
			$var+=($val-$ave)**2;
		}
		my $sterror=($var/($count-1))**0.5/$count**0.5;
		$ave=int($ave+0.5); $sterror=int($sterror*100+0.5)/100;
		print hand1 "\t$ave\t$sterror";
	}
	print hand1 "\n";
}
close hand1;

sub assign_26G {
	my ($iteration,$dis,$motif)=@_;
	my $seq =random_seq (); my %pos;

	my @pos; #@pos, G postion (+1), 1000 nts trimmed from 5' and 3'
	while ($seq =~ /G/g) {
		push @pos,pos($seq);
		$pos{pos($seq)}=1;
	}
	
	my $posno=$#pos+1; my %result; 
	for my $i (1..$g26no) {
		my $p=int(rand($posno));#randomly select an array element based on index;
		my $po=$pos[$p]; my $read=rand (10);
		
		#spread 26G downstream; spread $spreadstep steps
		$result{$po}+=$read;
		for my $j (1..$spreadstep) {
			$po+=23;
			last if $po > $genesize;

			my $seqdown=substr ($seq, $po-1);
			while ($seqdown =~ /G/g) {
				my $gp=pos($seqdown);
				last if $po+$gp-1> $genesize-25;
				my $ran=rand (1);
				if ($ran >$failrate) {
					$po+=($gp-1);
					$result{$po}+=$read;
					last;
				}
			}
		}
	}

	my $ccount;
	while ($seq =~ /G/g) {
		$ccount++;
	}
	
	#foreach my $lc (keys %result) {
		#print "$lc not exist\n" if !exists $pos{$lc};
	#}

	my $loci=scalar keys %result;
	#print "loci number $loci","\t",$ccount,"\n";

	foreach my $st (keys %result) {
		for my $i (-40..66) {
			if (exists $result{$i+$st-1}) {
				${$dis}{$i}{$iteration}+=$result{$i+$st-1};
			}
		}

		for my $ipo (-50..49) {
			my $uppo=$st+$ipo;
			if ($uppo>0 && $uppo <=length($seq)) {
				my $upnt=substr($seq,$uppo-1,1);
				${$motif}{$ipo}{$upnt}{$iteration}++;
			}
		}
	}
}
		
#generate randome sequence composed of 20% C/G and 30% A/T.	
sub random_seq	{
	my $s;
	for my $i (1..$genesize) {
		my $r=int(rand (100));
		if ($r <21) {
			$s.='C';
		}
		if ($r>20 && $r <42) {
			$s.='G';
		}
		if ($r> 41 && $r <71) {
			$s.='A';
		}
		if ($r >70) {
			$s.='T';
		}
	}
	return $s;
}
