use strict;
use warnings;
print "Sample name (no _all) separated by comma:";
my $sample=<>; chomp $sample;
print "26 (G, A, C, or T):";
my $nt26=<>; chomp $nt26;
my @sam=split /\s*\,+\s*/,$sample;
print "normalization standard:";
my $stan=<>; chomp $stan;
my $sizerange=30;
print "\n\n";
my %genome; &get_genome (\%genome);
my %alg; my %position;
&g26_stat (\%alg,\%position,\%genome);

foreach my $sa (@sam) {
	print "\t\t$sa:\n";
	my $ra=norm ($sa.'_all',$stan);
	my %result;
	&analyze_smRNA_distribution_flanking_26G ($sa.'_all',\$ra,\%genome,\%alg,\%position,\%result);

	open (hand0,">$sa\_all/$sa\_all_17_smRNA\_distribution_relative_to_ERGO126$nt26\_using_trans\_$stan.txt");
	foreach my $sr (keys %result) {
		for my $i (-66..39) {
			for my $j (17..$sizerange) {
				$result{$sr}{$i}{$j}+=0;
			}
		}
	}

	foreach my $sr (sort {$b cmp $a} keys %result) {
		print hand0 "\n$sr:\n";
		print hand0 'pos/size';
		for my $i (17..$sizerange) {
			print hand0 "\t$i";
		}
		print hand0 "\n";
		foreach my $p (sort {$b<=>$a} keys %{$result{$sr}}) {
			my $po=$p;
			$po++ if $po>=0;#-66 to -1 and 1 to 40
			print hand0 $po if $po % 10 ==0 || $po==-1 || $po==1; #print -60, -50, -1 and 1.....
			foreach my $si (sort {$a<=>$b} keys %{$result{$sr}{$p}}) {
				my $read=int($result{$sr}{$p}{$si}+0.5);
				print hand0 "\t$read";
			}
			print hand0 "\n";
		}
	}
	close hand0;
}

sub get_genome {
	my ($genome)=@_; my $chr;
	if (-e "ceWS215_gdna.fa") {
		open (hand1,"ceWS215_gdna.fa");
	}
	elsif (-e "ceWS215orsay_assemlong_gdna.fa") {
		open (hand1,"ceWS215orsay_assemlong_gdna.fa");
	}

	while (<hand1>) {
		$_ =~ s/\s+$//;
		if (/^>/) {
			$chr =$_;
			$chr =~ s/^>//;
			next;
		}
		${$genome}{$chr}.=$_;
	}
	close hand1;
}

sub g26_stat {
	my ($alg,$position,$ge)=@_; #26G target genes and positions
	open (hand1,"ERGO-1_targets.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my ($n1,$n2)=split /\t/;
		${$alg}{"$n1\t$n2"}=1;
	}
	close hand1; 

	my %po;#'key', genomic pos vs 'value', transcript pos
	my $gent=0; my $transnt=0; my $tranno=0; my $c=0; my %c;
	open (hand1,"ceWS215_cdna_pos_final.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//; my ($id,$po)=split /\|/;
		my @po=split /\s+/,$po; @po =sort {$a <=> $b} @po;
		my @id=split /\t/,$id; my @posi;
		next if !exists ${$alg}{"$id[1]\t$id[2]"};
		$tranno++;
		for (my $i=0; $i <$#po; $i+=2) {
			foreach my $j ($po[$i]..$po[$i+1]) {
				push @posi, $j;
				$transnt++;
				my $nt=substr (${$ge}{$id[5]},$j-1,1);
				if ($id[6] eq '-') {
					$nt =~ tr/ACGT/TGCA/;
				}	
				if (!exists $c{"$id[5]\t$id[6]\t$j"}) {
					if ($nt eq 'C') {
						$c++;
					}
					$gent++;
					$c{"$id[5]\t$id[6]\t$j"}=1;
				}	
			}
		}

		@posi=sort {$a <=> $b} @posi;
		if ($id[6] eq '-') {
			@posi=sort {$b <=> $a} @posi;
		}
		my $st=0;
		foreach my $i (@posi) {
			$st++;
			${$position}{"$id[1]\t$id[3]"}{$i}=$st;
		}
	}
	close hand1;
	print "average transcript size of ERGO-1 target mRNA ",$transnt/$tranno,"\n";
	print "genomic C rate for ERGO-1 26G target mRNA ", $c/$gent,"\n";
	print "genomic C sites for ERGO-1 26G target mRNA $c\n";
}

sub analyze_smRNA_distribution_flanking_26G {
	my ($sam,$r,$ge,$alg,$po,$res)=@_; my %an26G; my %all;
	my %loci;#match number based on transcripts
	open (hand1,"$sam/$sam\_17_cdna_normed.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		next if /^>/;
		$loci{$_}++;
	}
	close hand1;
			
	my $str; my $on=0; my $id; my %c;
	open (hand1,"$sam/$sam\_17_cdna_normed.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my @a1=split /\t/;
		if (/^>/) {
			$str=$a1[6];
			$id="$a1[1]\t$a1[3]";
			$on=0;
			if (exists ${$alg}{"$a1[1]\t$a1[2]"}) {
				$on=1;
			}
			next
		}

		#filters
		next if $on ==0;
		next if $a1[8] ne 'N';
		next if $a1[5]>5;
		next if length($a1[7]) >$sizerange;
		my $st=	${$po}{$id}{$a1[2]};
		$st=${$po}{$id}{$a1[2]+$a1[3]} if $a1[1] eq '-';

		if (length($a1[7]) ==26 || length($a1[7]) ==25) {
			if (substr($a1[7],0,1) eq $nt26 && $a1[1] ne $str) {
				$an26G{$id}{$st}+=$a1[4]/$a1[5]*$$r/$loci{$_};
				my $chst=$a1[2];
				if ($a1[1] eq '-') {
					$chst=$a1[2]+$a1[3];
				}
				$c{"$a1[0]\t$a1[1]\t$chst"}=1
			}
		}

		if ($a1[1] eq $str) {
			$all{$id}{'se'}{$st}{length($a1[7])}+=$a1[4]/$a1[5]*$$r/$loci{$_};
		}
		else {
			$all{$id}{'an'}{$st}{length($a1[7])}+=$a1[4]/$a1[5]*$$r/$loci{$_};
		}	
	}
	close hand1;
	my $cno=scalar keys %c; print "$sam has $cno ERGO-1 26G genomic loci\n\n\n";
	foreach my $trans (keys %an26G) {
		foreach my $st (keys %{$an26G{$trans}}) {
			my $r=$st+40; my $l=$st-65;
			foreach my $i ($l..$r) {
				foreach my $sr (keys %{$all{$trans}}) {
					if (exists $all{$trans}{$sr}{$i}) {
						foreach my $si (keys %{$all{$trans}{$sr}{$i}}) {
							${$res}{$sr}{$i-$st-1}{$si}+=$all{$trans}{$sr}{$i}{$si};
						}
					}
				}
			}
		}
	}
}

####normalization                                                                         
sub norm {
	my ($sam,$sd)=@_; my $ra=1;
	open (hand1,"$sam/sum\_$sam.txt") or die $!;
	my $geex=0; my $struc=0; my $cose=0;
	my $coan=0; my $coto=0; my $mir=0; my $u21=0;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		$geex=$1 if /^geex\t\S+\t\S+\t(\S+)/;
		$struc=$1 if /^struc\t(\S+)\t\S+\t\S+/;
		if (/^coding\t(\S+)\t(\S+)\t(\S+)/) {
			$cose=$1; $coan=$2; $coto=$3;
		}
		$mir=$1 if /^miRNA\t(\S+)\t\S+\t\S+/;
		$u21=$1 if /^u21RNA\t(\S+)\t\S+\t\S+/;
		last if /^star-miRNA/;
	}
	close hand1;

	my $nons=$geex-$struc;
	print "\t\t$sam\n\t\ttotal $geex\n\t\tsense struc $struc
		non-struc $nons\n\t\tsense coding $cose
		anti coding $coan\n\t\ttotal coding $coto
		sense miRNA $mir\n\t\tsense 21U-RNA $u21\n";
       
	if ($sd eq "mir") {
		$ra=1000000/$mir;
	}
	if ($sd eq "cose") {
		$ra=5000000/$cose;
	}
	if ($sd eq "coan") {
		$ra=2000000/$coan;
	}
	if ($sd eq "coto") {
		$ra=5000000/$coto;
	}
	if ($sd eq "nons") {
		$ra=5000000/$nons;
	}
	if ($sd eq "tot") {
		$ra=5000000/$geex;}
	if ($sd eq "21U") {
		$ra=500000/$u21;
	}
	print "\t\tNormalization ratio $ra\n\n";
	return $ra;
}
