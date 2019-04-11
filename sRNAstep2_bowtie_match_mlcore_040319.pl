use strict;
use warnings;
use Getopt::Std;

#current
print "sample name separated by comma:";
my $sample=<>; chomp $sample; my @sam=split /\s*\,+\s*/,$sample;
my $stan='nons';#mir, cose, coan, coto, tot, geex, 21U, nons
my $av='-v'; #-v for sequence only and -n for with quality score
my $mis=3; #mismatch
my $maxhit=400; #max hits
my $index='/home/wglabpred/Desktop/bowtie-0.12.7/indexes/ceWS215_all'; #ref index
#my $rec1=`export PATH=\$PATH:/home/wayf/Desktop/bowtie0127`; #bowtie path
my $mtbs=14;#mutation offset
my $min=17; #minimum size
#2 for gene analysis and 3 for transcript analysis
my $filum='yes'; # remove 21U or miRNA from protein coding

#############################################################
#add samples and run here

foreach my $sam (@sam) {
	my $rec1=`bowtie $av $mis -p 4 -f -B 1 -a --best --strata -m $maxhit $index fastaq/$sam/$sam\_$min\_uni.txt bowtie0127/$sam\_all\_$min\_all.temp`;

	unless (-e "$sam\_all") {
		mkdir "$sam\_all";
	}

	$rec1=`cp $0 $sam\_all`;
	#mkdir "$sam\_per";    
}

foreach my $sam (@sam) {
	&filter_major ("$sam\_all");
	foreach my $gr (qw(genome cdna intron rep celrep)) {
		&fix_coor ("$sam\_all",$gr);
	}                            
	&fix_cdna("$sam\_all");
	#&filter_perf ($sam);
                            
	&run_post ("$sam\_all");
	#&run_post ("$sam\_per");
}
#############################################################



# # # # # # # # # # # # # # all subs# # # # # # # # # # # # # # #

#### run post match
sub run_post {
	my ($sam)=@_;
	#log file, reopen it to avoid redundancy.
	open (hand111,">$sam/sum\_$sam.txt");
	print hand111 "Sample name: $sam
Sequence only (-v) or with quality score (-n): $av
Max seed mismatch: $mis\nMax hits: $maxhit
Reference index: $index\nMinimum size: $min
Gene or transcript: 2 for gene and 3 for transcript
Mutation base: $mtbs,int((length(read-$mtbs)/sqrt(mtbs))
Minimum read size: $min
Filter mi/piRNA for counting 22G: $filum
#####################################################\n";

	#sum before matching
	print hand111 sum_fasta($sam);

	# merge cdna and genome
	my %hitno; &merge($sam,\%hitno);

	#repeat-normalization
	foreach my $gr (qw(geex cdna intron rep celrep)) {
		&norm_rep ($sam,$gr,\%hitno);
	} 
	%hitno=();

	# count for each major category like cDNA,genome, repeat, etc
	print hand111 "\nReads matching each category of genomic loci:
	Type\tsense\tanti\tall\n";
	foreach  my $gr (qw(geex cdna intron rep)) {
		#folder,file,sort start,id str,fea str,read,hit
		my $sm=count_all($sam,$gr,2,6,1,4,5);
		print hand111 "$gr\t$sm\n";
	}
	my $sm=count_sim ($sam,'celrep',6,1,4,5);
	print hand111 "celrep\t$sm\n";

	# sum each category of cDNA - individually and wholly
	my %sss; &cdna_split($sam, \%sss,2);
	print hand111 "#####################################################\n";
	print hand111 "\nReads matching each category of cDNA:\nType\tsense\tanti\tall\n";
	my @cdna=(qw(coding ncRNA struc miRNA u21RNA pseudo pre-miRNA star-miRNA exex));
	foreach my $gr (@cdna) {
		&cdna_sum($sam,$gr,\%sss,2);
		my $sm=count_all($sam,"cdna\_$gr",2,6,1,4,5);
		print hand111 "$gr\t$sm\n";
	}
	%sss=(); &cdna_split($sam, \%sss,3);
	foreach my $gr (@cdna) {
		&cdna_sum($sam,$gr,\%sss,3);
	}
	%sss=();

	#sum celrep
	&count_celrep($sam);

	#sum the size and start of each major category
	foreach  my $gr (qw(geex cdna_coding cdna_u21RNA cdna_miRNA rep)) {
		print hand111 sum_loci($sam,$gr,4,5);
	}
	close hand111;

	#this generate 21U sum only covering 0 to 22 nt positions
	#this won't generate files containing RNA read information
	&get_ty1_ty2 ($sam,2); &get_ty1_ty2 ($sam,3);

	#make gff
	&gff($sam);
	opendir (dir1,$sam);
	while (my $file=readdir(dir1)) {
		next if -d "$sam/$file";
		unlink "$sam/$file" if $file =~ /\.std$/;
	}
	closedir dir1;
}

####split into all, per, and mut, and filter out non-specific
sub filter_major {
	my ($sam)=@_;
	my %tmp; my $id='?' x 10000; my $minsc=1000000000;
	open (hand1,"bowtie0127/$sam\_$min\_all.temp") or die $!;
	open (hand2,">$sam/$sam\_$min\_genome.temp");
	open (hand3,">$sam/$sam\_$min\_cdna.temp");
	open (hand4,">$sam/$sam\_$min\_intron.temp");
	open (hand5,">$sam/$sam\_$min\_rep.temp");
	open (hand6,">$sam/$sam\_$min\_celrep.temp");
	while (<hand1>) {
		$_ =~ s/\s+$/\n/;
		my @a1=split /\t/;

		#count the mutation number and filter
		my $sc=0;
		if (defined $a1[7]) {
			my @mt=$a1[7]=~ /\,/g;
			$sc=$#mt+2;
		}
		my $ba=sqrt (length($a1[4]));
		next if $sc > int((length($a1[4])-$mtbs)/$ba);  

		#find the smallest $sc for a read
		if ($a1[0] eq $id) {
			if ($sc>$minsc) {
				next;
			}
			elsif ($sc<$minsc) {
				$minsc=$sc;
				%tmp=();
				$tmp{$a1[2]}=$_;
			}
			else {
				$tmp{$a1[2]}.=$_;
			}
		}

		#print the best read matches and reset 
		else {
			foreach my $ty (keys %tmp) {
				if ($ty =~ /genome/ ) {
					print hand2 $tmp{$ty};
				} 
				elsif ($ty =~ /cdna/ ) {
					print hand3 $tmp{$ty};
				}
				if ($ty =~ /intron/ ) {
					print hand4 $tmp{$ty};
				} 
				if ($ty =~ /\|rep/ ) {
					print hand5 $tmp{$ty};
				}
				if ($ty =~ /celrep/ ) {
					print hand6 $tmp{$ty};
				}
			}       
			$id=$a1[0];$minsc=$sc;
			%tmp=(); $tmp{$a1[2]}=$_;
		}
	}

	foreach my $ty (keys %tmp) {
		if ($ty =~ /genome/ ) {
			print hand2 $tmp{$ty};
		} 
		elsif ($ty =~ /cdna/ ) {
			print hand3 $tmp{$ty};
		}
		if ($ty =~ /intron/ ) {
			print hand4 $tmp{$ty};
		} 
		if ($ty =~ /\|rep/ ) {
			print hand5 $tmp{$ty};
		}
		if ($ty =~ /celrep/ ) {
			print hand6 $tmp{$ty};
		}
	} 
	close hand1; close hand2; close hand3;
	close hand4; close hand5; close hand6;
	unlink "bowtie0127/$sam\_$min\_all.temp"; %tmp=();
}

####split all matches
sub fix_coor {
	my ($sam,$gr)=@_; my %id; &ref_id(\%id); my %mid; 
	my $oid='?' x 10000; 
	my $grt=$gr; $grt.='_step1' if $gr eq 'cdna';

	#reformat according to genome coordinates
	my $rec1=`sort -k3,3 $sam/$sam\_$min\_$gr.temp >$sam/$sam\_$min\_$gr.sorted`;
	unlink "$sam/$sam\_$min\_$gr.temp";

	open (hand1,"$sam/$sam\_$min\_$gr.sorted") or die $!;
	open (hand2,">$sam/$sam\_$min\_$grt.txt");
	while (<hand1>) {
		$_ =~ s/\s+$//; my @a1=split /\t/;
		my ($ch,$ty)=split /\|/,$a1[2]; $mid{$ch}=1;
		my $nam=$id{$gr}{$ch}; my @ch=split (/\t/,$nam); #id info
		$ch=$nam; $ch=$ch[5] if defined $ch[5];

		my $offset=0; # 0 for genome, cdna, celrep
		$offset=$ch[7]-1 if $gr eq 'intron' || $gr eq 'rep';

		my ($id,$re)=split (/_x/,$a1[0]);#read and id
		my $sr=$a1[1]; my $st=$a1[3]+$offset;
		my $en=length($a1[4])-1; my $rn=$a1[4];
		if ($a1[1] eq "-") {
			$rn=reverse $rn;
			$rn =~ tr/AGCT/TCGA/;
		}

		#fix mutation pos according to genome
		my $mut='N';
		if (defined $a1[7]) {
			my @mut=split /\,/,$a1[7];
			foreach my $m (@mut) {
				my ($po,$mt)=split (/\:/,$m);
				$mt=~ s/>//; my $mtpo=$po;
				$mtpo=$en-$po if $a1[1] eq '-';
				$mut.="\,".$mtpo."\:".$mt;
			}
			$mut=~ s/N\,//;
		}

		if ($nam ne $oid) {
			print hand2 ">$nam\n";
			$oid=$nam;
		}
		print hand2 "$ch\t$sr\t$st\t$en\t$re\t$id\t$rn\t$mut\n";
	}
	close hand1; unlink "$sam/$sam\_$min\_$gr.sorted"; 

	foreach my $id (keys %{$id{$gr}}) {
		next if exists $mid{$id};
		print hand2 ">$id{$gr}{$id}\n";
	}
	close hand2; %id=(); %mid=();
}

####fix cDNA coordinates and split ex-ex junctions
sub fix_cdna {
	my ($sam)=@_; my %expos;
	#genome coordinate intervals of cdna annotations
	open (hand1,'ceWS215_cdna_pos_final_mu_4.txt') or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my ($id,$pos)=split /\|/;
		@{$expos{'>'.$id}}=split (/ /,$pos);
	}
	close hand1;                                        

	#fix pos and exex
	my @expos;
	open (hand1,"$sam/$sam\_$min\_cdna_step1.txt") or die $!;
	open (hand2,">$sam/$sam\_$min\_cdna.txt");
	while (<hand1>) {
		$_ =~ s/\s+$//; my @a1=split /\t/; 
		if (/^>/) {
			@expos=();
			foreach my $in (@{$expos{$_}}) {
				my ($l,$r)=split (/\t/,$in);
				for my $p ($l..$r) {
					push (@expos,$p);
				}
		}
		@expos=sort {$a<=>$b} @expos;
		print hand2 $_,"\n"; next;
	}

	my $st=$expos[$a1[2]-1];
	my $en=$expos[$a1[2]+$a1[3]-1]-$st;

	my $mut='N';
	if ($a1[7] =~ /\:/) {
		my @mut=split (/\,/,$a1[7]);
		foreach my $m (@mut) {
			my ($po,$mt)=split (/\:/,$m);
			my $mtpo=$expos[$a1[2]+$po-1]-$st;
			$mut.="\,".$mtpo."\:".$mt;
		}
		$mut=~ s/N\,//;
	}

	#split ex-ex
	my @ex; my $end3=$a1[2]+$a1[3];
	my $dif=$expos[$end3-1]-$expos[$a1[2]-1]-$a1[3];
	if ($dif >0) {
		for my $p ($a1[2]..$end3) {
			push (@ex,$expos[$p-1]);
		}
		@ex=sort {$a<=>$b} @ex;
		$mut.="\t".$ex[0];
		for (my $i=0; $i<$#ex;$i++) {
			if ($ex[$i+1]-$ex[$i] >1) {
				my $r=$ex[$i]-$st;
				my $l=$ex[$i+1]-$st;
				$mut.='-'.$r.'_'.$l;
			}
		}
		my $last=$ex[-1]-$st;
		$mut.='-'.$last;
		}
		print hand2 "$a1[0]\t$a1[1]\t$st\t$en\t$a1[4]\t$a1[5]\t$a1[6]\t$mut\n";
    }
	close hand1; close hand2; %expos=(); @expos=();
	unlink "$sam/$sam\_$min\_cdna_step1.txt";
}

sub ref_id {
	my ($id)=@_;
	open (hand1,"ceWS215_all_id_type.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		my ($t,$i,$n)=split /\|/;
		${$id}{$t}{$i}=$n;
	}
	close hand1;
}

sub filter_perf {
	my ($sam)=@_;
	foreach my $ty(qw(genome cdna intron rep celrep)) {
		open (hand1,"$sam\_all/$sam\_all\_$min\_$ty.txt") or die $!;
		open (hand2,">$sam\_per/$sam\_per\_$min\_$ty.txt");
		while (<hand1>) {
			$_ =~ s/\s+$//; my @a1=split /\t/;
			print hand2 "$_\n" if /^>/;
			next unless $a1[7] eq 'N' && !/^>/;
			print hand2 $_,"\n";
		}       
		close hand1; close hand2;
	}
}

#### sum fasta
sub sum_fasta {
	my ($sam)=@_; my $id; my $read; my %size; my %sta; my $tot=0;
	$sam =~ s/\_(all|mut|per).*//;
	open (hand1,"fastaq/$sam/$sam\_$min\_uni.txt") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		if (/>/) {
			($id,$read)=split /_x/;
			next;
		}
		$size{length($_)}+=$read;
		$sta{substr($_,0,1)}+=$read;
	}
	close hand1;

	my $prt="\nSize before matching:\n";
	foreach my $size (sort {$a<=>$b} keys %size) {
		$tot+=$size{$size};
		$prt.="$size\t$size{$size}\n";
	}
	%size=();

	$prt.="\n1st nt before matching:\n";
	foreach my $nt (sort {$a cmp $b} keys %sta) {
		$prt.="$nt\t$sta{$nt}\n";
	}
	$prt.="\nTotal reads before matching: $tot\n";
	$prt.="#####################################################\n";
	%sta=(); return $prt;
}

#### merge and sum
sub merge {
	my ($sam,$hino)=@_; my %hit;
	open (hand1,"$sam/$sam\_$min\_genome.txt") or die $!;
	open (hand2,">$sam/$sam\_$min\_geex.txt");
	while (<hand1>) {
		my @a1=split /\t/; next if />/;
		print hand2 $_; ${$hino}{$a1[5]}++;
	}
	close hand1; unlink "$sam/$sam\_$min\_genome.txt";

	open (hand1,"$sam/$sam\_$min\_cdna.txt") or die $!;
	while (<hand1>) {
		my @a1=split /\t/;
		next unless !/>/ && $#a1 ==8 && !exists $hit{$_};
		print hand2 $_;
		${$hino}{$a1[5]}++; $hit{$_}=1;
	}
	close hand1; close hand2; %hit=();
}

#### normalization
sub norm_rep {
	my ($sam,$type,$hit) =@_;
	open (hand1,"$sam/$sam\_$min\_$type.txt") or die $!;
	open (hand2,">$sam/$sam\_$min\_$type\_normed.txt");
	if ($type eq 'geex') {
		open (hand3,">$sam/$sam\_$min\_$type\_normed.fff");
	}
	while (<hand1>) {
		if (/>/) {              
			print hand2 $_; next;
		}
		$_ =~ s/\s+$//; my @a1=split /\t/;
		my $ht =1; $ht=${$hit}{$a1[5]} if $type ne 'celrep' && $type ne 'est';
		print hand2 "$a1[0]\t$a1[1]\t$a1[2]\t$a1[3]\t$a1[4]\t$ht\t$a1[5]\t$a1[6]\t$a1[7]";
		my $lst="\n";
		my $en=$a1[2]+$a1[3]; my $se=$a1[2].'-'.$en;
		if (defined $a1[8]) {
			$lst="\t$a1[8]\n";
			my @ex=split /_|-/,$a1[8];
			$ex[0]=0; $se='';
			foreach my $p (@ex) {
				my $pp=$p+$a1[2];
				$se.=$pp.'-';
			}
			$se =~ s/-$//;
		}
		print hand2 $lst;
		next unless $type eq 'geex';
		print hand3 "$a1[0]\t$a1[1]\t$se\t$a1[4]\t$ht\n";
	}
	close hand3 if $type eq 'geex'; close hand1; close hand2;
	unlink "$sam/$sam\_$min\_$type.txt";
}

#### count file input
sub count_sim {  
	my ($di,$gr,$gs,$fs,$re,$ht)=@_;
	my $se=0;my $an=0;my %se=();my %an=();my $ssr='+';
	my $fi="$di\_$min\_$gr\_normed.txt";
	open (hand1, "$di/$fi") or die "no $fi\n";
	while (<hand1>) {
		my @a1=split /\t/;
		if (/>/) {
			$ssr=$a1[$gs];next;
		}    
		my $no=reformat_norm($a1[$re]/$a1[$ht]);    

		if ($a1[$fs] eq $ssr) {
			$se{$a1[6]}=$no; next;
		}
		$an{$a1[6]}=$no;
	}
	close hand1;

	foreach my $f (keys %se) {
		if (exists $an{$f}) {
			$se{$f}/=2;$an{$f}/=2;
		}     
		$se+=$se{$f};
	}

	foreach my $f (keys %an) {
		$an+=$an{$f};
	}
	%an=(); %se=();

	$an=int($an+0.5);$se=int($se+0.5);
	my $prt=$se+$an;
	$prt="$se\t$an\t$prt";
	return $prt;
}

#### count file input
sub count_all {  
	my ($di,$gr,$fst,$gs,$fs,$re,$ht)=@_; my $s=0; my $a=0;
	my $oid='?' x 10000; my $se=0; my $an=0; my $ssr='+';
	my $fi="$di\_$min\_$gr\_normed.txt";
	open (hand1, "$di/$fi") or die $!;
	open (hand2, ">$di/$fi.rs");
	while (<hand1>) {
		my @a1=split /\t/;
		if (/>/) {
			$ssr=$a1[$gs];next;
		} 
		if ($a1[$fs] eq $ssr) {
			print hand2 "+\t\|$_"; next;
		}
		print hand2 "-\t\|$_";
	}
	close hand1; close hand2;

	my $rec1=`sort -k$fst $di/$fi.rs >$di/$fi.std`;
	unlink "$di/$fi.rs";

	open (hand1,"$di/$fi.std");
	while (<hand1>) {
		my ($sr,$f)=split /\|/;
		my @f=split /\t/,$f;
		if ($oid eq $f) {
			if ($sr eq "+\t") {
				$s=$f[$re]/$f[$ht];
			}
			else {
				$a=$f[$re]/$f[$ht];
			}
		}

		else {
			if ($oid ne '?' x 10000) {
				if ($s>0 && $a>0) {
					$s/=2; $a/=2;
				}
				$se+=$s; $an+=$a;
			}
			$a=0; $s=0; $oid=$f;
			if ($sr eq "+\t") {
				$s=$f[$re]/$f[$ht];
			}
			else {
				$a=$f[$re]/$f[$ht];
			}
		}
	}

	if ($oid ne '?' x 10000) {
		if ($s>0 && $a>0) {
			$s/=2; $a/=2;
		}
		$se+=$s; $an+=$a;
	}
	close hand1;
	$se=int($se+0.5);$an=int($an+0.5);
	my $prt=$se+$an;
	$prt="$se\t$an\t$prt";
	return $prt;
}

#### split cDNA
sub cdna_split {
	my ($sam,$sss,$getr)=@_; my %ty=(); my %fil=();
	$ty{'u21RNA'}='u21RNA';$ty{'pre-miRNA'}='pre-miRNA';
	$ty{'pseudogene'}='pseudo'; $ty{'RNA_pseudogene'}='pseudo';
	$ty{'Coding_pseudogene'}='pseudo';$ty{'snRNA'}='struc';
	$ty{'scRNA'}='struc';$ty{'snoRNA'}='struc';$ty{'tRNA'}='struc';
	$ty{'rRNA'}='struc';$ty{'ncRNA'}='ncRNA'; $ty{'miRNA'}='miRNA';
	$ty{'star-miRNA'}='star-miRNA';$ty{'coding'}='coding';
	foreach my $cy (keys %ty) {
		open (hand1,">$sam/$sam\_$min\_cdna\_$ty{$cy}\_normed.txt");
		open (hand2,">$sam/$sam\_$min\_cdna\_$ty{$cy}\_normed.tmp");
		close hand1; close hand2;
	}

	# filter out miRNA and 21U-RNA from coding
	&filter_id("$sam/$sam\_$min\_cdna_normed.txt",6,\%fil) if $filum eq 'yes';

	my %all; my $id; my $ca; my $li=0;
	open (hand1,"$sam/$sam\_$min\_cdna\_normed.txt") or die $!;
	open (hand4,">$sam/$sam\_$min\_cdna_exex\_normed.txt");
	open (hand5,">$sam/$sam\_$min\_cdna_exex\_normed.tmp");
	while (<hand1>) {
		# print if too big
		if (/>/) {
			my @a1=split /\t/;
			if ($li >500) {
				foreach my $cy (keys %all) {
					open (hand2,">>$sam/$sam\_$min\_cdna\_$cy\_normed.txt");
					open (hand3,">>$sam/$sam\_$min\_cdna\_$cy\_normed.tmp");
					foreach my $id (keys %{$all{$cy}}) {
						print hand2 $id; my @id=split /\t/,$id;
						${$sss}{$cy}{"$id[1]\t$id[$getr]"}="0\t0\t0\n";
						if ($cy eq 'coding') {
							print hand4 $id;
							${$sss}{'exex'}{"$id[1]\t$id[$getr]"}="0\t0\t0\n";
						}
						foreach my $ff (keys %{$all{$cy}{$id}}) {
							my @ff=split /\t/,$ff;
							next if $ff eq 'empty';
							if (defined $ff[9] && $cy eq 'coding') {
								print hand4 $ff;
								print hand5 "$id[1]\t$id[$getr]\t$id[6]\t\|$ff";
							}
							next if $cy eq 'coding' && exists $fil{$ff[6]};
							print hand3 "$id[1]\t$id[$getr]\t$id[6]\t\|$ff";
							print hand2 $ff;
						}
					}
					close hand2; close hand3;
				}
				%all=(); $li=0;
			}

			$li++; $ca=$ty{$a1[4]};$id=$_;$all{$ca}{$id}{'empty'}=1;
		}
		else {
			$all{$ca}{$id}{$_}=1;
		}
	}

	foreach my $cy (keys %all) {
	open (hand2,">>$sam/$sam\_$min\_cdna\_$cy\_normed.txt");
	open (hand3,">>$sam/$sam\_$min\_cdna\_$cy\_normed.tmp");
	foreach my $id (keys %{$all{$cy}}) {
		print hand2 $id; my @id=split /\t/,$id;
		${$sss}{$cy}{"$id[1]\t$id[$getr]"}="0\t0\t0\n";
		if ($cy eq 'coding') {
			print hand4 $id;
			${$sss}{'exex'}{"$id[1]\t$id[$getr]"}="0\t0\t0\n";
		}
		foreach my $ff (keys %{$all{$cy}{$id}}) {
			my @ff=split /\t/,$ff;
			next if $ff eq 'empty';
			if (defined $ff[9] && $cy eq 'coding') {
				print hand4 $ff;
				print hand5 "$id[1]\t$id[$getr]\t$id[6]\t\|$ff";
                                      }
				next if $cy eq 'coding' && exists $fil{$ff[6]};
				print hand3 "$id[1]\t$id[$getr]\t$id[6]\t\|$ff";
				print hand2 $ff;
			}
		}
		close hand2; close hand3;
	}   
	%all=(); %ty=(); %fil=(); close hand1; close hand4; close hand5;
}

####sum each gene or transcript
sub cdna_sum {
	my ($sam,$ty,$sss,$getr)=@_; my %se; my %an; my $oid='?' x 10000;
	my %sum=%{${$sss}{$ty}};
	my $rec1=`sort -k1,2 $sam/$sam\_$min\_cdna\_$ty\_normed.tmp >$sam/$sam\_$min\_cdna\_$ty\_normed.sorted`;
	open (hand1,"$sam/$sam\_$min\_cdna\_$ty\_normed.sorted");
	while (<hand1>) {
		$_=~ s/\s+$//; my @a1=split /\t/;
		my ($g,$f)=split /\|/;
		if ($oid eq "$a1[0]\t$a1[1]") {
			if ($a1[2] eq $a1[4]) {$se{$f}=1;
			}
			else {
				$an{$f}=1;
			}
		}

		else {
			my $se=0; my $an=0; my $al=0;
			&count_single(\%se,\%an,\$se,\$an);
			$al=$se+$an;
			if ($oid ne "?" x 10000) {
				$sum{$oid}="$se\t$an\t$al\n";
			}

			%an=(); %se=();
			if ($a1[2] eq $a1[4]) {
				$se{$f}=1;
			}
			else {
				$an{$f}=1;
			} 
			$oid="$a1[0]\t$a1[1]";
		}
	}

	if ($oid ne "?" x 10000) {
		my $se=0; my $an=0; my $al=0;
		&count_single(\%se,\%an,\$se,\$an);
		$al=$se+$an;
		$sum{$oid}="$se\t$an\t$al\n";
	}
	close hand1; unlink "$sam/$sam\_$min\_cdna\_$ty\_normed.sorted";
	unlink "$sam/$sam\_$min\_cdna\_$ty\_normed.tmp";

	#sort and print sum for each gene
	open (hand1,">$sam/sum\_$sam\_$min\_cdna$getr\_$ty.txt");
	foreach my $na (sort {$a cmp $b} keys %sum) {
		print hand1 "$na\t$sum{$na}";
	}
	close hand1;
}

#### filter
sub filter_id {
	my ($file,$id,$fil)=@_; my $on=0;
	open (hand1,$file) or die $!;
	while (<hand1>) {
		my @a1=split /\t/;
		if (/^>/) {
			$on=0; $on=1 if $a1[4] eq 'miRNA' || $a1[4] eq 'u21RNA';
			next;
		}
		next if $on ==0; ${$fil}{$a1[$id]}=1;
	}
	close hand1; my $count=scalar keys %{$fil};
	print "miRNA/21U id used for filtering: $count\n";
}

####format_gff                                                    
sub gff {
	my ($sa)=@_;my $ra=norm($sa,$stan); my $j=0; my %ge; my $id;
	open (hand1,"ceWS215_gdna.fa") or die $!;
	while (<hand1>) {
		$_ =~ s/\s+$//;
		if (/^>/) {
			$id=$_; $id =~ s/^>//;
			next;
		}
		$ge{$id}.=$_;
	}
	close hand1;                
	my $rec1=`sort -k1,3 $sa/$sa\_$min\_geex\_normed.fff >$sa/$sa\_$min\_geex\_normed.sff`;
	#unlink "$sa/$sa\_$min\_geex\_normed.fff";

	open (hand0,"$sa/$sa\_$min\_geex\_normed.sff") or die $!;
	open (hand1, ">$sa/$sa\_$min\_geex\_$stan.gff");
	my $s=substr $sa,0,length($sa)-2; my $ol='?' x 10000;
	my $re=0; my $c; my $sr; my $ex;
	while (<hand0>) {
		$_ =~ s/\s+$//; my @a1=split /\t/; my $nl="$a1[0]\t$a1[1]\t$a1[2]";
		my $no=int($a1[3]/$a1[4]*1000+0.5)/1000;
		if ($nl eq $ol) {
			$re+=$no
		}
		else {
			if ($ol ne '?' x 10000) {
				my $ph=0; my $se; my @e=split /-/,$ex;
				foreach (my $p=0;$p<$#e;$p+=2) {
					$se.=substr $ge{$c},$e[$p]-1,$e[$p+1]-$e[$p]+1;
				} 
				if ($sr eq '-') {
					$se =reverse $se;
					$se=~ tr/AGCT/TCGA/;
				}

				$ph =1 if length($se) ==26 && $se =~ /^G/;
				$ph =2 if length($se) ==21 && $se =~ /^T/;

				$re*=$ra;
				if ($re >= 1) {
					$re=int($re*10+0.5)/10; $j++; my $id='R'.$s.$j;     

					if ($ex =~ /-\d+-/) {
						print hand1 "$c\t$s\tfr\t$e[0]\t$e[-1]\t$re\t$sr\t$ph\tID $id\n";
                                
						for (my $p=0; $p <$#e; $p+=2) {
							print hand1 "$c\t$s\tpr\t$e[$p]\t$e[$p+1]\t$re\t$sr\t$ph\tID $id\n";
						}
					}
					else {
						print hand1 "$c\t$s\tex\t$e[0]\t$e[-1]\t$re\t$sr\t$ph\tID $id\n";
					}
				}
			}
			$ol=$nl; $re=$no; $c=$a1[0]; $sr=$a1[1]; $ex=$a1[2];
		}
	}

	if ($ol ne '?' x 10000) {
		my $ph=0; my $se; my @e=split /-/,$ex;
		foreach (my $p=0;$p<$#e;$p+=2) {
			$se.=substr $ge{$c},$e[$p]-1,$e[$p+1]-$e[$p]+1;
		} 
		if ($sr eq '-') {
			$se =reverse $se;
			$se=~ tr/AGCT/TCGA/;}

			$ph =1 if length($se) ==26 && $se =~ /^G/;
			$ph =2 if length($se) ==21 && $se =~ /^T/;

			$re*=$ra;
			if ($re >= 1) {
				$re=int($re*10+0.5)/10; $j++; my $id='R'.$s.$j;     

				if ($ex =~ /-\d+-/) {
					print hand1 "$c\t$s\tfr\t$e[0]\t$e[-1]\t$re\t$sr\t$ph\tID $id\n";
                                
					for (my $p=0; $p <$#e; $p+=2) {
						print hand1 "$c\t$s\tpr\t$e[$p]\t$e[$p+1]\t$re\t$sr\t$ph\tID $id\n";
					}
				}
				else {
					print hand1 "$c\t$s\tex\t$e[0]\t$e[-1]\t$re\t$sr\t$ph\tID $id\n";
			}
		}
	}
	close hand1, close hand0; %ge=(); unlink "$sa/$sa\_$min\_geex\_normed.sff";
}

#### count celrep
sub count_celrep {
	my ($sam)=@_; my $oid; my $sr='+'; my $se; my %sum;
	open (hand1,"$sam/$sam\_$min\_celrep_normed.txt") or die $!;
	open (hand2,">$sam/$sam\_$min\_celrep_normed.sht");
	while (<hand1>) {
		my @a1=split /\t/;
		if (/^>/) {
			$sr=$a1[6];
			$oid="$a1[1]\t$a1[2]";
			$sum{$oid}="0\t0\t0\n";
			next;
		}

		if ($a1[1] eq $sr) {
			print hand2 "$oid\t+\t";
		}
		else {
			print hand2 "$oid\t-\t";
		}
		print hand2 "$a1[4]\t$a1[6]\n";
	}
	close hand1; close hand2;

	my $rec1=`sort -k1,2 $sam/$sam\_$min\_celrep_normed.sht >$sam/$sam\_$min\_celrep_normed.st`;
	unlink "$sam/$sam\_$min\_celrep_normed.sht"; my %se; my %an; $oid='?' x 10000;
	open (hand1,"$sam/$sam\_$min\_celrep\_normed.st");
	while (<hand1>) {
		$_=~ s/\s+$//; my @a1=split /\t/;
		if ($oid eq "$a1[0]\t$a1[1]") {
			if ($a1[2] eq '+') {
				$se{$a1[4]}=$a1[3];
			}
			else {
				$an{$a1[4]}=$a1[3];
			}
		}

		else {
			if ($oid ne "?" x 10000) {
				my $se=0; my $an=0; my $al=0;
				foreach my $rn (keys %se) {
					if (exists $an{$rn}) {
						$se+=$se{$rn}*0.5;
						$an{$rn}*=0.5;
					}
					$se+=$se{$rn};
				}

				foreach my $rn (keys %an) {
					$an+=$an{$rn};
				}

				$se=reformat_norm ($se);
				$an=reformat_norm ($an);
				$al=reformat_norm ($se+$an);
				$sum{$oid}="$se\t$an\t$al\n";
			}

			%an=(); %se=(); $oid="$a1[0]\t$a1[1]"; 
			if ($a1[2] eq '+') {
				$se{$a1[4]}=$a1[3];
			}
			else {
				$an{$a1[4]}=$a1[3];
			}
		}
	}

	if ($oid ne "?" x 10000) {
		my $se=0; my $an=0; my $al=0;
		foreach my $rn (keys %se) {
			if (exists $an{$rn}) {
				$se{$rn}*=0.5;
				$an{$rn}*=0.5;
			}
			$se+=$se{$rn};
		}

		foreach my $rn (keys %an) {
			$an+=$an{$rn};
		}

		$se=reformat_norm ($se);
		$an=reformat_norm ($an);
		$al=reformat_norm ($se+$an);
		$sum{$oid}="$se\t$an\t$al\n";
	}
	close hand1; unlink "$sam/$sam\_$min\_celrep\_normed.st";

	open (hand1,">$sam/sum\_$sam\_$min\_celrep.txt");
	foreach my $gene (sort {$a cmp $b} keys %sum) {
		print hand1 "$gene\t$sum{$gene}";
	}
	close hand1;
}

####count_single
sub count_single {
	my ($s,$a,$se,$an)=@_;
	foreach my $f (keys %{$s}) {
		my @a1=split /\t/,$f;
		if (exists ${$a}{$f}) {
			$$se+=$a1[4]/$a1[5]*0.5;
		}
		else {
			$$se+=$a1[4]/$a1[5];
		}
	}

	foreach my $f (keys %{$a}) {
		my @a1=split /\t/,$f;
		if (exists ${$s}{$f}) {
			$$an+=$a1[4]/$a1[5]*0.5;
		}
		else {
			$$an+=$a1[4]/$a1[5];
		}
	}

	$$se=reformat_norm($$se);
	$$an=reformat_norm($$an);
}

#### sum loci
sub sum_loci {
	my ($di,$gr,$re,$ht)=@_; my %ls;my %la; my %ns;my %na; 
	my $fi="$di\_$min\_$gr\_normed.txt"; my $oid='?' x 10000;
	my $os; my $s=0; my $r=0;
	open (hand1,"$di/$fi.std");
	while (<hand1>) {
		my ($sr,$f)=split /\|/;
		my @f=split /\t/,$f;
		if ($oid eq $f) {
			if ($sr eq "+\t") {
				$s=$f[$re]/$f[$ht];
			}
			else {
				$r=$f[$re]/$f[$ht];
			}
		}

		else {
			if ($oid ne '?' x 10000) {
				my $si=length($os);
				my $nt=substr($os,0,1);
				if ($s>0 && $r>0) {
					$s/=2; $r/=2;
				}
				$ls{$si}+=$s;$ns{$nt}+=$s; 
				$la{$si}+=$r;$na{$nt}+=$r;
			}
			$r=0; $s=0; $oid=$f; $os=$f[7];
			if ($sr eq "+\t") {
				$s=$f[$re]/$f[$ht];
			}
			else {
				$r=$f[$re]/$f[$ht];
			}
		}
	}

	if ($oid ne '?' x 10000) {
		my $si=length($os);
		my $nt=substr($os,0,1);
		if ($s>0 && $r>0) {
			$s/=2; $r/=2;
		}
		$ls{$si}+=$s;$ns{$nt}+=$s;
		$la{$si}+=$r;$na{$nt}+=$r;
	}
	close hand1; unlink "$di/$fi.std";

	my $prt="#####################################################\n";
	$prt.="\n$gr size matched:\ntype\tsense\tanti\tall\n";
	foreach my $l (sort {$a<=>$b} keys %ls) {
		my $al=int($ls{$l}+$la{$l}+0.5);
		my $se=int($ls{$l}+0.5);
		my $an=int($la{$l}+0.5);
		$prt.="$l\t$se\t$an\t$al\n";
	}

	$prt.="\n$gr 1st nt matched:\ntype\tsense\tanti\tall\n";
	foreach my $n (sort {$a cmp $b} keys %ns) {
		my $al=int($ns{$n}+$na{$n}+0.5);
		my $se=int($ns{$n}+0.5);
		my $an=int($na{$n}+0.5);
		$prt.="$n\t$se\t$an\t$al\n";
	}
	%ls=();%la=();%ns=();%na=(); return $prt;
}

sub reformat_norm {
	my ($no)=@_;
	$no =sprintf ("%.4f",$no);
	$no =~ s/\.?0+$//; return $no;
}

sub get_ty1_ty2 {
	my ($s,$getr)=@_; my %u; my %id; my %re;
	open (my $fh,"ce_type2_more_type1_21U_pos.txt") or die $!;
	while (<$fh>) {
		my ($id,$po)=split /\|/;
		my @id=split /\t/,$id;
		my $st=$id[7];
		if ($id[6] eq '-') {
			$st=$id[8];
		}
		$id{"$id[5]\t$id[6]\t$st"}="$id[1]\t$id[$getr]";
		$re{"$id[5]\t$id[6]\t$st"}=0;
	}
	close $fh;

	open ($fh,"ceWS215_cdna_pos_final_mu_4.txt");
	while (<$fh>) {
		my ($id,$po)=split /\|/;
		my @id=split /\t/,$id;
		next if $id[4] ne 'u21RNA';
		my $st=$id[7]+4;
		if ($id[6] eq '-') {
			$st=$id[8]-4;
		}
		$id{"$id[5]\t$id[6]\t$st"}="$id[1]\t$id[$getr]";
		$re{"$id[5]\t$id[6]\t$st"}=0;
	}
	close $fh;

	open ($fh, "$s/$s\_17_geex_normed.txt") or die $!;
	while (<$fh>) {
		my @a1=split /\t/;
		next if $a1[3] >21;
		my $st=$a1[2];
		if ($a1[1] eq '-') {
			$st=$a1[2]+$a1[3];
		}
		next if !exists $id{"$a1[0]\t$a1[1]\t$st"};
		$re{"$a1[0]\t$a1[1]\t$st"}+=$a1[4]/$a1[5];
	}
	close $fh;

	foreach (keys %re) {
		$u{$id{$_}}=int ($re{$_}*10000+0.5)/10000;
	}

	open ($fh,">$s/sum_$s\_$min\_cdna$getr\_all21U.txt");
	foreach my $u (sort {$a cmp $b} keys %u) {
		print $fh "$u\t$u{$u}\t0\t$u{$u}\n";
	}
	close $fh;
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
print "$sam\ntotal $geex\nsense struc $struc
non-struc $nons\nsense coding $cose
anti coding $coan\ntotal coding $coto
sense miRNA $mir\nsense 21U-RNA $u21\n";
       
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
	print "Normalization ratio $ra\n\n";
	return $ra;
}
