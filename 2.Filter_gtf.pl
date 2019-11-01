#! /usr/bin/perl -w
use warnings;
#use strict;
use Getopt::Std;
use POSIX;
use List::MoreUtils qw/uniq/;

my $usage =
"$0 <-a loc.fsim.txt -i merged.gtf> [options]

This script can filter merged gtf annotation file by junction reads

-a        loc.fsim.txt include the absolute path of the loc.fsim file form the bam files which produce gtf, the standardized coefficient is placed in the second column in tab interval,if not [default: 1]
-i        gtf file used to filtering
-o        output file                                   [default: ./merged.filter.gtf]
-m int    minimum filtered junction reads in one sample to support one exon-intron-exon       							    [default: 5]
-r int    minimun intron length                         [default: 30]
-x int    minimum exon length                           [default: 21]
-s int    minimum sample to support all the junction    [default: 1]
-c int    number of replication in each sample		[default: 2]
-n int    minimum exon number                           [default: 1]
-e int    minimum RPKM fot filter                       [default: 0]
";

my %options=();
getopts("a:i:o:m:r:x:s:n:e:c:",\%options);

if (!exists $options{a}||!exists $options{i}){
die "$usage";
}

if (!exists $options{o}){
	$options{o}="./merged.filter.gtf";
}

if (!exists $options{m}){
	$options{m}=5;
}

if (!exists $options{r}){
        $options{r}=30;
}

if (!exists $options{x}){
	$options{x}=21;
}

if (!exists $options{s}){
        $options{s}=1;
}

if (!exists $options{c}){
	$options{c}=2;
}

if (!exists $options{n}){
        $options{n}=1;
}

if (!exists $options{e}){
	$options{e}=0;
}

my $num=1;
open IN,"$options{a}"||die;
open OUT,">$options{o}"||die;
while (my $file=<IN>){
chomp $file;
my @fsim=split(/\s+/,$file);
if(!exists $fsim[1]){
	$fsim[1]=1;
}
	open DA,"$fsim[0]"||die;#reads_map_intron_loc_1,2,3...
	while (my $loc=<DA>){
		chomp $loc;
		my @line=split(/\s+/,$loc);
		if ($line[4] eq "I"){
			my $intron_loc=join("-",$line[0],$line[1],$line[2]);
			my $name= "junc_count_$num";
		 	   ${$name}{$intron_loc}=$line[3]*$fsim[1];
		}
		if ($line[4] eq "E"){
			my $exon_loc=join("-",$line[0],$line[1],$line[2]);
			my $chr_part=POSIX::ceil((($line[1]+$line[2])/2)/1000000);
			my $namein= "exon_count_$num";
			my $chr_exon= "chr_$num";
			my $new_exon_id=join("_",$line[0],$chr_part);
			   ${$namein}{$exon_loc}=$line[3]*$fsim[1];
			   push @{${$chr_exon}{$new_exon_id}},$exon_loc;
		}
	}
close DA;
$num+=1;
}
close IN;


#chr1A_part1     StringTie       exon    19823   21438   1000    .       .       gene_id "ABGENE.3"; transcript_id "ABGENE.3.1"; exon_number "1";
open GTF,"$options{i}"||die;#merged gtf
my @discard;
my (%trans_len,%trans_exon,%gene_all_trans,%exon_sta,%exon_end,%chr_trans,%trans_intron,%dis);
while (my $gtf=<GTF>){
chomp $gtf;
if ($gtf!~/^#/){
	my @tem=split(/\s+/,$gtf);
		if ($tem[2] eq "exon"){
			my @gene_trans=split(/"/,$gtf);
			my $exon_len=$tem[4]-$tem[3]+1;
			my $exon_id=join(";",$tem[0],$tem[3],$tem[4],$tem[6]);
			   $trans_len{$gene_trans[3]}+=$exon_len;
			   $trans_gene{$gene_trans[3]}=$gene_trans[1];
			push @{$trans_exon{$gene_trans[3]}},$exon_id;
				if ($tem[6] eq '.'){
					push @discard,$gene_trans[3];#no strand
					}
				if ($exon_len<$options{x}){
					push @discard,$gene_trans[3];#exon_len<=20
					}
				if (grep {$_ eq $gene_trans[3]} @{$gene_all_trans{$gene_trans[1]}}){
	}else{
				push @{$gene_all_trans{$gene_trans[1]}},$gene_trans[3];
	}
				if (grep {$_ eq $exon_id} @{$gene_all_exon{$gene_trans[1]}}){
	}else{
				push @{$gene_all_exon{$gene_trans[1]}},$exon_id;
	}
				push @{$exon_sta{$gene_trans[3]}},$tem[3];
				push @{$exon_end{$gene_trans[3]}},$tem[4];
				       $chr_trans{$gene_trans[3]}=$tem[0];
				push @{$exon_trans{$exon_id}},$gene_trans[3];
}
		if($tem[2] eq "transcript"){
			my @gene_trans=split(/"/,$gtf);
			my $trans_rpkm=$gene_trans[5];
			if($trans_rpkm<$options{e}){
				push @discard,$gene_trans[3];
		}
}
}
}
close GTF;


foreach $trans_fi (keys %trans_exon){
	my @all_exon=@{$trans_exon{$trans_fi}};
	if (@all_exon>=2){

	if  (grep {$_ eq $all_exon[0]} @{$gene_f_exon{$trans_gene{$trans_fi}}}){#排第一的外显子是否跨越别的内含子
	}else{
		push @{$gene_f_exon{$trans_gene{$trans_fi}}},$all_exon[0];
	}	

	if  (grep {$_ eq $all_exon[-1]} @{$gene_l_exon{$trans_gene{$trans_fi}}}){#排最后的外显子是否跨越别的内含子
	}else{
		push @{$gene_l_exon{$trans_gene{$trans_fi}}},$all_exon[-1];
	}

	for(my $m=0;$m<@all_exon-1;$m++){
		my @id=split(/;/,$all_exon[$m]);
		my @next_id=split(/;/,$all_exon[$m+1]);
		my $intron_start=$id[2]+1;
		my $intron_end=$next_id[1]-1;
		my $intron_len=$intron_end-$intron_start+1;
		if($intron_len<$options{r}){
			push @discard,$trans_fi;
		}
		my $intron_id=join(";",$id[0],$intron_start,$intron_end,$next_id[3]);
		push @{$trans_intron{$trans_fi}},$intron_id;
		
		if  (grep {$_ eq $intron_id} @{$gene_introns{$trans_gene{$trans_fi}}}){
	}else{
                push @{$gene_introns{$trans_gene{$trans_fi}}},$intron_id;
	}

	}
}
}



foreach (sort keys %gene_all_trans){
my @all_trans=@{$gene_all_trans{$_}};
my @all_gene_exon=@{$gene_all_exon{$_}};
if(exists $gene_introns{$_}){
@all_gene_intron=@{$gene_introns{$_}};
}
	for (my $i=0;$i<@all_trans;$i++){
		my @all=map {0}(0..$num-2);
		my @trans_exon_sta=@{$exon_sta{$all_trans[$i]}};
		my @trans_exon_end=@{$exon_end{$all_trans[$i]}};
		my $exon_num=@trans_exon_sta;
		my $intron_num=@trans_exon_sta-1;
		if (exists $trans_intron{$all_trans[$i]}){
		my @trans_all_intron=@{$trans_intron{$all_trans[$i]}};
		for (my $n=$i+1;$n<@all_trans;$n++){
			if(exists $trans_intron{$all_trans[$n]}){
			my @trans_all_intron_new=@{$trans_intron{$all_trans[$n]}};
			if (@trans_all_intron_new~~@trans_all_intron && $trans_len{$all_trans[$i]}<$trans_len{$all_trans[$n]}){
			push @discard,$all_trans[$i];
		}elsif(@trans_all_intron_new~~@trans_all_intron && $trans_len{$all_trans[$n]}<$trans_len{$all_trans[$i]}){
			push @discard,$all_trans[$n];
		}else{
		next;
	}
}
}
}

	if ($exon_num<$options{n}){
		push @discard,$all_trans[$i];#exon_num==1
	}

		for (my $j=1;$j<@trans_exon_sta;$j++){
			my $trans_intron_sta=$trans_exon_end[$j-1]+1;
			my $trans_intron_end=$trans_exon_sta[$j]-1;
			my $intron_bed=join("-",$chr_trans{$all_trans[$i]},$trans_intron_sta,$trans_intron_end);
		for (my $k=1;$k<=$num-1;$k++){
			my $ming="junc_count_$k";
		if ((exists ${$ming}{$intron_bed}) && (${$ming}{$intron_bed}>=$options{m})){
			$all[$k-1]+=1;
			}else{
			$all[$k-1]+=0;
		}
}
}



if ($intron_num>0){
my $new=0;
for(my $c=0;$c<@all;$c+=$options{c}){
	my $test=0;
	for(my $d=$c;$d<($c+$options{c});$d++){
		if($all[$d] == $intron_num){
			$test+=1;
		}

	}
	
	if($test == $options{c}){
	$new+=1;
}
}

if ($new<$options{s}){
push @discard,$all_trans[$i];
}
}

}

for (my $p=0;$p<@all_gene_exon;$p++){
	my @exon_id=split(/;/,$all_gene_exon[$p]);
	my $exon_chr_part=POSIX::ceil((($exon_id[1]+$exon_id[2])/2)/1000000);
	for(my $r=0;$r<@all_gene_intron;$r++){
	my @intron_id=split(/;/,$all_gene_intron[$r]);
	if($exon_id[1]<$intron_id[1] && $exon_id[2]>$intron_id[2]){
		my $gdexon=0;
		my @testexon=map {0}(0..$num-2);
		for (my $e=1;$e<=$num-1;$e++){
				undef my @mayexon;
				my $hua="chr_$e";
				my $tian="exon_count_$e";
				my $men1=0; my $men2=0;
				my $new_ir_exon_id=join("_",$exon_id[0],$exon_chr_part);
				my $new_ir_exon_id_u=join("_",$exon_id[0],$exon_chr_part-1);
				my $new_ir_exon_id_d=join("_",$exon_id[0],$exon_chr_part+1);
				if(exists ${$hua}{$new_ir_exon_id}){
				@mayexon=@{${$hua}{$new_ir_exon_id}};
				}
				if(exists ${$hua}{$new_ir_exon_id_u}){
				push @mayexon,@{${$hua}{$new_ir_exon_id_u}};
				}
				if(exists ${$hua}{$new_ir_exon_id_d}){
				push @mayexon,@{${$hua}{$new_ir_exon_id_d}};
				}
				foreach $exonin (@mayexon){
				@irexon=split(/-/,$exonin);
				if($irexon[1]<$intron_id[1] && $irexon[2]>$intron_id[1]){
			      	   $men1+=${$tian}{$exonin};
			}
				if($irexon[1]<$intron_id[2] && $irexon[2]>$intron_id[2]){
				   $men2+=${$tian}{$exonin};
			}
		}
				if($men1>=$options{m} && $men2>=$options{m}){
				   $testexon[$e-1]+=1;
			}
		}
				for(my $z=0;$z<@testexon;$z+=$options{c}){
					my $testexonnum=0;
					for(my $y=$z;$y<($z+$options{c});$y++){
						if($testexon[$y]==1){
						$testexonnum+=1;
					}
				}
				if($testexonnum==$options{c}){
					$gdexon+=1;
			}
			}
				if($gdexon<$options{s}){
					push @discard,@{$exon_trans{$all_gene_exon[$p]}};
				}
		}
}
}

if(exists $gene_f_exon{$_} && exists $gene_l_exon{$_} && exists $gene_introns{$_}){
my @gene_fe=@{$gene_f_exon{$_}};
my @gene_le=@{$gene_l_exon{$_}};
for (my $t=0;$t<@gene_fe;$t++){
	my @exon_id_f=split(/;/,$gene_fe[$t]);
	my $exon_chr_part=POSIX::ceil((($exon_id_f[1]+$exon_id_f[2])/2)/1000000);
	for(my $v=0;$v<@all_gene_intron;$v++){
		my @intron_id=split(/;/,$all_gene_intron[$v]);
		if($exon_id_f[1]<$intron_id[2] && $exon_id_f[2]>$intron_id[2]){
		my $fexon=0;
		my @testfexon=map {0}(0..$num-2);
		for (my $b=1;$b<=$num-1;$b++){
			undef my @f_exon;
			my $hua="chr_$b";
			my $tian="exon_count_$b";
			my $men_f=0;
			my $new_f_exon_id=join("_",$exon_id_f[0],$exon_chr_part);
			my $new_f_exon_id_u=join("_",$exon_id_f[0],$exon_chr_part-1);
			my $new_f_exon_id_d=join("_",$exon_id_f[0],$exon_chr_part+1);
			if(exists ${$hua}{$new_f_exon_id}){
			@f_exon=@{${$hua}{$new_f_exon_id}};
			}
			if(exists ${$hua}{$new_f_exon_id_u}){
			push @f_exon,@{${$hua}{$new_f_exon_id_u}};
			}
			if(exists ${$hua}{$new_f_exon_id_d}){
			push @f_exon,@{${$hua}{$new_f_exon_id_d}};
			}
			foreach $exonin_f (@f_exon){
			@f_exon_f=split(/-/,$exonin_f);
			if($f_exon_f[1]<$intron_id[2] && $f_exon_f[2]>$intron_id[2]){
			   $men_f+=${$tian}{$exonin_f};
		}
	}
			if($men_f>=$options{m}){
				   $testfexon[$b-1]+=1;
			}
		}
                                for(my $z=0;$z<@testfexon;$z+=$options{c}){
                                        my $testfexonnum=0;
                                        for(my $y=$z;$y<($z+$options{c});$y++){
                                                if($testfexon[$y]==1){
                                                $testfexonnum+=1;
                                        }
                                }
                                if($testfexonnum==$options{c}){
                                        $fexon+=1;
                        }
                        }
			if($fexon<$options{s}){
				push @discard,@{$exon_trans{$gene_fe[$t]}};
		}
	}
}
}

for (my $t=0;$t<@gene_le;$t++){
	my @exon_id_l=split(/;/,$gene_le[$t]);
	my $exon_chr_part=POSIX::ceil((($exon_id_l[1]+$exon_id_l[2])/2)/1000000);
	for(my $v=0;$v<@all_gene_intron;$v++){
		my @intron_id=split(/;/,$all_gene_intron[$v]);
		if($exon_id_l[1]<$intron_id[1] && $exon_id_l[2]>$intron_id[1]){
		my $lexon=0;
		my @testlexon=map {0}(0..$num-2);
		for (my $b=1;$b<=$num-1;$b++){
			undef my @l_exon;
			my $hua="chr_$b";
			my $tian="exon_count_$b";
			my $men_l=0;
			my $new_l_exon_id=join("_",$exon_id_l[0],$exon_chr_part);
			my $new_l_exon_id_u=join("_",$exon_id_l[0],$exon_chr_part-1);
			my $new_l_exon_id_d=join("_",$exon_id_l[0],$exon_chr_part+1);
			if(exists ${$hua}{$new_l_exon_id}){
			@l_exon=@{${$hua}{$new_l_exon_id}};
			}
			if(exists ${$hua}{$new_l_exon_id_u}){
			push @l_exon,@{${$hua}{$new_l_exon_id_u}};
			}
			if(exists ${$hua}{$new_l_exon_id_d}){
			push @l_exon,@{${$hua}{$new_l_exon_id_d}};
			}
			foreach $exonin_l (@l_exon){
			@l_exon_l=split(/-/,$exonin_l);
			if($l_exon_l[1]<$intron_id[1] && $l_exon_l[2]>$intron_id[1]){
			   $men_l+=${$tian}{$exonin_l};
		}
	}
			if($men_l>=$options{m}){
				   $testlexon[$b-1]+=1;
			}
		}
                                for(my $z=0;$z<@testlexon;$z+=$options{c}){
                                        my $testlexonnum=0;
                                        for(my $y=$z;$y<($z+$options{c});$y++){
                                                if($testlexon[$y]==1){
                                                $testlexonnum+=1;
                                        }
                                }
                                if($testlexonnum==$options{c}){
                                        $lexon+=1;
                        }
                        }

			if($lexon<$options{s}){
				push @discard,@{$exon_trans{$gene_le[$t]}};
		}
	}
}
}
}


}



my @discard_uniq = uniq @discard;
	foreach my $discard_uniq (@discard_uniq){
		   $dis{$discard_uniq}=1;
	}

open CE,"$options{i}"||die;#merge_gtf
while (<CE>){
chomp;
if($_!~/^#/){
	my @filter=split(/"/,$_);
		if (!exists $dis{$filter[3]}){
			print OUT "$_\n";
}	
}
}
close CE;
