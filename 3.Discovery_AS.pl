#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
my $usage =
"$0 <-i in.gtf> [options]

This script can discovery alternative splice from GTF file generate by cufflinks and so on

-i gtf    gtf file, must have gene_id and transcript_id
-o        output dir						 [default: ./DSAS]
-m int    bias of base length in ME/F/L/ALT3/ALT5		 [default: 10]
-l 	  name prefix for annotation AS				 [default: gtf file name]
";
my %options=();
getopts("i:o:m:l:",\%options);

if (!exists $options{i}){
die "$usage";
}

my $dir = $ENV{'PWD'};

if (!exists $options{o}){
	$options{o}="./DSAS/";
	}
if (-d "$dir/$options{o}"){
	}else{
	`mkdir $dir/$options{o}`;
}

if (!exists $options{m}){
	$options{m}=10;
}

if (!exists $options{l}){
	$options{l}=(split(/\//,$options{i}))[-1];
}

if ($options{o}!~/\/$/){
	$options{o}=$options{o}."/";
}

open GTF,"$dir/$options{i}"||die "can not open $options{i}\n";;
open EXON,">$dir/$options{o}./$options{l}.exon"||die;
open INTRON,">$dir/$options{o}./$options{l}.intron"||die;
open EIE,">$dir/$options{o}./$options{l}.exon_intron_exon"||die;
open IEI,">$dir/$options{o}./$options{l}.intron_exon_intron"||die;
open ES,">$dir/$options{o}./$options{l}.exon_skip"||die;
open IR,">$dir/$options{o}./$options{l}.intron_retention"||die;
open ALT5,">$dir/$options{o}./$options{l}.alt5"||die;
open ALT3,">$dir/$options{o}./$options{l}.alt3"||die;
open F,">$dir/$options{o}./$options{l}.first.exon"||die;
open L,">$dir/$options{o}./$options{l}.last.exon"||die;
open ME,">$dir/$options{o}./$options{l}.mutex_exon"||die;


print EXON "gene\ttranscript\tChr-exonstart-exonend\tstrand\n";
print INTRON "gene\ttranscript\tChr-intronstart-intronend\tstrand\n";
print EIE "gene\ttranscript\tChr-exon1start-exon1end-exon2start-exon2end\tstrand\n";
print IEI "gene\ttranscript\tChr-intron1start-exonstart-exonend-intron2end\tstrand\n";
print ES "gene\tskiped_exon\tskiped_exon_iei\thaveexon_trans\tintron_trans\tintron_trans_eie\tnothing\tstrand\n";
print ME "gene\tforward_exon\tforward_iei\tforward_trans\tback_trans\tback_exon\tback_iei\tstrand\n";
print F "gene\tback_e_e_all\tback_trans\tforward_all_trans\tforward_e_e_all\ttype\tnothing\tstrand\n";
print L "gene\tback_e_e_all\tback_trans\tforwarf_all_trans\tforward_e_e_all\ttype\tnothing\tstrand\n";
print IR "gene\tintron\teie\tisintron_trans\tisexon_trans\texon\tnothing\tstrand\n";
print ALT3 "gene\tshort_trans_intron\tshort_trans_eie\tshort_trans\tlong_trans\tlong_trans_intron\tlong_trans_eie\tstrand\n";
print ALT5 "gene\tshort_trans_intron\tshort_trans_eie\tshort_trans\tlong_transc\tlong_trans_intron\tlong_trans_eie\tstrand\n";
my(%chr,%trans_gene,%strand,%gene_trans,%gene_exon,%exon_trans,%trans_exon);
while (<GTF>) {
chomp;
next if ($_ eq "");
my @line=split(/\s+/,$_);
if ($line[2] eq "exon"){
my @tem=split(/"/,$_);
my $gene_id=(split(/"/,(split(/gene_id "/,$_))[1]))[0];
my $trans_id=(split(/"/,(split(/transcript_id "/,$_))[1]))[0];
my $exon=join("-",$line[0],$line[3],$line[4]);
$chr{$gene_id}=$line[0];
$trans_gene{$trans_id}=$gene_id;
$strand{$gene_id}=$line[6];
if (grep {$_ eq $trans_id} @{$gene_trans{$gene_id}}){
	}else{
		push @{$gene_trans{$gene_id}},$trans_id;
}
if (grep {$_ eq $exon} @{$gene_exon{$gene_id}}){
	}else{
		push @{$gene_exon{$gene_id}},$exon;
}
if (grep {$_ eq $trans_id} @{$exon_trans{$exon}}){
	}else{
		push @{$exon_trans{$exon}},$trans_id;
}
push @{$trans_exon{$trans_id}},$exon;
print EXON "$gene_id\t$trans_id\t$exon\t$line[6]\n";
}
}
close GTF;
close EXON;


my (%eie_all_id_gene,%gene_intron,%gene_eie_all,%intron_trans,%eie_all_trans,%iei_all_trans,%gene_iei_all,%iei_all_gene);
foreach my $key (sort keys %gene_trans){
	my @trans=@{$gene_trans{$key}};
	for (my $i=0;$i<@trans;$i++){	     
	   my @exon_id=@{$trans_exon{$trans[$i]}};
	   my $count_intron=@exon_id-1;
           for (my $j=0;$j<$count_intron;$j++){
	     my @exon_location=split(/-/,$exon_id[$j]);
	     my @next_exon_location=split(/-/,$exon_id[$j+1]);
	     my $intron_start=$exon_location[2]+1;
	     my $intron_end=$next_exon_location[1]-1;
	     my $intron_id=join("-",$chr{$key},$intron_start,$intron_end);
	     my $eie_all_id=join("-",$chr{$key},$exon_location[1],$exon_location[2],$next_exon_location[1],$next_exon_location[2]);
		$eie_all_id_gene{$eie_all_id}=$key;
	     if  (grep {$_ eq $intron_id} @{$gene_intron{$key}}){
}else{
		push @{$gene_intron{$key}},$intron_id;
}
	     if  (grep {$_ eq $eie_all_id} @{$gene_eie_all{$key}}){
}else{
		push @{$gene_eie_all{$key}},$eie_all_id;
}
	     if  (grep {$_ eq $trans[$i]} @{$intron_trans{$intron_id}}){
}else{
		push @{$intron_trans{$intron_id}},$trans[$i];
}
		push @{$eie_all_trans{$eie_all_id}},$trans[$i];
		print INTRON "$key\t$trans[$i]\t$intron_id\t$strand{$key}\n";
		print EIE "$key\t$trans[$i]\t$eie_all_id\t$strand{$key}\n";
		}
		if ($count_intron>1){
			for (my $k=0;$k<$count_intron-1;$k++){
			    my @exon_location=split(/-/,$exon_id[$k]);
			    my @next_exon_location=split(/-/,$exon_id[$k+1]);
			    my @third_exon_location=split(/-/,$exon_id[$k+2]);
			    my $iei_start=$exon_location[2]+1;
			    my $iei_end=$third_exon_location[1]-1;
			    my $iei_all_id=join("-",$chr{$key},$iei_start,$next_exon_location[1],$next_exon_location[2],$iei_end);
				$iei_all_gene{$iei_all_id}=$key;
				push @{$iei_all_trans{$iei_all_id}},$trans[$i];
				if (grep {$_ eq $iei_all_id} @{$gene_iei_all{$key}}){
}else{
				push @{$gene_iei_all{$key}},$iei_all_id;
}
				print IEI "$key\t$trans[$i]\t$iei_all_id\t$strand{$key}\n"
				}
			}
		}
}
close EIE;
close INTRON;
close IEI;


sub exon_skip {
	my ($key_iei,$exon_length,$exon_start,$exon_end,$exon_id,$iei_id,$trans_iei,$gene,$intron_start,$intron_end,$eie_start,$eie_end,$intron_trans,$compare_intron_s,$compare_intron_e);
	my (@exon,@gene_eie,@intron_location);
	foreach $key_iei (sort keys %iei_all_gene){
		@exon=split(/-/,$key_iei);
		$exon_length=$exon[3]-$exon[2]+1;
		$exon_start=$exon[2];
		$exon_end=$exon[3];
		$exon_id=join("-",$exon[0],$exon_start,$exon_end);
		$iei_id=join("-",$exon[0],$exon[1],$exon[4]);
		$trans_iei=join(",",@{$iei_all_trans{$key_iei}});
		if ($exon_length>=20){
			$gene=$iei_all_gene{$key_iei};
			@gene_eie=@{$gene_eie_all{$gene}};
			for (my $i=0;$i<@gene_eie;$i++){
				@intron_location=split(/-/,$gene_eie[$i]);
				$intron_start=$intron_location[2]+1;
				$intron_end=$intron_location[3]-1;
				$eie_start=$intron_location[1];
				$eie_end=$intron_location[4];
				$intron_trans=join(",",@{$eie_all_trans{$gene_eie[$i]}});
				if (($intron_start<$exon_start) && ($intron_end>$exon_end)){
					$compare_intron_s=$intron_start-$exon[1];
					$compare_intron_e=$intron_end-$exon[4];
						if(($compare_intron_s==0 && $eie_end>$exon[4]) || ($compare_intron_e==0 && $eie_start<$exon[1])){
						print ES "$gene\t$exon_id\t$iei_id\t$trans_iei\t$intron_trans\t$gene_eie[$i]\t.\t$strand{$gene}\n";
}
}
				}
			}
		}
	}
&exon_skip;
close ES;

sub mutex_exon {
my ($key_iei,$exon_length,$exon_id,$iei_id,$trans_iei,$gene,$back_iei_id,$back_exon_id,$compare_intron_s,$compare_intron_e,$back_trans);
my (@exon,@gene_all_iei,@iei);
	foreach $key_iei (sort keys %iei_all_gene){
		@exon=split(/-/,$key_iei);
		$exon_length=$exon[3]-$exon[2]+1;
                $exon_id=join("-",$exon[0],$exon[2],$exon[3]);
                $iei_id=join("-",$exon[0],$exon[1],$exon[4]);
                $trans_iei=join(",",@{$iei_all_trans{$key_iei}});
                if ($exon_length>=20){
                        $gene=$iei_all_gene{$key_iei};
                        @gene_all_iei=@{$gene_iei_all{$gene}};
			for(my $i=0;$i<@gene_all_iei;$i++){
			@iei=split(/-/,$gene_all_iei[$i]);
			$back_iei_id=join("-",$iei[0],$iei[1],$iei[4]);
			$back_exon_id=join("-",$iei[0],$iei[2],$iei[3]);
			$compare_intron_s=abs($iei[1]-$exon[1]);
			$compare_intron_e=abs($iei[4]-$exon[4]);
			if ($iei[1]<$exon[2] && $iei[2]>$exon[3] && $iei[3]<=$exon[4] && $compare_intron_s<=$options{m} && $compare_intron_e<=$options{m}){
			$back_trans=join(";",@{$iei_all_trans{$gene_all_iei[$i]}});	
				if ($strand{$gene} eq "+"){
				print ME "$gene\t$exon_id\t$iei_id\t$trans_iei\t$back_trans\t$back_exon_id\t$back_iei_id\t$strand{$gene}\n";
				}else{
				print ME "$gene\t$back_exon_id\t$back_iei_id\t$back_trans\t$trans_iei\t$exon_id\t$iei_id\t$strand{$gene}\n";
}
}
}
}
}
}
&mutex_exon;
close ME;



sub first_or_last_exon{
	my ($trans,$gene,$strand,$trans_all_exon_count,$right,$left,$all_trans_exon,$all_trans_exon_count,$one_eie_all,$other_eie_all,$compare_first_exon_s,$compare_first_exon_e,$last_one_eie_all,$other_last_eie_all,$compare_last_exon_e,$compare_last_exon_s);
	my (@trans_all_exon,@all_trans,@second_exon,@other_second_exon,@first_exon_loc,@other_first_exon_loc,@last_second_exon,@other_last_second_exon,@last_exon_loc,@other_last_exon_loc);
 foreach $trans(sort keys %trans_exon){
                $gene=$trans_gene{$trans};
                $strand=$strand{$gene};
                @trans_all_exon=@{$trans_exon{$trans}};
		$trans_all_exon_count=@{$trans_exon{$trans}};
                        @all_trans=@{$gene_trans{$gene}};
			if($trans_all_exon_count>1){
                        for (my $j=0;$j<@all_trans;$j++){
			$right=$left=0;
                        $all_trans_exon=join(",",@{$trans_exon{$all_trans[$j]}});
			$all_trans_exon_count=@{$trans_exon{$all_trans[$j]}};
			if($all_trans_exon_count>1){
                        for (my $i=0;$i<@trans_all_exon;$i++){
                if ($i==0 && $all_trans_exon=~/$trans_all_exon[$i]/){
                        $right+=1;
                }
                if ($i>0 && $i<@trans_all_exon && $all_trans_exon!~/$trans_all_exon[$i]/){
                        $right+=1
                }
		if ($i==(@trans_all_exon-1) && $all_trans_exon=~/$trans_all_exon[$i]/){
			$left+=1;
		}
		if ($i<@trans_all_exon-1 && $i>=0 && $all_trans_exon!~/$trans_all_exon[$i]/){
			$left+=1;
        }
}
		@second_exon=split(/-/,$trans_all_exon[1]);
                $one_eie_all=join("-",$trans_all_exon[0],$second_exon[1],$second_exon[2]);
                @other_second_exon=split(/-/,$trans_exon{$all_trans[$j]}[1]);
                $other_eie_all=join("-",$trans_exon{$all_trans[$j]}[0],$other_second_exon[1],$other_second_exon[2]);
                @first_exon_loc=split(/-/,$trans_all_exon[0]);
                @other_first_exon_loc=split(/-/,$trans_exon{$all_trans[$j]}[0]);
                $compare_first_exon_s=$first_exon_loc[1]-$other_first_exon_loc[1];
		$compare_first_exon_e=$first_exon_loc[2]-$other_first_exon_loc[2];
		@last_second_exon=split(/-/,$trans_all_exon[-1]);
                $last_one_eie_all=join("-",$trans_all_exon[-2],$last_second_exon[1],$last_second_exon[2]);
                @other_last_second_exon=split(/-/,$trans_exon{$all_trans[$j]}[-1]);
                $other_last_eie_all=join("-",$trans_exon{$all_trans[$j]}[-2],$other_last_second_exon[1],$other_last_second_exon[2]);
                @last_exon_loc=split(/-/,$trans_all_exon[-1]);
                @other_last_exon_loc=split(/-/,$trans_exon{$all_trans[$j]}[-1]);
                $compare_last_exon_e=$last_exon_loc[2]-$other_last_exon_loc[2];
		$compare_last_exon_s=$last_exon_loc[1]-$other_last_exon_loc[1];

		
                if ($compare_first_exon_s>$options{m} && $compare_first_exon_e!=0 && $right==0 && $strand eq "+" && $all_trans_exon_count==$trans_all_exon_count){
	print F "$gene\t$one_eie_all\t$trans\t$all_trans[$j]\t$other_eie_all\tF\t.\t$strand\n";
}elsif($compare_first_exon_s>$options{m} && $compare_first_exon_e!=0 && $right==0 && $strand eq "-" &&  $all_trans_exon_count==$trans_all_exon_count){
	print L "$gene\t$one_eie_all\t$trans\t$all_trans[$j]\t$other_eie_all\tL\t.\t$strand\n";
}else{
}
		if ($compare_last_exon_e<-$options{m} && $compare_last_exon_s!=0 && $left==0 && $strand eq "+" && $all_trans_exon_count==$trans_all_exon_count){
print L "$gene\t$last_one_eie_all\t$trans\t$all_trans[$j]\t$other_last_eie_all\tL\t.\t$strand\n";
}elsif($compare_last_exon_e<-$options{m} && $compare_last_exon_s!=0 && $left==0 && $strand eq "-" &&  $all_trans_exon_count==$trans_all_exon_count){
print F "$gene\t$last_one_eie_all\t$trans\t$all_trans[$j]\t$other_last_eie_all\tF\t.\t$strand\n";
}else{
}
}
}
}
}
}
&first_or_last_exon;
close F;
close L;




sub intron_retention {
	my ($key_eie,$eie_start,$eie_end,$gene,$intron_start,$intron_end,$intron_length,$eie_trans,$intron_id,$eie_id,$exon_start,$exon_end,$compare_exon_s,$compare_exon_e,$exon_trans);
	my (@eie_location,@all_gene_exon);
	foreach $key_eie (sort keys %eie_all_id_gene){
		@eie_location=split(/-/,$key_eie);
		$eie_start=$eie_location[1];
		$eie_end=$eie_location[4];
		$gene=$eie_all_id_gene{$key_eie};
		$intron_start=$eie_location[2]+1;
		$intron_end=$eie_location[3]-1;		
		$intron_length=$intron_end-$intron_start;
		$eie_trans=join(",",@{$eie_all_trans{$key_eie}});
		$intron_id=join("-",$eie_location[0],$intron_start,$intron_end);
		$eie_id=join("-",$eie_location[0],$eie_location[1],$eie_location[4]);
		if ($intron_length >= 50 ){
			@all_gene_exon=@{$gene_exon{$gene}};
			for (my $i=0;$i<@all_gene_exon;$i++){
			$exon_start=(split(/-/,$all_gene_exon[$i]))[1];
			$exon_end=(split(/-/,$all_gene_exon[$i]))[2];
				if($exon_start-$intron_start<=-20 && $exon_end-$intron_end>=20){
					$compare_exon_s=$exon_start-$eie_start;
					$compare_exon_e=$exon_end-$eie_end;
					$exon_trans=join(",",@{$exon_trans{$all_gene_exon[$i]}});
					if (($compare_exon_s<10 && $compare_exon_e==0) || ($compare_exon_s==0 && $compare_exon_e>-10)){
					print IR "$gene\t$intron_id\t$eie_id\t$eie_trans\t$exon_trans\t$all_gene_exon[$i]\t.\t$strand{$gene}\n";
}	
					}
				}
			}
		}
	}			
&intron_retention;
close IR;


sub alt{
	my ($eie_all_id,$strand,$eie_start,$eie_end,$intron_start,$intron_end,$eie_intron,$eie_location,$gene,$trans_eie_location_start,$trans_eie_location_end,$trans_intron_start,$trans_intron_end,$compare_eie_start,$compare_eie_end,$trans_intron,$trans_eie,$compare_intron_start,$compare_intron_end,$merge_1_trans,$merge_2_trans);
	my (@eie_location,@eie,@trans_eie_location);
	foreach $eie_all_id (sort keys %eie_all_id_gene){
		$strand=$strand{$eie_all_id_gene{$eie_all_id}};
		@eie_location=split(/-/,$eie_all_id);
		$eie_start=$eie_location[1];
		$eie_end=$eie_location[4];
		$intron_start=$eie_location[2]+1;
		$intron_end=$eie_location[3]-1;
		$eie_intron=join("-",$eie_location[0],$intron_start,$intron_end);
		$eie_location=join("-",$eie_location[0],$eie_start,$eie_end);
		$gene=$eie_all_id_gene{$eie_all_id};
		@eie=@{$gene_eie_all{$gene}};
				for (my $i=0;$i<@eie;$i++){
					@trans_eie_location=split(/-/,$eie[$i]);
					$trans_eie_location_start=$trans_eie_location[1];
					$trans_eie_location_end=$trans_eie_location[4];
					$trans_intron_start=$trans_eie_location[2]+1;
					$trans_intron_end=$trans_eie_location[3]-1;
					$compare_eie_start=$trans_eie_location_start-$eie_start;
					$compare_eie_end=$trans_eie_location_end-$eie_end;
					$trans_intron=join("-",$trans_eie_location[0],$trans_intron_start,$trans_intron_end);
					$trans_eie=join("-",$trans_eie_location[0],$trans_eie_location_start,$trans_eie_location_end);
					if (($compare_eie_start==0 && ($compare_eie_end>-$options{m} && $compare_eie_end<$options{m})) || ($compare_eie_end==0 && ($compare_eie_start>-$options{m} && $compare_eie_start<$options{m}))){
						
						$compare_intron_start=$trans_intron_start-$intron_start;
						$compare_intron_end=$trans_intron_end-$intron_end;
						$merge_1_trans=join(",",@{$eie_all_trans{$eie[$i]}});
						$merge_2_trans=join(",",@{$eie_all_trans{$eie_all_id}});
					if ($strand eq "+"){
					 	if ($compare_intron_start==0 && $compare_intron_end>=$options{m}){
						print ALT3 "$gene\t$trans_intron\t$trans_eie\t$merge_1_trans\t$merge_2_trans\t$eie_intron\t$eie_location\t$strand\n";
		}
						if ($compare_intron_end==0 && $compare_intron_start>=$options{m}){		
						print ALT5 "$gene\t$eie_intron\t$eie_location\t$merge_2_trans\t$merge_1_trans\t$trans_intron\t$trans_eie\t$strand\n";
		}
	}
					if ($strand eq "-"){
						if ($compare_intron_start==0 && $compare_intron_end>=$options{m}){
						print ALT5 "$gene\t$trans_intron\t$trans_eie\t$merge_1_trans\t$merge_2_trans\t$eie_intron\t$eie_location\t$strand\n";
		}
						if ($compare_intron_end==0 && $compare_intron_start>=$options{m}){
						print ALT3 "$gene\t$eie_intron\t$eie_location\t$merge_2_trans\t$merge_1_trans\t$trans_intron\t$trans_eie\t$strand\n";
		}
	}
}
}
}
}
&alt;
close ALT5;
close ALT3;
