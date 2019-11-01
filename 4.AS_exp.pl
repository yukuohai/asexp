#! /usr/bin/perl -w
use warnings;
#use strict;
use Getopt::Std;
use POSIX;
my $usage =
"$0 <-i Aligned_bam.loc.fsim.list -g merged.gtf -l readslength> [options]

This script can calculate alternative splicing event exp form junction reads file

-i        loc.fsim.txt include the absolute path of the loc.fsim file form the bam files which produce gtf, and the coefficient for Standard,if not [default: 1]
-g 	  gtf file and must have the AS file in the same path
-l        max reads length
-o        output dir                             [default: ./ASexp]

";
my %options=();
getopts("i:g:l:o:",\%options);

if(!exists $options{i} || !exists $options{g} || !exists $options{l}){
die "$usage";
}

if (!exists $options{o}){
	$options{o}="./ASexp";
}
my $dir = $ENV{'PWD'};

if (-d "$dir/$options{o}"){
        }else{
        `mkdir $dir/$options{o}`;
}

if ($options{o}!~/\/$/){
        $options{o}=$options{o}."/";
}

#chr1A_part1     100003442       100003560       1       E
my @all_id;
open IN,"$options{i}"||die;
while (my $file=<IN>){
chomp $file;
my @fsim=split(/\s+/,$file);
my $name=(split(/\//,$fsim[0]))[-1];
my $id=(split(/\./,$name))[0];
push @all_id,$id;
if(!exists $fsim[1]){
	$fsim[1]=1;
}
	open LOC,"$fsim[0]"||die;
	while (my $loc=<LOC>){
		chomp $loc;
		my @line=split(/\s+/,$loc);
		my $readsid=join("-",$line[0],$line[1],$line[2]);
		my $sample_exon_db="db_exon_all_$id";
		my $sample_intron_db="db_intron_all_$id";
		my $sample_exon_co="exon_id_$id";
		my $sample_intron_co="intron_id_$id";
		if ($line[4] eq "E"){
			my $chr_part=POSIX::ceil((($line[1]+$line[2])/2)/1000000);
			my $new_exon_id=join("_",$line[0],$chr_part);
			push @{${$sample_exon_db}{$line[0]}},$readsid;
			push @{${$sample_exon_db}{$new_exon_id}},$readsid;
				${$sample_exon_co}{$readsid}=$line[3]*$fsim[1];
			}
		if ($line[4] eq "I"){
			push @{${$sample_intron_db}{$line[0]}},$readsid;
				${$sample_intron_co}{$readsid}=$line[3]*$fsim[1];
			}
}
close LOC;
}
close IN;


sub alt3exp {
foreach my $sid(@all_id){
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open ALT3,"<$options{g}.alt3"||die;
readline ALT3;
open EXP,">$options{o}./$sid.exp"||die;
print EXP "#Alt3_gene\tshort_trans_intron\tshort_trans_eie\tshort_trans\tlong_trans\tlong_trans_intron\tlong_trans_eie\tstrand\tshort_count\tlong_count\tshort_rate\tlong_rate\n";
while (<ALT3>){
my ($db_alt3_short_count,$db_alt3_long_count,$db_alt3_short_rate,$db_alt3_long_rate);
chomp;
#DGENE.281       chr1D_part1-10986605-10986778   chr1D_part1-10986423-10987449   DGENE.281.1    DGENE.281.2      chr1D_part1-10986605-10986756   chr1D_part1-10986423-10987449   +
my @tem=split(/\s+/,$_);
$db_alt3_short_count=$db_alt3_long_count=0;
my $chr=(split(/-/,$tem[1]))[0];
foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
	if ($db_intron_all eq $tem[1]){
	$db_alt3_short_count+=${$sample_intron_co}{$db_intron_all};
	}
	if ($db_intron_all eq $tem[5]){
	$db_alt3_long_count+=${$sample_intron_co}{$db_intron_all};
	}
}

if ($db_alt3_short_count+$db_alt3_long_count>0){
$db_alt3_short_rate=$db_alt3_short_count/($db_alt3_short_count+$db_alt3_long_count);
$db_alt3_long_rate=$db_alt3_long_count/($db_alt3_short_count+$db_alt3_long_count);
}else{
$db_alt3_short_rate="NA";
$db_alt3_long_rate="NA";
}

print EXP "$_\t$db_alt3_short_count\t$db_alt3_long_count\t$db_alt3_short_rate\t$db_alt3_long_rate\n";
}
close ALT3;
}
}
&alt3exp;
close EXP;

sub alt5exp{
foreach my $sid(@all_id){
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open ALT5,"<$options{g}.alt5"||die;
readline ALT5;
open EXP,">>$options{o}./$sid.exp";
print EXP "#Alt5_gene\tshort_trans_intron\tshort_trans_eie\tshort_trans\tlong_trans\tlong_trans_intron\tlong_trans_eie\tstrand\tshort_count\tlong_count\tshort_rate\tlong_rate\n";
while (<ALT5>){
my ($db_alt5_short_count,$db_alt5_long_count,$db_alt5_short_rate,$db_alt5_long_rate);
chomp;
my @xiao=split(/\s+/,$_);
$db_alt5_short_count=$db_alt5_long_count=0;
my $chr=(split(/-/,$xiao[1]))[0];
foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
        if ($db_intron_all eq $xiao[1]){
        $db_alt5_short_count+=${$sample_intron_co}{$db_intron_all};
        }
        if ($db_intron_all eq $xiao[5]){
        $db_alt5_long_count+=${$sample_intron_co}{$db_intron_all};
        }
}

if ($db_alt5_short_count+$db_alt5_long_count>0){
$db_alt5_short_rate=$db_alt5_short_count/($db_alt5_short_count+$db_alt5_long_count);
$db_alt5_long_rate=$db_alt5_long_count/($db_alt5_short_count+$db_alt5_long_count);
}else{
$db_alt5_short_rate="NA";
$db_alt5_long_rate="NA";
}

print EXP "$_\t$db_alt5_short_count\t$db_alt5_long_count\t$db_alt5_short_rate\t$db_alt5_long_rate\n";
}
close ALT5;
}
}
&alt5exp;
close EXP;


sub esexp{
foreach my $sid(@all_id){
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open ES,"<$options{g}.exon_skip"||die;
readline ES;
open EXP,">>$options{o}./$sid.exp";
#DGENE.1355      chr1D_part1-104596742-104596809 chr1D_part1-104595949-104596929 DGENE.1355.4   DGENE.1355.1     chr1D_part1-104595874-104595948-104596930-104597022     .       +
print EXP "#ES_gene\tskiped_exon\tskiped_exon_iei\thave_exon_trans\tintron_trans\tintron_trans_eie\tnothing\tstrand\thaveexon_count\tnoexon_count\thaveexon_rate\tnoexon_rate\n";
while (<ES>){
my ($db_es_first_intron_count,$db_es_second_intron_count,$db_es_big_intron_count,$db_all,$db_have,$db_have_rate,$db_no_rate);
chomp;
my @line=split(/\s+/,$_);
my @len=split(/-/,$line[1]);
my $chr=$len[0];
my $first_intron_start=(split(/-/,$line[2]))[1];
my $first_intron_end=(split(/-/,$line[1]))[1]-1;
my $second_intron_start=(split(/-/,$line[1]))[2]+1;
my $second_intron_end=(split(/-/,$line[2]))[2];
my $first_intron_id=join("-",$chr,$first_intron_start,$first_intron_end);
my $second_intron_id=join("-",$chr,$second_intron_start,$second_intron_end);
my $big_intron_start=(split(/-/,$line[5]))[2]+1;
my $big_intron_end=(split(/-/,$line[5]))[3]-1;
my $big_intron_id=join("-",$chr,$big_intron_start,$big_intron_end);
$db_es_first_intron_count=$db_es_second_intron_count=$db_es_big_intron_count=0;
	foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
        if ($db_intron_all eq $first_intron_id){
        $db_es_first_intron_count+=${$sample_intron_co}{$db_intron_all};
        }
        if ($db_intron_all eq $second_intron_id){
        $db_es_second_intron_count+=${$sample_intron_co}{$db_intron_all};
        }
	if ($db_intron_all eq $big_intron_id){
	$db_es_big_intron_count+=${$sample_intron_co}{$db_intron_all};
	}
}	

$db_all=$db_es_first_intron_count+$db_es_second_intron_count;
$db_have=$db_all/2;
if ($db_have+$db_es_big_intron_count>0){
	$db_have_rate=$db_have/($db_have+$db_es_big_intron_count);
	$db_no_rate=$db_es_big_intron_count/($db_have+$db_es_big_intron_count);
}else{
	$db_have_rate="NA";
	$db_no_rate="NA";
}
print EXP "$_\t$db_have\t$db_es_big_intron_count\t$db_have_rate\t$db_no_rate\n";
}
close ES;
}
}
&esexp;
close EXP;


sub fexp{
#DGENE.10415     chr2D_part2-8018612-8019345-8020630-8021453     DGENE.10415.5   DGENE.10415.1  chr2D_part2-8018588-8019268-8020630-8021453      F       .       +
foreach my $sid(@all_id){
open F,"<$options{g}.first.exon"||die;
readline F;
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open EXP,">>$options{o}./$sid.exp";
print EXP "#FE_gene\tback_e_e_all\tback_trans\tforward_trans\tforward_e_e_all\ttype\tnothing\tstrand\tback_count\tforward_count\tback_rate\tforward_rate\n";
while (<F>){
my ($db_back_trans_count,$db_forward_trans_count,$db_back_rate,$db_forward_rate);
chomp;
my @tem=split(/\s+/,$_);
my $chr=(split(/-/,$tem[1]))[0];
my $intron1_start=(split(/-/,$tem[1]))[2]+1;
my $intron1_end=(split(/-/,$tem[1]))[3]-1;
my $intron2_start=(split(/-/,$tem[4]))[2]+1;
my $intron2_end=(split(/-/,$tem[4]))[3]-1;
my $intron1_id=join("-",$chr,$intron1_start,$intron1_end);
my $intron2_id=join("-",$chr,$intron2_start,$intron2_end);
$db_back_trans_count=$db_forward_trans_count=0;
	foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
	if ($db_intron_all eq $intron1_id){
	$db_back_trans_count+=${$sample_intron_co}{$db_intron_all};
	}
	if ($db_intron_all eq $intron2_id){
	$db_forward_trans_count+=${$sample_intron_co}{$db_intron_all};
	}
}
if ($db_back_trans_count+$db_forward_trans_count>0){
$db_back_rate=$db_back_trans_count/($db_back_trans_count+$db_forward_trans_count);
$db_forward_rate=$db_forward_trans_count/($db_back_trans_count+$db_forward_trans_count);
}else{
$db_back_rate="NA";
$db_forward_rate="NA";
}


print EXP "$_\t$db_back_trans_count\t$db_forward_trans_count\t$db_back_rate\t$db_forward_rate\n";
}
close F;
}
}
&fexp;
close EXP;

sub lexp{
foreach my $sid(@all_id){
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open L,"<$options{g}.last.exon"||die;
readline L;
open EXP,">>$options{o}./$sid.exp";
print EXP "#LE_gene\tback_e_e_all\tback_trans\tforward_trans\tforward_e_e_all\ttype\tnothing\tstrand\tback_count\tforward_count\tback_rate\tforward_rate\n";
while (<L>){
my($db_back_trans_count,$db_forward_trans_count,$db_back_rate,$db_forward_rate);
chomp;
my @tem=split(/\s+/,$_);
my $chr=(split(/-/,$tem[1]))[0];
my $intron1_start=(split(/-/,$tem[1]))[2]+1;
my $intron1_end=(split(/-/,$tem[1]))[3]-1;
my $intron2_start=(split(/-/,$tem[4]))[2]+1;
my $intron2_end=(split(/-/,$tem[4]))[3]-1;
my $intron1_id=join("-",$chr,$intron1_start,$intron1_end);
my $intron2_id=join("-",$chr,$intron2_start,$intron2_end);
$db_back_trans_count=$db_forward_trans_count=0;
        foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
        if ($db_intron_all eq $intron1_id){
        $db_back_trans_count+=${$sample_intron_co}{$db_intron_all};
        }
        if ($db_intron_all eq $intron2_id){
        $db_forward_trans_count+=${$sample_intron_co}{$db_intron_all};
        }
}
if ($db_back_trans_count+$db_forward_trans_count>0){
$db_back_rate=$db_back_trans_count/($db_back_trans_count+$db_forward_trans_count);
$db_forward_rate=$db_forward_trans_count/($db_back_trans_count+$db_forward_trans_count);
}else{
$db_back_rate="NA";
$db_forward_rate="NA";
}

print EXP "$_\t$db_back_trans_count\t$db_forward_trans_count\t$db_back_rate\t$db_forward_rate\n";
}
close L;
}
}
&lexp;
close EXP;


sub irexp {
foreach my $sid(@all_id){
my $sample_exon_db="db_exon_all_$sid";
my $sample_exon_co="exon_id_$sid";
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open IR,"<$options{g}.intron_retention"||die;
readline IR;
open EXP,">>$options{o}./$sid.exp";
#DGENE.1317      chr1D_part1-100118983-100119084 chr1D_part1-100118834-100119214 TraesCS1D01G107200.1    DGENE.1317.2    chr1D_part1-100118834-100119214 .       +
print EXP "#IR_gene\tintron\tintron_eie\tisintron_trans\tisexon_trans\texon\tnothing\tstrand\tisintron_count\tisexon_count\tisintron_rate\tisexon_rate\n";

while (<IR>){
my($length_intron,$db_intron_count,$db_exon_count,$db_exon_count_norm,$db_intron_rate,$db_exon_rate);
chomp;
undef my @chr_part_exon;
my @ir=split(/\s+/,$_);
my $chr=(split(/-/,$ir[1]))[0];
my $intron_start=(split(/-/,$ir[1]))[1];
my $intron_end=(split(/-/,$ir[1]))[2];
my $len=$intron_end-$intron_start+1;
my $chr_part=POSIX::ceil((($intron_start+$intron_end)/2)/1000000);
my $new_exon_id=join("_",$chr,$chr_part);
my $new_exon_id_u=join("_",$chr,$chr_part-1);
my $new_exon_id_d=join("_",$chr,$chr_part+1);
if(exists ${$sample_exon_db}{$new_exon_id}){
push @chr_part_exon,@{${$sample_exon_db}{$new_exon_id}};
}

if(exists ${$sample_exon_db}{$new_exon_id_u}){
push @chr_part_exon,@{${$sample_exon_db}{$new_exon_id_u}};
}

if(exists ${$sample_exon_db}{$new_exon_id_d}){
push @chr_part_exon,@{${$sample_exon_db}{$new_exon_id_d}};
}

if ($len>=$options{l}*2){
$length_intron=$options{l}*2;
}else{
$length_intron=$len;
}
$db_intron_count=$db_exon_count=0;
	foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
	if ($db_intron_all eq $ir[1]){
	$db_intron_count+=${$sample_intron_co}{$db_intron_all};
	}
}

	foreach my $db_exon_all (@chr_part_exon){
	my @db_exon=split(/-/,$db_exon_all);
		if (($db_exon[1]<$intron_start && $db_exon[2]>$intron_start)||($db_exon[1]<$intron_end && $db_exon[2]>$intron_end)){
			my $test_db_id1=join("-",$chr,$db_exon[2]+1,$intron_end);
			my $test_db_id2=join("-",$chr,$intron_start,$db_exon[1]-1);
			if (!exists ${$sample_intron_co}{$test_db_id1} && !exists ${$sample_intron_co}{$test_db_id2}){
			$db_exon_count+=${$sample_exon_co}{$db_exon_all};
		}
}
}

$db_exon_count_norm=$db_exon_count/(($options{l}*2+$length_intron)/($options{l}*2));
if ($db_intron_count+$db_exon_count_norm>0){
$db_intron_rate=$db_intron_count/($db_intron_count+$db_exon_count_norm);
$db_exon_rate=$db_exon_count_norm/($db_intron_count+$db_exon_count_norm);
}else{
$db_intron_rate="NA";
$db_exon_rate="NA";
}

print EXP "$_\t$db_intron_count\t$db_exon_count_norm\t$db_intron_rate\t$db_exon_rate\n";
}
close IR;
}
}
&irexp;
close EXP;



sub meexp{
foreach my $sid(@all_id){
my $sample_intron_db="db_intron_all_$sid";
my $sample_intron_co="intron_id_$sid";
open ME,"<$options{g}.mutex_exon"||die;
readline ME;
open EXP,">>$options{o}./$sid.exp";
#DGENE.1677      chr1D_part1-153772435-153772512 chr1D_part1-153772326-153773321 DGENE.1677.9   DGENE.1677.5;DGENE.1677.6;TraesCS1D01G131500.1   chr1D_part1-153772513-153772603 chr1D_part1-153772326-153773321 +
print EXP "#ME_gene\tforward_exon\tforward_iei\tforward_trans\tback_trans\tback_exon\tback_iei\tstrand\tforward_count\tback_count\tforward_rate\tback_rate\n";

while (<ME>){
my ($db_first_intron_count,$db_second_intron_count,$db_first_rate,$db_second_rate);
chomp;
my @mx=split(/\s+/,$_);
my @first_exon=split(/-/,$mx[1]);
my @first_iei=split(/-/,$mx[2]);
my @second_exon=split(/-/,$mx[5]);
my @second_iei=split(/-/,$mx[6]);
my $chr=$first_exon[0];
my $first_i=join("-",$chr,$first_iei[1],$first_exon[1]-1);
my $first_ii=join("-",$chr,$first_exon[2]+1,$first_iei[2]);
my $second_i=join("-",$chr,$second_iei[1],$second_exon[1]-1);
my $second_ii=join("-",$chr,$second_exon[2]+1,$second_iei[2]);
$db_first_intron_count=$db_second_intron_count=0;
	foreach my $db_intron_all (@{${$sample_intron_db}{$chr}}){
	if (($db_intron_all eq $first_i) || ($db_intron_all eq $first_ii)){
	$db_first_intron_count+=${$sample_intron_co}{$db_intron_all};
	}
	if (($db_intron_all eq $second_i) || ($db_intron_all eq $second_ii)){
	$db_second_intron_count+=${$sample_intron_co}{$db_intron_all};
	}
}

if ($db_first_intron_count+$db_second_intron_count>0){
$db_first_rate=$db_first_intron_count/($db_first_intron_count+$db_second_intron_count);
$db_second_rate=$db_second_intron_count/($db_first_intron_count+$db_second_intron_count);
}else{
$db_first_rate="NA";
$db_second_rate="NA";
}

print EXP "$_\t$db_first_intron_count\t$db_second_intron_count\t$db_first_rate\t$db_second_rate\n";
}
close ME;
}
}
&meexp;
close EXP;
