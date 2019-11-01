#! /usr/bin/perl -w
use warnings;
use strict;
use Getopt::Std;
my $usage =
"$0 <-i Aligned_bam.txt> [options]

This script can split junction reads and get mapping location from no-softclipping mapped bam file

-i     	  list.txt include the absolute path of the bam file
-o        output dir                                             [default: ./ReadsJ]
-m int    max distance to define concordant mapped               [default: 15000]
-d str    strict/loose                                           [default: strict]
-h int    minimum overhang of junction reads                     [default: 6]
-a int    maximum alignment times                                [default: 1]
-s int    maximum miss match bases for downstream analysis       [default: 2]
";
my $degree="Y|N";
my %options=();
getopts("i:o:m:d:h:a:s:",\%options);

if (!exists $options{i}){
die "$usage";
}
my $dir = $ENV{'PWD'};

if (!exists $options{o}){
        $options{o}="./ReadsJ/";
        }
if (-d "$dir/$options{o}"){
        }else{
        `mkdir $dir/$options{o}`;
}

if (!exists $options{m}){
        $options{m}=15000;
}

if ($options{o}!~/\/$/){
        $options{o}=$options{o}."/";
}

if (!exists $options{d}){
	$options{d}="strict";
}elsif($options{d} eq "loose"){
        $degree="Y|N|D|F";
}else{
}

if (!exists $options{h}){
	$options{h}=6;
}

if (!exists $options{a}){
        $options{a}=1;
}

if (!exists $options{s}){
        $options{s}=2;
}

open IN,"$options{i}"||die;
while  (my $sam=<IN>){
my %readsid;
chomp $sam;
next if ($sam eq "");
my $out_id=(split(/\//,$sam))[-1];
open OUT,">$options{o}./$out_id.loc"||die;
print OUT "Chr\tStart\tEnd\tReads\tStrand\tAligment_times\tF_overhang\tL_overhang\tMiss_match\tMate_type\tEorI\n";
open UNIQ,">$options{o}./$out_id.loc.fsim"||die;
foreach (`samtools view $sam`){
my (@line,@matchAM,@matchAI,@N,@number,@matchM,@matchI,@matchM2,@matchI2,@lastM,@lastI,@num,@num0,@num1,@num2,@num3,@matchM0,@matchI0,@matchM1,@matchI1,@matchM3,@matchI3,@matchM4,@matchI4);
my ($count_n,$match,$miss,$reads_num,$reads_strand,$type,$mate_pos,$con,$more,$next_more,$next_next_more,$next_next_next_more,$first,$last,$matchAI,$matchAM,$map_sta,$map_end,$map_sta1,$map_end1,$map_sta2,$map_end2,$map_sta3,$map_end3,$map_sta4,$map_end4,$map_sta5,$map_end5,$intron_sta,$intron_end,$intron_sta1,$intron_end1,$intron_sta2,$intron_end2,$intron_sta3,$intron_end3,$intron_sta4,$intron_end4,$matchM,$matchI,$matchM2,$matchI2,$lastM,$lastI,$matchM0,$matchI0,$matchM1,$matchI1,$matchM3,$matchI3,$matchM4,$matchI4,$fid,$fid_1,$fid_2,$fid_3);
chomp;
@line=split(/\s+/,$_);
if ($line[5] ne "*"){
$count_n=$line[5]=~tr/N/N/;

if ($_=~/NH\:i\:(\d+)/){
	$match=$1;
	}else{
}
if ($_=~/NM\:i\:(\d+)/){
	$miss=$1;
	}else{
}
$reads_num=(int($line[1]) & 192) >>6;
$reads_strand=(int($line[1]) & 32);
if ($reads_strand==32){
	$type="+";
	}else{
	$type="-";
}

$mate_pos=abs($line[8]);
if($line[6] eq "=" && $mate_pos<=$options{m} && $mate_pos!=0){
$con="Y";#concordant
}elsif($line[6] eq "*"){
$con="N";#mate unmap
}elsif($line[6] ne "=" && $line[6] ne "*"){
$con="D";#mate map different chr 
}else{
$con="F";#mate map too far
}

if ($count_n==0){
#HISEQ:790:CAKCCANXX:8:1310:18004:53914  385     chr1A_part1     5059    0       119M    chr1B_part1     80713951        0       TCAAAGATTAAGCCATGCATGTGCAAGTATGAACCAATTTGAACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTACGTGCTACTCGGATAACCGTAGTAA FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<BBFFFFFFBBFFFFFFFFFFBBFFFFB<FBFFFFF/FBFFFFFFFFFFFFFFFFFFBFFF/BBBFBFBFFFF<FFFFF/B AS:i:0  XN:i:0 XM:i:0   XO:i:0  XG:i:0  NM:i:0  MD:Z:119        YT:Z:UU NH:i:20 CC:Z:=  CP:i:28492      HI:i:0
	$more=0;
	@matchAM=($line[5]=~/(\d+)M/g);
		foreach $matchAM(@matchAM){
			$more+=$matchAM;
	}
	@matchAI=($line[5]=~/(\d+)D/g);
		foreach $matchAI(@matchAI){
			$more+=$matchAI;
	}
			$map_sta=$line[3];
			$map_end=$line[3]+$more-1;
	print OUT "$line[2]\t$map_sta\t$map_end\t$line[0]/$reads_num\t$type\t$match\t$miss\t$more\t$more\t$con\tE\n";
	if ($match<=$options{a} && $miss<=$options{s} && $more>=$options{h} && $con=~/$degree/){
		$fid=join("-",$line[2],$map_sta,$map_end,"E");
		$readsid{$fid}+=1;
	
}
}
if ($count_n==1){
#HISEQ:790:CAKCCANXX:8:1310:18004:53914  385     chr1A_part1     5059    0       119M    chr1B_part1     80713951        0       TCAAAGATTAAGCCATGCATGTGCAAGTATGAACCAATTTGAACTGTGAAACTGCGAATGGCTCATTAAATCAGTTATAGTTTGTTTGATGGTACGTGCTACTCGGATAACCGTAGTAA FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF<BBFFFFFFBBFFFFFFFFFFBBFFFFB<FBFFFFF/FBFFFFFFFFFFFFFFFFFFBFFF/BBBFBFBFFFF<FFFFF/B AS:i:0  XN:i:0 XM:i:0   XO:i:0  XG:i:0  NM:i:0  MD:Z:119        YT:Z:UU NH:i:20 CC:Z:=  CP:i:28492      HI:i:0
	$more=0;
	$next_more=0;
	@N=split(/N/,$line[5]);
	@num=split(/M|I|D/,$N[0]);
	@matchM=($N[0]=~/(\d+)M/g);
		foreach $matchM(@matchM){
			$more+=$matchM;
	}
#8M200        N     74M
	@matchI=($N[0]=~/(\d+)D/g);
		foreach $matchI(@matchI){
			$more+=$matchI;
	}
	@matchM2=($N[1]=~/(\d+)M/g);
		foreach $matchM2(@matchM2){
			$next_more+=$matchM2;
	}
	@matchI2=($N[1]=~/(\d+)D/g);
		foreach $matchI2(@matchI2){
			$next_more+=$matchI2;
	}
			$intron_sta=$line[3]+$more;
			$intron_end=$intron_sta+$num[-1]-1;
			$map_sta1=$line[3];
			$map_end1=$intron_sta-1;
			$map_sta2=$intron_end+1;
			$map_end2=$intron_end+$next_more;
	print OUT "$line[2]\t$intron_sta\t$intron_end\t$line[0]/$reads_num\t$type\t$match\t$miss\t$more\t$next_more\t$con\tI\n";
	print OUT "$line[2]\t$map_sta1\t$map_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$more\t$next_more\t$con\tE\n";
	print OUT "$line[2]\t$map_sta2\t$map_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$more\t$next_more\t$con\tE\n";
	if ($match<=$options{a} && $miss<=$options{s} && $more>=$options{h} && $next_more>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta,$intron_end,"I");
		$fid_2=join("-",$line[2],$map_sta1,$map_end1,"E");
		$fid_3=join("-",$line[2],$map_sta2,$map_end2,"E");
                $readsid{$fid_1}+=1;
		$readsid{$fid_2}+=1;
		$readsid{$fid_3}+=1;

}
}

if ($count_n==2){
undef @lastM;
undef @lastI;
	$more=0;
	$next_more=0;
	$first=0;
	$last=0;
	@N=split(/N/,$line[5]);
	@num0=split(/M|I|D/,$N[0]);
	@num1=split(/M|I|D/,$N[1]);
#8M200        N     74M20 N 8M1I2M
	@lastM=($N[2]=~/(\d+)M/g);
		foreach $lastM(@lastM){
			$last+=$lastM;
	}
	@lastI=($N[2]=~/(\d+)D/g);
		foreach $lastI(@lastI){
			$last+=$lastI;
	}
	@matchM0=($N[0]=~/(\d+)M/g);
		foreach $matchM0(@matchM0){
			$more+=$matchM0;
			$first+=$matchM0;
	}
	@matchI0=($N[0]=~/(\d+)D/g);
		foreach $matchI0(@matchI0){
		$more+=$matchI0;
		$first+=$matchI0;
	}
			$intron_sta1=$line[3]+$more;
			$intron_end1=$intron_sta1+$num0[-1]-1;
			$map_sta1=$line[3];
			$map_end1=$intron_sta1-1;
			$map_sta2=$intron_end1+1;
	print OUT "$line[2]\t$intron_sta1\t$intron_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta1\t$map_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
		        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta1,$intron_end1,"I");
		$fid_2=join("-",$line[2],$map_sta1,$map_end1,"E");
                $readsid{$fid_1}+=1;
		$readsid{$fid_2}+=1;
}

	@matchM1=($N[1]=~/(\d+)M/g);
		foreach $matchM1(@matchM1){
			$intron_end1+=$matchM1;
	}
	@matchI1=($N[1]=~/(\d+)D/g);
		foreach $matchI1(@matchI1){
			$intron_end1+=$matchI1;
	}
	@matchM2=($N[2]=~/(\d+)M/g);
		foreach $matchM2(@matchM2){
			$next_more+=$matchM2;
	}
	@matchI2=($N[2]=~/(\d+)D/g);
		foreach $matchI2(@matchI2){
			$next_more+=$matchI2;
	}
			$intron_sta2=$intron_end1+1;
			$intron_end2=$intron_sta2+$num1[-1]-1;
			$map_end2=$intron_sta2-1;
			$map_sta3=$intron_end2+1;
			$map_end3=$intron_end2+$next_more;
	print OUT "$line[2]\t$intron_sta2\t$intron_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta2\t$map_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
	print OUT "$line[2]\t$map_sta3\t$map_end3\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta2,$intron_end2,"I");
                $fid_2=join("-",$line[2],$map_sta2,$map_end2,"E");
                $fid_3=join("-",$line[2],$map_sta3,$map_end3,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;
                $readsid{$fid_3}+=1;

}
}

if ($count_n==3){
undef @lastM;
undef @lastI;
	$more=0;
	$next_more=0;
	$next_next_more=0;
	$first=0;
        $last=0;
	@N=split(/N/,$line[5]);
	@num0=split(/M|I|D/,$N[0]);
	@num1=split(/M|I|D/,$N[1]);
	@num2=split(/M|I|D/,$N[2]);
	        @lastM=($N[3]=~/(\d+)M/g);
                foreach $lastM(@lastM){
                        $last+=$lastM;
        }
        @lastI=($N[3]=~/(\d+)D/g);
                foreach $lastI(@lastI){
                        $last+=$lastI;
        }
	@matchM0=($N[0]=~/(\d+)M/g);
		foreach $matchM0(@matchM0){
			$more+=$matchM0;
			$first+=$matchM0;
	}
	@matchI0=($N[0]=~/(\d+)D/g);
		foreach $matchI0(@matchI0){
		$more+=$matchI0;
		$first+=$matchI0;
	}
			$intron_sta1=$line[3]+$more;
			$intron_end1=$intron_sta1+$num0[-1]-1;
			$map_sta1=$line[3];
			$map_end1=$intron_sta1-1;
			$map_sta2=$intron_end1+1;
	print OUT "$line[2]\t$intron_sta1\t$intron_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta1\t$map_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
	        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta1,$intron_end1,"I");
                $fid_2=join("-",$line[2],$map_sta1,$map_end1,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;

}
	@matchM1=($N[1]=~/(\d+)M/g);
		foreach $matchM1(@matchM1){
			$intron_end1+=$matchM1;
	}
	@matchI1=($N[1]=~/(\d+)D/g);
		foreach $matchI1(@matchI1){
			$intron_end1+=$matchI1;
	}
	@matchM2=($N[2]=~/(\d+)M/g);
		foreach $matchM2(@matchM2){
			$next_more+=$matchM2;
	}
	@matchI2=($N[2]=~/(\d+)D/g);
		foreach $matchI2(@matchI2){
			$next_more+=$matchI2;
	}
#          4M200        N     74M82      N        55M260     N       7M
			$intron_sta2=$intron_end1+1;
			$intron_end2=$intron_sta2+$num1[-1]-1;
			$map_end2=$intron_sta2-1;
			$map_sta3=$intron_end2+1;
			$map_end3=$intron_end2+$next_more;
	print OUT "$line[2]\t$intron_sta2\t$intron_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta2\t$map_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
	print OUT "$line[2]\t$map_sta3\t$map_end3\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta2,$intron_end2,"I");
                $fid_2=join("-",$line[2],$map_sta2,$map_end2,"E");
                $fid_3=join("-",$line[2],$map_sta3,$map_end3,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;
                $readsid{$fid_3}+=1;

}
	@matchM3=($N[3]=~/(\d+)M/g);
		foreach $matchM3(@matchM3){
			$next_next_more+=$matchM3;
	}
	@matchI3=($N[3]=~/(\d+)D/g);
		foreach $matchI3(@matchI3){
			$next_next_more+=$matchI3;
	}
			$intron_sta3=$map_end3+1;
			$intron_end3=$intron_sta3+$num2[-1]-1;
			$map_sta4=$intron_end3+1;
			$map_end4=$intron_end3+$next_next_more;
	print OUT "$line[2]\t$intron_sta3\t$intron_end3\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta4\t$map_end4\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta3,$intron_end3,"I");
                $fid_2=join("-",$line[2],$map_sta4,$map_end4,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;

}
}

if ($count_n==4){
#          4M200        N     74M82      N        55M260     N       7M5I30    N    5M   
undef @lastM;
undef @lastI;
	$more=0;
	$next_more=0;
	$next_next_more=0;
	$next_next_next_more=0;
	$first=0;
	$last=0;
	@N=split(/N/,$line[5]);
	@num0=split(/M|I|D/,$N[0]);
	@num1=split(/M|I|D/,$N[1]);
	@num2=split(/M|I|D/,$N[2]);
	@num3=split(/M|I|D/,$N[3]);
	@lastM=($N[4]=~/(\d+)M/g);
                foreach $lastM(@lastM){
                        $last+=$lastM;
        }
        @lastI=($N[4]=~/(\d+)D/g);
                foreach $lastI(@lastI){
                        $last+=$lastI;
        }

	@matchM0=($N[0]=~/(\d+)M/g);
		foreach $matchM0(@matchM0){
			$more+=$matchM0;
			$first+=$matchM0;
	}
	@matchI0=($N[0]=~/(\d+)D/g);
		foreach $matchI0(@matchI0){
		$more+=$matchI0;
		$first+=$matchI0;
	}
			$intron_sta1=$line[3]+$more;
			$intron_end1=$intron_sta1+$num0[-1]-1;
			$map_sta1=$line[3];
			$map_end1=$intron_sta1-1;
			$map_sta2=$intron_end1+1;
	print OUT "$line[2]\t$intron_sta1\t$intron_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta1\t$map_end1\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta1,$intron_end1,"I");
                $fid_2=join("-",$line[2],$map_sta1,$map_end1,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;

}
	@matchM1=($N[1]=~/(\d+)M/g);
		foreach $matchM1(@matchM1){
			$intron_end1+=$matchM1;
	}
	@matchI1=($N[1]=~/(\d+)D/g);
		foreach $matchI1(@matchI1){
			$intron_end1+=$matchI1;
	}
	@matchM2=($N[2]=~/(\d+)M/g);
		foreach $matchM2(@matchM2){
			$next_more+=$matchM2;
	}
	@matchI2=($N[2]=~/(\d+)D/g);
		foreach $matchI2(@matchI2){
			$next_more+=$matchI2;
	}
#          4M200        N     74M82      N        55M260     N       7M5I20    N   5M
			$intron_sta2=$intron_end1+1;
			$intron_end2=$intron_sta2+$num1[-1]-1;
			$map_end2=$intron_sta2-1;
			$map_sta3=$intron_end2+1;
			$map_end3=$intron_end2+$next_more;
	print OUT "$line[2]\t$intron_sta2\t$intron_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta2\t$map_end2\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
	print OUT "$line[2]\t$map_sta3\t$map_end3\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta2,$intron_end2,"I");
                $fid_2=join("-",$line[2],$map_sta2,$map_end2,"E");
                $fid_3=join("-",$line[2],$map_sta3,$map_end3,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;
                $readsid{$fid_3}+=1;

}
	@matchM3=($N[3]=~/(\d+)M/g);
		foreach $matchM3(@matchM3){
			$next_next_more+=$matchM3;
	}
	@matchI3=($N[3]=~/(\d+)D/g);
		foreach $matchI3(@matchI3){
			$next_next_more+=$matchI3;
	}
			$intron_sta3=$map_end3+1;
			$intron_end3=$intron_sta3+$num2[-1]-1;
			$map_sta4=$intron_end3+1;
			$map_end4=$intron_end3+$next_next_more;
	print OUT "$line[2]\t$intron_sta3\t$intron_end3\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta4\t$map_end4\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta3,$intron_end3,"I");
                $fid_2=join("-",$line[2],$map_sta4,$map_end4,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;

}
	@matchM4=($N[4]=~/(\d+)M/g);
		foreach $matchM4(@matchM4){
			$next_next_next_more+=$matchM4;
	}
	@matchI4=($N[4]=~/(\d+)D/g);
		foreach $matchI4(@matchI4){
			$next_next_next_more+=$matchI4;
	}
			$intron_sta4=$map_end4+1;
			$intron_end4=$intron_sta4+$num3[-1]-1;
			$map_sta5=$intron_end4+1;
			$map_end5=$intron_end4+$next_next_next_more;
	print OUT "$line[2]\t$intron_sta4\t$intron_end4\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tI\n";
	print OUT "$line[2]\t$map_sta5\t$map_end5\t$line[0]/$reads_num\t$type\t$match\t$miss\t$first\t$last\t$con\tE\n";
        if ($match<=$options{a} && $miss<=$options{s} && $first>=$options{h} && $last>=$options{h} && $con=~/$degree/){
                $fid_1=join("-",$line[2],$intron_sta4,$intron_end4,"I");
                $fid_2=join("-",$line[2],$map_sta5,$map_end5,"E");
                $readsid{$fid_1}+=1;
                $readsid{$fid_2}+=1;

}
}
}
}
`gzip $options{o}./$out_id.loc`;
foreach my $readskey (sort keys %readsid){
        my @filterloc=split(/-/,$readskey);
        print UNIQ"$filterloc[0]\t$filterloc[1]\t$filterloc[2]\t$readsid{$readskey}\t$filterloc[3]\n";
}
close UNIQ;
close OUT;
}
close IN;
