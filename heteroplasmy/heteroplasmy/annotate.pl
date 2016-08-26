#!/usr/bin/perl -w
use strict;
use zry qw/GetOptions Open basename dirname/;

my $usage = <<'USAGE';

Usage:
	perl annotate.pl 
    annotate heteroplasmy and homoplasmy site

	-in	FILE	input heteroplasmy file
	-out	FILE	output file
    -out2 FILE  output file for diease variants
	-h	HELP	help information

USAGE

my ( $in, $out, $out2, $pre, $help );

GetOptions(
    'in=s'    => \$in,
    'out=s'   => \$out,
    'out2=s'  => \$out2,
    'pre=s'   => \$pre,
);

die "$usage\n" if ( $help || !$in || !$out || !$out2);

#disease file
my $coding="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/coding.csv";
my $tRNA="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/tRNA.csv";

#region file
my $region="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/mtDNA.annotation";

#pscore file
my $pscore="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/pathogenic_score.csv";

#tRNA file
my $tcodon="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/tRNA.pos";

#CADD score
my $CADD2="/home/fs01/rz253/project/uk10k/twins/processed/MT.reads/200individuals/annotate/mt.cadd.1.3.tsv";

#my $poly="/home/fs01/rz253/project/1000genome/compare/het/het.5digit.cover3/data/polymorphic.site";

my @lcom=(66 .. 71,303 .. 309,514 .. 523,12418 .. 12425, 16184 .. 16193);
#map{print "$_\t"}@lcom;

my $C=Open($coding);
my $T=Open($tRNA);

$pre||='ann';

my %dis;
while(<$C>){
    chomp;
    my @line=split ',',$_;
    map{$_=~s/"//g}@line;
#    print @line;
    my ($pos,$name,$change,$type)=@line[0,3,4,5];
    my ($ori,$new)=split '-',$change;
    my $way;
    $type=~/syn/ ? $way='SY':$type=~/frame/?$way='FS':$type=~/non/?$way='noncoding':$way='NS';
    %{$dis{$name}}=('pos'=>$pos,'name'=>$name,'ori'=>$ori,'new'=>$new,'all'=>$_,'type'=>$way);
}
close $C;
#print $coding{"C3310T"}{'pos'};

my %trna;
while(<$T>){
    chomp;
    my @line=split ',',$_;
    map{$_=~s/"//g}@line;
#    print @line;
    my ($pos,$name,$change)=@line[0,3,4];
    my $ori = substr $change,0,1;
    my $new = substr $change,-1,1;
#    my ($ori,$new)=split '-',$change;
    %{$dis{$name}}=('pos'=>$pos,'name'=>$name,'ori'=>$ori,'new'=>$new,'all'=>$_,'type'=>'tRNA');
}
#print $trna{"T582C"}{'pos'};
close $T;

my %site;
my $R=Open($region);
while(<$R>){
    chomp;
    my ($pos,$ref,$ann)=(split)[0,1,2];
    $ann='NA' unless defined $ann;
    $site{$pos}=$ann;
}
close $R;

my $P=Open($pscore);
my %score;
while(<$P>){
    chomp;
    next if /^Position/;
    my ($pos,$change,$s1,$s2)=(split ',',$_)[0,1,-2,-1];
    my ($t1,$t2)=(split '>',$change);
    my $name=$t1.$pos.$t2;
    %{$score{$name}}=('pos'=>$pos,'s1'=>$s1,'s2'=>$s2,'ori'=>$t1,'new'=>$t2);
}
close $P;

my $A=Open($tcodon);
my %tpos;
while(<$A>){
    chomp;
    my ($type,@pos)=(split);
    map{$tpos{$_}=$type}@pos;
}
#print @anti;
close $A;

my $CAdd2=Open($CADD2);
my %cadd2;
while(<$CAdd2>){
    chomp;
    next if /^#/;
    my ($chr,$pos,$t1,$t2,$rawscore,$phred)=(split);
    my $name=$t1.$pos.$t2;
    %{$cadd2{$name}}=('pos'=>$pos,'rawscore'=>$rawscore,'phred'=>$phred,'ori'=>$t1,'new'=>$t2);
}
close $CAdd2;

my $IN=Open($in);
open O,">$out";
open O2,">$out2";

print O2 "change|tag|tRNA|rawscore|phred|detail\n";
foreach my $name (keys %dis){
    if (exists $cadd2{$name}{'rawscore'} or exists $score{$name} ){
        print O2 "$name|D|$dis{$name}{'type'}|$cadd2{$name}{'rawscore'}|$cadd2{$name}{'phred'}|$dis{$name}{'all'}\n";
    }
}



#print O "name\tpos\tref\tdnadepth\tdmajor\tdmajorf\tdminor\tdminorf\tdreff\trnadepth\t
print O "chr|pos|ref|cover|A|T|G|C|minorF|mle|tag|A1|A1f|A2|A2f|sampleID|type|lane|MajorA|Majorf|MinorA|Minorf|reff|altA|altf|redetail|retype|mutype|tRNA|pscore1|pscore2|rawscore|phred|change|transOrtranv|disease\n";
while(<$IN>){
    chomp;
    my ($chr,$pos,$ref,$cover,$A,$T,$G,$C,$minorF,$mle,$tag,$A1,$A1f,$A2,$A2f,$sampleID,$type,$lane,$MajorA,$Majorf,$MinorA,$Mionrf,$reff,$alt,$altf) = split ';',$_;
    next if ($pos eq "pos" or !defined $alt);
    map{print O "$_|"} ($chr,$pos,$ref,$cover,$A,$T,$G,$C,$minorF,$mle,$tag,$A1,$A1f,$A2,$A2f,$sampleID,$type,$lane,$MajorA,$Majorf,$MinorA,$Mionrf,$reff,$alt,$altf);
    my $name=$ref.$pos.$alt;

    ###position annotation
    my $re=$site{$pos};
    print O "$name|$re|";
    if ($re=~/Control/){
        print O "control_region|NA|NA|NA|";
    }elsif ($re=~/tRNA/){
        print O "tRNA|";
        if (exists $tpos{$pos}){
            print O "$tpos{$pos}|NA|NA|";
        }else{
            print O "NA|NA|NA|";
        }
    }elsif ($re=~/rRNA/){
        print O "rRNA|NA|NA|NA|";
    }elsif ($re=~/(ND)|(CO)|(ATP)|(CYB)/){
        print O "coding|";
        if (exists $score{$name}){
            print O "NS|$score{$name}{'s1'}|$score{$name}{'s2'}|";
        }else{
            print O "SY|NA|NA|";
        }
    }else{
        print O "NA|NA|NA|NA|";
    }

#CADD score
    print O "$cadd2{$name}{'rawscore'}|$cadd2{$name}{'phred'}|";

#base change 
    my @change=sort {$a cmp $b} ($ref,$alt);
    my $two=$change[0].$change[1];
    print O "$two|";

#transition or transverstion
    if ($two eq 'AG' or $two eq 'CT'){
        print O "transition|";
    }else{
        print O "transversion|";
    }


##disese
    if (exists $dis{$name}){
        print O "$dis{$name}{'all'}\n";
    }else{
        print O "NA\n";
    }

}
        

