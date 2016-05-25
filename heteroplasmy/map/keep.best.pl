#!/usr/bin/perl -w
use strict;
#use lib '/home/zry0510';
use zry qw/GetOptions Open make_path basename/;

my $usage=<<'USAGE';

Usage:
	perl keep.best.pl
    last update Dec 21st 2015, add mapping score filter ($map)

	-in	    FILE	input sequencing fasta file list
    -out	FILE	output dir
    -chr    STR     chromosom (only keep this chromosome)[chrM]
    -max    FLOAT   max mismatch rate [0.06]
	-h	    HELP	help information
                                      
USAGE


my ($in,$out,$chr,$max,$map,$help);
                                     
GetOptions(
	'in=s'	=> \$in,
	'out=s' => \$out,
	'chr=s' => \$chr,
    	'max=i' => \$max,
        'map=i' => \$map,
	);

die "$usage\n" if ($help || !$in || !$out);

#$chr||='';
$max||=0.06;
$map||=0;

my $IN=Open($in);
my $OUT=Open($out,">");

my $last='';
my $read1="";
my $read2="";
my $isread2="no";
my $discard="no";
my $score=0;
my $xscore=0;

while (<$IN>){
    chomp;
    if (/^@/){
        if (/HD/ or /chrM/ or /PG/){
            print $OUT "$_\n";
        }
            next;
    }

    my ($name,$flag,$c,$pos,$as,$match,$pair,$pair_p,$dis,$seq,$qua,@tag)=split;
    if ($name eq '' or $name ne $last){
        $read1=$_;
        $last=$name;
        $isread2="no";
        $score=0;
        $xscore=0;
        $discard='no';
    }else{
        $read2=$_;
        $isread2="yes";
        next if $discard eq "yes";
    }

    if ($c eq '*' or $as < $map){
        $discard="yes";
        next;
    }
    
    my $readlen=length($seq);

    if ($chr){
        if ($c ne $chr){
            $discard="yes";
            next;
        }
    }
    
    if ($_!~/NM:i:/){
        $discard="yes";
        next;
    }
    my ($miss) = /NM:i:(\d+)/;
    if ($miss > $max*$readlen){
        $discard="yes";
        next;
    }

    if ($_ !~ /XS:i:/){
        #print $OUT "$_\n";
        $discard="no";
    }else{
        my ($first) = /AS:i:(-?\d+|\d+)/;
        my ($second) = /XS:i:(-?\d+|\d+)/;
        $score=$score+$first;
        $xscore=$xscore+$second;
       # print "$first\t$second\t\t";
       # if ($first == $second){
       #     $discard="yes";
       #     next;
       # }else{
       #     $discard="no";
       # }
       # print "$first\t$second\t\t";
        #print $OUT "$_\n";
    }

    if ($isread2 eq "yes"){
        if ($score > $xscore){
            $discard="no";
        #    print "$score $xscore\n";
        }

        if ($discard eq "no"){
        print $OUT "$read1\n$read2\n";
        }
    }
}

close $IN;
close $OUT;
