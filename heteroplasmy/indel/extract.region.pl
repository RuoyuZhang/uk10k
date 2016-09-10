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


my ($in,$out,$chr,$max,$min,$help);
                                     
GetOptions(
	'in=s'	=> \$in,
	'out=s' => \$out,
	'chr=s' => \$chr,
    	'min=i' => \$min,
        'max=i' => \$max,
	);

die "$usage\n" if ($help || !$in || !$out);

$min||=240;
$max||=290;

my $IN=Open($in);
my $OUT=Open($out,">");

my $i=0;
while (<$IN>){
    chomp;
    if (/^@/){
        if (/HD/ or /chrM/ or /PG/){
            print $OUT "$_\n";
        }
            next;
    }

    my ($name,$flag,$c,$pos,$as,$match,$pair,$pair_p,$dis,$seq,$qua,@tag)=split;

    if ($pos <=$max and $pos >=$min and $c eq 'chrM' and $pair eq '='){
        print $OUT "$_\n";
    }


}

close $IN;
close $OUT;
