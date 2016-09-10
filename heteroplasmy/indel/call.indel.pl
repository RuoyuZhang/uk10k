#!/usr/bin/perl -w
use strict;
#use lib '/home/zry0510';
use zry qw/GetOptions Open/;
use File::Basename;

my $usage=<<'USAGE';

Usage:
	perl mapping.pl

	-in	FILE	input sequencing fasta file list
	-out	FILE	output dir
	-ref_mapping	FILE	reference for mapping
	-ref_mp	        FILE	reference for mpileup
	-r      FLAG    remove PCR duplication, 1 for remove, 0 for not [0]
    -uniq   FLAG    only keep best mapping [0]
	-d	INT	mpileup depth [100000]
	-h	HELP	help information

USAGE

my ($in,$out,$ref_mapping,$ref_mp,$uniq,$d,$r,$help);

GetOptions(
	'in=s'	=> \$in,
	'out=s' => \$out,
	'ref_mp=s' => \$ref_mp,
	'd=i'      => \$d,
	'r=i'	=> \$r,
    'uniq=i'=> \$uniq,
	);

die "$usage\n" if ($help || !$in || !$out);

$ref_mp||="/home/fs01/rz253/project/reference/rCRS/chrM.fa";
my $IN=Open($in);
#open (OUT,">$out");

my $outall = "$out/call.indel.sh";
open (OA,">$outall");

my $bowtie="bowtie2";
my $picard="/home/fs01/rz253/bin/picard/picard-tools-2.3.0/picard.jar";
my $gatk="/home/fs01/rz253/bin/gatk/";
my $extractperl="/home/fs01/rz253/project/uk10k/uk10k/heteroplasmy/indel/extract.region.pl";

my %sample;
while (<$IN>){
    chomp;
    my $samname=basename($_);
    my $sample=$samname;
    ( $sample)=$samname=~/(.*)\.temp\.sam/;
    my $diro="$out/$sample";
#  make_path($out/$sample);
    unless (-e $diro){
        `mkdir -p $out/$sample`;
    } 
    my $outfile="$out/$sample/$sample.indel.sh";
    open (OUT,">$outfile");
    print OUT "perl $extractperl -in $_ -out  $out/$sample/$sample.sam 2> $out/$sample/$sample.log\n";

    print OUT "java -jar $picard SamFormatConverter INPUT=$out/$sample/$sample.sam  OUTPUT=$out/$sample/$sample.bam  2>> $out/$sample/$sample.log\n";

    print OUT "java -jar $picard AddOrReplaceReadGroups INPUT=$out/$sample/$sample.bam OUTPUT=$out/$sample/$sample.gp.bam SORT_ORDER=coordinate RGID=$sample.gp RGLB=$sample.gb RGPL=illumina RGSM=$sample.gp RGPU=barcode 2>>$out/$sample/$sample.log\n";
    print OUT "java -jar $picard BuildBamIndex INPUT=$out/$sample/$sample.gp.bam 2>>$out/$sample/$sample.log\n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $out/$sample/$sample.gp.bam -R $ref_mp -o $out/$sample/$sample.indel.list 2>>$out/$sample/$sample.log\n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T IndelRealigner -l INFO -I $out/$sample/$sample.gp.bam -R $ref_mp -targetIntervals $out/$sample/$sample.indel.list -o $out/$sample/$sample.unsort.bam 2>>$out/$sample/$sample.log\n";
    print OUT "java -jar $picard SortSam INPUT=$out/$sample/$sample.unsort.bam OUTPUT=$out/$sample/$sample.sort.bam SORT_ORDER=coordinate 2>>$out/$sample/$sample.log \n";
    print OUT "java -jar $picard BuildBamIndex INPUT=$out/$sample/$sample.sort.bam 2>>$out/$sample/$sample.log\n";
#recalibration
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $ref_mp -I $out/$sample/$sample.sort.bam -o $out/$sample/$sample.indel.vcf -L $out/$sample/$sample.indel.list -glm INDEL 2>>$out/$sample/$sample.log\n";
    
    print OA "sh $out/$sample/$sample.indel.sh\n";
}

close OA;
close $IN;

