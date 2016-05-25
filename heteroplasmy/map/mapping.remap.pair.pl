#!/usr/bin/perl -w
use strict;
#use lib '/home/zry0510';
use zry qw/GetOptions Open/;

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
	'ref_mapping=s' => \$ref_mapping,
	'ref_mp=s' => \$ref_mp,
	'd=i'      => \$d,
	'r=i'	=> \$r,
    'uniq=i'=> \$uniq,
	);

die "$usage\n" if ($help || !$in || !$out);

$ref_mapping||="/home/fs01/rz253/project/reference/hg19.rCRS/hg19.rCRS";
$ref_mp||="/home/fs01/rz253/project/reference/rCRS/chrM.fa";
$d||=1000000;
$r||=1;
$uniq||=1;
#my $ref_mt="/home/zry0510/zry/genomes/bowtie2_index/rCRS/chrM.fa";

my $IN=Open($in);
#open (OUT,">$out");

my $outall = "$out/mapping.all.sh";
open (OA,">$outall");

my $bowtie="bowtie2";
my $picard="/home/fs01/rz253/bin/picard/picard-tools-2.3.0/picard.jar";
my $gatk="/home/fs01/rz253/bin/gatk/";
my $uniqperl="/home/fs01/rz253/project/uk10k/twins/heteroplasmy/map/keep.best.pl";

my %sample;
while (<$IN>){
  chomp;
  my ($read1,$read2,$dir)=split;
  my $sample=$read1;
  ($sample)=$read1=~/(.*)\.R1\.fastq/;
  my $diro="$out/$sample";
#  make_path($out/$sample);
  unless (-e $diro){
    `mkdir -p $out/$sample`;
  }
  my $outfile="$out/$sample/$sample.map.sh";
  open (OUT,">$outfile");
  print OUT "$bowtie -p 6 --local -x $ref_mapping -1 $dir/$read1 -2 $dir/$read2 > $out/$sample/$sample.sam 2>>  $out/$sample/$sample.mapping.log\n";

    if ($uniq != 0){
        print OUT "mv $out/$sample/$sample.sam $out/$sample/$sample.temp.sam\n";
        print OUT "perl $uniqperl -chr chrM -in $out/$sample/$sample.temp.sam -out  $out/$sample/$sample.sam 2>> $out/$sample/$sample.log\n";
    }

    print OUT "java -jar $picard SamFormatConverter INPUT=$out/$sample/$sample.sam  OUTPUT=$out/$sample/$sample.bam  2>> $out/$sample/$sample.log\n";

    print OUT "java -jar $picard AddOrReplaceReadGroups INPUT=$out/$sample/$sample.bam OUTPUT=$out/$sample/$sample.gp.bam SORT_ORDER=coordinate RGID=$sample.gp RGLB=$sample.gb RGPL=illumina RGSM=$sample.gp RGPU=barcode 2>>$out/$sample/$sample.log\n";
    print OUT "java -jar $picard BuildBamIndex INPUT=$out/$sample/$sample.gp.bam 2>>$out/$sample/$sample.log\n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $out/$sample/$sample.gp.bam -R $ref_mp -o $out/$sample/$sample.indel.list 2>>$out/$sample/$sample.log\n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T IndelRealigner -l INFO -I $out/$sample/$sample.gp.bam -R $ref_mp -targetIntervals $out/$sample/$sample.indel.list -o $out/$sample/$sample.unsort.bam 2>>$out/$sample/$sample.log\n";
    print OUT "java -jar $picard SortSam INPUT=$out/$sample/$sample.unsort.bam OUTPUT=$out/$sample/$sample.sort.bam SORT_ORDER=coordinate 2>>$out/$sample/$sample.log \n";
    print OUT "java -jar $picard BuildBamIndex INPUT=$out/$sample/$sample.sort.bam 2>>$out/$sample/$sample.log\n";
#recalibration
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $ref_mp -I $out/$sample/$sample.sort.bam -o $out/$sample/$sample.snp.vcf 2>>$out/$sample/$sample.log\n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -I  $out/$sample/$sample.sort.bam -R $ref_mp -knownSites $out/$sample/$sample.snp.vcf  -o $out/$sample/$sample.recalibration_report.grp 2>>$out/$sample/$sample.log \n";
    print OUT "java -Xmx4g -jar $gatk/GenomeAnalysisTK.jar -T PrintReads -R $ref_mp -I $out/$sample/$sample.sort.bam -BQSR $out/$sample/$sample.recalibration_report.grp -o $out/$sample/$sample.sort.recal.bam 2>>$out/$sample/$sample.log\n";

  if ($r==1){
      print OUT "java -jar $picard MarkDuplicates REMOVE_DUPLICATES=true INPUT=$out/$sample/$sample.sort.recal.bam OUTPUT=$out/$sample/$sample.sort.rmdup.bam M=$out/$sample/$sample.duplicate 2>>$out/$sample/$sample.log \n";
      print OUT "java -jar $picard BuildBamIndex INPUT=$out/$sample/$sample.sort.rmdup.bam OUTPUT=$out/$sample/$sample.sort.rmdup.bam.bai 2>>$out/$sample/$sample.log\n";
      print OUT "samtools mpileup -d $d -f $ref_mp $out/$sample/$sample.sort.rmdup.bam > $out/$sample/$sample.mp   2>>$out/$sample/$sample.log \n";
  }else{

  print OUT "java -jar $picard/BuildBamIndex.jar INPUT=$out/$sample/$sample.sort.bam OUTPUT=$out/$sample/$sample.sort.bam.bai 2>>$out/$sample/$sample.log\n";
  print OUT "samtools mpileup -d $d -f $ref_mp $out/$sample/$sample.sort.bam > $out/$sample/$sample.mp   2>>$out/$sample/$sample.log \n";
}

# `rm $out/$sample/$sample.sam 2>>$out/$sample/$sample.log`;
# `rm $out/$sample/$sample.bam $out/$sample/$sample.sort.bam 2>>$out/$sample/$sample.log &`;
 print OA "sh $out/$sample/$sample.map.sh\n";
}

close OA;
close $IN;

