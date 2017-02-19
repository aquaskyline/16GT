#!/usr/bin/perl -w

use File::Basename;
use lib dirname($0)."/lib";
use Inline C;
use Data::Dumper;

($ascInputFN, $sampleName, $reffn) = @ARGV;
die "perl $0 <txt Input> <Sample Name> <Referece FASTA>\n" unless $reffn;

my %ref=();
open my $fa,"$reffn" or die $!;
$/=">";$/=<$fa>;$/="\n";
while (<$fa>){
  chomp;
  (my $id=$_)=~s/\s+.*$//;

  $/=">";
  my $seq=<$fa>;
  chomp $seq;
  $seq=~s/\s+//g;
  $/="\n";
  $seq =~ tr/acgtn/ACGTN/;

  $ref{$id} = $seq;
}
close $fa;

open $ascInputFH, "gzip -dcf $ascInputFN |" or die $!;
$l=<$ascInputFH>;

# Input description
# [Chromosome] [Chromosome position] [Reference base character] [Pos-1 Reference base character] [pAA] [pAC] [pAG] [pAT] [pCC] [pCG] [pCT] [pGG] [pGT] [pTT] [pAX] [pCX] [pGX] [pTX] [pXX] [pXY] [Genotype] [Optimal Likelihood] [Sub-optimal Likelihood] [pD] [Random Forest Probability] [Random Forest Prediction] [INDEL mode] [INDEL pattern] [INDEL HQ count] [INDEL LQ count] [Total depth] [A depth] [C depth] [G depth] [T depth] [A+] [A-] [C+] [C-] [G+] [G-] [T+] [T-] [AA] [AC] [AG] [AT] [CC] [CG] [CT] [GG] [GT] [TT] [AX] [CX] [GX] [TX] [XX] [XY] [Strand Bias A] [Strand Bias C] [Strand Bias G] [Strand Bias T] [Left Strand Bias A] [Left Strand Bias C] [Left Strand Bias G] [Left Strand Bias T] [Right Strand Bias A] [Right Strand Bias C] [Right Strand Bias G] [Right Strand Bias T] [Base Quality Bias A] [Base Quality Bias C] [Base Quality Bias G] [Base Quality Bias T] [Left Base Quality Bias A] [Left Base Quality Bias C] [Left Base Quality Bias G] [Left Base Quality Bias T] [Right Base Quality Bias A] [Right Base Quality Bias C] [Right Base Quality Bias G] [Right Base Quality Bias T] [Read Position Bias] [Left GC ratio] [Right GC ratio] [Left Consecutive G Base Count] [Right Consecutive G Base Count] [Average Strand Count] [Left Read Quality] [Right Read Quality] [Polyrun Length] [Total Indel HQ Count] [Total Indel LQ Count]

%Rci = ("A"=>0, "C"=>1, "G"=>2, "T"=>3, "X"=>4, "Y"=>5);
%GTci = ("AA"=>0, "AC"=>1, "AG"=>2, "AT"=>3, "CC"=>4, "CG"=>5, "CT"=>6, "GG"=>7, "GT"=>8, "TT"=>9, "CA"=>1, "GC"=>5, "GA"=>2, "TG"=>8, "TC"=>6, "TA"=>3, "AX"=>10, "XA"=>10, "CX"=>11, "XC"=>11, "GX"=>12, "XG"=>12, "TX"=>13, "XT"=>13, "XX"=>14, "XY"=>15, "YX"=>15);
@GTiii = ([0,0],[0,1],[0,2],[0,3],[1,1],[1,2],[1,3],[2,2],[2,3],[3,3],[0,4],[1,4],[2,4],[3,4],[4,4],[4,5]);
$flt_max = fltMaxVal();
$qual_max = qualMaxVal();

print '##fileformat=VCFv4.1'."\n";
print '##FILTER=<ID=LowQual,Description="Low quality">'."\n";
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'."\n";
print '##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">'."\n";
print '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">'."\n";
print '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">'."\n";
print '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">'."\n";
print '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes, for each ALT allele, in the same order as listed">'."\n";
print '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, for each ALT allele, in the same order as listed">'."\n";
print '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">'."\n";
print '##INFO=<ID=FRAC,Number=1,Type=Float,Description="Alternative allele fraction">'."\n";
print '##INFO=<ID=BQFisherPhredRef,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test of ref vs. perfect base qualities">'."\n";
print '##INFO=<ID=BQFisherPhredAlt,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test of alt vs. perfect base qualities">'."\n";
print '##INFO=<ID=BQFisherPhredAlt2,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test of alt2 vs. perfect base qualities">'."\n";
print '##INFO=<ID=BestGT,Number=1,Type=String,Description="Best Genotype among the 16 candidate">'."\n";
print '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">'."\n";
print '##INFO=<ID=DP4,Number=4,Type=Integer,Description="Read depth of A,C,G,T">'."\n";
print '##INFO=<ID=DP8,Number=8,Type=Integer,Description="Read depth of A+,A-,C+,C-,G+,G-,T+,T-">'."\n";
print '##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Genotype with Indel at the genome position">'."\n";
print '##INFO=<ID=IndelMeta,Number=8,Type=String,Description="Details of the deepest and second deepest Indel of the position. Format: OptMode/OptPattern/OptHQ/OptLQ/SubOptMode/SubOptPattern/SubOptHQ/SubOptLQ">'."\n";
print '##INFO=<ID=ILDP,Number=1,Type=Integer,Description="Reads supporting the Indel patterns with average base quality<=13">'."\n";
print '##INFO=<ID=IHDP,Number=1,Type=Integer,Description="Reads supporting the Indel patterns with average base quality>13">'."\n";
print '##INFO=<ID=LeftG,Number=1,Type=Integer,Description="# of consecutive base G on the left of the Indel">'."\n";
print '##INFO=<ID=RightG,Number=1,Type=Integer,Description="# of consecutive base G on the right of the Indel">'."\n";
print '##INFO=<ID=HRun,Number=1,Type=Integer,Description="Largest Contiguous Homopolymer Run of Variant Allele In Either Direction">'."\n";
print '##INFO=<ID=QD,Number=1,Type=Float,Description="Variant (Confidence/Quality) by Depth">'."\n";
print '##INFO=<ID=SBFisherPhredRef,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias of ref allele">'."\n";
print '##INFO=<ID=SBFisherPhredAlt,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias of alt allele">'."\n";
print '##INFO=<ID=SBFisherPhredAlt2,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher\'s exact test to detect strand bias of alt2 allele">'."\n";
print '##INFO=<ID=ReadPosBias,Number=1,Type=Float,Description="Read position bias, [-1,1], 0:no bias, positive: bias into left 100bp, negative: bias into right 100bp">'."\n";
print '##INFO=<ID=DepthBias,Number=1,Type=Float,Description="Read depth balance between the left and right 100bp, [-1,1], 0: perfect balance, positive: deeper in left 100bp, negative: deeper in right 100bp">'."\n";
print '##INFO=<ID=SOR,Number=1,Type=Float,Description="Symmetric Odds Ratio of 2x2 contingency table to detect strand bias">'."\n";
print '##INFO=<ID=LeftGCPercent,Number=1,Type=Integer,Description="The percentage of GC content of the 100bp bases on the left">'."\n";
print '##INFO=<ID=RightGCPercent,Number=1,Type=Integer,Description="The percentage of GC content of the 100bp bases on the right">'."\n";
print '##INFO=<ID=LeftQuality,Number=1,Type=Integer,Description="The approximate average base quality of the 100bp bases on the left">'."\n";
print '##INFO=<ID=RightQuality,Number=1,Type=Integer,Description="The approximate average base quality of the 100bp bases on the right">'."\n";
print '##INFO=<ID=OptSuboptRatio,Number=1,Type=Float,Description="10*log10((Highest posterior probability of the 16 genotypes)/(Second highest posterior probability of the 16 genotypes))">'."\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sampleName\n";

while($l=<$ascInputFH>)
{
  @a=split /\s+/,$l;
  $chr = $a[0];
  $pos = $a[1];
  $refBase = $a[2];
  $preRefBase = $a[3];
  @probGT = @a[4..19];
  $altGT = $altBase = $a[20];
  $optProb = $a[21];
  $subOptProb = $a[22];
  $indelMode = $a[26];
  $indelPattern = $a[27];
  $indelHQCount = $a[28];
  $indelLQCount = $a[29];
  $indelMode2 = $a[30];
  $indelPattern2 = $a[31];
  $indelHQCount2 = $a[32];
  $indelLQCount2 = $a[33];
  $totalDepth = $a[34];
  @depth = @a[35..38];
  @strandCount = @a[39..46];
  @strandBias = @a[63..66];
  @strandBiasLeft = @a[67..70];
  @strandBiasRight = @a[71..74];
  @baseQualityBias = @a[75..78];
  @baseQualityBiasLeft = @a[79..82];
  @baseQualityBiasRight = @a[83..86];
  $ReadPosBias = $a[87];
  $LeftGCPercent = $a[88] * 100;
  $RightGCPercent = $a[89] * 100;
  $LeftGBaseCount = $a[90];
  $RightGBaseCount = $a[91];
  $DepthBalance = $a[92];
  $LeftQuality = $a[93] * 10;
  $RightQuality = $a[94] * 10;
  $polyRunCount = $a[95];
  $totalHQCount = $a[96];
  $totalLQCount = $a[97];

  next if ($optProb == 0);
  next if ($subOptProb == 0);
  if($altGT=~/^[ACGT]+$/) # is a SNP
  {
    # Positive and negative strand count aggregation
    $strandCountAggregate[0] = $strandCount[0] + $strandCount[1];
    $strandCountAggregate[1] = $strandCount[2] + $strandCount[3];
    $strandCountAggregate[2] = $strandCount[4] + $strandCount[5];
    $strandCountAggregate[3] = $strandCount[6] + $strandCount[7];
    $totalDepth += ($totalHQCount * 4);
    $DP = $strandCountAggregate[0] + $strandCountAggregate[1] + $strandCountAggregate[2] + $strandCountAggregate[3] + $totalHQCount + $totalLQCount;

    # Genotype likelihood normailzaiton
    for($i = 0; $i < scalar(@probGT); ++$i) { $probGT[$i] = $probGT[$i] / $optProb; }

    $gt = $GTci{$altBase};
    $altBase =~ s/$refBase//; $altBase =~ s/(.*)\1/$1/g; @altBase = split //,$altBase;
    $filter = "."; # TODO
    $AC = ""; $AF = ""; $AN = 2;
    $FRAC = $FSRef = $FSAlt = $FSAlt2 = $QD = $SBRef = $SBAlt = $SBAlt2 = $BaseQFisherRef = $BaseQFisherAlt = $BaseQFisherAlt2 = 0;
    $OptSuboptRatio = ($subOptProb==0)?($flt_max):($optProb/$subOptProb);
    $OptSuboptRatio = phred2($OptSuboptRatio);
    $GT = $AD = $DPF = $GQ = $PL = "";
    $INFO = "";
    $FORMAT = "";
    $rciRef = $rciAlt = $rciAlt0 = $rciAlt1 = 0;
    $qual = $qualMax = $qual_max;

    if($GTiii[$gt][0] eq $GTiii[$gt][1])
    {
      die "$.: scalar (\@altBase) != 1" if scalar (@altBase) != 1;
      $rciAlt = $Rci{$altBase};
      $rciRef = $Rci{$refBase};
      $AC = "2";
      $AF = 1.000;
      $FSAlt = phred($strandBias[$rciAlt]);
      $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
      $QD = $depth[$rciAlt] / $strandCountAggregate[$rciAlt] * 10;
      $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, $depth[$rciAlt], $strandCountAggregate[$rciAlt]));

      $SOR = 0.0;
      $SOR += (($strandCount[($rciRef*2)]+1) / ($strandCount[($rciRef*2+1)]+1)) * (($strandCount[($rciAlt*2+1)]+1) / ($strandCount[($rciAlt*2)]+1));
      $SOR += (($strandCount[($rciRef*2+1)]+1) / ($strandCount[($rciRef*2)]+1)) * (($strandCount[($rciAlt*2)]+1) / ($strandCount[($rciAlt*2+1)]+1));
      $refRatio = (cmin(($strandCount[($rciRef*2)]+1), ($strandCount[($rciRef*2+1)]+1))/cmax(($strandCount[($rciRef*2)]+1), ($strandCount[($rciRef*2+1)]+1)));
      $altRatio = (cmin(($strandCount[($rciAlt*2)]+1), ($strandCount[($rciAlt*2+1)]+1))/cmax(($strandCount[($rciAlt*2)]+1), ($strandCount[($rciAlt*2+1)]+1)));
      $SOR = $SOR * $refRatio / $altRatio;
      $SOR = clog($SOR);

      $GT = "1/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
      $DPF = $DP; $GQ = int($QD+0.499);
      $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase"}]), phred($probGT[$GTci{"$altBase$altBase"}]));
      $INFO = sprintf("AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;SOR=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $QD, $ReadPosBias, $DepthBalance, $SOR, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
    }
    else
    {
      if(scalar (@altBase) == 1)
      {
        $rciAlt = $Rci{$altBase};
        $rciRef = $Rci{$refBase};
        $AC = "1";
        $AF = 0.500;
        $FRAC = $strandCountAggregate[$rciAlt] / ($strandCountAggregate[$rciAlt]+$strandCountAggregate[$rciRef]);
        $BaseQFisherRef = phred($baseQualityBias[$rciRef]);
        $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
        $FSRef = phred($strandBias[$rciRef]);
        $FSAlt = phred($strandBias[$rciAlt]);
        $QD = ($depth[$rciAlt] + $depth[$rciRef]) / ($strandCountAggregate[$rciAlt] + $strandCountAggregate[$rciRef]) * 10;
        $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, ($depth[$rciAlt] + $depth[$rciRef]), ($strandCountAggregate[$rciAlt] + $strandCountAggregate[$rciRef])));

        $SOR = 0.0;
        $SOR += (($strandCount[($rciRef*2)]+1) / ($strandCount[($rciRef*2+1)]+1)) * (($strandCount[($rciAlt*2+1)]+1) / ($strandCount[($rciAlt*2)]+1));
        $SOR += (($strandCount[($rciRef*2+1)]+1) / ($strandCount[($rciRef*2)]+1)) * (($strandCount[($rciAlt*2)]+1) / ($strandCount[($rciAlt*2+1)]+1));
        $refRatio = (cmin(($strandCount[($rciRef*2)]+1), ($strandCount[($rciRef*2+1)]+1))/cmax(($strandCount[($rciRef*2)]+1), ($strandCount[($rciRef*2+1)]+1)));
        $altRatio = (cmin(($strandCount[($rciAlt*2)]+1), ($strandCount[($rciAlt*2+1)]+1))/cmax(($strandCount[($rciAlt*2)]+1), ($strandCount[($rciAlt*2+1)]+1)));
        $SOR = $SOR * $refRatio / $altRatio;
        $SOR = clog($SOR);

        $GT = "0/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
        $DPF = $DP; $GQ = int($QD+0.499);
        $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase"}]), phred($probGT[$GTci{"$altBase$altBase"}]));
        $INFO = sprintf("AC=%s;AF=%s;AN=%d;FRAC=%.3f;BQFisherPhredRef=%d;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredRef=%d;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;SOR=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $FRAC, $BaseQFisherRef, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSRef, $FSAlt, $QD, $ReadPosBias, $DepthBalance, $SOR, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
      }
      elsif(scalar (@altBase) == 2)
      {
        $rciAlt0 = $Rci{$altBase[0]};
        $rciAlt1 = $Rci{$altBase[1]};
        $AC = "1,1";
        $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
        $BaseQFisherAlt = phred($baseQualityBias[$rciAlt0]);
        $BaseQFisherAlt2 = phred($baseQualityBias[$rciAlt1]);
        $FSAlt = phred($strandBias[$rciAlt0]);
        $FSAlt2 = phred($strandBias[$rciAlt1]);

        $SOR = 0.0;
        $SOR += (($strandCount[($rciAlt0*2)]+1) / ($strandCount[($rciAlt0*2+1)]+1)) * (($strandCount[($rciAlt1*2+1)]+1) / ($strandCount[($rciAlt1*2)]+1));
        $SOR += (($strandCount[($rciAlt0*2+1)]+1) / ($strandCount[($rciAlt0*2)]+1)) * (($strandCount[($rciAlt1*2)]+1) / ($strandCount[($rciAlt1*2+1)]+1));
        $refRatio = (cmin(($strandCount[($rciAlt0*2)]+1), ($strandCount[($rciAlt0*2+1)]+1))/cmax(($strandCount[($rciAlt0*2)]+1), ($strandCount[($rciAlt0*2+1)]+1)));
        $altRatio = (cmin(($strandCount[($rciAlt1*2)]+1), ($strandCount[($rciAlt1*2+1)]+1))/cmax(($strandCount[($rciAlt1*2)]+1), ($strandCount[($rciAlt1*2+1)]+1)));
        $SOR = $SOR * $refRatio / $altRatio;
        $SOR = clog($SOR);

        $QD = ($depth[$rciAlt0] + $depth[$rciAlt1]) / ($strandCountAggregate[$rciAlt0] + $strandCountAggregate[$rciAlt1]) * 10;
        $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, ($depth[$rciAlt0] + $depth[$rciAlt1]), ($strandCountAggregate[$rciAlt0] + $strandCountAggregate[$rciAlt1])));
        $GT = "1/2"; $AD = "$strandCountAggregate[$Rci{$refBase}],$strandCountAggregate[$rciAlt0],$strandCountAggregate[$rciAlt1]";
        $DPF = $DP; $GQ = int($QD+0.499);
        $PL = sprintf("%d,%d,%d,%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase[0]"}]), phred($probGT[$GTci{"$refBase$altBase[1]"}]), phred($probGT[$GTci{"$altBase[0]$altBase[0]"}]), phred($probGT[$GTci{"$altBase[0]$altBase[1]"}]), phred($probGT[$GTci{"$altBase[1]$altBase[1]"}]));
        $INFO = sprintf("AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BQFisherPhredAlt2=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;SBFisherPhredAlt2=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;SOR=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $BaseQFisherAlt2, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $FSAlt2, $QD, $ReadPosBias, $DepthBalance, $SOR, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
      }
      else { die "$.: \@altBase != 1 and != 2 in Het.\n"; }
    }

    $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

    if ((substr($ref{$chr}, ($pos-1), length($refBase))) ne $refBase)
    {
      next;
    }

    printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refBase, (join ",",@altBase), $qual, $filter);
    printf ("%s\t", $INFO);
    printf ("GT:AD:DP:GQ:PL\t");
    printf ("%s\n", $FORMAT);
  }
  elsif($altGT=~/X/) # is a INDEL
  {
    # Positive and negative strand count aggregation
    $indelHQCount2 = 0 if $indelHQCount2 eq "*";
    $indelLQCount2 = 0 if $indelLQCount2 eq "*";
    $strandCountAggregate[0] = $strandCount[0] + $strandCount[1];
    $strandCountAggregate[1] = $strandCount[2] + $strandCount[3];
    $strandCountAggregate[2] = $strandCount[4] + $strandCount[5];
    $strandCountAggregate[3] = $strandCount[6] + $strandCount[7];
    $strandCountAggregate[4] = $indelHQCount + $indelLQCount;
    $strandCountAggregate[5] = $indelHQCount2 + $indelLQCount2;
    $DP = &accumulate(@strandCountAggregate[0..3], $totalHQCount, $totalLQCount);
    $totalDepth += ($totalHQCount * 4);
    $depth[4] += $indelHQCount * 4; $depth[5] += $indelHQCount2 * 4;

    # Unify the Strand Bias and Base Quality Bias by selecting the minimum from both the values from left and right of the Indel position
    $strandBias[4] = $strandBias[5] = 1;
    foreach(@strandBiasLeft, @strandBiasRight)
    { $strandBias[4] = $strandBias[5] = $_ if $_ < $strandBias[4]; }
    $baseQualityBias[4] = $baseQualityBias[5] = 1;
    foreach(@baseQualityBiasLeft, @baseQualityBiasRight)
    { $baseQualityBias[4] = $baseQualityBias[5] = $_ if $_ < $baseQualityBias[4]; }

    # Genotype likelihood normailzaiton
    for($i = 0; $i < scalar(@probGT); ++$i) { $probGT[$i] = $probGT[$i] / $optProb; }

    $gt = $GTci{$altBase};
    $altBase =~ s/$refBase//; $altBase =~ s/(.*)\1/$1/g; @altBase = split //,$altBase;
    $filter = "."; # TODO
    $AC = ""; $AF = ""; $AN = 2;
    $FRAC = $FSRef = $FSAlt = $FSAlt2 = $QD = $SBRef = $SBAlt = $SBAlt2 = $BaseQFisherRef = $BaseQFisherAlt = $BaseQFisherAlt2 = 0;
    $OptSuboptRatio = ($subOptProb==0)?($flt_max):($optProb/$subOptProb);
    $OptSuboptRatio = phred2($OptSuboptRatio);
    $GT = $AD = $DPF = $GQ = $PL = "";
    $INFO = "";
    $FORMAT = "";
    $rciRef = $rciAlt = $rciAlt0 = $rciAlt1 = 0;
    $qual = $qualMax = $qual_max;
    $refPattern = $altPattern = "";

    if($GTiii[$gt][0] eq $GTiii[$gt][1]) #XX
    {
      die "$.: scalar (\@altBase) != 1" if scalar (@altBase) != 1;
      $rciAlt = $Rci{$altBase};
      $rciRef = $Rci{$refBase};
      $DP -= $strandCountAggregate[$rciRef];
      $AC = "2";
      $AF = 1.000;
      $FSAlt = phred($strandBias[$rciAlt]);
      $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
      $QD = $depth[$rciAlt] / $strandCountAggregate[$rciAlt] * 10;
      $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, $depth[$rciAlt], $strandCountAggregate[$rciAlt]));
      $GT = "1/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
      $DPF = $DP; $GQ = int($QD+0.499);
      $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase"}]), phred($probGT[$GTci{"$altBase$altBase"}]));
      $INFO = sprintf("INDEL;AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $QD, $ReadPosBias, $DepthBalance, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
      if($indelMode=~/I$/)
      {
        $refPattern = "$preRefBase";
        $altPattern = "$preRefBase$indelPattern";
      }
      elsif($indelMode=~/D$/)
      {
        $refPattern = "$preRefBase$indelPattern";
        $altPattern = "$preRefBase";
      }
      else
      { die "$.: Unknown Indel Mode: $indelMode\n"; }

      $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

      if ((substr($ref{$chr}, ($pos-1), length($refPattern))) ne $refPattern)
      {
        next;
      }

      printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refPattern, $altPattern, $qual, $filter);
      printf ("%s\t", $INFO);
      printf ("GT:AD:DP:GQ:PL\t");
      printf ("%s\n", $FORMAT);
    }
    else
    {
      if(scalar (@altBase) == 1) #ref{A,C,G,T}X
      {
        $rciAlt = $Rci{$altBase};
        $rciRef = $Rci{$refBase};
        $DP -= int($strandCountAggregate[$rciRef]/2);
        $AC = "1";
        $AF = 0.500;
        $FRAC = $strandCountAggregate[$rciAlt] / ($strandCountAggregate[$rciAlt] + $strandCountAggregate[$rciRef]);
        $FSRef = phred($strandBias[$rciRef]);
        $BaseQFisherRef = phred($baseQualityBias[$rciRef]);
        $FSAlt = phred($strandBias[$rciAlt]);
        $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
        $QD = ($depth[$rciAlt] + $depth[$rciRef]) / ($strandCountAggregate[$rciAlt] + $strandCountAggregate[$rciRef]) * 10;
        $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, ($depth[$rciAlt] + $depth[$rciRef]), ($strandCountAggregate[$rciAlt] + $strandCountAggregate[$rciRef])));
        $GT = "0/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
        $DPF = $DP; $GQ = int($QD+0.499);
        $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase"}]), phred($probGT[$GTci{"$altBase$altBase"}]));
        $INFO = sprintf("INDEL;AC=%s;AF=%s;AN=%d;FRAC=%.3f;BQFisherPhredRef=%d;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredRef=%d;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $FRAC, $BaseQFisherRef, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSRef, $FSAlt, $QD, $ReadPosBias, $DepthBalance, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
        if($indelMode=~/I$/)
        {
          $refPattern = "$preRefBase";
          $altPattern = "$preRefBase$indelPattern";
        }
        elsif($indelMode=~/D$/)
        {
          $refPattern = "$preRefBase$indelPattern";
          $altPattern = "$preRefBase";
        }
        else
        { die "$.: Unknown Indel Mode: $indelMode\n"; }

        $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

        if ((substr($ref{$chr}, ($pos-1), length($refPattern))) ne $refPattern)
        {
          next;
        }

        printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refPattern, $altPattern, $qual, $filter);
        printf ("%s\t", $INFO);
        printf ("GT:AD:DP:GQ:PL\t");
        printf ("%s\n", $FORMAT);
      }
      elsif(scalar (@altBase) == 2) #non-ref{A,C,G,T}X or XY
      {
        if($altGT=~/^XY$/)
        {
          $rciAlt0 = $Rci{$altBase[0]};
          $rciAlt1 = $Rci{$altBase[1]};
          $rciRef = $Rci{$refBase};
          $DP -= $strandCountAggregate[$rciRef];
          $AC = "1,1";
          $BaseQFisherAlt = phred($baseQualityBias[$rciAlt0]);
          $BaseQFisherAlt2 = phred($baseQualityBias[$rciAlt1]);
          $FSAlt = phred($strandBias[$rciAlt0]);
          $FSAlt2 = phred($strandBias[$rciAlt1]);
          $QD = ($depth[$rciAlt0] + $depth[$rciAlt1]) / ($strandCountAggregate[$rciAlt0] + $strandCountAggregate[$rciAlt1]) * 10;
          $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, ($depth[$rciAlt0] + $depth[$rciAlt1]), ($strandCountAggregate[$rciAlt0] + $strandCountAggregate[$rciAlt1])));
          $GT = "1/2"; $AD = "$strandCountAggregate[$Rci{$refBase}],$strandCountAggregate[$rciAlt0],$strandCountAggregate[$rciAlt1]";
          $DPF = $DP; $GQ = int($QD+0.499);
          $refY = $YY = 0;
          $PL = sprintf("%d,%d,%d,%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase[0]"}]), phred($refY), phred($probGT[$GTci{"$altBase[0]$altBase[0]"}]), phred($probGT[$GTci{"$altBase[0]$altBase[1]"}]), phred($YY));
          if($altBase[0]=~/[ACGT]/ && $altGT!~/^XY$/)
          {
            $AF = 0.500;
            if($indelMode=~/D$/)
            {
              $refPattern = "$preRefBase$indelPattern";
              $altPattern = "$preRefBase,$altBase[0]";
            }
            elsif($indelMode=~/I$/)
            {
              $refPattern = "$preRefBase";
              $altPattern = "$preRefBase$indelPattern,$altBase[0]";
            }
            else
            { die "$.: Unknown Indel Mode: $indelMode\n"; }
          }
          elsif($altBase[0]!~/[ACGT]/ && $altGT=~/^XY$/)
          {
            if($indelMode=~/D$/ && $indelMode2=~/D$/)
            {
              if(length($indelPattern) > length($indelPattern2))
              {
                $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
                $refPattern = "$preRefBase$indelPattern";
                die "$.: Unable to substract $indelPattern2 from $indelPattern.\n" unless(($tmpIndelPattern = $indelPattern) =~ s/$indelPattern2//);
                $altPattern = "$preRefBase,$preRefBase$tmpIndelPattern";
              }
              elsif(length($indelPattern2) > length($indelPattern))
              {
                $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
                $refPattern = "$preRefBase$indelPattern2";
                die "$.: Unable to substract $indelPattern from $indelPattern2.\n" unless(($tmpIndelPattern2 = $indelPattern2) =~ s/$indelPattern//);
                $altPattern = "$preRefBase,$preRefBase$tmpIndelPattern2";
              }
              else
              { die "$.: Equal length of deletion: $indelMode $indelPattern $indelMode2 $indelPattern2\n"; }
            }
            elsif($indelMode=~/I$/ && $indelMode2=~/I$/)
            {
              $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
              $refPattern = "$preRefBase";
              $altPattern = "$preRefBase$indelPattern,$preRefBase$indelPattern2";
            }
            elsif($indelMode=~/I$/ && $indelMode2=~/D$/)
            {
              $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
              $refPattern = "$preRefBase$indelPattern2";
              $altPattern = "$preRefBase,$preRefBase$indelPattern2$indelPattern";
            }
            elsif($indelMode=~/D$/ && $indelMode2=~/I$/)
            {
              $AF = sprintf("%.3f,%.3f", 0.500, 0.500);
              $refPattern = "$preRefBase$indelPattern";
              $altPattern = "$preRefBase,$preRefBase$indelPattern$indelPattern2";
            }
            else
            { die "$.: Unknown Indel Mode: $indelMode\n"; }
          }
          else
          { die "$.: Unknown alternative genotype: $altGT\n"; }

          $INFO = sprintf("INDEL;AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BQFisherPhredAlt2=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;SBFisherPhredAlt2=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $BaseQFisherAlt2, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $FSAlt2, $QD, $ReadPosBias, $DepthBalance, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
          $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

          if ((substr($ref{$chr}, ($pos-1), length($refPattern))) ne $refPattern)
          {
            next;
          }

          printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refPattern, $altPattern, $qual, $filter);
          printf ("%s\t", $INFO);
          printf ("GT:AD:DP:GQ:PL\t");
          printf ("%s\n", $FORMAT);
        }
        elsif($altGT=~/^[ACGT]X$/)
        {
          $rciAlt = $Rci{$altBase[1]};
          $rciRef = $Rci{$refBase};
          $DP -= $strandCountAggregate[$rciRef];
          $AC = "1";
          $AF = 0.500;
          $FSAlt = phred($strandBias[$rciAlt]);
          $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
          $QD = ($depth[$rciAlt]) / ($strandCountAggregate[$rciAlt]) * 10;
          $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, $depth[$rciAlt], ($strandCountAggregate[$rciAlt])));
          $GT = "0/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
          $DPF = $DP; $GQ = int($QD+0.499);
          $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase[1]"}]), phred($probGT[$GTci{"$altBase[1]$altBase[1]"}]));
          $INFO = sprintf("INDEL;AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $QD, $ReadPosBias, $DepthBalance, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
          if($indelMode=~/I$/)
          {
            $refPattern = "$preRefBase";
            $altPattern = "$preRefBase$indelPattern";
          }
          elsif($indelMode=~/D$/)
          {
            $refPattern = "$preRefBase$indelPattern";
            $altPattern = "$preRefBase";
          }
          else
          { die "$.: Unknown Indel Mode: $indelMode\n"; }
          $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

          if ((substr($ref{$chr}, ($pos-1), length($refPattern))) ne $refPattern)
          {
            goto NEXTALLELE;
          }

          printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refPattern, $altPattern, $qual, $filter);
          printf ("%s\t", $INFO);
          printf ("GT:AD:DP:GQ:PL\t");
          printf ("%s\n", $FORMAT);

          NEXTALLELE:
          $rciAlt = $Rci{$altBase[0]};
          $rciRef = $Rci{$refBase};
          $AC = "1";
          $AF = 0.500;
          $BaseQFisherAlt = phred($baseQualityBias[$rciAlt]);
          $FSAlt = phred($strandBias[$rciAlt]);
          $QD = ($depth[$rciAlt]) / ($strandCountAggregate[$rciAlt]) * 10;
          $qual = ($subOptProb == 0)?($qualMax):(qual($subOptProb, $optProb, $depth[$rciAlt], ($strandCountAggregate[$rciAlt])));
          $GT = "0/1"; $AD = "$strandCountAggregate[$rciRef],$strandCountAggregate[$rciAlt]";
          $DPF = $DP; $GQ = int($QD+0.499);
          $PL = sprintf("%d,%d,%d", phred($probGT[$GTci{"$refBase$refBase"}]), phred($probGT[$GTci{"$refBase$altBase[0]"}]), phred($probGT[$GTci{"$altBase[0]$altBase[0]"}]));
          $INFO = sprintf("AC=%s;AF=%s;AN=%d;FRAC=1.000;BQFisherPhredAlt=%d;BestGT=%s;DP=%d;DP4=%s;DP8=%s;IndelMeta=%s;SBFisherPhredAlt=%d;QD=%.3f;ReadPosBias=%.5f;DepthBias=%.5f;LeftGCPercent=%d;RightGCPercent=%d;LeftQuality=%d;RightQuality=%d;OptSuboptRatio=%d;LeftG=%d;RightG=%d;HRun=%d;IHDP=%d;ILDP=%d", $AC, $AF, $AN, $BaseQFisherAlt, $altGT, $DP, join(",", @strandCountAggregate[0..3]), join(",", @strandCount[0..7]), join(",", $indelMode, $indelPattern, $indelHQCount, $indelLQCount, $indelMode2, $indelPattern2, $indelHQCount2, $indelLQCount2), $FSAlt, $QD, $ReadPosBias, $DepthBalance, $LeftGCPercent, $RightGCPercent, $LeftQuality, $RightQuality, $OptSuboptRatio, $LeftGBaseCount, $RightGBaseCount, $polyRunCount, $totalHQCount, $totalLQCount);
          $FORMAT = sprintf("%s:%s:%d:%d:%s", $GT, $AD, $DPF, $GQ, $PL);

          if ((substr($ref{$chr}, ($pos-1), length($refBase))) ne $refBase)
          {
            next;
          }

          printf ("%s\t%d\t.\t%s\t%s\t%.1f\t%s\t", $chr, $pos, $refBase, $altBase[0], $qual, $filter);
          printf ("%s\t", $INFO);
          printf ("GT:AD:DP:GQ:PL\t");
          printf ("%s\n", $FORMAT);
        }
        else
        { die "$.: Problematic alternative GT: $altGT\n"; }
      }
      else { die "$.: \@altBase != 1 and != 2 in Het.\n"; }
    }
  }
  else
  { die "$.: Not SNP and Indel for entry: $_"; }
}

0;

sub accumulate {
  $rt = 0;
  $rt += $_ foreach (@_);
  return $rt;
}

__END__
__C__
#include <float.h>
#include <math.h>
#include <stdio.h>

double clog(double n)
{
  return log(n);
}

double cmax(double m, double n)
{
  return m>n?m:n;
}

double cmin(double m, double n)
{
  return m<n?m:n;
}

int phred(double num)
{
  if(num < 0 || num > 1.01)
  {
     fprintf(stderr, "Input for phred < 0 or > 1: %f\n", num);
    exit(EXIT_FAILURE);
  }
  if(num>1.){num=1.;}
  return (int)((num==0.)?(10*DBL_MAX_10_EXP):(-10*log10(num)));
}

int phred2(double num)
{
  if(num < 0)
  {
    fprintf(stderr, "Input for phred2 < 0: %f \n", num);
    exit(EXIT_FAILURE);
  }
  return (int)((num==0.)?(0):(10*log10(num)));
}

double qual(double subOpt, double opt, double qd_nom, double qd_den)
{
  double rt = (-4.343 * log(subOpt/opt)) * (1 - (log2(5.f-qd_nom/qd_den)/log2(5.f))) ;
  if(rt != rt)
  {
    fprintf(stderr, "%e, %e, %e, %f, %f, %f\n", rt, subOpt, opt, qd_nom, qd_den, qd_nom/qd_den);
  }
  return rt;
}

double qualMaxVal()
{
  return (4.343 * log(DBL_MAX));
}

double fltMaxVal()
{
  return DBL_MAX;
}

double SNPHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
  {
    printf("FATAL ERROR - SNP-HWE: Current genotype configuration (%d  %d %d ) includes a negative count", obs_hets, obs_hom1, obs_hom2);
    exit(EXIT_FAILURE);
  }

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;

  double * het_probs = (double *) malloc((size_t) (rare_copies + 1) * sizeof(double));
  if (het_probs == NULL)
  {
    printf("FATAL ERROR - SNP-HWE: Unable to allocate array for heterozygote probabilities" );
    exit(EXIT_FAILURE);
  }
  memset((void*)het_probs, '\0', (size_t) (rare_copies + 1) * sizeof(double));

  int i;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);

  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;

  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
  {
    het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0) / (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
    sum += het_probs[curr_hets - 2];

    /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
    curr_homr++;
    curr_homc++;
  }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
  {
    het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc /((curr_hets + 2.0) * (curr_hets + 1.0));
    sum += het_probs[curr_hets + 2];

    /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
    curr_homr--;
    curr_homc--;
  }

  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
  double p_hi = het_probs[obs_hets];
  for (i = obs_hets + 1; i <= rare_copies; i++)
    p_hi += het_probs[i];

  double p_lo = het_probs[obs_hets];
  for (i = obs_hets - 1; i >= 0; i--)
    p_lo += het_probs[i];

  double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++)
  {
    if (het_probs[i] > het_probs[obs_hets])
    continue;
    p_hwe += het_probs[i];
  }

  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  free(het_probs);

  return p_hwe;
}

