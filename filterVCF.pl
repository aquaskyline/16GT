#!/usr/bin/perl

use strict;
use warnings;

my $defMinQual = 10;
my $defMinFracSNP = 0.3;
my $defMinFracINDEL = 0.15;
my $defMaxSORSNP = 3;
my $defMaxFSINDEL = 30;
my $defMinDP = 3;
my $defMaxDPMultiplier = 5;
my $defMinStrandCount = 1;
my $defMinQDSNP = 24;
my $defMinQDINDEL = 20;

if(not defined $ARGV[1])
{
  print STDERR "perl $0 <16GT VCF Input> <dbSNP VCF Input> [Min QUAL, $defMinQual] [Min FRAC SNP, $defMinFracSNP] [Min FRAC INDEL, $defMinFracINDEL] [Max SOR SNP, $defMaxSORSNP] [Max FS INDEL, $defMaxFSINDEL] [Min DP, $defMinDP] [Max DP Multiplier, meanDP+FLOAT*sqrt(meanDP), $defMaxDPMultiplier] [Min Strand Count, $defMinStrandCount] [Min QD SNP, $defMinQDSNP] [Min QD INDEL, $defMinQDINDEL]\n"; exit;
}

my $fh;
unless(-e $ARGV[0]) { die "Failed to open $ARGV[0], exiting...\n"; }
unless(-e $ARGV[1]) { die "Failed to open $ARGV[1], exiting...\n"; }

my %dbSNP = ();
open $fh, "gzip -dcf $ARGV[1] |" or die $!;
print STDERR "Loading dbSNP $ARGV[1]...\n";
while(my $l = <$fh>)
{
  next if $l =~ /^#/;
  my @F = split /\t/, $l;
  my @alts = split /,/, $F[4];
  foreach my $alt (@alts)
  {
    my $key = "$F[0]-$F[1]-$F[3]-$alt";
    $dbSNP{$key} = $F[2];
  }
}
close $fh;
print STDERR "Loaded ".scalar(keys %dbSNP)." dbSNP entries from $ARGV[1]\n";

my $isListInput = 0;
my @files2bFiltered = ();
if($ARGV[0]=~/\.list$/)
{
  $isListInput = 1;
  open $fh, $ARGV[0] or die $!;
  @files2bFiltered = <$fh>; chomp @files2bFiltered;
  close $fh;
}
elsif($ARGV[0]=~/\.vcf$/ || $ARGV[0]=~/.vcf.gz/)
{ push @files2bFiltered, $ARGV[0]; }
else
{ print STDERR "Input is neither a VCF file ended with .vcf or .vcf.gz, nor a list file ended with .list\n"; }

foreach my $fn (@files2bFiltered)
{
  my $ofh;
  if($isListInput == 1)
  {
    my $ofn = $fn;
    if($ofn=~s/\.vcf$//) {}
    elsif($ofn=~s/\.vcf.gz$//) {}
    $ofn .= ".filtered.vcf";
    open $ofh, ">$ofn" or die "$!";
    select $ofh;
  }

  open $fh, "gzip -dcf $fn |" or die $!;
  print STDERR "Loading VCF $fn...\n";
  my @lines = <$fh>;
  close $fh;
  print STDERR "Loaded ".scalar(@lines)." lines from $fn\n";

  my $meanDP = 0;
  {
    print STDERR "First pass (determine mean depth) ...\n";
    my @dpCount = ();
    foreach my $l (@lines)
    {
      next if ($l=~/^#/);
      if($l=~/DP=(\d+?);/) { ++$dpCount[$1]; }
    }
    my $max = 0;
    for(my $i=0; $i<@dpCount; ++$i)
    {
      if(defined $dpCount[$i])
      {
        if($dpCount[$i] >= $max)
        {
          $meanDP = $i;
          $max = $dpCount[$i];
        }
      }
    }
    print STDERR "Mean DP: $meanDP\n";
  }

  my $minQual = defined $ARGV[2] ? $ARGV[2] : $defMinQual;
  my $minFracSNP = defined $ARGV[3] ? $ARGV[3] : $defMinFracSNP;
  my $minFracINDEL = defined $ARGV[4] ? $ARGV[4] : $defMinFracINDEL;
  my $maxSORSNP = defined $ARGV[5] ? $ARGV[5] : $defMaxSORSNP;
  my $maxFSINDEL = defined $ARGV[6] ? $ARGV[6] : $defMaxFSINDEL;
  my $minDP = defined $ARGV[7] ? $ARGV[7] : $defMinDP;
  my $maxDPMultiplier = defined $ARGV[8] ? $ARGV[8] : $defMaxDPMultiplier;
  my $maxDP = $meanDP + $maxDPMultiplier * sqrt($meanDP);
  $maxDP = int($maxDP + 0.51);
  my $minStrandCount = defined $ARGV[9] ? $ARGV[9] : $defMinStrandCount;
  my $minQDSNP = defined $ARGV[10] ? $ARGV[10] : $defMinQDSNP;
  my $minQDINDEL = defined $ARGV[11] ? $ARGV[11] : $defMinQDINDEL;

  print STDERR "Parameters:\n";
  print STDERR "Min Qual: $minQual\n";
  print STDERR "Min FRAC SNP: $minFracSNP\n";
  print STDERR "Min FRAC INDEL: $minFracINDEL\n";
  print STDERR "Max SOR SNP: $maxSORSNP\n";
  print STDERR "Max FS INDEL: $maxFSINDEL\n";
  print STDERR "Min DP: $minDP\n";
  print STDERR "Max DP Multiplier: $maxDPMultiplier\n";
  print STDERR "Max DP: $maxDP\n";
  print STDERR "Min Strand Count: $minStrandCount\n";
  print STDERR "Min QD SNP: $minQDSNP\n";
  print STDERR "Min QD INDEL: $minQDINDEL\n";

  print STDERR "Second pass (add filter tags) ...\n";
  my $prev = 0;
  my %typeSummary = ();
  foreach my $l (@lines)
  {
    if($l=~/^#/)
    {
      if($l=~/^#CHROM/)
      {
        print "##FILTER=<ID=QUALlow,Description=\"Low Quality: QUAL<$minQual\">\n";
        print "##FILTER=<ID=DPhigh,Description=\"High read depth: DP>$maxDP\">\n";
        print "##FILTER=<ID=DPlow,Description=\"Low read depth: DP<$minDP\">\n";
        print "##FILTER=<ID=SORhigh,Description=\"Large SOR at SNPs: SOR>$maxSORSNP\">\n";
        print "##FILTER=<ID=FShigh,Description=\"Large Fisher-Strand bias at INDELs: SBFisherPhredAlt>$maxFSINDEL\">\n";
        print "##FILTER=<ID=FRAClow,Description=\"Low fraction of non-reference reads: FRAC<$minFracSNP at SNPs or FRAC<$minFracINDEL at INDELs\">\n";
        print "##FILTER=<ID=SClow,Description=\"Low double-strand support at SNPs: Lowest Strand Count Support of Alt<$minStrandCount\">\n";
        print "##FILTER=<ID=QDlow,Description=\"Low QD: QD<$minQDSNP at SNPs or QD<$minQDINDEL at INDELs\">\n";
        print $l;
      }
      else
      {
        print $l;
      }
    }
    else{
      my @a = split /\s+/, $l;
      my @filterTag = ();
      if($a[5] < $minQual)
      { push @filterTag, "QUALlow"; }
      if($a[7] =~ /DP=(\d+?);/)
      {
        if($1 > $maxDP) { push @filterTag, "DPhigh"; }
        if($1 < $minDP) { push @filterTag, "DPlow"; }
      } else {die "no DP tag: $a[7]";}
      if($a[7] =~ /INDEL;/)
      {
        if($a[7] =~ /SBFisherPhredAlt=(\d+?);/)
        {
          if($1 > $maxFSINDEL) { push @filterTag, "FShigh"; }
        } #else { die "no SBFisherPhredAlt tag: $a[7]"; }
        if($a[7] =~ /FRAC=([-+]?[0-9]*\.?[0-9]+);/)
        {
          if($1 < $minFracINDEL) { push @filterTag, "FRAClow"; }
        } else { die "no FRAC tag: $a[7]"; }
        if($a[7] =~ /QD=([-+]?[0-9]*\.?[0-9]+);/)
        {
          if($1 < $minQDINDEL) { push @filterTag, "QDlow"; }
        } else { die "no QD tag: $a[7]"; }
      }
      else
      {
        if($a[7] =~ /SOR=([-+]?[0-9]*\.?[0-9]+);/)
        {
          if($1 > $maxSORSNP) { push @filterTag, "SORhigh"; }
        }
        elsif($a[7] =~ /SBFisherPhredAlt=(\d+?);/)
        {
          if($1 > $maxFSINDEL) { push @filterTag, "FShigh"; }
        } #else { die "no SOR and SBFisherPhredAlt tag: $a[7]"; }
        if($a[7] =~ /FRAC=([-+]?[0-9]*\.?[0-9]+);/)
        {
          if($1 < $minFracSNP) { push @filterTag, "FRAClow"; }
        } else { die "no FRAC tag: $a[7]"; }
        if($a[7] =~ /QD=([-+]?[0-9]*\.?[0-9]+);/)
        {
          if($1 < $minQDSNP) { push @filterTag, "QDlow"; }
        } else { die "no QD tag: $a[7]"; }
        if($a[7] =~ /DP8=(\S+?);/)
        {
          my @sc = split /,/, $1;
          my $flag = 0;
          foreach my $scidx (split ",", $a[4])
          {
            $scidx =~ tr/ACGT/0246/;
            if($sc[$scidx] < $minStrandCount || $sc[$scidx+1] < $minStrandCount) { $flag = 1; }
          }
          if($flag == 1) { push @filterTag, "SClow"; }
        } else { die "no DP8 tag: $a[7]"; }
      }

      my @alts = split /,/, $a[4];
      my %rsID = ();
      foreach my $alt (@alts)
      {
        my $key = "$a[0]-$a[1]-$a[3]-$alt";
        if(defined $dbSNP{$key})
        {
          $rsID{$dbSNP{$key}} = 1;
        }
      }
      if(scalar(keys %rsID) > 0)
      {
        $a[2] = join ",", (keys %rsID);
      }

      if(scalar(@filterTag) == 0)
      { print join "\t", @a; print "\n"; $typeSummary{"."}++; }
      elsif($a[2] ne ".")
      { print join "\t", @a; print "\n"; $typeSummary{"dbSNP"}++; }
      else
      { $a[6] = join ",", @filterTag; print join "\t", @a; print "\n"; $typeSummary{"$a[6]"}++; }

      $prev = $a[1];
    }
  }

#foreach(sort keys %typeSummary)
#{ print STDERR "$_: $typeSummary{$_}\n"};

  print STDERR "Filtering $fn done.\n";
  if($isListInput == 1)
  {
    select STDOUT;
  }
}

0;


