#!/usr/bin/perl

use warnings;
use strict;

my $vcfFN = shift @ARGV;
my $bamFN = shift @ARGV;
my $fastaFN = shift @ARGV;
my $readLength = 101;

if(not defined $fastaFN) { print STDERR "perl $0 <vcfFN> <bamFN> <fastaFN>\n"; exit; }
if(!-e $vcfFN) { print STDERR "perl $0 <vcfFN> <bamFN> <fastaFN>\n"; exit; }
if(!-e $bamFN) { print STDERR "perl $0 <vcfFN> <bamFN> <fastaFN>\n"; exit; }
if(!-e $fastaFN) { print STDERR "perl $0 <vcfFN> <bamFN> <fastaFN>\n"; exit; }

open my $vcfFH, "$vcfFN" or die $!;

while(my $l = <$vcfFH>)
{
  next if $l=~/^#/;
  chomp $l;
  my @F = split /\t/, $l;
  my @alts = split /,/, $F[4];
  foreach my $alt (@alts)
  {
    my $ref = $F[3];
    my $cigarMarker = "";
    my $posCount = 0; my $negCount = 0; my $snpMingle = 0; my $indelMingle = 0;
    my $len = 0;
    if(length($ref) > length($alt)) #deletion
    {
      $len = length($ref) - length($alt);
      $cigarMarker = "${len}D";
    }
    elsif(length($ref) < length($alt)) #insertion
    {
      $len = length($alt) - length($ref);
      $cigarMarker = "${len}I";
    }
    else
    { print STDERR "Should never reach here. Ref/Alt: $ref/$alt\n";}

    my $cord1 = $F[1] - $readLength; my $cord2 = $F[1] + $readLength;
    my $range = "$F[0]:$cord1-$cord2";
    my @alignments = `samtools view $bamFN $range`;
    if(scalar(@alignments) == 0)
    { print STDERR "Command error samtools view $bamFN $range\n"; exit;}
    my @fastas = `samtools faidx $fastaFN $range`;
    if(scalar(@fastas) <= 1)
    { print STDERR "Command error samtools tabix $fastaFN $range\n"; exit;}
    chomp(@fastas);
    my @fasta = ();
    for(my $i = 1; $i < @fastas; ++$i) { push @fasta, (split //, $fastas[$i]); }

    foreach my $aln (@alignments)
    {
      chomp $aln;
      my @A = split /\t/, $aln; print STDERR "aln parseFail\n" if (scalar(@A) == 0);
      if($A[5] =~ /[IDMSH]$cigarMarker/)
      {
        my $genomeOffset = 0;
        my $readOffset = 0;
        my $errFlag = 1;
        my $prevCigar = "";
        #print STDERR "Cigar: $A[5], find: $cigarMarker\n"; ##
        while($A[5] =~ /(\d+)([IDMSH])/g)
        {
          #print STDERR "Parsed part: $1$2\n"; ##
          if("$1$2" eq "$cigarMarker")
          { $errFlag = 0; last; }
          if("$2" eq "M" || "$2" eq "D") { $genomeOffset += $1; $readOffset += $1; }
          elsif("$2" eq "I" || "$2" eq "S") { $readOffset += $1; }
          elsif("$2" eq "H") {}
          else { print STDERR "Should never reach here! Parse Part.\n"; exit; }
          $prevCigar .= "$1$2";
        }
        if($errFlag == 1) { print STDERR "Cigar string parsing failure: $aln\n"; exit; }
        if(($A[3] + $genomeOffset - 1) == $F[1])
        {
          if(($A[1] & 16) == 0) { ++$posCount; }
          elsif(($A[1] & 16) != 0) { ++$negCount; }
          else{ print STDERR "Cannot decide strand: $aln\n"; exit;}

          my $leftSearchBp = 10; my $rightSearchBp = 10;
          if($prevCigar =~ /[ID](\d+)M/)
          {
            if($1 <= 10)
            {
              ++$indelMingle;
              $leftSearchBp = $1;
              #print STDERR "Indel Prev \t$aln\n"; ##
            }
          }
          my $postCigar = $A[5];
          unless($postCigar =~ s/^$prevCigar$cigarMarker//)
          {print STDERR "Cannot remove $prevCigar$cigarMarker from $A[5]\n"; exit;}
          if($postCigar =~ /(\d+)M\d+[ID]/)
          {
            if($1 <= 10)
            {
              ++$indelMingle;
              $rightSearchBp = $1;
              #print STDERR "Indel Post \t$aln\n"; ##
            }
          }

          my $fastaOffset = $A[3] - $cord1 + $genomeOffset;
          my @readBases = split //, $A[9];
          my $mismatchFlag = 0;
          for(my $i = 1; $i < $leftSearchBp; ++ $i)
          {
            last if ($readOffset-$i < 0);
            last if ($fastaOffset-$i < 0);
            last if ($readOffset-$i >= @readBases);
            if("$fasta[$fastaOffset-$i]" ne "$readBases[$readOffset-$i]")
            {
              $mismatchFlag = 1;
              #print STDERR "SNP Prev ".($readOffset-$i)."\t$aln\n"; ##
            }
          }
          if($cigarMarker =~ /D$/)
          {
            for(my $i = 1; $i < $rightSearchBp; ++ $i)
            {
              last if (($readOffset+$i) >= @readBases);
              last if (($fastaOffset+$len+$i) >= @fasta);
              if("$fasta[$fastaOffset+$len+$i]" ne "$readBases[$readOffset+$i]")
              {
                $mismatchFlag = 1;
                #print STDERR "SNP Post D".($readOffset+$i)."\t$aln\n"; ##
              }
            }
          }
          elsif($cigarMarker =~ /I$/)
          {
            for(my $i = 1; $i < $rightSearchBp; ++ $i)
            {
              last if (($readOffset+$len+$i) >= @readBases);
              last if (($fastaOffset+$i) >= @fasta);
              if("$fasta[$fastaOffset+$i]" ne "$readBases[$readOffset+$len+$i]")
              {
                $mismatchFlag = 1;
                #print STDERR "SNP Post I ".($readOffset+$len+$i)."\t$aln\n"; ##
              }
            }
          }
          else { print STDERR "Should never reach here!\n"; exit; }
          if($mismatchFlag == 1) { ++$snpMingle; }

        }
      }
    }

    print "$F[0]\t$F[1]\t$posCount\t$negCount\t$snpMingle\t$indelMingle\n";
  }
}

