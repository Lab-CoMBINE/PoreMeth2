#!/usr/bin/perl
#############################################
#Author: Alberto magi
#email: albertomagi@gmail.com
#Created Time: 08 May 2017   
#University of Florence
#############################################
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use IO::File;

my ($Chr,$readname,$CpGPos,$PathOut,$Label,$CpGGroupCount,$MinReads,$beta,@CpGVec,@ChrVec,$fout1,@data,@MethSupp,$MethLike,%ReadHash,%CovHash,%MethHash,$MethStatus,$lenArr,@counter,@keyReadHash,@arrSize,$string2compare,@index);
#$SamIn = $ARGV[0];

my @stringRef=('000','001','011','010','100','101','111','110');
# Open the SAM file
#open (SAM, "$SamIn") || die "Major problem: cannot open $SamIn for reading: $!";

my $CpGPosOld=0;
my $NumCpG=3;
$MinReads=2;
my $ChrOld="chrS";
my $lcount=1;
while (<STDIN>){
  if (! /^@/){
    @data = split /\t/, $_;
    $readname = $data[3];
    $MethLike = $data[4];
    $CpGPos = $data[1];
    $Chr = $data[0];
    $CpGGroupCount = $data[7];
    if ($MethLike> .75){
      $MethStatus=1;
    }
    if ($MethLike< .75){
      $MethStatus=0;
    }
    if ($Chr ne $ChrOld){
      $lcount=1;
      %ReadHash = ();
      %MethHash = ();
      %CovHash = ();
      @ChrVec = ();
      @CpGVec = ();

    }
    if ($CpGPos!=$CpGPosOld){
      if ($lcount>$NumCpG){
        @counter = (0) x 8;
        my @keyReadHash = keys %ReadHash;
        for my $keyReadHash (@keyReadHash){
          my @subkeyS = keys %{$ReadHash{$keyReadHash}};
          my $arrSize = @subkeyS;
          #print "$arrSize\n";
          if ($arrSize==$NumCpG){
            @MethSupp=();
            foreach my $subkeyS (@subkeyS) {
              push @MethSupp, @{$ReadHash{$keyReadHash}{$subkeyS}};
              #print "$ReadHash{$keyReadHash}{$subkeyS}\n";
            }
            $string2compare = join("", @MethSupp);
            #print "$string2compare\n";
            @index = grep { $stringRef[$_] eq $string2compare } 0..$#stringRef;
            my $ind_size = @index;
            if ($ind_size!=0)
            {
              #print "@stringRef\n";
              $counter[$index[0]] += 1;
            }
          }
          #print "@CpGVec\n";
          delete($ReadHash{$keyReadHash}{$CpGVec[0]});
          @subkeyS = keys %{$ReadHash{$keyReadHash}};
          $arrSize = @subkeyS;
          if ($arrSize==0){
            delete($ReadHash{$keyReadHash});
          }
        }
        my $sumcounter = 0;
        foreach my $var (@counter) {
          $sumcounter = $sumcounter + $var;
        }
        if ($sumcounter>$MinReads && $CovHash{$CpGVec[0]}>$MinReads){
          my $counterNorm = 0;
          foreach my $var (@counter) {
            if ($var!=0){
              $counterNorm = $counterNorm + ($var/$sumcounter)*(log($var/$sumcounter)/log(2));
            }
          }
          $counterNorm= - $counterNorm * (1/($NumCpG));
          $beta=$MethHash{$CpGVec[0]}/$CovHash{$CpGVec[0]};
          print "$ChrVec[0]";
          print "\t";
          print "$CpGVec[0]";
          print "\t";
          print "$counterNorm";
          print "\t";
          print "$sumcounter";
          print "\t";
          print "$beta";
          print "\t";
          print "$CovHash{$CpGVec[0]}\n"
        }
        delete($MethHash{$CpGVec[0]});
        delete($CovHash{$CpGVec[0]});
        shift @CpGVec;
        shift @ChrVec;
      }
      push @ChrVec, $Chr;
      push @CpGVec, $CpGPos;
      $lcount += 1;
      #print "$lcount\n";
    }
    @{$ReadHash{$readname}{$CpGPos}}=$MethStatus;
    $MethHash{$CpGPos} += $MethStatus;
    $CovHash{$CpGPos} += 1;
    $CpGPosOld=$CpGPos;
    $ChrOld=$Chr;  
    }
}

