#!/usr/bin/perl
$| = 1;

use strict;
use warnings;

use IO::Socket::INET;

my $directory = 'AllelesMatrix';
my $alleleFile = $directory.'/Illumina.Alleles.Allelelist.v4.txt';
my $matrixFile = $directory.'/Illumina.Alleles.matrix.v4';
my $chipFile = $directory.'/Illumina.Alleles.Chips.v4.txt';

my $socket;
my $client_socket;
my $peer_address;
my $peer_port;
my $port = '5001';
my $host = '127.0.0.1';
my $state = 0;
my $chiplist;
my $bimfile;
my %alleles;

my %snpChipMapping ; #list of numeric chip ids, indexed by snp name
my @chipNames;
my %alleleMapping;

open LOG, ">>ChipAlignmentServer.log" or die $!;

print LOG "Reading data\n";
getData();
print LOG "Read data\n";
print "Server Ready for queries\n";
#set up new socket
$socket = new IO::Socket::INET 
 (
 LocalHost => $host,
 LocalPort => $port,
 Proto     => 'tcp',
 Listen    => 5,
 ReuseAddr => 1
 )
or die "ERROR in Socket Creation : $!\n";
print LOG "Ready, Waiting for connections on $host $port\n";
print "Ready, Waiting for connections on $host $port\n";
close LOG;

while(1)
 {
 $client_socket = $socket->accept();
 
 $peer_address = $client_socket->peerhost();
 $peer_port = $client_socket->peerport();
 #print LOG "Accepted client connection from : $peer_address, $peer_port\n";
 
 #$bimfile = <$client_socket>;
 my $dataline = <$client_socket>;
 if ($dataline)
  {
  #my $summary = getAlignment($bimfile, $chiplist);
  my $summary = getAlignment($dataline);
  print $client_socket "$summary\n";
  }
 }
 
$socket->close();
 
sub getData
 {
 my %allAlleles;
 print LOG "Reading $alleleFile\n"; 
 open IN, "$alleleFile" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  $allAlleles{$temp[1]} = $temp[0];
  }
 close IN;
 print LOG "Done\n";
  
 print LOG "Reading $matrixFile\n";
 open IN, "$matrixFile" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  if ($temp[1] > 0 and $temp[1] != 2)
   {
   if ($. % 1000000 == 0)
    {
    print LOG "$.\n";
    }
   #$alleleMapping{$temp[0]} = $temp[1];
   getAllAlleleCombinations($temp[0], $allAlleles{$temp[1]});
   $snpChipMapping{$temp[0]} = $temp[2];
   }
  }
 close IN;
 print LOG "\nDone\n";
 
 print LOG "Reading $chipFile\n";
 open IN, "$chipFile" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  $chipNames[$temp[1]] = $temp[0];
  }
 close IN;
 print LOG "Done\n";
  
 }

sub getAllAlleleCombinations
 {
 my $snp = $_[0];
 my $alleleString = $_[1];
 
 $alleleString =~ /(\w)(\w)(\d)(\d)/;
 
 my $top1 = $1;
 my $top2 = $2;
 my $illm = $3;
 my $src = $4;
 my $topalleles = $top1.$top2;
 
 my $bot1 = $top1;
 my $bot2 = $top2;
 $bot1 =~ tr/ACGT/TGCA/;
 $bot2 =~ tr/ACGT/TGCA/;
 my $botalleles = $bot2.$bot1;
 
 if ($topalleles eq 'AT' or $topalleles eq 'CG')
  {
  $alleles{$snp}{$topalleles} = '1111';
  $alleles{$snp}{$top1} = '1111';
  $alleles{$snp}{$top2} = '1111';
  }
 else
  {
  $alleles{$snp}{$topalleles} = '10'.$illm.$src;
  $alleles{$snp}{$top1} = '10'.$illm.$src;
  $alleles{$snp}{$top2} = '10'.$illm.$src;
  
  $illm =~ tr/01/10/;
  $src =~ tr/01/10/;
 
  $alleles{$snp}{$botalleles} = '01'.$illm.$src;
  $alleles{$snp}{$bot1} = '01'.$illm.$src;
  $alleles{$snp}{$bot2} = '01'.$illm.$src;
  }
 #print LOG "$snp\t$topalleles\t$top1\t$top2\t$alleles{$snp}{$topalleles}\t$botalleles\t$bot1\t$bot2\t$alleles{$snp}{$botalleles}\n";
 }

sub getAlignment
 {
 my $dataline = $_[0];
 #my $bimfile = $_[0];
 #my $chipstring = $_[1];
 my $summary;
 my @result;
 open LOG, ">>ChipAlignmentServer.log";
 
 #my @chiplist = split(/\t/, $chipstring);
 #print "Reading $bimfile\n";
 #my $startTime = time();
 #my %bimAlleles = getBimAlleles($bimfile); #sorted alleles indexed by snp name
 my %bimAlleles = getBimAlleles($dataline);
 my @bimSNPs = keys %bimAlleles;
 my $bimTotal = @bimSNPs;
 #my $endTime = time();
 #my $runTime = $endTime - $startTime;
 #print "Done $runTime seconds\n";
 #iterate through SNPs on chip
 foreach my $snp (@bimSNPs)
  {  
  if ($alleles{$snp}{$bimAlleles{$snp}})
   {
   my @counts = split(//, $alleles{$snp}{$bimAlleles{$snp}});
   $result[0] += $counts[0]; #TOP
   $result[1] += $counts[1]; #BOT
   $result[2] += $counts[2]; #Illmn
   $result[3] += $counts[3]; #SRC
   $result[4]++;             #total matches for chip
   }
  }
 $summary = createSummary(\@result, $bimTotal);
 close LOG;
 return $summary;
 }

sub createSummary
 {
 my @result = @{$_[0]};
 my $bimTotal = $_[1];
 my $summary;
 
 #for (my $i = 1; $i <= $#chipNames; $i++)
 # {
  if ($result[4])
   {
   my $strand = 'TOP';
   my $strandpct = 0;
   
   my $toppct = $result[0]/$result[4]*100;
   my $botpct = $result[1]/$result[4]*100;
   my $illpct = $result[2]/$result[4]*100;
   my $srcpct = $result[3]/$result[4]*100;
   
   if ($toppct > $strandpct)
    {
    $strand = 'TOP';
    $strandpct = $toppct;
    }
   
   if ($botpct > $strandpct)
    {
    $strand = 'BOT';
    $strandpct = $botpct;
    }
   
   if ($illpct > $strandpct)
    {
    $strand = 'ILMN';
    $strandpct = $illpct;
    }
   
   if ($srcpct > $strandpct)
    {
    $strand = 'SRC';
    $strandpct = $srcpct;
    }
   
   #print LOG "$chipNames[$i]\t$result[$i][0] ($toppct)\t$result[$i][1] ($botpct)\t$result[$i][2] ($illpct)\t$result[$i][3] ($srcpct)\t$result[$i][4]\t$bimTotal\n";
   
   #$summary .= "$chipNames[$i]\t$toppct\t$botpct\t$illpct\t$srcpct|";
   $summary = "$strand\t$strandpct|";
   }
  else
   {
   #$summary .= "$chipNames[$i]\t0\t0\t0\t0|";
   #$summary = "|";
   }
  #}
 return $summary;
 }

sub getBimAlleles
 {
 my $data = $_[0];
 my %alleles;
 my @entries = split(/\|/, $data);
 foreach my $entry (@entries)
  {
  my @temp = split(/\t/, $entry);
  $alleles{$temp[0]} = $temp[1];
  }
 #print "Read $#entries SNPs\n"; 
 return %alleles; 
 }

sub getChipList
 {
 my $file = $_[0];
 my %chiplist;
 open IN, "$file" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  $chiplist{$temp[1]} = $temp[0];
  }
 close IN;
 return %chiplist;
 }
