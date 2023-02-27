#!/usr/bin/perl
$| = 1;

use strict;
use warnings;
use IO::Socket::INET;

my $snpDir = 'SNPLists';
my $file = $snpDir.'/chipSNPlist.txt';
my $countsFile = $snpDir.'/chipCountlist.txt';
my $mappingFile = $snpDir.'/chipMapping.txt';

my $socket;
my $client_socket;
my $peer_address;
my $peer_port;
my $port = '5000';
my $host = '127.0.0.1';

#read data on SNP/chip 
#print "Reading data on SNPs per chip\n";
my %data = getData($file);
#print "Done\n";

#read data on SNP counts per chip
#print "Reading SNP counts per chip\n";
my %chipCounts = getData($countsFile);
#print "Done\n";

#read mapping
#print "Reading Mapping\n";
my %mapping = getData($mappingFile);
#print "Done\n";

# start server and listen for requests to check
#print "Starting server on $port $host\n";

$socket = new IO::Socket::INET 
 (
 LocalHost => $host,
 LocalPort => $port,
 Proto     => 'tcp',
 Listen    => 5,
 ReuseAddr => 1
 ) or die "ERROR in Socket Creation : $!\n";

print "Ready for queries\n";

while(1)
 {
 $client_socket = $socket->accept();
 
 $peer_address = $client_socket->peerhost();
 $peer_port = $client_socket->peerport();
# print "Accepted client connection from : $peer_address, $peer_port\n";
 
 my $snplist = <$client_socket>;
 #print "$snplist\n";
 if ($snplist)
  {
  my $summary = checkChip($snplist);
  #print "$summary\n";
  #$summary =~ s/\n/\|/g;
  print $client_socket "$summary\n";
  }
 }
 
$socket->close();
 
sub getData
 {
 my $file = $_[0];
 my %data;
 open IN, "$file" or die $!;
 while (<IN>)
  {
  chomp;
  my @temp = split/\t/;
  if ($temp[0])
   {
   $data{$temp[0]} = $temp[1];
   }
  }
 close IN;
 return %data;
 }
 
sub checkChip 
 {
 my $bimcount = 0;
 open LOG, ">>ChipFinderServer.log";
 
 my $snplist = $_[0];
 my @snps = split(/\s+/, $snplist);
 my $summary;
 my @count;
 my $bestMatch = 0;
 my $bestChipMatch = 0;
 my $bestChipName;
 
 $bimcount = @snps;
 print LOG "Received $bimcount SNPs\n";
 print "Searching ...";
 
 foreach my $snp (@snps)
  {
  if ($data{$snp})
   {
   my @chips = split(/\s+/, $data{$snp});
   foreach my $chip (@chips)
    {
    $count[$chip]++;
    }
   }
  }
  
 #my @list = keys %count;
 my $total = @count;
 print LOG "\nFound $total possible chips\n";
 my @chips = keys %mapping;
 
 foreach my $chip (@chips)
  {
  
  if ($count[$chip])
   {
   my $pct = sprintf("%0.2f", $count[$chip]/$bimcount*100);
   my $pct2 = sprintf("%0.2f", $bimcount/$chipCounts{$mapping{$chip}}*100);
   my $pct3 = sprintf("%0.2f", $count[$chip]/$chipCounts{$mapping{$chip}}*100);
   print LOG "$chip\t$mapping{$chip}\t$count[$chip]\t$bimcount\t$chipCounts{$mapping{$chip}}\t$pct\t$pct2\t$pct3\n";
   
   if ($pct > $bestMatch)
    {
    $bestMatch = $pct;
    $bestChipMatch = $pct2;
    $bestChipName = $mapping{$chip};
    }
   
   if ($pct > 95 and $pct2 > 50 and $pct2 <= 100)
    { 
    $summary .= "$mapping{$chip}\t$pct\t$count[$chip]\t$bimcount\t$chipCounts{$mapping{$chip}}\t$pct2\t$pct3|";
    } 
   }
  }  
 if (!$summary) 
  {
  $summary = "No Chips Found above threshold, overlap is $bestChipMatch% to best match: $bestChipName\t$bestMatch\t0\t$bimcount\t0\t$bestChipMatch\t0";
  }
 #print "$summary\n"; 
 close LOG;
 return $summary;
 }
