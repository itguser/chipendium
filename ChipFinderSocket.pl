#!/usr/bin/perl

use strict;
use warnings;
use IO::Socket::INET;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use Socket qw( :crlf );

$| = 1;
#$/ = "\x{00}";

#mapping for numeric alleles
my %alleleMapping = ('1', 'A', '2', 'C', '3', 'G', '4', 'T');

my $outputDir = 'Output';

#server connection information
my $host = '127.0.0.1';
my $finderPort = '5000';
my $alignmentPort = '5001';

#server set up information
my $port = '5050';
my $client_socket;
my $peer_address;
my $peer_port;

#data strings for SNP names and SNP names/alleles
my $finderSNPlist;
my $allelesList;

#global flags on data available
my $snpListOnly = 0; #flag for the bim vs SNP list, default bim file
my $bimFilePresent = 0;
my $logflag = 0; #flag to if any log messages printed

#set up new socket
my $socket = new IO::Socket::INET 
 (
 LocalHost => $host,
 LocalPort => $port,
 Proto     => 'tcp',
 Listen    => 5,
 ReuseAddr => 1
 )
or die "ERROR in Socket Creation : $!\n";
print "Ready, Waiting for connections on $host $port\n";

while (1)
 {
 $client_socket = $socket->accept();
  
 $peer_address = $client_socket->peerhost();
 $peer_port = $client_socket->peerport();
  
 my $file = <$client_socket>;
 if ($file)
  {
  chomp $file;
  my $summary = processData($file);
  #$client_socket->send($summary);
  $client_socket->autoflush(0);
  print $client_socket "$summary\n";
  $client_socket->flush();
  $finderSNPlist = '';
  $allelesList = '';
  }
 }

###############################################################################

sub processData
 {
 open LOG, ">>chipFinder.log";

 my $file = $_[0];
 $file =~ /.*\/(.*)$/;
 my %summaryData;

 $summaryData{'FILENAME'} = $1;
 $summaryData{'STATUS'} = 'OKAY';
 $summaryData{'MSG'} = 'OKAY';
 $summaryData{'TYPE'} = 'TXT';
 $summaryData{'MARKER_COUNT'} = '0';
 $summaryData{'MOST_LIKELY_COUNT'} = '0';
 $summaryData{'STRAND_MATCH'} = ' ';
 $summaryData{'STRAND_PROBABILITY'} = '0';
 $summaryData{'MOST_LIKELY'} = '';
 $summaryData{'CHIP_MATCHES'} = '';
  
 print LOG "$file\n";
 
 #validate file (can be bim for SNP list only, SNP list can't be aligned to strand)
 my $fileCheckFlag = checkFile($file);

 if ($fileCheckFlag)
  {
  print LOG "Found $fileCheckFlag SNPs\n";
  $summaryData{'MARKER_COUNT'} = $fileCheckFlag + 1;
  
  my $finderSocket = openSocket($host, $finderPort);
  print $finderSocket "$finderSNPlist\n";
  my $summary = <$finderSocket>;
  $finderSocket->close();
  $summary =~ s/\|/\n/g;
  chomp $summary;
  my @lines = split(/\n/, $summary);
  #$summary is "$mapping{$chip}\t$pct\t$count[$chip]\t$bimcount\t$chipCounts{$mapping{$chip}}\t$pct2\t$pct3\n";

  foreach my $l (@lines)
   {
   print LOG "$l\n";
   my @data = split(/\t/, $l);
   $summaryData{'CHIP_MATCHES'} .= "$data[0]\t$data[2]\t$data[3]\t$data[4]\n";
   }

  if ($summary =~ /^No Chips Found/)
   {
   $summaryData{'STATUS'} = 'ERROR';
   $summaryData{'MSG'} = $summary;
   $summaryData{'MOST_LIKELY_COUNT'} = 0;
   $summaryData{'MOST_LIKELY'} = '';
   $summaryData{'CHIP_MATCHES'} = '';
   }
  else
   {
   #generate list of the chips
   my @listOfMostLikelyChips = mostLikelyChips($summary);
   my $countOfMostLikelyChips = @listOfMostLikelyChips;
   $summaryData{'MOST_LIKELY_COUNT'} = $countOfMostLikelyChips;
   foreach my $mostLikely (@listOfMostLikelyChips)
    {
    $summaryData{'MOST_LIKELY'} .= "$mostLikely\n";
    }
   my %mostLikelyChips;
   $mostLikelyChips{$_}++ for (@listOfMostLikelyChips);
 
   #if (!$snpListOnly) #flags bim/SNP list
   if ($bimFilePresent) #flag bim file present
    {
    $summaryData{'TYPE'} = 'BIM';
    #list of all chips to report strand alignment
    my %chips = getChips($summary);
       
    my $alignSocket = openSocket($host, $alignmentPort);
    print $alignSocket "$allelesList\n";
    my $alignment = <$alignSocket>;
    $alignment =~ s/\|//;
    chomp $alignment;
    my @entries = split(/\t/, $alignment);
    $summaryData{'STRAND_MATCH'} = $entries[0];
    $summaryData{'STRAND_PROBABILITY'} = $entries[1];
    }
   else
    {
    $summaryData{'MSG'} = 'SNP list uploaded, cannot check strand alignment';
    #print "\nSNP list only uploaded cannot check strand\n";
    }
   } 
  }
 else
  {
  $summaryData{'STATUS'} = 'ERROR';
  $summaryData{'MSG'} = 'Unrecognised file format';
  }
 print LOG "\nFin\n";
 
 my $outputSummary = createSummary(\%summaryData);
 close LOG;
 return $outputSummary;
 }
 
sub createSummary
 {
 my %data = %{$_[0]};
 
 my $summary = "FILENAME\t$data{'FILENAME'}\nSTATUS\t$data{'STATUS'}\nMSG\t$data{'MSG'}\nTYPE\t$data{'TYPE'}\nMARKER_COUNT\t$data{'MARKER_COUNT'}\nMOST_LIKELY_COUNT\t$data{'MOST_LIKELY_COUNT'}\nSTRAND_MATCH\t$data{'STRAND_MATCH'}\nSTRAND_PROBABILITY\t$data{'STRAND_PROBABILITY'}\n#MOST_LIKELY\n$data{'MOST_LIKELY'}#CHIP_MATCHES\n$data{'CHIP_MATCHES'}##";
 print LOG "$summary\n"; 
 return $summary;
 }


sub getChips
 {
 my $summary = $_[0];
 my @lines = split(/\n/, $summary);
 my %chips;
 foreach my $line (@lines)
  {
  my @entries = split(/\t/, $line);
  #push (@chips, $entries[0]);
  $chips{$entries[0]}++;
  }
 return %chips; 
 }

sub mostLikelyChips
 {
 my $summary = $_[0];
 my $mostLikely = 'None';
 my $match = 10;
 my $match2 = 10;
 my @chipList;
 my @lines = split(/\n/, $summary);
 #$summary is "$mapping{$chip}\t$pct\t$count[$chip]\t$bimcount\t$chipCounts{$mapping{$chip}}\t$pct2\t$pct3\n";
 
 foreach my $line (@lines)
  {
  my @entries = split(/\t/, $line);
  my $pct  = $entries[1];
  my $pct2 = $entries[5]; 
  my $pct3 = $entries[6];
  
  my $composite = $pct + $pct3;
  if ($composite > $match)
   {
   $match = $composite;
   $mostLikely = $entries[0]; 
   }
  elsif ($composite == $match) 
   {
   $mostLikely = $mostLikely."\t$entries[0]";    
   }
  }   
 @chipList = split(/\t/, $mostLikely);
 return @chipList;
 }

sub outputFile
 {
 my $file = $_[0];
 $file =~ /.*\/(.*)\.(bim|txt)/;
 my $stem = $1;
 my $outfile = $stem.'.ChipCheck.txt';
 return $outfile;
 }

sub openSocket
 {
 my $host = $_[0];
 my $port = $_[1];
 my $socket = new IO::Socket::INET
  (
  PeerHost => $host,
  PeerPort => $port,
  Proto => 'tcp',
  )or die "ERROR in Socket Creation : $!\n";
 #print "Connected to Server $host $port\n";
 return $socket;
 }
 
sub sortAlleles
 {
 my @alleles = @{$_[0]};
 my $sortedAlleles ;
 my @sa = sort {uc($a) cmp uc($b)} @alleles;
 for (my $i = 0; $i <= $#sa; $i++)
  {
  if ($sa[$i] =~ m/A|C|G|T/)
   {
   $sortedAlleles .= $sa[$i];
   }
  elsif ($sa[$i] =~ m/1|2|3|4/)
   {
   $sortedAlleles .= $alleleMapping{$sa[$i]};
   }
  }
 return $sortedAlleles;
 }

sub checkFile
 {
 my $file = $_[0];
 my $check = 0;
 if (-e $file)
  {
  my $count = 0;
  my $zipped = isGzipped($file);
  my $z;
  my $check = 0;
 
  if ($zipped)
   {
   $z = new IO::Uncompress::Gunzip "$file" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
   $check = 1;
   }
  elsif ($zipped == -1)
   {
   }
  else 
   {
   open $z, "$file" or die $!;
   $check = 1;
   }
  
  while (<$z>)
   {
   chomp;
   my @temp = split/\s+/;
   if ($#temp == 5) #bim file
    {
    $bimFilePresent++;
    if (length($temp[1]) < 100)
     {
     if ($temp[4] eq '0' and $temp[5] eq '0')
      {
      }
     elsif ($temp[4] eq 'I' or $temp[5] eq 'I')
      {
      }
     elsif ($temp[4] eq 'D' or $temp[5] eq 'D')
      {
      }
     else
      {
      $finderSNPlist .= " $temp[1]";
      my @alleles = ($temp[4],$temp[5]);
      my $sortedAlleles = sortAlleles(\@alleles);
      if ($sortedAlleles)
       {
       if ($allelesList)
        {
        $allelesList .= "|$temp[1]\t$sortedAlleles";
        }
       else
        {
        $allelesList .= "$temp[1]\t$sortedAlleles";
        }
       $count++;
       }
      }
     }
    else
     {
     print LOG "ERROR: rejecting line $. as entry too long\n";
     $logflag = 1;
     }
    }
   elsif ($#temp == 0) #SNP list
    {
    #$snpListOnly = 1;
    if (length($temp[0]) < 100)
     {
     $finderSNPlist .= " $temp[0]";
     $count++;
     }
    else
     {
     print LOG "ERROR: rejecting line $. as entry too long\n";
     $logflag = 1;
     }
    }
   else #not recognised
    {
    #print LOG "$file has a format not recognised\n";
    $logflag = 1;
    }
   }
  close $z;
  print LOG "Read $count SNPs from $file\n"; 
  #print "Read $count SNPs from $file\n";
  if ($check)
   {
   return $count;
   }
  else
   {
   return $check;
   }
  }
 else
  {
  return $check;
  }
 }
 

sub isGzipped
 {
 my $file = $_[0];
 my $gzipped = 0;
 if ($file =~ /.*\.gz$/)
  {
  $gzipped = 1;
  print LOG "$file is gzipped\n";
  }
 elsif ($file =~ /.*\.zip$/)
  {
  $gzipped = -1;
  print LOG "ERROR: $file file is zipped\n";
  }

 return $gzipped;
 }
