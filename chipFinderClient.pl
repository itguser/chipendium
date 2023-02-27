#!/usr/bin/perl

use strict;
use warnings;
use IO::Socket::INET;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);

$| = 1;

#mapping for numeric alleles
my %alleleMapping = ('1', 'A', '2', 'C', '3', 'G', '4', 'T');

my $outputDir = 'Output';

#server connection information
my $host = '127.0.0.1';
my $finderPort = '5000';
my $alignmentPort = '5001';

#data strings for SNP names and SNP names/alleles
my $finderSNPlist;
my $allelesList;

#global flags on data available
my $snpListOnly = 0; #flag for the bim vs SNP list #removed as inadvertantlety flagging as not bim files
my $bimFilePresent = 0;
my $logflag = 0; #flag to if any log messages printed

open LOG, ">>chipFinder.log" or die $!;
my $file = $ARGV[0];

#validate file (can be bim for SNP list only, SNP list can't be aligned to strand)
my $fileCheckFlag = checkFile($file);

if ($fileCheckFlag)
 {
 my $outfile = outputFile($file);
 $outfile = $outputDir.'/'.$outfile;
 my $finderSocket = openSocket($host, $finderPort);
 
 print $finderSocket "$finderSNPlist\n";
 
 my $summary = <$finderSocket>;
  
 open OUT, ">$outfile" or die $!;
 
 print "\nAll genotype arrays matching > 95%\n";
 print OUT "\nAll genotype arrays matching > 95%\n";
 
 $summary =~ s/\|/\n/g;
 my @lines = split(/\n/, $summary);
 foreach my $l (@lines)
  {
  my @data = split(/\t/, $l);
  print "$data[0]\t$data[1]\%\n";
  print OUT "$data[0]\t$data[1]\%\n";
  }

 print OUT "$summary\n";

 $finderSocket->close();

 if ($summary =~ /^No Chips Found/)
  {
  print "Number of most likely chips: 0\n";
  print OUT "Number of most likely chips: 0\n";
  }
 else
  {
  #generate list of the chips
  my @listOfMostLikelyChips = mostLikelyChips($summary);
  my $countOfMostLikelyChips = @listOfMostLikelyChips;
  
  print "\nMost Likely Genotype Arrays\n";
  print OUT "\nMost Likely Genotype Arrays\n";
  print "This is based on overall match and size of the bim file relative to the annotation file and so may not be the highest match listed above\n"; 
  print OUT "This is based on overall match and size of the bim file relative to the annotation file and so may not be the highest match listed above\n";
  print "\nNumber of most likely arrays: $countOfMostLikelyChips\nThese are:\n\n";
  print OUT "\nNumber of most likely chips: $countOfMostLikelyChips\nThese are:\n\n";
  
  foreach my $mostLikely (@listOfMostLikelyChips)
   {
   print "$mostLikely\n";
   print OUT "$mostLikely\n";
   }
  my %mostLikelyChips;
  $mostLikelyChips{$_}++ for (@listOfMostLikelyChips);
 
  if ($bimFilePresent) #flags bim/SNP list
   {
   
   #list of all chips to report strand alignment
   print "\n\nChecking strand alignment for $file\n";
   print OUT "\n\nChecking strand alignment for $file\n";
   
   my %chips = getChips($summary);
   
   #my @allChips = keys %chips;
   #my $allChipsCount = @allChips;
   #print "Have $allChipsCount Chips\n";
   
   my $alignSocket = openSocket($host, $alignmentPort);
   print $alignSocket "$allelesList\n";
   my $alignment = <$alignSocket>;
   $alignment =~ s/\|//;
   chomp $alignment;
   print "\nData are aligned to $alignment%\n";
   print OUT "\n$alignment\n";
   
   #my @alignments = split(/\|/, $alignment);
   #print "Received Alignments for $#alignments chips\n";
   #foreach my $a (@alignments)
   # {
   # #print "$a\n";
   # my @temp = split(/\t/, $a);
   # if ($chips{$temp[0]})
   #  {
   #  print "$a";
   #  print OUT "$a";
   #  if ($mostLikelyChips{$temp[0]})
   #   {
   #   print "\t*\n";
   #   print OUT "\t*\n";
   #   }
   #  else
   #   {
   #   print "\t\n";
   #   print OUT "\t\n";
   #   }
   #  }
   # }
   }
  else
   {
   print "\nSNP list only uploaded cannot check strand\n";
   }
  } 
 close OUT; 
 }
print "\nFin\n";
close LOG;

###############################################################################

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
  #if ($pct > $match and $pct3 > $match2)
  # {
  # $match = $pct;
  # $match2 = $pct3;
  # $mostLikely = $entries[0];
  # }
  #elsif ($pct > $match and $pct3 == $match2)
  # {
  # $match = $pct;
  # $mostLikely = $entries[0];
  # } 
  #elsif ($pct == $match and $pct3 > $match2) 
  # {
  # $match2 = $pct3;
  # $mostLikely = $entries[0];
  # }
  #elsif ($pct == $match and $pct3 == $match2) 
  # {
  # $mostLikely = $mostLikely."\t$entries[0]"; 
  # }
  }   
 #$summary .= "\n\nMost likely chip(s):\n$mostLikely\n\n";
 #print "$mostLikely\n";
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
 my $count = 0;
 my $zipped = isGzipped($file);
 my $z;
 my $check = 1;
 
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
   if (length($temp[0]) < 80)
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
   print LOG "$file has a format not recognised\n";
   $logflag = 1;
   }
  }
 close $z;
 print LOG "Read $count SNPs from $file\n"; 
 print "Read $count SNPs from $file\n";
 return $check;
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
