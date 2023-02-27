#!/usr/bin/perl -w

use strict;
use warnings;use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
my $dir = './Annotation'; #path to folder of genotype chip csvs
my $outdir = 'SNPLists';
if (!-e $outdir)
 {
 mkdir $outdir;
 }

my @annotation = get_annotation_files($dir);
my %data;
my $outfile = $outdir.'/chipSNPlist.txt';
my $outcounts = $outdir.'/chipCountlist.txt';
my $mappingfile = $outdir.'/chipMapping.txt';
my $mapnumber = 0;

open O, ">$outcounts" or die $!;
open M, ">$mappingfile" or die $!;

foreach my $file (@annotation)
 {
 print "$file\n";
 $file =~ /(.*)\.csv.*/;
 my $result = createlist($dir, $file, $outfile);
 }
open OUT, ">$outfile" or die $!;

my @snplist = keys %data;
foreach my $snp (@snplist)
 {
 print OUT "$snp\t$data{$snp}\n";
 }

close OUT;

sub createlist
 {
 my $dir = $_[0];
 my $file = $_[1];
 my $outfile = $_[2];
 
 my $z;
 my $dataflag = 0;
 my $namecolumn = 1;
 my $total = 0;
 my $match = 0;
 my $pct = 0;
 my $filepath = $dir.'/'.$file;
 my $count = 0;
 my $chip;
 
 if ($file =~ /(.*)\.csv\.gz/)
  {
  $z = new IO::Uncompress::Gunzip "$filepath" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  $chip = $1;
  }
 elsif ($file =~ /(.*)\.csv/) 
  {
  open $z, "$filepath";
  $chip = $1;
  }
 else
  {
  return -1;  
  }

 $mapnumber++;
 print M "$mapnumber\t$chip\n"; 
 
 while (<$z>)
  {
  if ($dataflag)
   {
   if (/^\[Controls\].*/ or /^\[Gentrain Request\]/) 
    {
    print "Reached Control section, line $.\n";
    $dataflag = 0;
    }
   else
    {
    $count++;
    my @temp = split/\,/;
    #print OUT "$temp[$namecolumn]\n";
    if ($data{$temp[$namecolumn]})
     {
     $data{$temp[$namecolumn]} = $data{$temp[$namecolumn]}." $mapnumber";
     }
    else
     {
     $data{$temp[$namecolumn]} = $mapnumber;
     }
    }
   }
  if (/IlmnID/ || /Ilmn ID/) 
   {
   my @titles = split/\,/; 
   $dataflag = 1;
   for (my $i = 0; $i <= $#titles; $i++) 
    {
    if ($titles[$i] eq 'Name')
     {
     $namecolumn = $i;
     }
    }
   }
  }
 print O "$chip\t$count\n";
 close $z;
 }

sub get_annotation_files
 {
 my $directory = $_[0];
 my @chips;

 opendir(DIR, "$directory") or die "Not a directory $directory";
 my @all_files =  grep !/^\.\.?\z/, readdir DIR;
 
 foreach my $file (@all_files)
  {
  if ($file =~ /.*\.csv$/)
   {
   push (@chips, $file);
   }
  elsif ($file =~ /.*\.csv\.gz$/)
   {
   push (@chips, $file);
   }
  }
 closedir DIR;
 return @chips;
 }
