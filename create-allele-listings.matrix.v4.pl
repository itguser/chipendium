$| = 1;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use strict;
use warnings;

my $outdir = 'AllelesMatrix';
my $directory = './Annotation';

if (!-e $outdir)
 {
 mkdir $outdir;
 }

my @chiplist = get_chip_list($directory);
my $outfile = 'Illumina.Alleles.matrix.v4';
my $outpath = $outdir.'/'.$outfile;
my $snppath = $outdir.'/Illumina.Alleles.Allelelist.v4.txt';

my $chiplistfile = $outdir.'/Illumina.Alleles.Chips.v4.txt';
open CHIPMAP, ">$chiplistfile" or die $!;

my %data;
my %allchips;
my %allSNPs;

my $chipcount = @chiplist;
my $counter = 0;
my $allAllelesCounter = 0;
my %allAllelesList;

foreach my $chip (@chiplist)
 {
 $counter++;
 $chip =~ /(.*)\.csv.*/;
 my $chipname = $1;
 print CHIPMAP "$chipname\t$counter\n";
 my $snpcount = 0;
 my $flag = 0;
 my $chippath = $directory.'/'.$chip;
 print "Processing $chip ($counter of $chipcount)\n";
 my %headers = get_headers($chippath);
 my $z;
 
 if ($chip =~ /.*\.csv$/)
  {
  open $z, "$chippath" or die $!;
  }
 elsif ($chip =~ /.*\.csv.gz$/)
  {
  $z = new IO::Uncompress::Gunzip "$chippath" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
  
 while (<$z>)
  {
  s/\"//g;
   
  if (/^\[Controls\].*/)
   {
   $flag = 0;
   }
  
  if ($flag)
   {
   $snpcount++;
   my $topalleles = 0;
   my $botalleles = 0;
   my $srcalleles = 0;
   my $illmalleles = 0;
   my $alleles = 'NN';
   
   chomp;
   my @temp = split/\,/;
   if ($headers{'SNP'})
    {
    if ($temp[$headers{'SNP'}])
     {
     $temp[$headers{'SNP'}] =~ /^.*?\[(.*?)\/(.*?)\].*?$/;
     $alleles = $1.$2;
     }
    
    if ($headers{'IlmnStrand'} and $headers{'SourceStrand'})
     {
     if (uc $temp[$headers{'IlmnStrand'}] eq 'TOP')
      {
      $topalleles = $alleles;
      $illmalleles = '1';
      }
     elsif (uc $temp[$headers{'IlmnStrand'}] eq 'BOT')
      {
      $alleles =~ tr/ACGT/TGCA/;
      $topalleles = $alleles;
      $illmalleles = '0';
      }
     else
      {
       
      }
       
     if (uc $temp[$headers{'SourceStrand'}] eq 'TOP')
      {
      $srcalleles = '1';
      }
     elsif (uc $temp[$headers{'SourceStrand'}] eq 'BOT')
      {
      $srcalleles = '0'; 
      }
     else
      {
       
      }
      
     } 
    else
     {
     $topalleles = $botalleles = $srcalleles = 0;
     }
    }
   else
    {
    $topalleles = $botalleles = $srcalleles = $illmalleles = 0;
    }
    
   my $allAlleles = $topalleles.$illmalleles.$srcalleles;
   
   if (!$allAllelesList{$allAlleles}) #create a list of all the possible allele combinations and number them
    {
    $allAllelesList{$allAlleles} = $allAllelesCounter;
    $allAllelesCounter++;
    }
   
   if ($data{$temp[$headers{'Name'}]})
    {
    if ($data{$temp[$headers{'Name'}]} eq $allAllelesList{$allAlleles}) 
     {
     if($allSNPs{$temp[$headers{'Name'}]} )
      {
      $allSNPs{$temp[$headers{'Name'}]} .= ",$counter";
      }
     else
      {
      $allSNPs{$temp[$headers{'Name'}]} = $counter;
      }
     }
    else # not the same drop SNP from the comparisons
     {
     $data{$temp[$headers{'Name'}]} = -1;
     }
    }
   else
    {
    $data{$temp[$headers{'Name'}]} = $allAllelesList{$allAlleles};
    $allSNPs{$temp[$headers{'Name'}]} = $counter;
    }
   }
  
  if (/IlmnID/ || /Ilmn ID/) #reads header from Illumina annnotation file
   {
   $flag = 1;
   }
  
  }
 close $z;
 print "Read $snpcount SNPs\n";

 }

my @snps = keys %data;

open OUT, ">$outpath" or die $!;
 
for (my $i = 0; $i <= $#snps; $i++)
 {
 print OUT "$snps[$i]\t$data{$snps[$i]}\t$allSNPs{$snps[$i]}";
 print OUT "\n"; 
 } 

close OUT;

my @allelesList = keys %allAllelesList;
open ALLELE, ">$snppath" or die $!;
foreach my $allele (@allelesList)
 {
 print ALLELE "$allele\t$allAllelesList{$allele}\n";
 }


sub get_headers
 {
 my $annotation = $_[0];
 my %header;
 my $z;
 
 if ($annotation =~ /.*\.csv$/)
  {
  open $z, "$annotation" or die;
  }
 elsif ($annotation =~ /.*\.csv.gz$/)
  {
  $z = new IO::Uncompress::Gunzip "$annotation" or die "IO::Uncompress::Gunzip failed: $GunzipError\n"; 
  }
  
 while (<$z>)
  {
  chomp;
  s/\"//g;
  if (/IlmnID/ || /Ilmn ID/) #reads header from Illumina annnotation file
   {
   my @titles = split/\,/; 
   for (my $i = 0; $i <= $#titles; $i++) 
    {
    $header{$titles[$i]} = $i;
    }
   }
  } 
 close $z;
 #print "IllmnStrand\t$header{'IlmnStrand'}\nSourceStrand\t$header{'SourceStrand'}\nSNP\t$header{'SNP'}\n";
 return %header;
 }

sub get_chip_list
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
   print "Adding file $file to processing list\n";
   }
  elsif ($file =~ /.*\.csv\.gz$/)
   {
   push (@chips, $file);
   print "Adding file $file to processing list\n";
   }
  }
 return @chips;
 }
