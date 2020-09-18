#!/usr/bin/perl -w
#
# extract_Annotation.pl
#
# extract a simplified version of annotation information for a subset of positions (from an Ann.vcf file)
# 
# April 23, 2015

use strict;

# Get the name of the annotated vcf file, a text file with positions (chr, pos) to extract, and the name of an output file from the command line
my ($USAGE) = "$0 <input.ann.vcf> <snpPos.txt> <output.txt>
\t\tinput.ann.vcf = An annotated VCF file, created by a program like snpEff
\t\tsnpPos.txt = A tab-delimited file with the list of positions to extract (need 2 columns: CHR, POS)
\t\t\tNote: THIS FILE MUST BE SORTED BY POSITION!!!
\t\toutput.txt = The name of the output file to create\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}

my ($vcffile, $posfile, $outfile) = @ARGV;

# Open the output file for printing
open (OUT, ">$outfile") || die "\nUnable to open the file $outfile!\n";

# Open the vcf file, and get the header information from the first few lines
open (VCF, $vcffile) || die "\nUnable to open the file $vcffile!\n";

HEADERLOOP:while(<VCF>) {
  chomp $_;
  if ($_ =~ /^\#\#/) {
    if ($_ =~ /ID\=ANN/) {
      $_ =~ s/^[^\"]*\"//g;
      $_ =~ s/Functional\sannotations\:\s{1,}\'//g;
      $_ =~ s/\'\s{1,}\">//g;
      $_ =~ s/\s{1,}//g;
      my @headings = split(/\|/, $_);
      print OUT "CHR\t", "POS\t", "REF\t", "ALT";
      foreach my $ann_heading (@headings) {
	print OUT "\t", $ann_heading;
      }
      print OUT "\n";
    } else {
      next HEADERLOOP;
    }
  } elsif ($_ =~ /^\#CHROM/) {
    last HEADERLOOP;
  }
}

# Open the file with the position information
open (POS, $posfile) || die "\nUnable to open the file $posfile!\n";

# Read through the position file one line at a time (the file has to be sorted!)
# For each position, start looking through the vcf file until that position is found
# Then, record the relevant info, and go to the next position
POSLOOP:while (<POS>) {
  my $posline = $_;
  chomp $posline;
  my ($chr, $pos) = split(/\s{1,}/, $_);
 
 VCFLOOP:while (<VCF>) {
    my $vcfline = $_;
    chomp $vcfline;
    my @info = split(/\t/, $vcfline);
    if (($info[0] == $chr) && ($info[1] == $pos)) {
      print OUT $info[0], "\t", $info[1], "\t", $info[3], "\t", $info[4];
      $info[7] =~ s/^ANN\=//ig;
      my @ann_info = split(/\|/, $info[7]);
      foreach my $field (@ann_info) {
	print OUT "\t", $field;
      }
      print OUT "\n";
      last VCFLOOP;
    } elsif ($info[0] > $chr) {
      next POSLOOP;
    } else {
      next VCFLOOP;
    }
  }
}
close(POS);
close(VCF);

close(OUT);
exit;
