#!/usr/bin/perl -w
#
# SigSNP_geneScan_GEMMA.pl
#
# Find any genes near SNPs found to be significant in a GWAS
#
# November 17, 2014

use strict;
use Data::Dumper;

# From the command line, get the following information:
# 1. The name of the (merged) para.txt file created from the R script to process GEMMA files
# 2. The name of the gff3 format file for the genome
# 3. The name of the gene info file (tab-delimited)
# 4. The significance threshold for significant SNPs
# 5. The distance in kb to search for genes
# 6. The name of a tab-delimited output file to create

my ($USAGE) = "$0 <input.param.txt> <input.gene> <gene.info> <p-value_threshold> <Max_distance (kb)> <output.txt>
\tinput.param.txt = The name of the GWAS results param.txt file from GEMMA
\tinput.gene = The name of the gff3 format file for the genome
\tgene.info = The name of the gene info file (tab-delimited)
\tp-value_threshold = The significance threshold for significant SNPs
\tMax_distance (kb) = The distance in kb to search for genes
\toutput.txt = The name of a tab-delimited output file to create\n";

unless (@ARGV) {
  print $USAGE;
  exit;
}

my ($snpfile, $gfffile, $infofile, $p_cutoff, $kb_dist, $output) = @ARGV;
my $max_dist = $kb_dist * 1000;

my @sig_snps = ();

# Read through the GWAS results file first and find the SNPs above the p-value cut-off
open (SNPS, $snpfile) || die "\nUnable to open the file $snpfile!\n";
while (<SNPS>) {
  chomp $_;
  if ($_ =~ /^CHR/i) {
    next;
  } else {  
    my @info = split(/\s{1,}/, $_);
  
    # First, check if the SNP is significant enough to save
    if ($info[5] >= $p_cutoff) {
      push (@sig_snps, $_);
    } else {
      next;
    }
  }
}
close(SNPS);

#Group SNPs within a certain distance together

my %snp_groups = ();

my $prev_chr = 0;
my $prev_pos = 0;
my @snps = ();
my $snp_group = 1;

foreach my $line (@sig_snps) {
  
  my @info = split(/\s{1,}/, $line);
  
  # First, check if the SNP is significant enough to save
  if ($info[5] >= $p_cutoff) {
    
    # If the SNP is significant, then check if it should be grouped with the previous SNP
    if (($prev_chr == $info[0]) && (($prev_pos + $max_dist) >= $info[2])) {

      # If the SNP groups with the previous SNP, then save it in the array
      # and update the prev_pos value
      push (@snps, $line);
      $prev_pos = $info[2];
    } 
    
    # If this is the first significant SNP to be saved, then update the values
    # and save the SNP in the array
    elsif ($prev_chr == 0) {
      $prev_chr = $info[0];
      $prev_pos = $info[2];
      push (@snps, $line);
    }

    # Otherwise, find the maximum and minimum values of the search window for the group of snps
    # Save this info. in the hash, and update the snp group number
    else {
      my $first_snp = $snps[0];
      my @first_info = split(/\s{1,}/, $first_snp);
      my $saved_chr = $first_info[0];
      my $min_pos = $first_info[2] - $max_dist;
      my $last_snp = $snps[(scalar @snps) - 1];
      my @last_info = split(/\s{1,}/, $last_snp);
      my $max_pos = $last_info[2] + $max_dist;
      my @temp = ();
      push (@temp, $saved_chr, $first_info[2], $min_pos, $max_pos, scalar @snps);
      @{$snp_groups{$snp_group}} = @temp;
      $snp_group++;

      $prev_chr = $info[0];
      $prev_pos = $info[2];
      @snps = ();
      push (@snps, $line);
    }
  } else {
    next;
  }
}

# Process the last SNP group
my $first_snp = $snps[0];
my @first_info = split(/\s{1,}/, $first_snp);
my $saved_chr = $first_info[0];
my $min_pos = $first_info[2] - $max_dist;
my $last_snp = $snps[(scalar @snps) - 1];
my @last_info = split(/\s{1,}/, $last_snp);
my $max_pos = $last_info[2] + $max_dist;
my @temp = ();
push (@temp, $saved_chr, $first_info[2], $min_pos, $max_pos, scalar @snps);
@{$snp_groups{$snp_group}} = @temp;

# For each SNP group, search the gff3 (gene) file
# for any features that overlap the range
# Find the feature in the info. file,
# Calculate the distance from the actual SNP, and print the info. to the output

# Open the output file for writing
open (OUT, ">$output") || die "\nUnable to open the file $output!\n";
print OUT "SNP_GROUP\t", "CHR\t", "SNP_POS\t", "NUM_SNPs\t", "GENE_ID\t", "GENE_Start\t", "GENE_End\t", "GeneDesc\t", "PantherID\t", "PfamID\t", "GO\t", "EC\t", "KeggID\t", "KogID\n";

# Sort the gene groups (which should sort snps by position), and then search one group at a time
my @groups = sort {$a <=> $b} (keys %snp_groups);

foreach my $group (@groups) {
  my ($c, $pos, $min, $max, $num_snp) = @{$snp_groups{$group}};
  my $check = 0;

  # Open the gene file to start sorting through it
  open (GENE, $gfffile) || die "\nUnable to open the file $gfffile!\n";
  
  FILELOOP:while (<GENE>) {
    chomp $_;
    if ($_ =~ /^\#/) {
      next;
    } 
    my @info = split(/\s{1,}/, $_);
    if ($info[2] =~ /mRNA/) {
      $info[0] =~ s/^Chr0//g;
      $info[0] =~ s/^Chr//g;
      $info[0] =~ s/super_/1/g;

      if ($info[0] < $c) {
	next FILELOOP;
      } elsif ($info[0] == $c) {
	if ((($min <= $info[3]) && ($max >= $info[3]) && ($max <= $info[4])) || (($min >= $info[3]) && ($max <= $info[4])) || (($min >= $info[3]) && ($min <= $info[4]) && ($max >= $info[4])) || (($min <= $info[3]) && ($max >= $info[4]))) {
	  my $geneID = (split(/;/, $info[8]))[0];
	  $geneID =~ s/ID=//g;
	  $geneID =~ s/\.v*\.1//ig;
	  
	  # Look up the gene info in the annotation info file
	  my %geneInfo = ();
	  open (INFO, $infofile) || die "\nUnable to open the file $infofile!\n";
	  while (<INFO>) {
	    chomp $_;
	    if ($_ =~ /^SI/) {
	      next;
	    }
	    my @ann_info = split(/\t/, $_);
	    if ($ann_info[2] =~ /$geneID/) {
	      if (exists $geneInfo{$ann_info[5]}) {
		my $temp = $geneInfo{$ann_info[5]};
		my @array = split(/,/, $temp);
		push (@array, $ann_info[4]);
		my $string = join(",", @array);
		$geneInfo{$ann_info[5]} = $string;
	      } else {
		$geneInfo{$ann_info[5]} = $ann_info[4];
	      }
	    }
	  }
	  close(INFO);

	  # Now print out all of the information for the gene
	  my $upstream = $pos - $info[3];
	  my $downstream = $pos - $info[4];
	  print OUT $group, "\t", $c, "\t", $pos, "\t", $num_snp, "\t", $geneID, "\t", $info[3], "\t", $info[4], "\t";
	  if (exists $geneInfo{"GeneDesc"}) {
	    print OUT $geneInfo{"GeneDesc"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"PantherID"}) {
	    print OUT $geneInfo{"PantherID"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"PfamID"}) {
	    print OUT $geneInfo{"PfamID"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"GO"}) {
	    print OUT $geneInfo{"GO"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"EC"}) {
	    print OUT $geneInfo{"EC"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"KeggID"}) {
	    print OUT $geneInfo{"KeggID"}, "\t",
	  } else {
	    print OUT "NA\t";
	  }
	  if (exists $geneInfo{"KogID"}) {
	    print OUT $geneInfo{"KogID"}, "\n",
	  } else {
	    print OUT "NA\n";
	  }
	  $check++;
	} else {
	  next FILELOOP;
	}
      } else {
	last FILELOOP;
      }
    } else {
      next FILELOOP;
    }
  }
  unless ($check > 0) {
    print OUT $group, "\t", $c, "\t", $pos, "\t", $num_snp, "\t", "No genes\t", "NA\t", "NA\t","NA\t","NA\t","NA\t","NA\t","NA\t", "NA\n";
  }
  close(GENE);
}

close(OUT);
exit;
      

