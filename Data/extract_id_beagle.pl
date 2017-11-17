#!/usr/bin/perl
#script to get the desired data of a set of individuals from the beagle output (to run saber+)
#Inputs: phased data file, ID file, chr #
#30 may 2012

use List::MoreUtils qw(indexes);
use strict; use warnings;

open(IN, "<$ARGV[0]") or die "error reading $ARGV[0]"; #name of phased data file
open(IDS, "<$ARGV[1]") or die "error reading $ARGV[1]"; #name of file with IDs to extract

my $popname = $ARGV[1];
my $chr=$ARGV[2];
my @subs=split(/\//,$popname);
my $popname2=$subs[-1];
$popname2 =~ s/\.txt//;
my $outfile=$popname2.".".$chr;
open(OUT, ">$outfile") or die "error creating $outfile";

my @sampleIDs=();
while (<IDS>) {
	chomp;
	push(@sampleIDs,$_);
}

my $long=scalar(@sampleIDs);

for (my $i=0;$i<$long;$i++){
	if ($i==$long-1){
		print OUT "$sampleIDs[$i] $sampleIDs[$i]\n";
	}
	else {	
		print OUT "$sampleIDs[$i] $sampleIDs[$i] ";
	}
}

my $nSNPs=-1;
my @indexes=();
my $nChrom=0;

#go through each line of the phased data file and convert it to the desired format
while (my $line=<IN>) {
	chomp($line);
	$nSNPs++; #SNPs start on line 2 so nSNPs=1 at line 2
	my $nBases = 0;
	my @curBases = ();
	my @split_line = split(' ', $line);
	#print "@split_line\n";
	#find columns corresponding to the provided GRC numbers
	if ($nSNPs==0) {
		shift(@split_line); # split_line is now just GRC numbers
		foreach my $sID (@sampleIDs) {
			my @locs = indexes {$_ eq $sID } @split_line;
			push(@indexes, @locs);
		}
		$nChrom=@indexes;
		#print "@indexes \n";
		#print "number of chromosomes= $nChrom \n";
	}
	
	elsif ($nSNPs > 0) {
		#remove the first element that is the rsNumber
		shift(@split_line);
		@split_line=@split_line[@indexes]; #now split_line is just the columns for the desired samples
		#print "@split_line\n";
		my $long2=scalar(@split_line);
		for(my $j=0;$j<$long2;$j++){
			if ($j==$long2-1){
				print OUT "$split_line[$j]\n";
			}
			else{
				print OUT "$split_line[$j] ";
			}
		}
	}
}

close IN;
close IDS;
close OUT;
print "Done!\n";



