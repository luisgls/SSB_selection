#!/usr/bin/perl
#
#

use strict;
use warnings;
use Data::Dumper;

open FH, $ARGV[0];

my $nonsil=0;
my $totnonsil=0;
my $syn=0;
my $totsyn=0;
my $globaldnds;

while (my $line = <FH>) {
	chomp $line;
#header
	next if $line =~ "EnsemblID";
#EnsemblID       nonsilent_variant       synonymous_variant      nonsyn_sites    syn_sites       nonsynsites_SSB192      synsites_SSB192 dnds_SSB192     pval_SSB192     pval_SSB192.adj Hugo_symbol
	my @tmp=split("\t",$line);
	$nonsil  += $tmp[1];
	$totnonsil += $tmp[5];
	$syn += $tmp[2];
	$totsyn += $tmp[6];
}

if ( $totnonsil > 0 && $totsyn >0){
	$globaldnds = ($nonsil/$totnonsil)/($syn/$totsyn);
}else{
	print STDERR "check total number of sites\n";
	$globaldnds = ($nonsil/($totnonsil+1))/($syn/($totsyn+1));	 
}

print $globaldnds."\n";
