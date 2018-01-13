#!/usr/bin/perl
#
#
use Data::Dumper;

open FH,$ARGV[0];
my %hugo=getHugo($ARGV[1]);

#print Dumper(%hugo);

while (my $line = <FH>){
	chomp $line;
	my @tmp=split("\t",$line);
	my $hugo="no_symbol";	
	if($line =~ m/ensembl_id|EnsemblID|Ensembl_ID|Ensembl_id/){
		print $line."\tHugo_symbol"."\n";
		next;
	}

	if(exists $hugo{$tmp[0]}{gene}){
		$hugo=$hugo{$tmp[0]}{gene};
	}
	print $line."\t".$hugo."\n";

}


sub getHugo{
	open FH2, $_[0];
	my %hash;
	while (my $line2 = <FH2>){
		chomp $line2;
		my @tmp2=split("\t",$line2);
		$hash{$tmp2[1]}{gene}=$tmp2[0];	
	}
	return %hash;
}
