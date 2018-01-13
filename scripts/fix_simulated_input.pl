#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $file=$ARGV[0];
open FH, $file;

while (my $line = <FH>) {
    chomp $line;
    my @tmp=split("\t",$line);
    
    #next if ($tmp[1] =~ m/MT/);
    #next if ($tmp[1] =~ m/GL/);
    
    if ($tmp[13] =~ m/STRAND=-1/) {
        $tmp[0] =~ tr/TCGA/AGCT/; 
        $tmp[2] =~ tr/TCGA/AGCT/;
    }
    
print join("\t",@tmp)."\n";    
}