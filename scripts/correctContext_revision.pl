#!/usr/local/bin/perl

use strict;
#use warnings;
use Data::Dumper;

my %exp = GetExpFreq($ARGV[0]);

my %obs = GetObsFreq($ARGV[1]);
#print Dumper(%obs);

open FH, $ARGV[1];

while (my $line = <FH>) {
    chomp $line;
    next if ($line =~ /Ensemblid/);
    
    #print $line."\n";
    my @tmp=split("\t",$line);
    
    my $FCna=$exp{$tmp[1]}{freqObs}/$obs{$tmp[1]}{freqNaExp};   
    my $FCns=$exp{$tmp[1]}{freqObs}/$obs{$tmp[1]}{freqNsExp};   
    
    my $newNa=$tmp[2]*$FCna;
    my $newNs=$tmp[3]*$FCns;
    
    print join ("\t",$tmp[0],$tmp[1],$tmp[2],$tmp[3],$newNa,$newNs
                #,$tmp[5]/$tmp[6],$newNa/$newNs
                )."\n";
}


sub GetObsFreq{
    my $file=shift;
    my %freq;
    open( FQ, $file );
    while (my $line = <FQ>){
        chomp $line;
        next if ($line =~ /Ensemblid/);
        my @tmp = split("\t",$line);
        $freq{$tmp[1]}{count_Na} += $tmp[2] ;
        $freq{$tmp[1]}{count_Ns} += $tmp[3] ;
    }
    
    my $totalNs   = 0;
    my $totalNa   = 0;
    foreach my $changes (keys(%freq)) {
        $totalNa += $freq{$changes}{count_Na};
        $totalNs += $freq{$changes}{count_Ns};
    }


    foreach my $changes2 (keys(%freq)) {
         $freq{$changes2}{freqNaExp}= ($freq{$changes2}{count_Na}+$freq{$changes2}{count_Ns})/($totalNa + $totalNs);
         $freq{$changes2}{freqNsExp}= ($freq{$changes2}{count_Na}+$freq{$changes2}{count_Ns})/($totalNa + $totalNs);
    }
    
    return %freq;
}

    
sub GetExpFreq{
    my $file2=shift;
    my %freq2;
    open( FQ2, $file2 );
    while (my $line2 = <FQ2>){
        chomp $line2;
        my @tmp2 = split("\t",$line2);
        
        $freq2{$tmp2[0]}{count_muts} += $tmp2[1];
    }
    my $totalmuts = 0;
    foreach my $changes2 (keys(%freq2)) {
            #print $changes."\t".$freq{$changes}{count_muts}."\n";
    $totalmuts += $freq2{$changes2}{count_muts};
    }
    
    foreach my $changes3 (keys(%freq2)) {
        $freq2{$changes3}{freqObs}= $freq2{$changes3}{count_muts}/$totalmuts;
    }
    
    return %freq2;
}    
