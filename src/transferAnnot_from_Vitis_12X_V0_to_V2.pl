#!/usr/bin/perl

#@author: camille.rustenholz@colmar.inra.fr
#Small changes: nacer.mohellibi@versailles.inra.fr
#@License: GPL-3+

#use strict;
use warnings;
use Getopt::Long;
use FileHandle;

my $posi;
my $annot;
my $output = "results.gff3";

my $result = &GetOptions("pos|p=s" => \$posi,
			 "annot|a=s" => \$annot,
			 "output|o=s" => \$output);

unless (defined($posi) && defined($annot))
	{
	print STDERR sprintf("\n");
	print STDERR sprintf("Usage: perl transferAnnot_from_Vitis_12X_V0_to_V2.pl --pos|-p <position file> --annot|-a <annot file in gff format> --output|-o <Output file in gff format>\n\n");
	exit;
}

#First step, read position file and store informations in the hash %list

open $posi, $posi or die $!; 

my %list;
my @line;
my $ori;

while (<$posi>){
	next if($_ =~ "#.*");
	chomp;
	@line = split ("\t",$_);
	if ($line[6] eq "\*"){
		$ori = $line[2]."+";
	} else {
		$ori = $line[2].$line[6];
	}
	push @{$list{$line[1]}{$line[0]}}, $ori, $line[3], $line[4], $line[5], $line[7], $line[8], $line[9], $line[10], $line[11], $line[12], $line[13];

}

close $posi;
open $annot, $annot or die $!;
open $output, ">$output";

# Second step, read annotation file, find on which scaffold they are located and transfer coordinates

my $start;
my $stop;

while (<$annot>){
	chomp;
	@line = split ("\t",$_);

    #END ITERATE OVER CHROMOSOMES
	foreach my $chr (keys %list){
		#for the concerned chromosome
		if ($chr eq $line[0]){
            #iterate over scaffolds
			foreach my $scaff (keys %{$list{$chr}}){
				my $scaffold_start_v0 = $list{$chr}{$scaff}[1];
				my $scaffold_end_v0 = $list{$chr}{$scaff}[2];
				my $feature_start_vo = $line[3];
				my $feature_end_vo = $line[4];
				my $ori_scaff = $list{$chr}{$scaff}[0];
				my $scaffold_start_v2 = $list{$chr}{$scaff}[4];
				my $previous_scaf_v0 = $list{$chr}{$scaff}[6];
				my $next_scaf_v0 = $list{$chr}{$scaff}[7];
				my $previous_scaf_v2 = $list{$chr}{$scaff}[8];
                my $next_scaf_v2 = $list{$chr}{$scaff}[9];
                my $stop_next_scaf_vo = $list{$chr}{$scaff}[10];
                
				#Features totally in a scaffold
				if ($scaffold_start_v0 <=$feature_start_vo && $scaffold_end_v0 >=$feature_end_vo){
					$len = $feature_end_vo - $feature_start_vo;

					if (($ori_scaff eq "++") || ($ori_scaff eq "--") || ($ori_scaff eq "+*")){
						
						$start = $feature_start_vo - $scaffold_start_v0 + $scaffold_start_v2;
					}
					else {
						
						$start = $scaffold_start_v2 + $scaffold_end_v0 - $feature_end_vo;
						
						if ($line[6] eq "+"){
							$line[6] = "-";
						}
						else {
							$line[6] = "+";
						}
					}
					$stop = $start + $len;
					print {$output} $list{$chr}{$scaff}[3], "\t",join("\t",@line[1..2]), "\t", $start, "\t", $stop, "\t", join("\t",@line[5..8]), "\n";
					last;
				}#End Feature in a scaffold
				#Features having the same next scaffold in both assemblies...
				elsif( ($next_scaf_v0 eq $next_scaf_v2) && ($next_scaf_v0 ne "--")  && ($next_scaf_v0 ne "----")){
					#...and that overlap two successive scaffolds
					if( ($scaffold_start_v0 <=$feature_start_vo) && ($scaffold_end_v0 >= $feature_start_vo) && ($feature_end_vo <= $stop_next_scaf_vo) ){
						#if they have the same orientation
						if (($ori_scaff eq "++") || ($ori_scaff eq "--") || ($ori_scaff eq "+*")){
							$len = $feature_end_vo - $feature_start_vo;
							$start = $feature_start_vo - $scaffold_start_v0 + $scaffold_start_v2;
							$stop = $start + $len;
	                        print {$output} $list{$chr}{$scaff}[3], "\t",join("\t",@line[1..2]), "\t", $start, "\t", $stop, "\t", join("\t",@line[5..8]), "\n";
	                        last;
						}
						
					}
				}
			}
			last;
		}#END ITERATE OVER CHROMOSOMES
	}
}

close $annot;
close $output;
exit;
