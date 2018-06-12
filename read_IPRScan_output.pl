# June 11th 2018
# Read TSV output file from InterProScan and print out information for each protein in 1 row
# ProteinId	Pfam	Interpro	TMHMM	Pathway	GO

#! /usr/bin/perl -w
use strict;

my $filein_iprscantsv="/home/mnguyen/Research/For_Marcos/Dicsqu464_1/IPRScan/Dicsqu464_1_GeneCatalog_proteins_20151220.aa.fasta.tsv";
my $fileout="/home/mnguyen/Research/For_Marcos/Dicsqu464_1/IPRScan/Disc464_1_IPRScan_1prot_1row_11June2018.txt";
my $filein_GOID_GOterms="/home/mnguyen/Research/For_Marcos/Dicsqu464_1/Dics464_1_GOID_GOterms.tsv";

################################################################################################################
open(GO,"<$filein_GOID_GOterms") || die "Cannot open file $filein_GOID_GOterms";
my %hash_go_id_term;
my %hash_go_aspect;
while (<GO>)
{
	$_=~s/\s*$//;
	if ($_=~/^GO Term/){next;}
	my @cols=split(/\t/,$_);
	my $goid=$cols[0];
	my $goterm=$cols[1];
	my $aspect=$cols[2];
	$goid=~s/\s*//g;
	$goterm=~s/\s*$//;
	$aspect=~s/\s*//g;
	$hash_go_id_term{$goid}=$goterm;
	$hash_go_aspect{$goid}=$aspect;
}
close(GO);

################################################################################################################


################################################################################################################
open(In,"<$filein_iprscantsv") || die "Cannot open file $filein_iprscantsv";
open(Out,">$fileout") || die "Cannot open file $fileout";
my %hash_pfam;
my %hash_ipr;
my %hash_tmhmm;
my %hash_pathways;
my %hash_go_f;
my %hash_go_p;
my %hash_go_c;
my %hash_all_protids;
while (<In>)
{
	$_=~s/\s*$//;
	#print "\ntest\n";
	my @cols=split(/\t/,$_);
	my $protid=$cols[0];
	my $analysis=$cols[3];
	my $accession=$cols[4];
	my $desc=$cols[5];
	my $start=$cols[6];
	my $end=$cols[7];
	my $evalue=$cols[8];
	my $ipr_id=$cols[11];
	my $ipr_desc=$cols[12];
	my $go_id=$cols[13];
	my $pathway=$cols[14];
	
	$hash_all_protids{$protid}++;
	
	$protid=~s/\s*//g;
	$analysis=~s/\s*//g;
	if ($analysis eq "TMHMM")
	{
		my $tmhmm_location=$start."..".$end;
		if ($hash_tmhmm{$protid}){$hash_tmhmm{$protid}=$hash_tmhmm{$protid}.";".$tmhmm_location;}
		else{$hash_tmhmm{$protid}=$tmhmm_location;}
	}elsif($analysis eq "Pfam")
	{
		my $new_pfam=$accession."(".$start."-".$end.",".$evalue."):".$desc;
		if ($hash_pfam{$protid}){$hash_pfam{$protid}=$hash_pfam{$protid}." | ".$new_pfam;}
		else{$hash_pfam{$protid}=$new_pfam;}
	}
	
	#######################################################################################################
	# GO
	$go_id=~s/\s*//g;
	if ($go_id)
	{
		my @go_ids=split(/\|/,$go_id);
		foreach my $go (@go_ids)
		{
			my $aspect=$hash_go_aspect{$go};
			my $go_term=$hash_go_id_term{$go};
			my $new_go=$go.":".$go_term;
			if ($aspect)
			{
				if ($aspect eq "molecular_function")
				{
					if ($hash_go_f{$protid}){$hash_go_f{$protid}=$hash_go_f{$protid}."|".$new_go;}
					else{$hash_go_f{$protid}=$new_go;}
				}elsif ($aspect eq "biological_process")
				{
					if ($hash_go_p{$protid}){$hash_go_p{$protid}=$hash_go_p{$protid}."|".$new_go;}
					else{$hash_go_p{$protid}=$new_go;}
				}else
				{
					if ($hash_go_c{$protid}){$hash_go_c{$protid}=$hash_go_c{$protid}."|".$new_go;}
					else{$hash_go_c{$protid}=$new_go;}
				}
			}
		}
	}
	#######################################################################################################
	$pathway=~s/\s*//g;
	if ($pathway)
	{
		if ($hash_pathways{$protid}){$hash_pathways{$protid}=$hash_pathways{$protid}."|".$pathway;}
		else{$hash_pathways{$protid}=$pathway;}
	}
	$ipr_id=~s/\s*//g;
	if ($ipr_id)
	{
		my $new_ipr=$ipr_id.":".$ipr_desc;
		if ($hash_ipr{$protid}){$hash_ipr{$protid}=$hash_ipr{$protid}."|".$new_ipr;}
		else{$hash_ipr{$protid}=$new_ipr;}
	}
}

########################################################################################################################################################

print Out "#ProteinID	Pfam(start-end,evalue)\tInterpro\tTMHMM locations\tPathways\tGO_Function\tGO_Process\tGO_Cellular_Component\n";
my @all_proids=keys(%hash_all_protids);
foreach my $each_protid (@all_proids)
{
	my $pfam_domains=$hash_pfam{$each_protid};
	my $interpro=$hash_ipr{$each_protid};
	my $tmhmm_locations=$hash_tmhmm{$each_protid};
	my $all_pathways=$hash_pathways{$each_protid};
	my $go_function=$hash_go_f{$each_protid};
	my $go_process=$hash_go_p{$each_protid};
	my $go_component=$hash_go_c{$each_protid};
	
	my $interpro_nr=&Remove_duplicates($interpro);
	my $pathways_nr=&Remove_duplicates($all_pathways);
	my $go_function_nr=&Remove_duplicates($go_function);
	my $go_process_nr=&Remove_duplicates($go_process);
	my $go_component_nr=&Remove_duplicates($go_component);
	
	print Out "$each_protid\t$pfam_domains\t$interpro_nr\t$tmhmm_locations\t$pathways_nr\t$go_function_nr\t$go_process_nr\t$go_component_nr\n";
}
close(In);
close(Out);
################################################################################################################


################################################################################################################
sub Remove_duplicates
{
	my $x=$_[0];
	my @array_temp=split(/\|/,$x);
	my %hash_temp;
	foreach my $temp (@array_temp){$hash_temp{$temp}++;}
	my @array_temp_nr=keys(%hash_temp);
	my $y=join(" | ",@array_temp_nr);
}
################################################################################################################