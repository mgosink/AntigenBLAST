#!/usr/bin/perl

################################################################
#
#  Author(s) - Mark M. Gosink, Ph.D.
#  Company -   Pfizer Inc
#
#  Creation Date - March 16, 2017
#
#  Function - Script for parsing NCBI's GenBank protein sequences
#		and NCBI Gene gene_info file and formating to make
#		BLASTable databases
#	Dependancies:
#		Required: NCBI's BLAST executables
#			Note: Files can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
#			Note: Appropriate executables must be installed; path to makeblastdb can be set below 
#		Required: Organism-specific protein files in GenBank format and gz compressed
#			Note: Files can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/genomes/<organism>/protein/protein.gbk.gz
#		Required: gene_info file from NCBI Gene database
#			Note: File can be downloaded at ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
#			Note: Tool assumes file is uncompressed after downloading
#		Optional: gene_info file from NCBI Gene database
#			Note: File can be downloaded at ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz
#			Note: Tool assumes file is uncompressed after downloading
#
################################################################

#	Default organisms for proteome databases
#		Note: edit paths as appropriate
%Organisms = (
	'9606'	=>	{
		'NAME'	=>	'Human',
		'DATA_DIR'	=>	'./data/human/'
	},
	'9541'	=>	{
		'NAME'	=>	'Cynomolgus',
		'DATA_DIR'	=>	'./data/cyno/'
	},
	'10090'	=>	{
		'NAME'	=>	'Mouse',
		'DATA_DIR'	=>	'./data/mouse/'
	},
	'10116'	=>	{
		'NAME'	=>	'Rat',
		'DATA_DIR'	=>	'./data/rat/'
	},
	'9615'	=>	{
		'NAME'	=>	'Dog',
		'DATA_DIR'	=>	'./data/dog/'
	},
	'9986'	=>	{
		'NAME'	=>	'Rabbit',
		'DATA_DIR'	=>	'./data/rabbit/'
	}
);

#$makeblastcmd = '/path/to/makeblastdb/executable/' . 'makeblastdb';
$makeblastcmd = '' . 'makeblastdb';

#$geneinfofile = '/path/to/geneinfo/file/' . 'gene_info';
$geneinfofile = '' . 'gene_info';

#$genegofile = '/path/to/gene2go/file/' . 'gene2go';
$genegofile = '' . 'gene2go';

#	Build a hash mapping GeneIDs to gene symbols
%GeneID_2_Sym = ();
open(FILE, "$geneinfofile");
while ($line = <FILE>) {
	if ($line =~ /^#/) { next } 	#	skip comment lines
	$line =~ s/[\n\r\f]+//;
	my @Vals = split(/\t/, $line);
	my $tax_id = $Vals[0];
	if (not exists $Organisms{$tax_id})	{ next }	#	skip entries for organisms we are not capturing
	my $geneid = $Vals[1];
	my $symbol = $Vals[2];
	if (($symbol eq '') || ($symbol eq '-')) { next }
	$GeneID_2_Sym{$geneid} = $symbol;
}
close(FILE);

#	Parse the Gene databases Gene Ontology mapping file to capture the Cellular Component entries
open (FILE, "$genegofile");
while ($line = <FILE>) {
	$line =~ s/[\n\r\f]+//;
	my @Vals = split(/\t/, $line);
	my $tax_id = $Vals[0];
	my $go_cat = $Vals[7];
	if ($go_cat ne 'Component') { next }	#	just need cellular component data
	if (exists $Organisms{$tax_id})	{
		my $species = $Organisms{$tax_id}{NAME};
		my $data_dir = $Organisms{$tax_id}{DATA_DIR};
		my $go_file = $data_dir . $species . '.GO_Component';
		open (OUTFILE, ">>$go_file");
		print OUTFILE "$line\n";
		close (OUTFILE);
	}
}
close(FILE);

#	Foreach chosen species, generate a fasta file and then make a BLASTable form
foreach $tax_id (keys(%Organisms)) {
	my $data_dir = $Organisms{$tax_id}{DATA_DIR};
	my $protein_file = $data_dir . 'protein.gbk.gz';
	my $species = $Organisms{$tax_id}{NAME};
	my $fasta_file = $data_dir . $species . '.tfa';
	open(OUTFILE, ">$fasta_file");
	my %Proteins = ();
	my $in_seq_flag = 'F';
	my $sequence = "";
	open(FILE, "gunzip -c $protein_file |") or die "gunzip $protein_file: $!";
	while ($line = <FILE>) {
		$line =~ s/[\n\f\r]+//;
		if ($line =~ /^LOCUS\s+(\S+)/) { $current_prot = $1 }
		if ($line =~ /^DEFINITION\s+(.+)/) { $current_def = $1 }
		if ($line =~ /^\s+\/db_xref="GeneID:(\d+)/) { $geneid = $1 }
		if ($line =~ /^ORIGIN/) {
			$in_seq_flag = 'T';
			$sequence = "";
			next;
		}
		if ($line =~ /^\/\//) {
			$in_seq_flag = 'F';
			my $symbol = $GeneID_2_Sym{$geneid};
			print OUTFILE ">$current_prot  GENEID:$geneid SYM:$symbol $current_def\n";
			print OUTFILE "$sequence\n";
			next;
		}
		if ($in_seq_flag eq 'T') {
			my $res = $line;
			$res =~ s/[\s\d]+//g;
			$sequence .= "$res\n";
		}
	}
	close(FILE);
	close(OUTFILE);
	
	my $cmd = $makeblastcmd . ' -dbtype prot -in ' . $fasta_file . " -title '$species Protein Sequences' -out $data_dir$species";
	if (system("$cmd") == 0) { print "Successfully created the BLASTable $species database!\n" }
	else { die "Something went wrong making the BLASTable $species database\n\tCMD: $cmd\n" }
}

exit;
