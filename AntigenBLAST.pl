#!/usr/bin/perl

################################################################
#
#  Author(s) - Mark M. Gosink, Ph.D.
#  Company -   Pfizer Inc
#
#  Function - Compares all possible n-mer peptides within an input
#		protein for matches in the selected proteome; output is
#		formatted into a line graph with hyperlinks to table of hits
#	Dependancies:
#		Required: NCBI's BLAST executables
#			Note: Files can be downloaded at ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST
#			Note: Appropriate executables must be installed; path to blastp can be set below
#		Required: BLASTable databases created using the included Make_Proteome_Sets.pl
#			Note: Database files are assumed to be in current directory location
#
################################################################

$infile = $ARGV[0];
$outfile_prefix = $infile;
$outfile_prefix =~ s/\.[^\.]+$//;

#	Set the minimally acceptable percent identity for full-length BLAST search hits
$protein_ident_thres = 50.0;

#	Set the minimally acceptable percent identity for peptide BLAST search hits
#		Note: Mismatches can be allowed by setting the appropriate percentage threshold
$peptide_ident_thres = 100.0;

#	Set the sliding window size for peptide fragments
@Pep_Sizes = (6, 8);

#	base link to NCBI's gene & protein databases
$gene_base_link = 'http://www.ncbi.nlm.nih.gov/gene/';
$protein_base_link = 'https://www.ncbi.nlm.nih.gov/protein/';

#$blastpcmd = '/path/to/blastp/executable/' . 'blastp';
$blastpcmd = '' . 'blastp';

$species = uc($ARGV[1]);
if ($species eq 'RAT') {
	$data_dir = './data/rat/';
	$db = $data_dir . 'Rat';
	$db_fasta = $data_dir . 'Rat.tfa';
	$go_file = $data_dir . 'Rat.GO_Component';
	$title = $infile . ' versus Rat RefSeq Proteins';
}
elsif ($species eq 'RABBIT') {
	$data_dir = './data/rabbit/';
	$db = $data_dir . 'Rabbit';
	$db_fasta = $data_dir . 'Rabbit.tfa';
	$go_file = $data_dir . '';		#	Rabbit genes have not been annotated with GO as of Feb 2018
	$title = $infile . ' versus Rabbit RefSeq Proteins';
}
elsif ($species eq 'MOUSE') {
	$data_dir = './data/mouse/';
	$db = $data_dir . 'Mouse';
	$db_fasta = $data_dir . 'Mouse.tfa';
	$go_file = $data_dir . 'Mouse.GO_Component';
	$title = $infile . ' versus Mouse RefSeq Proteins';
}
elsif ($species eq 'DOG') {
	$data_dir = './data/dog/';
	$db = $data_dir . 'Dog';
	$db_fasta = $data_dir . 'Dog.tfa';
	$go_file = $data_dir . 'Dog.GO_Component';
	$title = $infile . ' versus Dog RefSeq Proteins';
}
elsif ($species eq 'CYNO') {
	$data_dir = './data/cyno/';
	$db = $data_dir . 'Cynomolgus';
	$db_fasta = $data_dir . 'Cynomolgus.tfa';
	$go_file = $data_dir . '';		#	Cynomolgus genes have not been annotated with GO as of Feb 2018
	$title = $infile . ' versus Cynomolgus RefSeq Proteins';
}
else {
	$data_dir = './data/human/';
	$db = 'Human';
	$db_fasta = 'Human.tfa';
	$go_file = 'Human.GO_Component';
	$title = $infile . ' versus Human RefSeq Proteins';
	$species = 'HUMAN';
}

foreach $peptide_size (@Pep_Sizes) {
	my $outfile = $outfile_prefix . '_' . $peptide_size . 'aa.html';
	#	open the fasta protein seq file and generate a multifasta n-mer file
	open(FILE, $infile);
	@Sequence = ();
	while ($line = <FILE>) {
		if ($line =~ /^>/) { next }
		$line =~ s/[\n\f\r]+//;
		$line =~ s/[^A-Za-z]//g;
		my @Res =split(//, $line);
		push(@Sequence, @Res);
	}
	close(FILE);
	%Frag_Seqs = ();
	
	#	Unique file names for temporary files so multiple instances can run simultaneously
	$process_id = $$;
	$frags_tfa_file = $process_id . '_frags.tfa';
	$frags_blast_res_file = $process_id . '_frags.blast_results';
	$prot_blast_res_file = $process_id . '_whole.blast_results';
	
	open (OUTFILE, ">$frags_tfa_file");
	$count = 1;
	for ($idx = 0; $idx < ($#Sequence-$peptide_size); $idx++) {
		my $frag_name = 'frag_' . $count;
		my $frag_region = ($idx+1) . ".." . ($idx+$peptide_size);
		my $frag_seq = join('', @Sequence[$idx..($idx+($peptide_size-1))]);
		print OUTFILE ">$frag_name  $frag_region\n";
		print OUTFILE "$frag_seq\n";
		$Frag_Seqs{$frag_region} = $frag_seq;
		$count++;
	}
	close(OUTFILE);
	
	#	BLAST the original full length fasta seq
	my $cmd = $blastpcmd . ' -db ' . $db . ' -query ' . $infile . ' -seg no -evalue 1e-5 -outfmt 7 -out ' . $prot_blast_res_file;
	if (system("$cmd") != 0) { die "Something went wrong executing the BLAST command on the full length sequence.\n\tCMD: $cmd\n" }
	
	#	BLAST the generated n-mer peptides
	$cmd = $blastpcmd . ' -db ' . $db . ' -query ' . $frags_tfa_file . ' -seg no -num_alignments 4000 -word_size 2 -evalue 20000 -outfmt 7 -out ' . $frags_blast_res_file;
	if (system("$cmd") != 0) { die "Something went wrong executing the BLAST command on the peptide fragments.\n\tCMD: $cmd\n" }
	
	#	pull in symbol, Gene ID and protein name from sequence file
	%Gene_Info = ();
	%Prot_2_Gene = ();
	open (FILE, "$db_fasta");
	while ($line = <FILE>) {
		if ($line !~ /^>/) { next }
		if ($line =~ /^>(\S+)/) { $prot_id = $1 }
	
		if ($line =~ /GENEID:(\d+)\s+SYM:(\S+)\s+(.+)/) {
			my $gene_id = $1;
			$Gene_Info{$gene_id}{SYMBOL} = $2;
			$Gene_Info{$gene_id}{DEF} = $2;
			$Prot_2_Gene{$prot_id} = $gene_id;
		}
	}
	close(FILE);

	#	pull in the GO component annotations for each gene
	open (FILE, "$go_file");
	while ($line = <FILE>) {
		$line =~ s/[\n\r\f]+//;
		my @Vals = split(/\t/, $line);
		my $gene_id = $Vals[1];
		my $go = $Vals[5];
		$Gene_Info{$gene_id}{GO}{$go}++;
	}
	close(FILE);
	
	#	parse the BLAST results of the full-length protein search
	$found_first_hit = 'F';
	$input_gene_id = '';
	%HasFullLengthHits = ();
	%FullLenHitInfo = ();
	open (FILE, "$prot_blast_res_file");
	while ($line = <FILE>) {
		$line =~ s/[\n\r\f]+//;
		if ($line !~ /^#/) {
			my @Vals = split(/\t/, $line);
			my $hit = $Vals[1];
			my $identity = $Vals[2];
			my $aln_length = $Vals[3];
			my $mismatches = $Vals[4];
			my $gap_opens = $Vals[5];
			#	If the first hit in the BLAST result is a 100% match then it is presumed to be the same
			#	protein as the query and will be ignored for output
			if (($found_first_hit eq 'F') && ($identity == 100)) {
				$input_gene_id = $Prot_2_Gene{$hit};
				$found_first_hit = 'T';
				next;
			}
			if ($identity < $protein_ident_thres) { next }
			my $gene_id = $Prot_2_Gene{$hit};
			if ($gene_id eq $input_gene_id) { next }	#	skip reporting on protein if it is the input protein
			$HasFullLengthHits{$gene_id} = "T";
			$FullLenHitInfo{$gene_id}{IDENT} = $identity;
			$FullLenHitInfo{$gene_id}{LEN} = $aln_length;
			$FullLenHitInfo{$gene_id}{MISSMATCH} = $mismatches;
			$FullLenHitInfo{$gene_id}{GAPS} = $gap_opens;
		}
	}
	close(FILE);
	
	#	parse the BLAST results of the individual mer seq searchs
	$peptide_max_gaps = $peptide_size - int($peptide_size * ($peptide_ident_thres/100));
	$peptide_max_mismatch = $peptide_max_gaps;
	%PepHits = ();
	%PepHitInfo = ();
	@AllRegions = ();
	%RegionGeneHitCounts = ();
	$current_pos = 0;
	open (FILE, "$frags_blast_res_file");
	while ($line = <FILE>) {
		$line =~ s/[\n\r\f]+//;
		if ($line =~ /^# Query: frag_(\d+)  (\d+\.\.\d+)/) {
			$current_pos = $1;
			$current_region = $2;
			push(@AllRegions, $current_region);
		}
		# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
		if ($line =~ /^frag_\d+/) {
			my @Vals = split(/\t/, $line);
			my $hit = $Vals[1];
			my $identity = $Vals[2];
			my $aln_length = $Vals[3];
			my $mismatches = $Vals[4];
			my $gap_opens = $Vals[5];
			my $hit_region = $Vals[1] . ' (' . $Vals[8] . '..' . $Vals[9] . ')';
			my $hit_link = $protein_base_link . $Vals[1];
			if (($identity < $peptide_ident_thres)
				 || ($mismatches > $peptide_max_mismatch)
				 || ($aln_length < $peptide_size)
				 || ($gap_opens > $peptide_max_gaps)) { next }
			my $gene_id = $Prot_2_Gene{$hit};
			if ($gene_id == $input_gene_id) { next }	#	skip reporting on protein if it is the input protein
			$PepHits{$current_region}{$gene_id} = "T";
#		Uncomment the next line to add back links to hit protein sequences (also see below)
#			$PepHitInfo{$current_region}{$gene_id}{$hit_region}{LINK} = $hit_link;
			$PepHitInfo{$current_region}{$gene_id}{$hit_region}{IDENT} = $identity;
			$PepHitInfo{$current_region}{$gene_id}{$hit_region}{LEN} = $aln_length;
			$PepHitInfo{$current_region}{$gene_id}{$hit_region}{MISSMATCH} = $mismatches;
			$PepHitInfo{$current_region}{$gene_id}{$hit_region}{GAPS} = $gap_opens;
			$RegionGeneHitCounts{$current_pos}{$gene_id} = "T";
		}
	}
	close(FILE);
	
	open (OUTFILE, ">$outfile");
	print OUTFILE "
		<HTML>
		<HEAD>
		<TITLE>$title</TITLE>
		<style>
			.google-visualization-tooltip {
				position:relative;
			}
		</style>
		<script type='text/javascript' src='https://www.google.com/jsapi'></script>
		<script type='text/javascript'>
			google.load('visualization', '1', {packages:['corechart']});
			google.setOnLoadCallback(drawChart);
			function drawChart() {
				var data = new google.visualization.arrayToDataTable([
					['Position', 'Hits', {type:'string', role:'tooltip'}, 'link' ],
	";

	my $array_vals = '';
	for ($idx = 1; $idx <= $current_pos; $idx++) {
		my @Symbols = ();
		foreach $gene_id (keys(%{$RegionGeneHitCounts{$idx}})) {
			my $symbol = $Gene_Info{$gene_id}{SYMBOL};
			if ($symbol eq '') { $symbol = 'GeneID:' . $gene_id}
			push(@Symbols, $symbol);
		}
		@Symbols = sort {$a cmp $b} @Symbols;
		my $num_hits = $#Symbols + 1;
		my $tooltip = "Pos:$idx Num Hits:$num_hits\\n";
		my $symbols = "";
		for ($idx_2 = 0; $idx_2 <= $#Symbols; $idx_2++) {
			if (int(($idx_2+1)/6) == (($idx_2+1)/6)) { $symbols .= "$Symbols[$idx_2];\\n" }
			else { $symbols .= "$Symbols[$idx_2];" }
		}
		$symbols =~ s/;$//;
		$symbols =~ s/;\\n$//;
		$tooltip = $tooltip . $symbols;
		$array_vals .=  "\t\t\t['$idx', $num_hits, '$tooltip', '#POS_$idx' ],\n";
	}
	$array_vals =~ s/,$//;
	print OUTFILE  $array_vals;
	print OUTFILE "
				]);
			var view = new google.visualization.DataView(data);
			view.setColumns([0, 1, 2]);
			var options = {
				title: 'Peptide Matches',
				tooltip: { isHtml: true },
				legend: { position: 'none' }
			};
			var chart = new google.visualization.LineChart(document.getElementById('chart_div'));
			chart.draw(view, options);
			var selectHandler = function(e) {
				self.location = data.getValue(chart.getSelection()[0]['row'], 3 );
			}
			google.visualization.events.addListener(chart, 'select', selectHandler);
		}
		</script>
	";
	print OUTFILE "</HEAD>\n";

	print OUTFILE "<BODY><CENTER><H1>Homology analysis of all possible $peptide_size-mer antigens from '$infile'.</H1></CENTER>\n";
	if ($input_gene_id ne '') {
		my $def = $Gene_Info{$input_gene_id}{DEF};
		my $symbol = $Gene_Info{$input_gene_id}{SYMBOL};
		if ($symbol eq '') { $symbol = 'GeneID:' . $input_gene_id}
		print OUTFILE "The sequence in '$infile' displays 100% identity with '$symbol'.<BR>\n";
		print OUTFILE "  $symbol ($def) is presumed to be the input and hits will be ignored.<BR>\n<BR>\n";
	}
	print OUTFILE "$species proteins with similarity > $protein_ident_thres% against the full-length query are:<BR>\n";
	foreach $gene_id (sort {$FullLenHitInfo{$b}{IDENT} <=> $FullLenHitInfo{$a}{IDENT}} ((keys(%HasFullLengthHits)))) {
		my $symbol = $Gene_Info{$gene_id}{SYMBOL};
		if ($symbol eq '') { $symbol = 'GeneID:' . $gene_id}
		my $go = join("; ", keys(%{$Gene_Info{$gene_id}{GO}}));
		my $info = ' <FONT SIZE = "-2">' . $FullLenHitInfo{$gene_id}{IDENT} . '% Identity' .
						' over ' . $FullLenHitInfo{$gene_id}{LEN} . 'aa' .
						' with ' . $FullLenHitInfo{$gene_id}{MISSMATCH} . ' mismatches' .
						' and ' . $FullLenHitInfo{$gene_id}{GAPS} . ' gap-opens(s).</FONT>';
		print OUTFILE "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;$symbol",
				" (<A HREF='$gene_base_link$gene_id'>$gene_id</A> ",
				" - $go) $info<BR>\n";
	}
	
	
	print OUTFILE "<HR>\n<H2>Peptide Comparisons</H2>";
	print OUTFILE "<div id='chart_div' style='width: auto; height: 350pt;'></div>\n";
	print OUTFILE "<H3>List of Peptide Matches Sorted By Number of Matches</H3>\n";
	%Symbol_2_Loc = ();
	print OUTFILE "$species proteins with similarity &gt; $peptide_ident_thres% against one or more query $peptide_size-mers are:<BR>\n<BR>\n";
	for ($idx = 0; $idx<=$#AllRegions; $idx++) {
		my $current_region = $AllRegions[$idx];
		if (defined $PepHits{$current_region}) {
			my $frag_seq = $Frag_Seqs{$current_region};
			foreach $gene_id (keys(%{$PepHits{$current_region}})) {
				my $symbol = $Gene_Info{$gene_id}{SYMBOL};
				if ($symbol eq '') { $symbol = 'GeneID:' . $gene_id} 
				my $go = join("; ", keys(%{$Gene_Info{$gene_id}{GO}}));
				$Symbol_2_Loc{$symbol}{LOC} .= "\t$current_region";
				$Symbol_2_Loc{$symbol}{MATCH_NUM}++;
				$Symbol_2_Loc{$symbol}{ENTREZ} = $gene_id;
				$Symbol_2_Loc{$symbol}{GO} = $go;
			}
			
		}
	}
	
	@SymList = sort by_match (keys(%Symbol_2_Loc));
	foreach $symbol (@SymList) {
		print OUTFILE "$symbol",
				" (<A HREF='$gene_base_link$Symbol_2_Loc{$symbol}{ENTREZ}'>$Symbol_2_Loc{$symbol}{ENTREZ}</A>)<BR>\n",
				"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;GO:$Symbol_2_Loc{$symbol}{GO}<BR>\n",
				"&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;MATCHES:$Symbol_2_Loc{$symbol}{MATCH_NUM} => $Symbol_2_Loc{$symbol}{LOC}<BR>\n";
	}
	print OUTFILE "<BR>\n<CENTER>-----------------</CENTER><BR><H3>List of Peptide Matches Sorted By Start Position</H3>\n";
	
	%Symbol_2_Loc = ();
	print OUTFILE "$species proteins with similarity &gt; $peptide_ident_thres% against one or more query $peptide_size-mers are:<BR>\n<BR>\n";
	for ($idx = 0; $idx<=$#AllRegions; $idx++) {
		my $current_region = $AllRegions[$idx];
		if (defined $PepHits{$current_region}) {
			my $frag_seq = $Frag_Seqs{$current_region};
			my $start_loc = $current_region;
			$start_loc =~ s/\.\.\d+$//;
			print OUTFILE "<a name='POS_$start_loc'>$current_region</a> ($frag_seq)<UL style='margin:0;'>\n";
			foreach $gene_id (keys(%{$PepHits{$current_region}})) {
				my $symbol = $Gene_Info{$gene_id}{SYMBOL};
				if ($symbol eq '') { $symbol = 'GeneID:' . $gene_id} 
				my $go = join("; ", keys(%{$Gene_Info{$gene_id}{GO}}));
				$Symbol_2_Loc{$symbol}{LOC} .= "\t$current_region";
				$Symbol_2_Loc{$symbol}{MATCH_NUM}++;
				$Symbol_2_Loc{$symbol}{ENTREZ} = $gene_id;
				$Symbol_2_Loc{$symbol}{GO} = $go;
				my $info = '  <FONT SIZE = "-2" COLOR="Gray">';
				foreach $hit_region (keys($PepHitInfo{$current_region}{$gene_id})) {
#		Uncomment the next line and comment out the line following that to add back links to hit protein sequences (also see above)
#					$info .= "<A HREF='$PepHitInfo{$current_region}{$gene_id}{$hit_region}{LINK}'>" . $hit_region . '</A>' .
					$info .= $hit_region .
									' -> ' . $PepHitInfo{$current_region}{$gene_id}{$hit_region}{IDENT} . '% Identity' .
									' over ' . $PepHitInfo{$current_region}{$gene_id}{$hit_region}{LEN} . 'aa' .
									' with ' . $PepHitInfo{$current_region}{$gene_id}{$hit_region}{MISSMATCH} . ' mismatches' .
									' and ' . $PepHitInfo{$current_region}{$gene_id}{$hit_region}{GAPS} . ' gap-open(s). | ';
				}
				$info =~ s/ \| $//;
				$info .= '</FONT>';
				print OUTFILE "<LI>$symbol",
						" (<A HREF='$gene_base_link$gene_id'>$gene_id</A>",
						" - $go)$info</LI>";
			}
			print OUTFILE "</UL><BR>\n";
		}
	}
	
	print OUTFILE "</BODY>\n",
			"</HTML>\n";
	close(OUTFILE)
}

#	cleanup temporary files
unlink("$frags_tfa_file");
unlink("$frags_blast_res_file");
unlink("$prot_blast_res_file");

exit;

sub by_match {
	if ($Symbol_2_Loc{$a}{MATCH_NUM} > $Symbol_2_Loc{$b}{MATCH_NUM}) { return -1 }
	elsif ($Symbol_2_Loc{$a}{MATCH_NUM} < $Symbol_2_Loc{$b}{MATCH_NUM}) { return 1 }
	else { return $a cmp $b }
};
