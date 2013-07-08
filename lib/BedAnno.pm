package BedAnno;

use 5.010001;
use strict;
use warnings;
use Carp;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use BedAnno ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.

our @EXPORT = qw(
    fetchseq get_codon parse_var
);

our $VERSION = '0.11';

=head1 NAME

BedAnno - Perl module for annotating variation depend on the BED +1 format database.

=head1 SYNOPSIS

  use BedAnno;
  my $beda = BedAnno->new( db => "in.bed.gz", regbed => 'in.region.bed', codon => 'in.coding.fa', trans => \%trans);
  my $anno = $beda->anno( 'chr20', 1234567, 'TG', 'TGGG' ); 

=head1 DESCRIPTION

This module need bgzipped BED+1 format database, combined with tabix index,
tabix is require for fast extract anno infomation randomly.
The fasta sequences are indexed by samtools faidx, sequence splicing also requires
samtools.

=head2 EXPORT

fetchseq()  - used to batched fetch genomic sequences (flanks), 
	      this depend on samtools installed, and a faidx(ed) 
	      fasta file.
get_codon() - used to get codon information through cds position.
parse_var() - parse variation information stand-alone, giving
	      the guess information and possible coordinates.

=head1 Methods
=cut

our %codon3 = (
    TTT=>"Phe",	CTT=>"Leu", ATT=>"Ile",	GTT=>"Val",
    TTC=>"Phe", CTC=>"Leu", ATC=>"Ile", GTC=>"Val",
    TCT=>"Ser", CCT=>"Pro", ACT=>"Thr", GCT=>"Ala",
    TCC=>"Ser", CCC=>"Pro", ACC=>"Thr", GCC=>"Ala",
    TAT=>"Tyr", CAT=>"His", AAT=>"Asn", GAT=>"Asp",
    TAC=>"Tyr", CAC=>"His", AAC=>"Asn", GAC=>"Asp",
    TGT=>"Cys", CGT=>"Arg", AGT=>"Ser", GGT=>"Gly",
    TGC=>"Cys", CGC=>"Arg", AGC=>"Ser", GGC=>"Gly",
    TTA=>"Leu", CTA=>"Leu", ATA=>"Ile", GTA=>"Val",
    TCA=>"Ser", CTG=>"Leu", ACA=>"Thr", GTG=>"Val",
    TAA=>"*",   CCA=>"Pro", AAA=>"Lys", GCA=>"Ala",
    TGA=>"*",   CCG=>"Pro", AGA=>"Arg", GCG=>"Ala",
    TTG=>"Leu", CAA=>"Gln", ATG=>"Met", GAA=>"Glu",
    TCG=>"Ser", CAG=>"Gln", ACG=>"Thr", GAG=>"Glu",
    TAG=>"*",   CGA=>"Arg", AAG=>"Lys", GGA=>"Gly",
    TGG=>"Trp", CGG=>"Arg", AGG=>"Arg", GGG=>"Gly",

    TCN=>"Ser", CCN=>"Pro", ACN=>"Thr", GTN=>"Val",
		CTN=>"Leu",		GCN=>"Ala",
		CGN=>"Arg",		GGN=>"Gly"
);

our %Polar = (
    Ala => "NP", Asn => "P0", Arg => "P+", Asp => "P-",
    Ile => "NP", Cys => "P0", His => "P+", Glu => "P-",
    Leu => "NP", Gln => "P0", Lys => "P+",
    Met => "NP", Gly => "P0",
    Phe => "NP", Ser => "P0",
    Pro => "NP", Thr => "P0",
    Trp => "NP", Tyr => "P0",
    Val => "NP",

    '*' => '.'
);

=head2 new

    About   : Creat a new annotation entry, automatically call load_anno method.
    Usage   : my $beda = BedAnno->new( db => 'in.bed.gz', regbed => 'in.region.bed', codon => 'in.coding.fa[.gz]');
    Args    : The args db, codon are necessary arg, this codon is cut from hg19 the following are optional args.
	      "genes"  => {"ABC" => 1, "DEF" => 1} / "genes.list"		    limit the genes to be annotate
	      "trans"  => {"NM_0012.1" => 1, "NM_0034.2" => 1} /  "trans.list"	    limit the transcripts to be annotate
	      "region" => "chr20:1234567-1234568"				    give the target region in region string
	      "regbed" => "in.region.bed"					    give the target region in bed format region files.
	      "onlyPr" => 1							    to limit the anno range to the primary transcript for a gene.
	      "mmap"   => 1							    to output all needed multiple mapping record of transcript.
	      when use "genes", "trans", "onlyPr" together, it will extract transcript of "genes" with primary tag, besides all the "trans",
	      when use "genes" only, will extract all transcripts of genes without multiple mapping.
	      when use "trans" only, will extract only these transcripts.
	      when use "mmap", multiple mapping record will be involved.
	      when use "onlyPr", record without primary tag and not involved in the transcript list, will be skipped.
    Returns : return annoed blocks which are merged depend on the genes, trans, and region restriction.

=cut
sub new {
    my ($class, @args) = @_;
    my $self = {@args};
    bless $self, ref($class) || $class;

    if (   ( !exists $$self{db} )
        or ( !-e $$self{db} )
        or ( !exists $$self{codon} )
        or ( !-e $$self{codon} ) )
    {
        $self->throw("Error filename or too few args, need db, codon file to be available at least.");
    }

    $$self{codonseq} = readfa($$self{codon});
    my %open_args;
    if (exists $$self{region}) {
	$open_args{region} = $$self{region};
    }
    if (exists $$self{regbed}) {
	$open_args{regbed} = $$self{regbed};
    }
    if (exists $$self{genes}) {
	if (ref($$self{genes}) eq 'HASH') {
	    $open_args{genes} = $$self{genes};
	}
	else {
	    open (GENE, $$self{genes}) or $self->throw("$$self{genes} : $!");
	    my %genes = map { s/\s+//g; $_ => 1 } <GENE>;
	    close GENE;
	    $open_args{genes} = \%genes;
	}
    }
    if (exists $$self{trans}) {
	if (ref($$self{trans}) eq 'HASH') {
	    $open_args{trans} = $$self{trans};
	}
	else {
	    open (TRAN, $$self{trans}) or $self->throw("$$self{trans} : $!");
	    my %trans = map { s/\s+//g; $_ => 1 } <TRAN>;
	    close TRAN;
	    $open_args{trans} = \%trans;
	}
    }
    if (exists $$self{onlyPr}) {
	$open_args{onlyPr} = $$self{onlyPr};
    }
    if (exists $$self{mmap}) {
	$open_args{mmap} = $$self{mmap};
    }

    return $self->load_anno(%open_args);
}

sub readfa {
    my $file = shift;
    open (FAS, "zcat -f $file |") or confess "$file: $!";
    local $/ = ">";
    my %seqs = ();
    while (<FAS>) {
	s/[\s>]+$//g;
	next if (/^\s*$/);
	my $hd = $1 if (s/^(\S+)[^\n]*\n//);
	s/\s+//g;
	$seqs{$hd} = $_;
    }
    close FAS;
    return \%seqs;
}

=head2 load_anno

    About   : load all needed annotation infomation into memory for multi-proc annotation
    Usage   : $beda->load_anno( region => "chr20:1234567-1234568", trans => \%trans );
	     or: my $rlocal_db = $beda->load_anno( region => "chr20:1234567-1234568", local => 1 );
    Args    : using %args to override new methods args:
	      region, regbed, genes, trans
	      the tag local can return a local annotate db without overwrite the class entry's.
    Returns : BedA entry with full merged annotation blocks, or a localized merged anno db.

=cut
sub load_anno {
    my ( $self, %args ) = @_;
    my $cmd = "<$$self{db}";
    my $tabix_args = qq['$$self{db}'];
    if ( -e $$self{db} && $$self{db} =~ /\.gz/i ) {
        if ( exists $args{region} and defined $args{region} ) {
	    $tabix_args .= qq[ '$args{region}'];
            $cmd = "tabix $tabix_args |";
        }
	elsif ( exists $args{regbed} and defined $args{regbed} and -e $args{regbed} ) {
	    $tabix_args = '-B '.$tabix_args.qq[ '$args{regbed}'];
            $cmd = "tabix $tabix_args |";
	}
        else {
            $cmd = "zcat -f '$$self{db}' |";
        }
    }
    open( ANNO, $cmd ) or $self->throw("$cmd: $!");

    # trans filter is always be ahead of genes
    my $prTag = (exists $args{onlyPr}) ? 1 : 0;
    my $mmapTag = (exists $args{mmap}) ? 1 : 0;
    my $geneTag = (exists $args{genes}) ? 1 : 0;
    my $tranTag = (exists $args{trans}) ? 1 : 0;
    my $geneList = $args{genes} if ($geneTag);
    my $tranList = $args{trans} if ($tranTag);
    my %pureTran = map { s/\-\d+$//; $_ => 1 } keys %$tranList if ($tranTag);

    my %annodb = ();
    while (<ANNO>) {
	s/\s+$//;
	my ($chr, $start, $stop, $annostr) = split(/\t/);
	my @annos = split(/; /, $annostr);
	
	#		 0-based	1-based
	my %ent = (sta => $start, sto => $stop);
	foreach my $anno_ent (@annos) {
	    my @cont = split(/\|/, $anno_ent);
	    # no multiple mapping transcripts if !$mmapTag
            next
              if (  $cont[0] =~ /\-\d+$/
                and !$mmapTag
                and ( !$tranTag or !exists( $$tranList{ $cont[0] } ) ) );
            next
              if (  $prTag
                and ( $cont[10] ne 'Y' )
                and ( !$tranTag or !exists( $$tranList{ $cont[0] } ) ) );
	    my $ori_cont = $cont[0];
            $cont[0] =~ s/\-\d+$//;
            next
              if (
                    $geneTag
                and !exists( $$geneList{ $cont[1] } )
                and (  !$tranTag
                    or !exists( $pureTran{ $cont[0] } ) )
              );
            next
              if (
                $tranTag
                and ( !exists( $pureTran{ $cont[0] } )
                    or ( !$mmapTag and !exists $$tranList{$ori_cont} ) )
                and ( !$geneTag or !exists( $$geneList{ $cont[1] } ) )
              );

	    # NM_152486.2|SAMD11|+|IC7|IVS8|949|950|869|870|829|Y|0
            my $ofst = $1
              if ( $anno_ent =~ s/\|(\d+)$// )
              or $self->throw("db format error: [$anno_ent]");
            $ent{annos}{$anno_ent} = $ofst;
	}

	next if (!exists $ent{annos});
	push (@{$annodb{$chr}}, {%ent});
    }
    close ANNO;

    if (exists $args{local} and $args{local}) {
	return $self->region_merge( db => \%annodb ) if ($prTag or $mmapTag or $geneTag or $tranTag);
	return \%annodb;
    }
    else {
	$self->throw("Error: empty available-db!") if (0 == keys %annodb);

	$$self{annodb} = \%annodb;
	return $self->region_merge() if ($prTag or $mmapTag or $geneTag or $tranTag);
	return $self;
    }
}

=head2 write_using

    About   : write the current using database information to files
    Usage   : $beda->write_using( $file, $type );
    Args    : file gives the filename to output, and type is one of the following:
		g:  gene symbol list
		t:  transcript acc.ver list
		b:  standard bed format of only exon region
		a:  BED +1 format, exon region with '+1' annotation, 
		    oneline per one exon for one transcript, this
		    format allow redundancy exists, and not sorted by
		    chromosomal coordinates, but by the transcript acc.ver
		    and the exon number, this file is mainly for 
		    transcript originated statistics for only exon.

=cut
sub write_using {
    my ($self, $file, $type) = @_;
    open (F, ">", $file) or $self->throw("$file: $!");
    my (%genes, %trans, @beds, %anno) = ();
    foreach my $chr (sort keys %{$$self{annodb}}) {
	foreach my $bedent (@{$$self{annodb}{$chr}}) {
	    my $exOpt = 0;
	    foreach my $annoblk (keys %{$$bedent{annos}}) {
		# NM_152486.2|SAMD11|+|IC7|IVS8|949|950|869|870|829|Y
		my @info = split(/\|/, $annoblk);

		given ($type) {
		    when ('g') {
			$genes{$info[1]} = 1;
		    }
		    when ('t') {
			$trans{$info[0]} = 1;
		    }
		    when ('b') {
			$exOpt = 1 if ($info[4] =~ /^EX/);
		    }
		    when ('a') {
			if ( $info[4] =~ /^EX/) {
			    if (   !exists( $anno{ $info[0] } )
				or !exists( $anno{ $info[0] }{ $info[4] } ) )
			    {
				$anno{ $info[0] }{ $info[4] }{chr} = $chr;
				$anno{ $info[0] }{ $info[4] }{sta} = $$bedent{sta};
				$anno{ $info[0] }{ $info[4] }{sto} = $$bedent{sto};
				$anno{ $info[0] }{ $info[4] }{blk} = $annoblk;
			    }
			    elsif ($$bedent{sto} > $anno{ $info[0] }{ $info[4] }{sto}) {
				$anno{ $info[0] }{ $info[4] }{sto} = $$bedent{sto};
			    }
			    else {
				$self->throw("Error: bad db, non-departed beds, or no-sort!");
			    }
			}
		    }
		    default { $self->throw("type: $type not regconized."); }
		}
	    }
	    push (@beds, [ $chr, $$bedent{sta}, $$bedent{sto} ]) if ($type eq 'b' and $exOpt);
	}
    }

    given ($type) {
	when ('g') {
	    say F (sort keys %genes);
	}
	when ('t') {
	    say F (sort keys %trans);
	}
	when ('b') { # merge bed regions due to pre-sorted db.
	    if (0 == @beds) {
		$self->warn("no region in db.\n");
		return;
	    }
	    my ($pre_chr, $pre_sta, $pre_sto) = @{$beds[0]};
	    for (my $i = 1; $i < @beds; $i++) {
		my ($cur_chr, $cur_sta, $cur_sto) = @{$beds[$i]};
		if (($cur_chr ne $pre_chr) or ($cur_sta > $pre_sto)) {
		    say F join("\t", $pre_chr, $pre_sta, $pre_sto);
		    ($pre_chr, $pre_sta, $pre_sto) = ($cur_chr, $cur_sta, $cur_sto);
		}
		elsif ($cur_sta == $pre_sto) {
		    $pre_sto = $cur_sto;
		}
		else {
		    $self->throw("Error: bad db, non-departed beds, or no-sort!");
		}
	    }
	    say F join("\t", $pre_chr, $pre_sta, $pre_sto);
	}
	when ('a') {
	    foreach my $nm (sort keys %anno) {
		foreach my $ex (sort exsort keys %{$anno{$nm}}) {
		    say F join("\t", $anno{$nm}{$ex}{chr}, $anno{$nm}{$ex}{sta},
			$anno{$nm}{$ex}{sto}, $anno{$nm}{$ex}{blk});
		}
	    }
	}
	default { $self->throw("type: $type not regconized."); }
    }
    close F;
    return;
}

sub exsort {
    my ($sym, $anum, $bnum);
    if ($a =~ /^EX([\+\-\*]?)(\d+)$/) {
	$sym = $1;
	$anum = $2;
    }
    if ($b =~ /^EX[\+\-\*]?(\d+)$/) {
	$bnum = $1;
    }
    confess "ExIn number format error. [$a, $b]" if (!defined $anum or !defined $bnum);
    if (!defined $sym or $sym !~ /\-/) {
	$anum <=> $bnum;
    }
    else {
	$bnum <=> $anum;
    }
}

=head2 region_merge

    About   : merge consecutive same-entries regions
    Usage   : $beda->region_merge();
	     or:
		my $rAnnodb = $beda->region_merge( db => \%anndb );
    Args    : using 'db' tag to locally merge regions,
	      without updating the db of class entry.

=cut
sub region_merge {
    my ($self, %args) = @_;

    my $radb = (exists $args{db}) ? $args{db} : $$self{annodb};

    foreach my $chr (keys %$radb) {
	my @annoents = @{$$radb{$chr}};
	my $oricount = scalar (@annoents);
	next if ($oricount <= 1);

	# merge from bottom to top
	my $curid = $oricount - 1;
	ANNOENT:while ($curid > 0) {
	    my $latter = $annoents[$curid];
	    my $former = $annoents[$curid - 1];
	    $curid --;
	    
	    next if ($$latter{sta} != $$former{sto});
	    my @former_ann = keys %{$$former{annos}};
	    my $formern = scalar @former_ann;
	    next if ($formern != (scalar keys %{$$latter{annos}}));
	    foreach my $ann (@former_ann) {
		next ANNOENT if (!exists $$latter{annos}{$ann});
	    }

	    # splice the latter one and correct the former one
	    $$former{sto} = $$latter{sto};
	    splice (@annoents, ($curid + 1), 1);
	}

	$$radb{$chr} = [@annoents] if ($oricount > (scalar @annoents));
    }
    if (exists $args{db}) {
	return $radb;
    }
    else {
	return $self;
    }
}

=head2 anno

    About   : Annotate single short variation by annotation db.
    Usage   : my $rAnnoRst = $beda->anno( 'chr20', 1234567, 'TG', 'TGGG' );
    Args    : chromosome id, chromosome position, reference, alternative.
    Returns : a hash ref of annotation informations, see pairanno().

=cut
sub anno {
    my ($self, $chr, $pos, $ref, $alt) = @_;
    my $var = parse_var($chr, $pos, $ref, $alt);
    return $self->varanno($var);
}

=head2 varanno

    About   : generate all the needed annotation for a var entry
    Usage   : my $rAnnoRst = $beda->varanno($var, genes=>\%genes, trans=>\%trans);
    Args    : see parse_var(), and load_anno()
    Returns : see pairanno()

=cut
sub varanno {
    my ($self, $var) = @_;

    if (!exists $$var{sel}) {	# no position annotation
	my $pos1 = get_anno_pos($var);
	if ( $pos1 >  0 ) {
	    my $hit_db = 0;
	    foreach my $rbed (@{$$self{annodb}{$$var{chr}}}) {
		next if ($$rbed{sto} < $pos1);
		last if ($$rbed{sta} >= $pos1);
		if (!exists $$rbed{detail}) {
		    foreach my $annoblk (keys %{$$rbed{annos}}) {
			my $offset1 = $$rbed{annos}{$annoblk};
			my ($tid, $ranno) = parse_annoent($annoblk);
			$$rbed{detail}{$tid} = $ranno;
			$$rbed{detail}{$tid}{offset} = $offset1;
		    }
		}
		$hit_db = 1;
		my $ofst1 = $pos1 - $$rbed{sta} - 1;
		my %oneline_cPos;
		foreach my $t (sort keys %{$$rbed{detail}}) {
		    $oneline_cPos{$pos1}{$t} = in_reg_cPos($$rbed{detail}{$t}, $ofst1);
		}
		$var = $self->select_position($var, \%oneline_cPos);
		last;
	    }
	    if ($hit_db == 0) {
		return { var => $var, info => {} };
	    }
	}
	else {
	    $var = $self->select_position($var, {});
	}
    }

    my %ret_anno = ( var => $var, info => {} );
    if (exists $$var{sel} and exists $$var{sel}{std}) {
	my %infos = ();
	foreach my $sel (@{$$var{sel}{std}}) {
	    my ($tid, $ranno) = $self->pairanno($var, $sel);
	    $infos{$tid} = $ranno;
	}
	$ret_anno{info} = \%infos;
    }
    return \%ret_anno;
}

=head2 pairanno

    About   : get infomations from select_position(), (std grp) besides function,
	      codon-change, flank_region position, polar and HGVS strings:
		  for coding gene region gives c. HGVS and p. HGVS
		  for non-coding RNA gives r. HGVS, 
		  for non-available gives '.'.
	      function of mutation, follow the definitions mainly from dbSNP:
		  for intergenic or non-target region: using '.'
		  utr-3, utr-5, ncRNA, 
		  altstart, missense, nonsense, coding-synon, 
		  intron, splice-3, splice-5,
		  cds-indel, frameshift, stop-gain, stop-loss 
		  unknown for 'c.?'
		  abnormal-inseq-stop for non-end stop codon in refseq.
		  abnormal-intron for too short intron (<4bp) due to mapping deduced.
	      Sub-region and Exon/Intron number infomation
    Usage   : my $h_annos = $beda->pairanno($var, $sel);
    Args    : a var entry and a selected anno
    Returns : ( $tid, { c => c.HGVS, p => p.HGVS,  cc => condon-change,
			r => region, exin => exin, func => func,  polar => pol, 
			flanks => { l => left_region, r => right_region, strd => [+-] }
	              } )

=cut
sub pairanno {
    my ($self, $var, $sel) = @_;

    # annotated by paired (sta, sto) cPos
    my ($tid, $cL, $cR) = @$sel;
    return ( $tid, {} ) if (!exists $$cL{cpos});

    my ( $ref, $alt, $pos, $rl, $al ) =
      ( $$var{ref}, $$var{alt}, $$var{pos}, $$var{reflen}, $$var{altlen} );
    if (exists $$var{'+'} and exists $$var{'+'}{p}) { # for diff strand position cases
	$ref = $$var{$$cL{strd}}{r};
	$alt = $$var{$$cL{strd}}{a};
	$pos = $$var{$$cL{strd}}{p};
	$rl  = $$var{$$cL{strd}}{rl};
	$al  = $$var{$$cL{strd}}{al};
    }

    if ($rl > 1 or $al > 1) { # get the real diff base in ref
	$ref = substr($ref, 1); 
	$alt = substr($alt, 1);
	$rl --;
	$al --;
	$pos ++;
    }

    # get flank region string ready for fetchseq
    my $flanks = get_flanks($$var{chr}, $pos, $rl);
    $$flanks{strd} = $$cL{strd};

    my $strandopt = ($$cL{strd} eq '-') ? 0 : 1;
    my $t_ref = ($strandopt) ? uc($ref) : rev_comp(uc($ref));
    my $t_alt = ($strandopt) ? uc($alt) : rev_comp(uc($alt));

    my ($cHGVS, $pHGVS, $cc, $reg, $exin, $func, $polar, $bc_cHGVS) = ('.') x 8;

    if ($$cL{cpos} =~ /\?/) {
	return (
	    $tid,
	    {
		c      => '.',
		p      => '.',
		cc     => '.',
		r      => $$cL{reg},
		exin   => $$cL{exin},
		func   => 'unknown',
		polar  => '.',
		bc     => $bc_cHGVS,
		flanks => $flanks
	    }
	);
    }

    my $query_tid = $tid;
    $query_tid =~ s/\-\d+$//;
    if (!defined $cR) {  # snv or 1bp deletion or 1bp repeat
	$cR = $cL;
    }
    if (!$strandopt) { # swap LR is '-'
	my $tmp = $cL;
	$cL = $cR;
	$cR = $tmp;
    }
    # get cPos without c./r. 
    my $cP_L = substr($$cL{cpos}, 2);
    my $cP_R = substr($$cR{cpos}, 2);
    my $pre_L = substr($$cL{cpos}, 0, 2);
    my $pre_R = substr($$cR{cpos}, 0, 2);

    my $coding_max = int(length($$self{codonseq}{$query_tid})/3) if (exists $$self{codonseq}{$query_tid});

    given ($$var{guess}) {

	when ('snv') { # only one position
	    $cHGVS = $$cL{cpos}.$t_ref.'>'.$t_alt; # c.123A>G
	    $reg = $$cL{reg};
	    $exin = $$cL{exin};
	    ($func, $cc, $pHGVS, $polar) = $self->parse_cPos($tid, $$cL{cpos}, $t_alt);
	}
	when ('ins') {

	    if ($pre_L eq $pre_R) {
		$cHGVS = $$cL{cpos}.'_'.$cP_R.'ins'.$t_alt;
	    }
	    else {
		$cHGVS = $$cL{cpos}.'_'.$$cR{cpos}.'ins'.$t_alt;
		return (
		    $tid,
		    {
			c      => $cHGVS,
			p      => $pHGVS,
			cc     => $cc,
			r      => $reg,
			exin   => $exin,
			func   => $func,
			polar  => $polar,
			bc     => $bc_cHGVS,
			flanks => $flanks
		    }
		);
	    }

	    if ($$cL{reg} eq $$cR{reg}) { # non-edge insertion

		$reg  = $$cL{reg};
		$exin = $$cL{exin};

		given ($reg) {
		    when (/^C/) { # coding region
			( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_L, $t_alt, 1);
		    }
		    when (/^I/) {
			if ($$cL{bd} =~ /0/ and $$cR{bd} =~ /0/) {
			    $func = 'intron';
			}
			elsif ($$cL{bd} =~ /b/) {
			    $func = 'splice-5';
			}
			elsif ($$cR{bd} =~ /B/) {
			    $func = 'splice-3';
			}
			else {
			    $func = 'abnormal-intron';
			}
		    }
		    when (/^5/) {
			$func = 'utr-5';
		    }
		    when (/^3/) {
			$func = 'utr-3';
		    }
		    when (/^R/) {
			$func = 'ncRNA';
		    }
		    default { $self->throw("unexpected region string [$reg]"); }
		}
	    }
	    else { # insertion on the edge of region
		$reg = $$cL{reg}.'-'.$$cR{reg};
		$exin = $$cL{exin}.'-'.$$cR{exin};

		my $indicator = substr($$cL{reg},0,1).substr($$cR{reg},0,1);
		given ($indicator) {
		    when ('IC') { # 5' coding exon insertion
			( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_R, $t_alt, -1);
		    }
		    when ('CI') { # 3' conding exon insertion
			( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_L, $t_alt, 1);
		    }
		    when (/5/) {
			$func = 'utr-5';
		    }
		    when (/3/) {
			$func = 'utr-3';
		    }
		    when (/R/) {
			$func = 'ncRNA';
		    }
		    default { $self->throw("unexpected indicator [$indicator]"); }
		}
	    }
	}
	when ([ 'del', 'delins' ]) {

            if ( $pre_L eq $pre_R ) {
                if ( $$cL{cpos} eq $$cR{cpos} ) {
                    $cHGVS = $$cL{cpos} . 'del' . $t_ref;
                }
                else {
                    $cHGVS = $$cL{cpos} . '_' . $cP_R . 'del' . $t_ref;
                }
		
		$cHGVS .= 'ins' . $t_alt if (/delins/);

            }
            else { # $pre_L ne $pre_R
                $cHGVS = $$cL{cpos} . '_' . $$cR{cpos} . 'del' . $t_ref;
		$cHGVS .= 'ins' . $t_alt if (/delins/);

		if (($pre_L eq '--') and ($pre_R eq '++')) { # whole transcript is missing
		    $pHGVS = 'p.0?';
		    $func  = 'knockout'; # which may not be exists in variation calling result.
		}
		elsif ($pre_L eq '--') { # the 5' end of transcript is missing
		    given ($$cR{reg}) {
			when (/^5/) {
			    $func = 'utr-5';
			}
			when (/^C/) {
			    $pHGVS = 'p.0?';
			    $func  = 'altstart';
			}
			when (/^3/) {
			    $pHGVS = 'p.0?';
			    $func  = 'knockout';
			}
			when (/^R/) {
			    $func  = 'ncRNA';
			}
			when (/^I/) {
			    if ($pre_R eq 'r.') {
				$func = 'ncRNA';
			    }
			    elsif ($cP_R =~ /^\-|^1\-\d/) {
				$func = 'utr-5';
			    }
			    elsif ($cP_R =~ /^\d+/) {
				$pHGVS = 'p.0?';
				$func  = 'altstart';
			    }
			    elsif ($cP_R =~ /^\*/) {
				$pHGVS = 'p.0?';
				$func  = 'knockout';
			    }
			    else {
				$self->throw("unexpected cP_R [$cP_R]");
			    }
			}
			default {
			    $self->throw("unexpected region [$$cR{reg}]");
			}
		    }
		}
		else { # $pre_R eq '++' the 3'end of transcript is missing.
		    given ($$cL{reg}) {
			when (/^5/) {
			    $pHGVS = 'p.0?';
			    $func = 'knockout';
			}
			when (/^C/) {
			    my $pP_L = aaNum($cP_L);
			    my $aa_L = get_aa($$self{codonseq}, $query_tid, $pP_L);
			    if ($aa_L ne '') {
				$pHGVS = 'p.'.$aa_L.$pP_L.'_*'.$coding_max.'del';
				$func  = 'stop-loss';
			    }
			    else {
				$func  = 'unknown';
			    }
			}
			when (/^3/) {
			    $func  = 'utr-3';
			}
			when (/^R/) {
			    $func  = 'ncRNA';
			}
			when (/^I/) {
			    if ($pre_L eq 'r.') {
				$func = 'ncRNA';
			    }
			    elsif ($cP_L =~ /^\-/) {
				$pHGVS = 'p.0?';
				$func  = 'knockout';
			    }
			    elsif ($cP_L =~ /^\*/) {
				$func = 'utr-3';
			    }
			    elsif ($cP_L =~ /^(\d+)\+/) {
				if (exists $$self{codonseq}{$query_tid}) {
				    if ($1 == length($$self{codonseq}{$query_tid})) {
					$func = 'utr-3';
				    }
				    else {
					my $pP_L = aaNum($1+1);
					if ($pP_L < $coding_max) {
					    my $aa_L = get_aa($$self{codonseq}, $query_tid, $pP_L);
					    $pHGVS = 'p.'.$aa_L.$pP_L.'_*'.$coding_max.'del';
					}
					else {
					    $pHGVS = 'p.'.$coding_max.'del*';
					}
					$func  = 'stop-loss';
				    }
				}
				else {
				    $func = 'unknown';
				}
			    }
			    elsif ($cP_L =~ /^(\d+)\-/) {
				if (exists $$self{codonseq}{$query_tid}) {
				    my $pP_L = aaNum($1);
				    if ($pP_L < $coding_max) {
					my $aa_L = get_aa($$self{codonseq}, $query_tid, $pP_L);
					$pHGVS = 'p.'.$aa_L.$pP_L.'_*'.$coding_max.'del';
				    }
				    else {
					$pHGVS = 'p.'.$coding_max.'del*';
				    }
				    $func  = 'stop-loss';
				}
				else {
				    $func = 'unknown';
				}
			    }
			    else {
				$self->throw("unexpected cP_L [$cP_L]");
			    }
			}
			default {
			    $self->throw("unexpected region [$$cR{reg}]");
			}
		    }
		}

		return (
		    $tid,
		    {
			c      => $cHGVS,
			p      => $pHGVS,
			cc     => $cc,
			r      => $reg,
			exin   => $exin,
			func   => $func,
			polar  => $polar,
			bc     => $bc_cHGVS,
			flanks => $flanks
		    }
		);
		
            }

	    my $indicator = substr($$cL{reg},0,1).substr($$cR{reg},0,1);

            if ( $$cL{reg} eq $$cR{reg} ) {
                $reg  = $$cL{reg};
                $exin = $$cL{exin};

		given ($indicator) {
		    when ('CC') {
			($pHGVS, $func) = $self->get_aaDelInfo($query_tid, $cP_L, $cP_R, $t_alt);
		    }
		    when ('II') {
			if ($$cL{bd} =~ /1/) {
			    $func = 'abnormal-intron';
			}
			elsif ($$cL{bd} =~ /0/ and $$cR{bd} =~ /0/) {
			    $func = 'intron';
			}
			elsif ($$cL{bd} =~ /b/ and $$cR{bd} =~ /B/) {
			    $func = 'splice';
			}
			elsif ($$cL{bd} =~ /b/) {
			    $func = 'splice-5';
			}
			elsif ($$cR{bd} =~ /B/) {
			    $func = 'splice-3';
			}
		    }
		    when ('55') {
			$func = 'utr-5';
		    }
		    when ('33') {
			$func = 'utr-3';
		    }
		    when ('RR') {
			$func = 'ncRNA';
		    }
		    default {
			$self->throw("unexpected region [$$cL{reg}]");
		    }
		}
            }
            else { # multiple region del/delins
                $reg  = $$cL{reg} . '-' . $$cR{reg};
                $exin = $$cL{exin} . '-' . $$cR{exin};
		given ($indicator) {
		    when ('CC') {
			($pHGVS, $func) = $self->get_aaDelInfo($query_tid, $$cL{cpos}, $$cR{cpos}, $t_alt);
		    }
		    when ('55') {
			$func = 'utr-5';
		    }
		    when ('33') {
			$func = 'utr-3';
		    }
		    when ('RR') {
			$func = 'ncRNA';
		    }
		    when ('II') {
			if ($$cL{bd} =~ /1/) {
			    $func = 'abnormal-intron';
			}
			elsif ($$cL{bd} =~ /0/ and $$cR{bd} =~ /0/) {
			    if ($pre_L eq 'r.') {
				$func = 'ncRNA';
			    }
			    else {
				my ($cds_L, $cds_R);
				if ($cP_L =~ /^(\d+)\-/) {
				    $cds_L = $1;
				}
				elsif ($cP_L =~ /^(\d+)\+/) {
				    $cds_L = $1 + 1;
				}

				if ($cP_R =~ /^(\d+)\+/) {
				    $cds_R = $1;
				}
				elsif ($cP_R =~ /^(\d+)\-/) {
				    $cds_R = $1 - 1;
				}

				($pHGVS, $func) = $self->get_aaDelInfo($query_tid, $cds_L, $cds_R);
			    }
			}
			else { # ($$cL{bd} =~ /b/i or $$cR{bd} =~ /b/i) the splice cases is too complex to parse
			    # give splice to indicate either or both of side are at splice site.
			    $func = 'splice';
			}
		    }
		    when ('5C') {
			$pHGVS = 'p.0?';
			$func  = 'altstart';
		    }
		    when ('C3') {
			my $pP_L = aaNum($cP_L);
			my $aa_L = get_aa($$self{codonseq}, $query_tid, $pP_L);
			if ($aa_L eq '') {
			    $func = 'unknown';
			}
			else {
			    if ($pP_L == $coding_max) {
				$pHGVS = 'p.'.$coding_max.'del*';
			    }
			    else {
				$pHGVS = 'p.'.$aa_L.$pP_L.'_*'.$coding_max.'del';
			    }
			    $func  = 'stop-loss';
			}
		    }
		    when ('53') {
			$pHGVS = 'p.0?';
			$func  = 'knockout';
		    }
		    default { # for other complex cases
			if ($$cR{exin} =~ /C1/) {
			    $pHGVS = 'p.0?';
			    $func = 'splice-3';
			}
			elsif ($$cL{exin} =~ /C\d+E/) {
			    my $pP_L = aaNum($cP_L);
			    my $aa_L = get_aa($$self{codonseq}, $query_tid, $pP_L);
			    if ($aa_L ne '') {
				if ($pP_L == $coding_max) {
				    $pHGVS = 'p.'.$coding_max.'del*';
				}
				else {
				    $pHGVS = 'p.'.$aa_L.$pP_L.'_*'.$coding_max.'del';
				}
			    }
			    $func = 'splice-5';
			}
			else {
			    $func = 'splice';
			}
		    }
		}
            }
	}
	# when meet repeats, the input selection must come from
	# the standard group: (repeat element cds position), 
	# the backward compatible cds position is used to calculate
	# the bcHGVS, pHGVS and give the function prediction as well.
	# the standard cds position is only for the cHGVS.
	when ('rep') {
	    my $t_rep = ($$cL{strd} eq '+') ? $$var{rep} : rev_comp($$var{rep});
	    if ($$var{ref_cn} == 1 and $$var{alt_cn} == 2) { # duplication
		if ($$var{replen} > 1) {
		    if ($pre_L ne $pre_R) {
			$cHGVS = $$cL{cpos}.'_'.$$cR{cpos}.'dup'.$t_rep;
		    }
		    else {
			$cHGVS = $$cL{cpos}.'_'.$cP_R.'dup'.$t_rep;
		    }
		}
		else {
		    $cHGVS = $$cL{cpos}.'dup'.$t_rep;
		}
	    }
	    else { # simple repeats
		$cHGVS = $$cL{cpos}.$t_rep.'['.$$var{ref_cn}.'>'.$$var{alt_cn}.']';
	    }

	    # if backward compatible selection exists, the pHGVS and function
	    # will be assigned to the variation.
	    if (exists $$var{sel}{bc}) {
		my $sel_for_bc;
		foreach my $bc_sel (@{$$var{sel}{bc}}) {
		    my ($bc_tid, $bc_cL, $bc_cR) = @{$bc_sel};
		    next if ($bc_tid ne $tid);
		    $sel_for_bc = $bc_sel;
		    last;
		}
		if (defined $sel_for_bc) {
		    my %bc_var; # create a backward compatible var
		    my $rbc_var_info = $$var{$$cL{strd}};
		    $bc_var{chr}   = $$var{chr};
		    $bc_var{guess} = guess_type($$rbc_var_info{brl}, $$rbc_var_info{bal});
		    $bc_var{pos}   = $$rbc_var_info{bp};
		    $bc_var{ref}   = $$rbc_var_info{br};
		    $bc_var{alt}   = $$rbc_var_info{ba};
		    $bc_var{reflen} = $$rbc_var_info{brl};
		    $bc_var{altlen} = $$rbc_var_info{bal};

                    my ( $ret_tid, $ret_rbc_annos ) =
                      $self->pairanno( \%bc_var, $sel_for_bc );
                    $pHGVS         = $$ret_rbc_annos{p};
                    $func          = $$ret_rbc_annos{func};
                    $reg           = $$ret_rbc_annos{r};
                    $exin          = $$ret_rbc_annos{exin};
                    $bc_cHGVS      = $$ret_rbc_annos{c};
                    $flanks        = $$ret_rbc_annos{flanks};
                    $$flanks{strd} = $$cL{strd};
		}
	    }
	}
	default { $self->throw("unexpected type of variation [$$var{guess}]."); }
    }

    return (
        $tid,
        {
            c      => $cHGVS,
            p      => $pHGVS,
            cc     => $cc,
            r      => $reg,
            exin   => $exin,
            func   => $func,
            polar  => $polar,
	    bc     => $bc_cHGVS,
            flanks => $flanks
        }
    );
}


# get aa sequence from nucl string
sub get_aaseq {
    my $new_codon = shift;
    my $new_aa = "";
    for (my $i = 0; $i < length($new_codon); $i+=3) {
	my $add_codon = substr($new_codon, $i, 3);
	my $add_aa = (exists $codon3{$add_codon}) ? $codon3{$add_codon} : "";
	if ($add_aa ne '') {
	    $new_aa .= $add_aa;
	}
	last if ($add_aa eq '' or $add_aa eq '*');
    }
    return $new_aa;
}


=head2 parse_cPos

    About   : parse the cds position of snv variation
    Usage   : my ($func, $cc, $pHGVS, $polar) = $beda->parse_cPos($query_tid, $cpos, $t_alt);
    Args    : query_tid is transcript id with out tail [-n], cpos is cds pos in hgvs c./r. format.
	      t_alt are strand specific alt char.
    Returns : the function code, codon-change, 'p.' format HGVS string, polar-change. see pairanno()

=cut
sub parse_cPos {
    my ($self, $query_tid, $cpos, $t_alt) = @_;
    my ($func, $cc, $pHGVS, $polar) = ('.') x 4;
    given ($cpos) {
	when (/^c\.(\d+)$/) { # coding region
	    my $cP = $1;
	    my ($codon, $aa, $pol, $frame) = get_codon($$self{codonseq}, $query_tid, $cP);
	    my $new_codon = $codon;
	    if ($aa ne '') {
		substr($new_codon, $frame, 1, $t_alt);
		my $new_aa = (exists $codon3{$new_codon}) ? $codon3{$new_codon} : "";
		my $new_pol = (exists $Polar{$new_aa}) ? $Polar{$new_aa} : "";
		$func = $self->get_aafunc($aa, $new_aa, $query_tid, $cP);
		$pHGVS = get_aaHGVS($aa, $new_aa, $func, $cP);
		$cc = $codon.'=>'.$new_codon;
		$polar = ($pol eq $new_pol) ? '.' : $pol.'=>'.$new_pol;
		if ($aa eq '*' and $cP + 3 <= length($$self{codonseq}{$query_tid})) {
		    $func = 'abnormal-inseq-stop';
		}
	    }
	}
	when (/^[cr]\.[\-\*]?\d+[\-\+]\d+/) { # intron region
	    if ($$cpos{bd} =~ /0/) {
		$func = 'intron';
	    }
	    elsif ($$cpos{bd} =~ /b/) {
		$func = 'splice-5';
	    }
	    elsif ($$cpos{bd} =~ /B/) {
		$func = 'splice-3';
	    }
	    else {
		$func = 'abnormal-intron';
	    }
	}
	when (/^c\.-\d+$/) { # 5' UTR
	    $func = 'utr-5';
	}
	when (/^c\.\*\d+$/) { # 3' UTR
	    $func = 'utr-3';
	}
	when (/^r\.\d+$/) { # ncRNA
	    $func = 'ncRNA';
	}
	default { $func = 'unknown'; }
    }
    return ($func, $cc, $pHGVS, $polar);
}

=head2 get_aaDelInfo
    
    About   : get pHGVS and function code from coding deletion or delins.
    Usage   : my ($pHGVS, $func) = $beda->get_aaDelInfo($qid, $cP_L, $cP_R, $insBase);
    Args    : query_tid see get_aaInsInfo(), cP_L and cP_R is the coding region deletion start and end.
	      insBase is the inserted bases if delins.
	      currently pHGVS prediction won't extend outside the cds region, so it will mis-predict
	      the represent of stop-codon and will issue a "stop-loss".
    Returns : the protein level HGVS and the function code.

=cut
sub get_aaDelInfo {
    my ($self, $query_tid, $cP_L, $cP_R, $t_alt) = @_;
    $t_alt ||= '';
    my ($pHGVS, $func) = ( '.', '.' );
    my ($codon_L, $aa_L, $pol_L, $frame_L) = get_codon($$self{codonseq}, $query_tid, $cP_L);
    if ($aa_L eq '') { # no cds seq
	return ( $pHGVS, 'unknown' );
    }
    my $pP_L = aaNum($cP_L);
    my $remain_L = substr($codon_L, 0, $frame_L);
    my ($codon_R, $aa_R, $pol_R, $frame_R, $pP_R);
    if ($cP_L == $cP_R) { # 1bp deletion
	($codon_R, $aa_R, $pol_R, $frame_R) = ($codon_L, $aa_L, $pol_L, $frame_L);
	$pP_R = $pP_L;
    }
    else {
	($codon_R, $aa_R, $pol_R, $frame_R) = get_codon($$self{codonseq}, $query_tid, $cP_R);
	$pP_R = aaNum($cP_R);
    }

    my $todel_aa = "";
    my $todel_aalen = $pP_R - $pP_L + 1;
    for (my $i = $pP_L; $i <= $pP_R; $i++) {
	$todel_aa .= get_aa($$self{codonseq}, $query_tid, $i);
    }

    my $remain_R = substr($codon_R, ($frame_R + 1));
    my $fsOpt = (($cP_R - $cP_L + 1 - length($t_alt)) % 3 == 0) ? 0 : 1;

    my $total_remain   = $remain_L.$t_alt.$remain_R;
    my $tail_flank_len = (3 - (length($total_remain) % 3)) % 3;

    my $add_tail = "";
    if ($tail_flank_len > 0 and $cP_R < length($$self{codonseq}{$query_tid})-$tail_flank_len) {
	$add_tail = substr($$self{codonseq}{$query_tid}, $cP_R, $tail_flank_len);
    }
    $total_remain .= $add_tail;

    my $ins_aa = get_aaseq($total_remain);
    my $ins_aalen = (length($ins_aa) % 3 == 0) ? (length($ins_aa)/3) : (int(length($ins_aa)/3)+1);
    given ($ins_aa) {
	when (/^$/) { #frameshit or cds deletion without insertion 
	    if ($pP_L == 1) {

		$pHGVS = 'p.0?';
		$func  = 'altstart';
		return ( $pHGVS, $func );
		
	    }

	    if ($pP_L == $pP_R) {
                $pHGVS =
                  ($fsOpt)
                  ? 'p.' . $aa_L . $pP_L . 'fs*?'
                  : 'p.' . $aa_L . $pP_L . 'del';
	    }
	    else {
                $pHGVS =
                  ($fsOpt)
                  ? 'p.' . $aa_L . $pP_L . '_' . $aa_R . $pP_R . 'del'
                  : 'p.' . $aa_L . $pP_L . '_' . $aa_R . $pP_R . 'delfs*?';
	    }

            $func =
              ($fsOpt) ? 'frameshift'
              : ( ( $cP_R + 3 > length( $$self{codonseq}{$query_tid} ) )
                ? 'stop-loss'
                : 'cds-indel' );

	}
	when (/\*/) { # stop gain after delins
	    if ($aa_R eq '*') { # cds-indel
		my $local_ins_aa = $`; # $` prematch
		my $local_ins_aalen = $ins_aalen - 1;

		if ($cP_L == $cP_R) {
		    if ($local_ins_aa eq '') {

			$pHGVS = 'p.=';
			$func  = 'coding-synon';

		    }
		    else { # delins caused insertion
			my $aa_pre = get_aa($$self{codonseq}, $query_tid, ($pP_L - 1));

			$pHGVS = 'p.'. $aa_pre. ($pP_L - 1) . '_' . $aa_R . $pP_R . 'ins' . $local_ins_aa;
			$func  = 'cds-indel';
		    }
		}
		else { # for inconvinience of realign the protein del and ins, delins may be false predicting, except for head/tail align.
		    my $local_todel = $todel_aa;
		    $local_todel =~ s/\*$//;
		    my $local_todel_len = $todel_aalen - 1;

		    if ($local_todel eq $local_ins_aa) {

			$pHGVS = 'p.=';
			$func  = 'coding-synon';
			return ( $pHGVS, $func );

		    }

		    $func = 'cds-indel';

		    my ($local_pP_L, $local_pP_R, $local_aa_L, $local_aa_R, $local_toins_aa);
		    if ($local_todel =~ /^$local_ins_aa/) { # deletion is longer
			$local_pP_L = $pP_L + $local_ins_aalen;
			$local_pP_R = $pP_R - 1;
			$local_aa_L = get_aa($$self{codonseq}, $query_tid, $local_pP_L);
			if ($local_pP_L == $local_pP_R) {
			    $pHGVS = 'p.'.$local_aa_L.$local_pP_L.'del';
			}
			else {
			    $local_aa_R = get_aa($$self{codonseq}, $query_tid, $local_pP_R);
			    $pHGVS = 'p.'.$local_aa_L.$local_pP_L.'_'.$local_aa_R.$local_pP_R.'del';
			}
		    }
		    elsif ($local_todel =~ /$local_ins_aa$/) {
			$local_pP_R = $pP_R - 1 - $local_ins_aalen;
			$local_aa_R = get_aa($$self{codonseq}, $query_tid, $local_pP_R);
			if ($local_pP_R == $pP_L) {
			    $pHGVS = 'p.'.$aa_L.$pP_L.'del';
			}
			else {
			    $pHGVS = 'p.'.$aa_L.$pP_L.'_'.$local_aa_R.$local_pP_R.'del';
			}
		    }
		    elsif ($local_ins_aa =~ /^$local_todel/) { # insertion is longer
			$local_pP_L = $pP_R - 1;
			$local_toins_aa = $'; # $' post match
			$local_aa_L = get_aa($$self{codonseq}, $query_tid, $local_pP_L);
			
			$pHGVS = 'p.'.$local_aa_L.$local_pP_L.'_'.$aa_R.$pP_R.'ins'.$local_toins_aa;
		    }
		    elsif ($local_ins_aa =~ /$local_todel$/) {
			$local_pP_R = $pP_L + 1;
			$local_toins_aa = $`; # $` prematch
			$local_aa_R = get_aa($$self{codonseq}, $query_tid, $local_pP_R);

			$pHGVS = 'p.'.$aa_L.$pP_L.'_'.$local_aa_R.$local_pP_R.'ins'.$local_toins_aa;
		    }
		    else { # delins
			$local_pP_R = $pP_R - 1;
			$local_aa_R = get_aa($$self{codonseq}, $query_tid, $local_pP_R);

			$pHGVS = 'p.'.$aa_L.$pP_L.'_'.$local_aa_R.$local_pP_R.'delins'.$ins_aa;
		    }
		}
	    }
	    else { # $aa_R ne '*'
		$func  = 'stop-gain';
		if ($pP_L == $pP_R) {
		    if (length($ins_aa) < 4) {
			$pHGVS = 'p.'.$aa_L.$pP_L.$ins_aa;
		    }
		    else {
			$pHGVS = 'p.'.$aa_L.$pP_L.'delins'.$ins_aa;
		    }
		}
		else {
		    $pHGVS = 'p.'.$aa_L.$pP_L.'_'.$aa_R.$pP_R.'delins'.$ins_aa;
		}
	    }
	}
	default {
	    if ($todel_aa eq $ins_aa) { # frameshift or delins 

		$pHGVS = ($fsOpt) ? 'p.' . $aa_R . $pP_R . 'fs*?' : 'p.=';
		$func  = ($fsOpt) ? 'frameshift' : 'coding-synon';

	    }
	    elsif ($pP_L == 1 and substr($todel_aa, 0, 3) ne substr($ins_aa, 0, 3)) {

		$pHGVS = 'p.0?';
		$func  = 'altstart';

	    }
	    else {
		if ($pP_L == $pP_R) {
		    if (length($ins_aa) < 4) {
			$pHGVS = 'p.' . $aa_L . $pP_L . $ins_aa;
		    }
		    else {
			$pHGVS = 'p.' . $aa_L . $pP_L . 'delins' . $ins_aa;
		    }
		}
		else {
		    $pHGVS = 'p.' . $aa_L . $pP_L . '_' . $aa_R . $pP_R . 'delins' . $ins_aa;
		}
		$pHGVS .= 'fs*?' if ($fsOpt);

		$func =
		  ($fsOpt) ? 'frameshift'
		  : ( ( $cP_R + 3 > length( $$self{codonseq}{$query_tid} ) )
		    ? 'stop-loss'
		    : 'cds-indel' );

	    }
	}
    }
    return ( $pHGVS, $func );
}

=head2 get_aaInsInfo
    About   : get pHGVS and function code from coding insertion
    Usage   : my ($pHGVS, $func) = $beda->get_aaInsInfo($query_tid, $cP, $t_alt, $tp);
    Args    : query_tid is the tid for query refseq codon, [no -N tail]
	      cP is the anchored cds position
	      t_alt is the strand specific insertion string
	      tp gives the insert position of t_alt:
	      -1 for inserting before cP
	       1 for inserting after cP
    Returns : see pairanno()

=cut
sub get_aaInsInfo {
    my ($self, $query_tid, $cP, $t_alt, $tp) = @_;
    my ($pHGVS, $func) = ('.', '.');

    my ($codon, $aa, $pol, $frame) = get_codon($$self{codonseq}, $query_tid, $cP);
    if ($aa eq '') {
	return ( $pHGVS, 'unknown' );
    }

    my $pP = aaNum($cP);
    my $pP_pre = $pP - 1;
    my $pP_nex = $pP + 1;

    return ( $pHGVS, 'utr-5') if ( $tp < 0 and $frame == 0 and $pP_pre < 0 ); # insert before init codon
    return ( $pHGVS, 'utr-3') if ( $tp > 0 and $frame == 2 and ( $cP + 3 ) > length( $$self{codonseq}{$query_tid} ) ); # insert after stop codon

    my ($aa_pre, $aa_nex);
    $aa_pre = get_aa($$self{codonseq}, $query_tid, $pP_pre);
    $aa_pre = '?' if ($aa_pre eq '');
    $aa_nex = get_aa($$self{codonseq}, $query_tid, $pP_nex);
    $aa_nex = '?' if ($aa_nex eq '');

    if (($tp < 0 and $frame == 0) or ($tp > 0 and $frame == 2)) {
	my $ret_aa = get_aaseq($t_alt);

	# insert before anchor aa, transfer to insert after the previous aa.
	if ($tp < 0 and $frame == 0) { 
	    $aa = $aa_pre;
	    $aa_nex = $aa;
	    $pP = $pP_pre;
	    $pP_nex = $pP;
	}

	if ($ret_aa =~ /\*/) {

	    $pHGVS = 'p.'. $aa . $pP . '_' . $aa_nex. $pP_nex . 'ins' . $` .'*';    # $` prematch
	    $func  = 'stop-gain';

	}
	elsif ( length($t_alt) % 3 != 0 ) {

	    $pHGVS = 'p.'. $aa_nex. $pP_nex. 'fs*?';
	    $func  = 'frameshift';

	}
	else {

	    $pHGVS = 'p.'. $aa . $pP . '_' . $aa_nex. $pP_nex . 'ins' . $ret_aa;
	    $func  = 'cds-indel';

	}
    }
    else { # insert into anchor aa : p.Arg21delinsLeuTyr*, p.Arg21fs*?
	if ($tp < 0) {
	    $frame -= 1;  # transfer to insert after the previous frame
	}

	my $new_codon = $codon;
	substr($new_codon, ($frame+1), 0, $t_alt); # get alternate codon seq
	my $fsOpt = (length($new_codon) % 3 == 0) ? 0 : 1;

	my $new_aa = get_aaseq($new_codon);

	if ($pP == 1 and ($fsOpt or $new_aa !~ /$aa/)) {

	    $pHGVS = 'p.0?';
	    $func  = 'altstart';

	}
	elsif ($pP == 1 and $new_aa =~ /$aa$/) { # non-frameshift with inition code tail

	    $pHGVS = 'p.=';
	    $func  = 'coding-synon';

	}
	elsif ($pP == 1 and $new_aa =~ /$aa/) { # non-frameshift with insert a inition code
	    my $ins_aa = $';			    # $' is the post match, actually the insertion

            $pHGVS = 'p.' . $aa . $pP . '_' . $aa_nex . $pP_nex . 'ins' . $ins_aa;
	    $func  = 'cds-indel';

	}
	elsif ($new_aa !~ /\*/ and $fsOpt) { # non-stop-gain  with frameshift

	    $pHGVS = 'p.'.$aa.$pP.'fs*?';
	    $func = 'frameshift';

	}
	else {

	    if ($new_aa !~ /\*/) { # non-frameshift, non-stop-gain

		if ($new_aa !~ /^$aa/ and $new_aa !~ /$aa$/) {

		    $pHGVS = 'p.' . $aa . $pP . 'delins' . $new_aa;
		    $func  = ($aa eq '*') ? 'stop-loss' : 'cds-indel';

		}
		elsif ($new_aa =~ /^$aa/ ) {

		    $pHGVS = 'p.' . $aa . $pP . '_' . $aa_nex . $pP_nex . 'ins' . $'; # $' is the actually post insertion
		    $func  = 'cds-indel';

		}
		else {

		    $new_aa =~ /$aa$/;
		    $pHGVS  = 'p.' . $aa_pre . $pP_pre . '_' . $aa . $pP . 'ins' . $`; # $` is the actually prefix insertion
		    $func   = 'cds-indel';

		}

	    }
	    else { # $new_aa =~ /\*/ gain stop codon

		if ($new_aa =~ /^\*/ and $aa eq '*') {

		    $pHGVS = 'p.=';
		    $func  = 'coding-synon';

		}  
		elsif ($aa eq '*') {

		    $pHGVS = 'p.' . $aa_pre . $pP_pre . '_' . $aa . $pP . 'ins' . $`;  # $` prematch
		    $func  = 'cds-indel';

		}
		else {

		    $new_aa =~ /\*/;
		    $pHGVS = 'p.' . $aa . $pP . 'delins' . $`;			       # $` prematch
		    $func  = 'stop-gain';

		}
	    }
	}
    }

    if ( ( $cP + 3 ) <= length( $$self{codonseq}{$query_tid} )
        and $aa eq '*' )
    {
        $func = 'abnormal-inseq-stop';
    }

    return ( $pHGVS, $func );
}


# get pHGVS for aa change.
sub get_aaHGVS {
    my ($aa, $new_aa, $func, $cP) = @_;
    return '.' if ($aa eq '' or $new_aa eq '' or $func eq '');
    if ($func eq 'altstart') {
	return 'p.0?';
    }
    elsif ($func eq 'coding-synon') {
	return 'p.='; # HGVS suggest to not to give this kind of changing in pHGVS.
    }
    else {
	my $codon_num = aaNum($cP);
	return 'p.'.$aa.$codon_num.$new_aa;
    }
}

# get aa char by aa position
sub get_aa {
    my ($rseq, $qtid, $aa_num) = @_;
    my $start = ($aa_num - 1) * 3;
    if ($start >= (length($$rseq{$qtid}) - 2)) {
	return "";
    }
    else {
	my $codon = uc(substr($$rseq{$qtid}, $start, 3));
	return ((!exists $codon3{$codon}) ? "" : $codon3{$codon});
    }
}

sub aaNum {
    my $cP = shift;
    return (($cP % 3 == 0) ? ($cP / 3) : (int($cP/3)+1));
}

# get function code from aa change
sub get_aafunc {
    my ($self, $aa, $new_aa, $query_tid, $cP) = @_;
    my $func = '.';
    return $func if ($aa eq '' or $new_aa eq '');
    if ($aa eq $new_aa) {
	$func = 'coding-synon';
    }
    elsif ($cP < 4) {
	$func = 'altstart';
    }
    elsif ($aa eq '*') {
	if ($cP + 3 > (length($$self{codonseq}{$query_tid}))) {
	    $func = 'stop-loss';
	}
	else {
	    $func = 'abnormal-inseq-stop';
	}
    }
    elsif ($new_aa eq '*') {
	$func = 'nonsense';
    }
    else {
	$func = 'missense';
    }
    return $func;
}

=head2 get_codon

    About   : Get codon and its frame info
    Usage   : my ($codon, $aa, $pol, $frame) = get_codon($rseqs, $tid, $cpos);
    Args    : a hash ref of coding sequences, transcript id, and the cds postition.
    Returns : codon string, amino-acid, polar and the cpos's frame.

=cut
sub get_codon {
    my ($rseq, $tid, $cpos) = @_;
    my $frame = ($cpos - 1) % 3; # 0, 1, 2
    if (!exists $$rseq{$tid}) {
	carp "no coding sequence for $tid";
	return ("", "", "", $frame);
    }
    my $codon = uc(substr($$rseq{$tid}, ($cpos - $frame - 1), 3));
    my $aa = (exists $codon3{$codon}) ? $codon3{$codon} : "";
    my $pol = (exists $Polar{$aa}) ? $Polar{$aa} : "";
    return ($codon, $aa, $pol, $frame);
}

=head2 get_flanks

    About   : get region string for flanks, ready to fetchseq in batch mode.
    Usage   : my $flank = get_flanks($chr, $sta, $len);
    Returns : a hash ref like this:
		{ l => 'chr20:123-124', r => 'chr20:125-126' } 

=cut
sub get_flanks {
    my ($chr, $sta, $len) = @_;
    my %flank = (
	    l => $chr . ':' . ( $sta - 2 ) . '-' . ( $sta - 1 ),
	    r => $chr . ':' . ( $sta + $len ) . '-' . ( $sta + $len + 1 )
    );
    return \%flank;
}


=head2 select_position

    About   : Select the position that should be annotated on and get pairs by transcript ids
    Usage   : $var = $beda->select_position($var, \%cPos);
    Args    : var entry and cPos hash ref.
    Returns : var with standard {sel}{std}, and possible backward compitable {sel}{bc}

=cut
sub select_position {
    my ($self, $var, $rcPos) = @_;

    $$var{sel} = {};
    if (!exists $$var{'+'} or (!exists $$var{'+'}{p} and $$var{guess} ne 'rep')) {
    # simple variation or same pos(strand) non-repeat variation which caused by complex delins or multiple samples
	my $anno_sels = [];
	given ($$var{guess}) {
	    when ('snv') { # only select pos for the case: c.123A>G
		if (exists $$rcPos{$$var{pos}}) {
		    foreach my $tid (sort keys %{$$rcPos{$$var{pos}}}) {
			push (@$anno_sels, [ $tid, $$rcPos{$$var{pos}}{$tid} ]);
		    }
		}
	    }
	    when ('ins') { # the pos and pos+1 are selected: c.123_124insTG
		$anno_sels = $self->get_cover($$var{chr}, $$var{pos}, ($$var{pos}+1));
	    }
	    when (['del', 'delins']) { # the pos+1 and pos+reflen-1 are selected
		if ($$var{reflen} == 2) { # 1 bp deletion : c.123delT, c.123delTinsGAC
		    if (exists $$rcPos{($$var{pos}+1)}) {
			foreach my $tid (sort keys %{$$rcPos{($$var{pos}+1)}}) {
			    push (@$anno_sels, [ $tid, $$rcPos{($$var{pos}+1)}{$tid} ]);
			}
		    }
		}
		else { # multiple bases deletion: c.124_125delTG
                    $anno_sels = $self->get_cover(
                        $$var{chr},
                        ( $$var{pos} + 1 ),
                        ( $$var{pos} + $$var{reflen} - 1 )
                    );
		}
	    }
	    default { $self->throw("Error: unexpected guess for normal case: $$var{guess}"); }
	}
	$$var{sel}{std} = $anno_sels;
    }
    elsif (exists $$var{'+'}{p}) { # strand different pos
	my ($f_annos, $r_annos) = ([], []);
	my @sel_std = ();
	given ($$var{guess}) { # modified repeat caused ambigous 
	    when ('ins') {
                $f_annos = $self->get_cover( $$var{chr}, $$var{'+'}{p},
                    ( $$var{'+'}{p} + 1 ) );
                $r_annos = $self->get_cover( $$var{chr}, $$var{'-'}{p},
                    ( $$var{'-'}{p} + 1 ) );
	    }
	    when (['del', 'delins']) {
		if ($$var{'+'}{rl} == 2) { # 1bp deletion or delins : c.123delT  c.123delAinsGT
                    my $localrcp = $self->get_cPos( $$var{chr},
                        [ ( $$var{'+'}{p} + 1 ), ( $$var{'-'}{p} + 1 ) ]
                    );
                    if ( exists $$localrcp{ ( $$var{'+'}{p} + 1 ) } ) {
                        foreach my $tid (
                            sort
                            keys %{ $$localrcp{ ( $$var{'+'}{p} + 1 ) } } )
                        {
                            push(
                                @$f_annos,
                                [
                                    $tid,
                                    $$localrcp{ ( $$var{'+'}{p} + 1 ) }{$tid}
                                ]
                              );
                        }
                    }
                    if ( exists $$localrcp{ ( $$var{'-'}{p} + 1 ) } ) {
                        foreach my $tid (
                            sort
                            keys %{ $$localrcp{ ( $$var{'-'}{p} + 1 ) } } )
                        {
                            push(
                                @$r_annos,
                                [
                                    $tid,
                                    $$localrcp{ ( $$var{'-'}{p} + 1 ) }{$tid}
                                ]
                              );
                        }
                    }
		}
		else {
                    $f_annos = $self->get_cover(
                        $$var{chr},
                        ( $$var{'+'}{p} + 1 ),
                        ( $$var{'+'}{p} + $$var{'+'}{rl} - 1 )
                    );
                    $r_annos = $self->get_cover(
                        $$var{chr},
                        ( $$var{'-'}{p} + 1 ),
                        ( $$var{'-'}{p} + $$var{'-'}{rl} - 1 )
                    );
		}
	    }
	    default { $self->throw("Error: unexpected guess for different fr pos $$var{guess}"); }
	}
	
	foreach my $f_cpp (@$f_annos) {
	    push (@sel_std, $f_cpp) if ($$f_cpp[1]{strd} eq '+');
	}
	foreach my $r_cpp (@$r_annos) {
	    push (@sel_std, $r_cpp) if ($$r_cpp[1]{strd} eq '-');
	}

	$$var{sel}{std} = [@sel_std];
    }
    else {  # repeat
	my $whole_ref_annos = $self->get_cover( $$var{chr}, ($$var{pos}+1), ($$var{pos}+$$var{reflen}-1) );

	if ($$var{ref_cn} == 1 and $$var{alt_cn} == 2 and $$var{replen} > 1) { # for c.123_125dupGCA
	    my ($left_annos, $right_annos);
	    $left_annos = $self->get_cover(
		$$var{chr},
		($$var{pos} + 1),
		( $$var{pos} + $$var{replen} )
	    );
	    $right_annos = $self->get_cover(
		$$var{chr},
		( $$var{pos} + $$var{reflen} - $$var{replen} ),
		( $$var{pos} + $$var{reflen} - 1 )
	    );

	    my @sel_std = ();
	    foreach my $forw (@$left_annos) {
		push (@sel_std, $forw) if ($$forw[1]{strd} eq '+');
	    }
	    foreach my $rev (@$right_annos) {
		push (@sel_std, $rev)  if ($$rev[1]{strd} eq '-');
	    }
	    $$var{sel}{std} = [@sel_std];
	}
	else { # c.172AGT[3>5], c.123dupA
	    my @std_single = ();
	    foreach my $w_cpp (@$whole_ref_annos) {
		my ($tid, $rcpL, $rcpR) = @$w_cpp;
		if ($$rcpL{strd} eq '+') { # forward strand tid
		    push (@std_single, [ $tid, $rcpL ]);
		}
		else {
		    push (@std_single, [ $tid, $rcpR ]);
		}
	    }
	    $$var{sel}{std} = [@std_single];
	}
	
	# for convinience of plotting protein give the backward compatible mutation name
	if ($$var{replen} == 1 and ($$var{ref_cn} - $$var{alt_cn} == 1)) { # actually 1bp deletion: c.123delT
	    my @bc_single_del = ();
	    foreach my $w_cpp (@$whole_ref_annos) {
		my ($tid, $rcpL, $rcpR) = @$w_cpp;
		if ($$rcpL{strd} eq '+') { # forward strand tid
		    push (@bc_single_del,  [ $tid, $rcpR ]);
		}
		else { # reverse strand tid
		    push (@bc_single_del,  [ $tid, $rcpL ]);
		}
	    }
	    $$var{sel}{bc}  = [@bc_single_del];
	}
	else {
	    my ($bcf_annos, $bcr_annos); 
	    my @bc_annos = ();
	    if ($$var{ref_cn} > $$var{alt_cn}) { # deletion: c.123_125delTGT
		$bcf_annos = $self->get_cover(
		    $$var{chr},
		    ( $$var{'+'}{bp} + 1 ),
		    ( $$var{'+'}{bp} + $$var{'+'}{brl} - 1 )
		  );
		$bcr_annos = $self->get_cover(
		    $$var{chr},
		    ( $$var{'-'}{bp} + 1 ),
		    ( $$var{'-'}{bp} + $$var{'-'}{brl} - 1 )
		  );
	    }
	    else { # insertion: c.234_235insTTTT
		$bcf_annos = $self->get_cover(
		    $$var{chr},
		    ( $$var{'+'}{bp} + $$var{'+'}{brl} - 1 ),
		    ( $$var{'+'}{bp} + $$var{'+'}{brl} )
		);
		$bcr_annos = $self->get_cover(
		    $$var{chr},
		    $$var{'-'}{bp}, ($$var{'-'}{bp}+1)
		);
	    }
	    foreach my $bc_cppf (@$bcf_annos) {
		push (@bc_annos, $bc_cppf) if ($$bc_cppf[1]{strd} eq '+');
	    }
	    foreach my $bc_cppr (@$bcr_annos) {
		push (@bc_annos, $bc_cppr) if ($$bc_cppr[1]{strd} eq '-');
	    }
	    $$var{sel}{bc}  = [@bc_annos];
	}
    }
    return $var;
}


=head2 batch_anno

    About   : The fastest way to annotate multiple snv and 1bp deletion variations,
	      indel and other types also can be annotated, but no faster than annotated
	      one by one.
    Usage   : $beda = BedAnno->new('in.bed.gz', 'in.chr.fa', 'in.codon.fa'); 
	      @all_annos = ();
	      foreach chr {
		my @vars = (); 
		foreach pos, ref, alt {
		  $var = parse_var($chr, $pos, $ref, $alt);
		  push (@vars, $var);
		}
		my $rChrAnnos = $beda->batch_anno( chr, \@vars); 
		push (@all_annos, @$rChrAnnos);
	      }
    Args    : an array ref of vars and an array ref of cPos
    Returns : an array ref of annos, see anno().

=cut
sub batch_anno {
    my ($self, $chr, $rVars) = @_;
    my @all_annoRst = ();

    my %all_pos = ();
    foreach my $v (@$rVars) {
	my $single_pos = get_anno_pos($v);
	$all_pos{$single_pos} = 1 if ($single_pos > 0);
    }
    my @poses = sort {$a<=>$b} keys %all_pos;
    
    my $cPos = $self->get_cPos( $chr, \@poses );
    for (my $i = 0; $i < @$rVars; $i ++) {
	my $var = $self->select_position($$rVars[$i], $cPos);
	push (@all_annoRst, $self->varanno($var));
    }
    return \@all_annoRst;
}


=head2 get_anno_pos

    About   : get snv and variation with 1bp deletion positions.
    Usage   : my @toAnno = get_anno_pos($var);

=cut
sub get_anno_pos {
    my $var = shift;
    if (exists $$var{'+'}) {	# for multiple samples caused pseudo-complex
	if ($$var{'+'}{brl} == 1 and $$var{'+'}{bal} == 1) { # back compatible snv
	    return $$var{'+'}{bp};
	}
	elsif ($$var{'+'}{brl} == 2) { # 1bp deletion or delins
	    return ($$var{'+'}{bp} + 1);
	}
    }
    elsif ($$var{guess} eq 'snv') {
	return $$var{pos};
    }
    elsif ($$var{guess} eq 'del' and $$var{reflen} == 2) {
	return ($$var{pos} + 1);
    }
    return 0;
}

=head2 get_cPos

    About   : get cPos for all needed positions on the same chromosome, 
	      and also parse the anno bedent by the way.
    Usage   : my $rcPos = $beda->get_cPos( $chr, \@pos, local => 1 );
    Args    : chromosome positions should be sorted from small to large, 
	      if 'local' tag used, a localized bed will be load without
	      override the original db entry, other restriction will
	      also be available to use, just like genes, trans...
    Returns : an hash ref of cPos, cPos format:
	      border-indicater
		0 for non-border
		b for 5'end border
		B for 3'end border
		1 for too short block

	       { $pos1 => { 
			    $tid => {
				      gsym => $gsym, 
				      reg  => $region, 
				      exin => $exin,
				      cpos => $cpos, 
				      strd => [+-], 
				      bd   => [bB01]
				    }, ... 
			   }, ... 
	       }

=cut
sub get_cPos {
    my ($self, $chr, $rpos, %args) = @_;

    my $rAnnos;
    if (exists $args{local} and $args{local}) {
	$rAnnos = $self->load_anno( %args, region => $chr.":".$$rpos[0]."-".$$rpos[-1] );
    }
    else {
	$rAnnos = $$self{annodb};
    }
    return {} if (!exists $$rAnnos{$chr});

    my $rbed = $$rAnnos{$chr};

    my %rmdup = map {$_ => 1} @$rpos;
    my @sortp = sort {$a <=> $b} keys %rmdup;
    my %cPos = ();
    my $k = 0;
    foreach my $pos (@sortp) { # one by one from small to large
	while ($k < @$rbed) {
	    if ($$rbed[$k]{sto} < $pos) {
		$k ++;
		next;
	    }
	    if ($$rbed[$k]{sta} >= $pos) {
		$k = ($k - 1 < 0) ? 0 : ($k - 1);
		last;
	    }
	    if (!exists $$rbed[$k]{detail}) {
		foreach my $annoblk (keys %{$$rbed[$k]{annos}}) {
		    my $offset1 = $$rbed[$k]{annos}{$annoblk};
		    my ($tid, $ranno) = parse_annoent($annoblk);
		    $$rbed[$k]{detail}{$tid} = $ranno;
		    $$rbed[$k]{detail}{$tid}{offset} = $offset1;
		}
	    }

	    my $offset = $pos - $$rbed[$k]{sta} - 1;
	    foreach my $t (sort keys %{$$rbed[$k]{detail}}) {
		my $rh = $$rbed[$k]{detail}{$t};
		$cPos{$pos}{$t} = in_reg_cPos($rh, $offset);
	    }
	    last;
	}
    }
    return \%cPos;
}

sub in_reg_cPos {
    my ($rh, $left_offset) = @_;
    my $total_left_offset  = $left_offset + $$rh{offset};
    my $total_right_offset = $$rh{wlen} - $total_left_offset - 1;
    my $border;
    if ( $total_left_offset < 2 and $total_right_offset < 2 ) {
	$border = 1;
    }
    elsif (( $total_left_offset < 2 and $$rh{strand} eq '+' )
	or ( $total_right_offset < 2 and $$rh{strand} eq '-' ) )
    {
	$border = 'b';
    }
    elsif (( $total_left_offset < 2 and $$rh{strand} eq '-' )
	or ( $total_right_offset < 2 and $$rh{strand} eq '+' ) )
    {
	$border = 'B';
    }
    else {
	$border = 0;
    }
    
    if ($$rh{nsta} eq '') { # for incomplete cds cases
	return { gsym => $$rh{gsym}, reg => $$rh{blka}, exin => $$rh{exin}, cpos => 'c.?', strd => $$rh{strand}, bd => $border };
    }

    my $cpos;
    my $strandopt = ($$rh{strand} eq '+') ? 1 : 0;
    given ($$rh{blka}) {
	when (/^C/ or /^5/) { # cds or utr-5
	    if ($strandopt) {
		$cpos = 'c.'.($$rh{csta} + $total_left_offset);
	    }
	    else {
		$cpos = 'c.'.($$rh{csto} + $total_right_offset);
	    }
	}
	when (/^I/) { # intron
	    if ($$rh{csta} ne '') { # coding rna
		if ($total_left_offset < $total_right_offset) {
		    my $opt = ($strandopt) ? '+' : '-';
		    $cpos = 'c.'.$$rh{csta}.$opt.($total_left_offset + 1);
		}
		else {
		    my $opt = ($strandopt) ? '-' : '+';
		    $cpos = 'c.'.$$rh{csto}.$opt.($total_right_offset + 1);
		}
	    }
	    else { # ncRNA
		if ($total_left_offset < $total_right_offset) {
		    my $opt = ($strandopt) ? '+' : '-';
		    $cpos = 'r.'.$$rh{nsta}.$opt.($total_left_offset + 1);
		}
		else {
		    my $opt = ($strandopt) ? '-' : '+';
		    $cpos = 'r.'.$$rh{nsto}.$opt.($total_right_offset + 1);
		}
	    }
	}
	when (/^3/) { # utr-3
	    if ($strandopt) {
		my $p = $$rh{csta}; $p =~ s/^\*//;
		$cpos = 'c.*'.($p + $total_left_offset);
	    }
	    else {
		my $p = $$rh{csto}; $p =~ s/\*//;
		$cpos = 'c.*'.($p + $total_right_offset);
	    }
	}
	when (/^R/) { # ncRNA
	    if ($strandopt) {
		$cpos = 'r.'.($$rh{nsta} + $total_left_offset);
	    }
	    else {
		$cpos = 'r.'.($$rh{nsto} + $total_right_offset);
	    }
	}
	default { confess "unrecognized block attribution!"; }
    }
    return { gsym => $$rh{gsym}, reg => $$rh{blka}, exin => $$rh{exin}, cpos => $cpos, strd => $$rh{strand}, bd => $border };
}

=head2 get_cover

    About   : get covered region infos for deletions, and also for insertion pos pairs.
    Usage   : my $rcover = $beda->get_cover($chr, $start, $stop, local => 1);
    Args    : use 'local' tag to call tabix to get bed annotation locally.
    Returns : an array ref of covered region cPos pairs, cPos with "--" or "++" indicate 
	      outside of transcripts.
		[ [ $tid, $left_cPos, $right_cPos ], ... ]

=cut
sub get_cover {
    my ($self, $chr, $start, $stop, %args) = @_;

    my $rAnnos;
    if (exists $args{local} and $args{local}) {
	# also allow other open_args to restrict the return bed ents, override the ori class args.
	$rAnnos = $self->load_anno( %args, region => $chr.":".$start."-".$stop );
    }
    else {
	$rAnnos = $$self{annodb};
    }

    return [] if (!exists $$rAnnos{$chr});
    my @covbeds = ();
    my %min_left = ();
    my %min_right = ();
    my %left_most = ();
    my %right_most = ();
    my $k = 0;
    foreach my $bedent (@{$$rAnnos{$chr}}) {
	next if ($start > $$bedent{sto});
	last if ($stop  <= $$bedent{sta});
	if (!exists $$bedent{detail}) {
	    foreach my $annoblk (keys %{$$bedent{annos}}) {
		my $offset1 = $$bedent{annos}{$annoblk};
		my ($tid, $ranno) = parse_annoent($annoblk);
		$$bedent{detail}{$tid} = $ranno;
		$$bedent{detail}{$tid}{offset} = $offset1;
	    }
	}
	push (@covbeds, $bedent);
	foreach my $tid (sort keys %{$$bedent{detail}}) {
	    if (!exists $min_left{$tid}) {
		$min_left{$tid} = $$bedent{sta} + 1 - $start;
		$left_most{$tid} = $k;
	    }
	    if (!exists $min_right{$tid}) {
		$min_right{$tid} = $stop - $$bedent{sto};
		$right_most{$tid} = $k;
	    }
	    elsif ($min_right{$tid} > ($stop - $$bedent{sto})) {
		$min_right{$tid} = $stop - $$bedent{sto};
		$right_most{$tid} = $k;
	    }
	}
	$k ++;
    }
    return [] if (0 == @covbeds);

    my @ret_pairs;
    foreach my $tid (sort keys %left_most) {
	my ($left_cPos, $right_cPos);
        my ( $left_bed, $right_bed ) =
          ( $covbeds[ $left_most{$tid} ], $covbeds[ $right_most{$tid} ] );

	if ($min_left{$tid} > 0) {  # extend off the left of tid
	    my $left_cp;
	    my $bd = 0;
	    if ($$left_bed{detail}{$tid}{strand} eq '+') {
		$bd = "b" if ($min_left{$tid} < 2);
		$left_cp = "--".$min_left{$tid};
	    }
	    else {
		$bd = "B" if ($min_left{$tid} < 2);
		$left_cp = "++".$min_left{$tid};
	    }
            $left_cPos = {
                gsym => $$left_bed{detail}{$tid}{gsym},
                reg  => $$left_bed{detail}{$tid}{blka},
                exin => $$left_bed{detail}{$tid}{exin},
                strd => $$left_bed{detail}{$tid}{strand},
                cpos => $left_cp,
                bd   => $bd
              };
	}
	else {
	    $left_cPos = in_reg_cPos($$left_bed{detail}{$tid}, -$min_left{$tid});
	}
	if ($min_right{$tid} > 0) { # extend off the right of tid
	    my $right_cp;
	    my $rbd = 0;
	    if ($$right_bed{detail}{$tid}{strand} eq '+') {
		$rbd = "B" if ($min_right{$tid} < 2);
		$right_cp = "++".$min_right{$tid};
	    }
	    else {
		$rbd = "b" if ($min_right{$tid} < 2);
		$right_cp = "--".$min_right{$tid};
	    }
            $right_cPos = {
                gsym => $$right_bed{detail}{$tid}{gsym},
                reg  => $$right_bed{detail}{$tid}{blka},
                exin => $$right_bed{detail}{$tid}{exin},
                strd => $$right_bed{detail}{$tid}{strand},
                cpos => $right_cp,
                bd   => $rbd
              };
	}
	else {
            my $converse_left =
              ( $$right_bed{sto} - $$right_bed{sta} + $min_right{$tid} - 1 );
            $right_cPos =
              in_reg_cPos( $$right_bed{detail}{$tid}, $converse_left );
	}
	push (@ret_pairs, [ $tid, $left_cPos, $right_cPos ]);
    }
    return \@ret_pairs;
}


# just parse annoents when mutation hits.
sub parse_annoent {
    my $annoent = shift;
    my %annoinfo = ();
    # NM_152486.2|SAMD11|+|IC7|IVS8|949|950|869|870|Y
    my @infos = split(/\|/, $annoent);
    my $tid = shift @infos;
    confess "Error format of anno ents [$annoent]" if (10 != @infos);
    my @tags  = qw(gsym strand blka exin nsta nsto csta csto wlen pr);
    @annoinfo{@tags} = @infos;
    return ($tid, \%annoinfo);
}

=head2 parse_var

    About   : parse the variation directly by the ref and alt string
    Usage   : my $var = parse_var( $chr, $pos, $ref, $alt );
    Returns : a hash ref of variation
	      chr, pos, ref, alt, reflen, altlen, guess

	      # if fpos is the same with rpos, (pos, ref, alt, reflen, altlen) will be updated

	      # for different fpos and rpos cases only, usually complex case
	      '+' => {
		  p - forward strand offset
		  r - forward strand ref string
		  a - forward strand alt string
		  rl
		  al
		}
	      '-' => {
		  same to '+'
		}

	      # for repeat type only
	      rep    - repeat element
	      replen - repeat elementlen
	      ref_cn - copy number in reference
	      alt_cn - copy number in alternative

	      # backward compatible to previous results if parsed complex
	      '+' => {
		  bp - forward strand offset
		  br - forward strand ref string
		  ba - forward strand alt string
		  brl
		  bal
		}
	      '-' => {
		  same to '+'
		}

=cut
sub parse_var {
    my ($chr, $pos, $ref, $alt) = @_;
    my $ref_len = length ($ref);
    my $alt_len = length ($alt);
    my %var = (
        chr    => $chr,
        pos    => $pos,
        ref    => $ref,
        alt    => $alt,
        reflen => $ref_len,
        altlen => $alt_len
    );
    if ($ref_len > 1 and $alt_len > 1) {
        my $complex = parse_complex( $ref, $ref_len, $alt, $alt_len );
	my $ori_pos = $var{pos};

        $var{guess} = $$complex{guess};

        $var{'+'}{bp} = $ori_pos + $$complex{'+'}{bcOffst};
        $var{'+'}{br} =
          substr( $ref, $$complex{'+'}{bcOffst}, $$complex{bcRlen} );
        $var{'+'}{ba} =
          substr( $alt, $$complex{'+'}{bcOffst}, $$complex{bcAlen} );
        $var{'+'}{brl} = $$complex{bcRlen};
        $var{'+'}{bal} = $$complex{bcAlen};

        $var{'-'}{bp}    = $ori_pos + $$complex{'-'}{bcOffst};
        $var{'-'}{br} =
          substr( $ref, $$complex{'-'}{bcOffst}, $$complex{bcRlen} );
        $var{'-'}{ba} =
          substr( $alt, $$complex{'-'}{bcOffst}, $$complex{bcAlen} );
        $var{'-'}{brl} = $$complex{bcRlen};
        $var{'-'}{bal} = $$complex{bcAlen};

        if ( $$complex{'+'}{offset} == $$complex{'-'}{offset} ) {
            $var{pos} += $$complex{'+'}{offset};
            $var{reflen} = $$complex{newRlen};
            $var{altlen} = $$complex{newAlen};
            $var{ref} =
              substr( $ref, $$complex{'+'}{offset}, $$complex{newRlen} );
            $var{alt} =
              substr( $alt, $$complex{'+'}{offset}, $$complex{newAlen} );
        }
        else {
            $var{'+'}{p} = $ori_pos + $$complex{'+'}{offset};
            $var{'+'}{r} =
              substr( $ref, $$complex{'+'}{offset}, $$complex{newRlen} );
            $var{'+'}{a} =
              substr( $alt, $$complex{'+'}{offset}, $$complex{newAlen} );
            $var{'+'}{rl} = $$complex{newRlen};
            $var{'+'}{al} = $$complex{newAlen};

            $var{'-'}{p}    = $ori_pos + $$complex{'-'}{Offset};
            $var{'-'}{r} =
              substr( $ref, $$complex{'-'}{offset}, $$complex{newRlen} );
            $var{'-'}{a} =
              substr( $alt, $$complex{'-'}{offset}, $$complex{newAlen} );
            $var{'-'}{rl} = $$complex{newRlen};
            $var{'-'}{al} = $$complex{newAlen};
        }
        if ( $$complex{guess} eq 'rep' ) {
            $var{rep}    = $$complex{rep};
            $var{replen} = $$complex{replen};
            $var{ref_cn} = $$complex{refcn};
            $var{alt_cn} = $$complex{altcn};
        }
    }
    else {
	$var{guess} = guess_type( $ref_len, $alt_len );
    }

    return \%var;
}

=head2 parse_complex
    
    About   : parse complex ref and alt and guess the standard variation depond on the strand
	      if the strand is '+', the leading coordinate string will be trimmed, and for '-', 
	      the tail will be trimmed.
    Usage   : my $rguess = parse_complex( $ref, $len_ref, $alt, $len_alt );
    Returns : a hash ref which contain the following information 
	      {
	       guess  => $guess,	   (set( 'ref' / 'snv' / 'ins' / 'del' / 'delins' / 'rep' ))
	       
	       newRlen => $new_ref_len,
	       newAlen => $new_alt_len,

	       bcRlen  => $bcRlen,
	       bcAlen  => $bcAlen,
	       # for backward compatible when dealing complex delins (usually rep with modifications)

	       bcGuess => $bcGuess,

	       '+' => {
		   offset  => $offset,
		   bcOffst => $bcOffst,
		},

	       '-' => {
		   the same with '+'
	        },

	       # for repeat only
	       rep    => $repeat,
	       replen => $repeatlen,
	       ref_cn => $ref_copy_num,
	       alt_cn => $alt_copy_num
	     }
=cut
sub parse_complex {
    my ($ref, $len_ref, $alt, $len_alt) = @_;

    return { guess => "ref" } if ($ref eq $alt);

    # trim the first leading string for ref and alt
    my $cref = substr($ref, 1);
    my $calt = substr($alt, 1);
    my $reflead = substr($ref, 0, 1);
    my $altlead = substr($alt, 0, 1);

    my $rc_ref = count_content($cref);
    my $rc_alt = count_content($calt);
    my @diff = map { $$rc_ref[$_] - $$rc_alt[$_] } (0 .. 4);

    my $get_rst = get_internal( $ref, $len_ref, $alt, $len_alt );
    my $guess = guess_type( $$get_rst{r}, $$get_rst{a} );

    my %rst = ();
    @rst{ ("bcGuess", "bcRlen", "bcAlen") } = ($guess, $$get_rst{r}, $$get_rst{a});
    $rst{'+'}{bcOffst} = $$get_rst{'+'};
    $rst{'-'}{bcOffst} = $$get_rst{'-'};

    my $sign_coord =
      check_sign( \@diff ); # check if the sign of all base diff are consistent.
    if ( ( $sign_coord == 0 ) or ( $diff[0] == 0 ) or ( $reflead ne $altlead ) )
    {
      # for uncompatible lead char, same length and diff-sign incoordinate cases
        @rst{ ( "guess", "newRlen", "newAlen" ) } =
          ( $guess, $$get_rst{r}, $$get_rst{a} );
	$rst{'+'}{offset} = $$get_rst{'+'};
	$rst{'-'}{offset} = $$get_rst{'-'};
    }
    elsif ($sign_coord) { # possible short tandom repeat variation
	my @absdiff = map {abs} @diff;
	my ($larger, $smaller, $llen, $slen);
        if ( $len_ref > $len_alt ) {
            $larger  = $cref;
            $llen    = $len_ref - 1;
            $smaller = $calt;
            $slen    = $len_alt - 1;
        }
        else {
            $larger  = $calt;
            $llen    = $len_alt - 1;
            $smaller = $cref;
            $slen    = $len_ref - 1;
        }

	my %has = ();
	for (my $rep = $llen; $rep > 0; $rep--) {
	    while ($larger =~ /([ACGT]+)(?:\1){$rep}/g) {
		next if (exists $has{$1});

                my $rep_el = $1;
                my $lofs   = length($`);    # $` is the prematched string

                my $cn = check_div( $rep_el, \@absdiff );
		if ($cn and check_insrep($larger, $smaller, $rep_el, $cn)) {
		    my $lenrep = length($rep_el);
                    @rst{
                        qw(guess rep replen)
                    } = ( 'rep', $rep_el, $lenrep );

		    $rst{'+'}{offset} = $lofs;
		    $rst{'-'}{offset} = $lofs;

		    $rep += 1; # add the first copy of element

		    if ($llen == $len_ref - 1) { # ref is longer
                        @rst{
                            qw(refcn altcn newRlen newAlen)
                          } = (
                            $rep, ( $rep - $cn ),
                            ( $lenrep * $rep + 1 ), ( $lenrep * ( $rep - $cn ) + 1 )
                          );
		    }
		    else { # alt is longer
                        @rst{
                            qw(altcn refcn newAlen newRlen)
                          } = (
                            $rep, ( $rep - $cn ),
                            ( $lenrep * $rep + 1 ), ( $lenrep * ( $rep - $cn ) + 1 )
                          );
		    }
		    return \%rst;
		}
		$has{$rep_el} = 1;
	    }
	}

	# for non-tandom repeat case, sparsed case which will lead to difference between samtools mpileup and GATK.
	if (check_trim_tail($larger, $llen, $smaller, $slen)) { 
	    # check if the tail of smaller is uniq occurence in larger, if uniq, this may caused by multiple samples.
	    @rst{ ( "guess", "newRlen", "newAlen" ) } =
	      ( $guess, $$get_rst{r}, $$get_rst{a} );
	    $rst{'+'}{offset} = $$get_rst{'+'};
	    $rst{'-'}{offset} = $$get_rst{'-'};
	}
	else { # try get backward compatible variation and keep the original delins 
	    @rst{ ("guess", "newRlen", "newAlen") } = ('delins', $len_ref, $len_alt);
	    $rst{'+'}{offset} = 0;
	    $rst{'-'}{offset} = 0;
	}
    }
    return \%rst;
}


# deal with the GATK multiple-sample cases
# if the smaller cstr (string after the first base) only exists at the end of 
# the larger cstr, 
sub check_trim_tail {
    my ($larger, $llen, $smaller, $slen) = @_;
    if (($llen - $slen) == index($larger, $smaller) and (substr($larger, -$slen, $slen) eq $smaller)) {
	return 1;
    }
    return 0;
}

sub guess_type {
    my ($reflen, $altlen) = @_;
    my $guess;
    if ($reflen == $altlen) {
	$guess = ($reflen == 1) ? 'snv' : 'delins';
    }
    elsif ($reflen > 1 and $altlen > 1) {
	$guess = 'delins';
    }
    elsif ($reflen == 1) {
	$guess = 'ins';
    }
    else {
	$guess = 'del';
    }
    return $guess;
}


# check whether the smaller with inserted $cn copies repeat elements
# is the same with larger one
sub check_insrep {
    my ($larger, $smaller, $repeat, $cn) = @_;
    my $ind = index($smaller, $repeat);
    if ($ind > -1) {
	my $temp = $smaller;
	substr($temp, $ind, 0, ($repeat x $cn));
	if ($larger eq $temp) {
	    return 1;
	}
    }
    return 0;
}


=head2 get_internal

    About   : recalculate ref alt for delins and mutiple sample caused ref-alt pair.
	      depend on different strand of the to-annotated gene or transcript, 
	      the offset may be different for the same ref and alt,
	      because of the 3'end nearest annotation rules.
    Usage   : my $rephase = get_internal( $ref, $reflen, $alt, $altlen );
    Returns : a hash ref of : 
		{
		    '+' => $f_lofs,
		    '-' => $r_lofs,
		    'r' => $new_ref_len, 
		    'a' => $new_alt_len
		}

=cut
sub get_internal {
    my ($ref, $reflen, $alt, $altlen) = @_;
    my $shorter = ($reflen < $altlen) ? $reflen : $altlen;
    my ($lgo, $loff, $rgo, $roff) = (1, 0, 1, 0);
    for (my $i = 0; $i < $shorter; $i++) {
	if ($lgo and substr($ref, $i, 1) eq substr($alt, $i, 1)) {
	    $loff ++;
	}
	else {
	    $lgo = 0;
	}

	if ($rgo and substr($ref,-($i+1),1) eq substr($alt,-($i+1),1)) {
	    $roff ++;
	}
	else {
	    $rgo = 0;
	}

	last if ($lgo == 0 and $rgo == 0);
    }
    my ($lofs, $new_ref_len, $new_alt_len);
    if ($shorter >= $loff + $roff) {
	$new_ref_len = $reflen - $loff - $roff + 1;
	$new_alt_len = $altlen - $loff - $roff + 1;
	$lofs = ($loff - 1 < 0) ? 0 : ($loff - 1);
	return {
	    '+' => $lofs,
	    '-' => $lofs,
	    'r' => $new_ref_len,
	    'a' => $new_alt_len
	};
    }
    else {
	my $trim_len = $shorter - 1;
	$new_ref_len = $reflen - $trim_len;
	$new_alt_len = $altlen - $trim_len;

	my ($f_lofs, $r_lofs);
	$f_lofs = ($loff - 1 < 0) ? 0 : ($loff - 1);
	$r_lofs = ($shorter - $roff - 1 < 0) ? 0 : ($shorter - $roff - 1 < 0);

	return {
	    '+' => $f_lofs,
	    '-' => $r_lofs,
	    'r' => $new_ref_len,
	    'a' => $new_alt_len
	};
    }
}

# check sign for diff-array
sub check_sign {
    my $rc = shift;
    return 0
      if ( ( $$rc[0] > 0 and ( grep { $_ < 0 } @$rc ) )
        or ( $$rc[0] < 0 and ( grep { $_ > 0 } @$rc ) ) );
    return 1;
}

# check length and content consistent-divisability for string and diff-array
sub check_div {
    my ($s, $rc) = @_;
    my $rcs = count_content($s);
    return 0 unless ($$rc[0] % $$rcs[0] == 0);
    my $div = $$rc[0] / $$rcs[0];
    for (1 .. 4) {
        return 0 if ($$rc[$_] != $$rcs[$_] * $div);
    }
    return $div;
}

sub count_content {
    my $s = uc(shift);
    my $l = ($s =~ tr/ACGT/1234/);
    my @count = ($l, 0, 0, 0, 0);
    while ($s=~/(\d)/g) {
	$count[$1] ++;
    }
    return \@count;
}


=head2 fetchseq

    About   : get sequence from fasta db using samtools faidx
    Usage   : my $rhash = fetchseq('db.fasta', \@regions);
    Args    : a list of region in format: (chr1:123-456, chr1:789-1000, chr2:234-567, or NM_01130.1:345 )
    Returns : a hash ref of { region => seq }

=cut
sub fetchseq {
    my ($fasta, $rRegs) = @_;

    my %seqs = ();
    local $/ = ">";
    while (0 < @$rRegs) {
	my @set = splice (@$rRegs, 0, 1000);
	my $regions = join(" ", @set);
	open (FASTA, "samtools faidx $fasta $regions |") or confess "samtools faidx: $!";
	my @all_seq = <FASTA>;
	close FASTA;
	shift @all_seq;
	foreach my $rec (@all_seq) {
	    chomp $rec;
	    my $hd = $1 if ($rec =~ s/^(\S+)//);
	    $rec =~ s/[\s>]+//g;
	    $seqs{$hd} = uc($rec);
	}
    }
    return \%seqs;
}

sub rev_comp {
    my $Seq = shift;
    $Seq = reverse($Seq);
    $Seq =~ tr/ATCG/TAGC/; # only deal with 'A', 'T', 'G', 'C'
    return $Seq;
}

sub throw {
    my ($self,@msg) = @_;
    confess @msg,"\n";
}

sub warn {
    my ($self, @msg) = @_;
    carp @msg,"\n";
}


1;
__END__

=head1 SEE ALSO

Mention other useful documentation such as the documentation of
related modules or operating system documentation (such as man pages
in UNIX), or any relevant external documentation such as RFCs or
standards.

If you have a mailing list set up for your module, mention it here.

If you have a web site set up for your module, mention it here.

=head1 AUTHOR

liutao, E<lt>liutao@genomics.cn<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by liutao

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
