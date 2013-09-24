package BedAnno;

use strict;
use warnings;
use threads::shared; # if threads used before BedAnno, then this will be threaded
use Thread::Queue;
use Scalar::Util qw(refaddr);
use Carp;

require Exporter;

our @ISA = qw(Exporter);

our @EXPORT_OK = qw(
    parse_var fetchseq Code2Pep C3toC1 C1toC3
);

our %EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

our $VERSION = '0.32';

=head1 NAME

BedAnno - Perl module for annotating variation depend on the BED +1 format database.

=head2 VERSION v0.32

From version 0.32 BedAnno will change to support CG's variant shell list
and use ncbi annotation release 104 as the annotation database
with reformatted database format, and won't give any individual
annotation, so the individual_anno() is no longer available.

=head1 SYNOPSIS

  use threads;
  use BedAnno;
  my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', short => 1 );
  my $anno = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );

=head1 DESCRIPTION

This module need bgzipped BED+1 format database, combined with tabix index,
tabix is require for fast extract anno infomation randomly.
The fasta sequences are indexed by samtools faidx, sequence splicing also requires
samtools. If allele frequency information, prediction information or cytoBand 
information is needed, then extra dependencies will be required

=head2 EXPORT_OK for all

    parse_var()		- parse variation information stand-alone, giving
			  the guess information and possible coordinates.
    fetchseq()		- used to batched fetch genomic sequences (flanks), 
			  this depend on samtools installed, and a faidx(ed) 
			  fasta file.
    Code2Pep()		- encode 3 bp nucleotides to peptide code, the format
			  of return code, depend on the option you used in 
			  inition method, if "short" exists, then C1 will be
			  the output format, default C3 format.
    C3toC1()		- encode 3 chars peptide code to 1 char peptide code.
    C1toC3()		- encode 1 chars peptide code to 3 char peptide code.

=head1 Methods

=cut

our (%C3, %C1, %Code2Pep, %SO2Name, %Name2SO, %Polar, %C3toC1, %C1toC3) :shared;

%C3 = (
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
		CGN=>"Arg",		GGN=>"Gly",
    
    # inseq stop codon
    UAA=>"X",   UAG=>"X",

    # selenocysteine
    UGA=>"Sec"
);


%C1 = (
    TTT=>"F", CTT=>"L", ATT=>"I", GTT=>"V",
    TTC=>"F", CTC=>"L", ATC=>"I", GTC=>"V",
    TCT=>"S", CCT=>"P", ACT=>"T", GCT=>"A",
    TCC=>"S", CCC=>"P", ACC=>"T", GCC=>"A",
    TAT=>"Y", CAT=>"H", AAT=>"N", GAT=>"D",
    TAC=>"Y", CAC=>"H", AAC=>"N", GAC=>"D",
    TGT=>"C", CGT=>"R", AGT=>"S", GGT=>"G",
    TGC=>"C", CGC=>"R", AGC=>"S", GGC=>"G",
    TTA=>"L", CTA=>"L", ATA=>"I", GTA=>"V",
    TCA=>"S", CTG=>"L", ACA=>"T", GTG=>"V",
    TAA=>"*", CCA=>"P", AAA=>"K", GCA=>"A",
    TGA=>"*", CCG=>"P", AGA=>"R", GCG=>"A",
    TTG=>"L", CAA=>"Q", ATG=>"M", GAA=>"E",
    TCG=>"S", CAG=>"Q", ACG=>"T", GAG=>"E",
    TAG=>"*", CGA=>"R", AAG=>"K", GGA=>"G",
    TGG=>"W", CGG=>"R", AGG=>"R", GGG=>"G",

    TCN=>"S", CCN=>"P", ACN=>"T", GTN=>"V",
	      CTN=>"L",		  GCN=>"A",
	      CGN=>"R",		  GGN=>"G",
    
    UAA=>"X",   UAG=>"X",

    UGA=>"U"
);

%C3toC1 = ();
%C1toC3 = ();
foreach my $threebase (sort keys %C3) {
    $C3toC1{$C3{$threebase}} = $C1{$threebase};
    $C1toC3{$C1{$threebase}} = $C3{$threebase};
}

%Polar = (
    Ala => "NP", Asn => "P0", Arg => "P+", Asp => "P-",
    Ile => "NP", Cys => "P0", His => "P+", Glu => "P-",
    Leu => "NP", Gln => "P0", Lys => "P+",
    Met => "NP", Gly => "P0",
    Phe => "NP", Ser => "P0",
    Pro => "NP", Thr => "P0",
    Trp => "NP", Tyr => "P0",
    Val => "NP",
    Sec => "NP",

    'X' => '.',  '*' => '.'
);

%Code2Pep = %C3;
@Polar{@C1{(sort keys %C1)}} = @Polar{@C3{(sort keys %C1)}};

%SO2Name = (

    # variant type
    "SO:0000159" => 'del',
    "SO:1000032" => 'delins',
    "SO:0001483" => 'snv',
    "SO:0000667" => 'ins',
    "ref"        => 'ref',
    "no-call"    => 'no-call',

    # Gene Parts
    "SO:0000316"      => 'CDS',
    "SO:0000204"      => 'five_prime_UTR',
    "SO:0000205"      => 'three_prime_UTR',
    "SO:0000655"      => 'ncRNA',
    "SO:0000191"      => 'interior_intron',
    "SO:0000448"      => 'three_prime_UTR_intron',
    "SO:0000447"      => 'five_prime_UTR_intron',
    "SO:0000163"      => 'five_prime_cis_splice_site',
    "SO:0000164"      => 'three_prime_cis_splice_site',
    "SO:0000167"      => 'promoter',
    "span"            => 'span',
    "abnormal_intron" => 'abnormal_intron',
    "annotation-fail" => 'annotation-fail',

    # Function Parts
    "SO:0001819" => 'synonymous_variant',
    "SO:0001583" => 'missense_variant',
    "SO:0001587" => 'stop_gained',
    "SO:0001578" => 'stop_lost',
    "SO:0001822" => 'inframe_deletion',
    "SO:0001821" => 'inframe_insertion',
    "SO:0001589" => 'frameshift_variant',
    "SO:0001582" => 'initiator_codon_variant',
    "SO:0001893" => 'transcript_ablation',
    "SO:0001567" => 'stop_retained_variant',

    # for non-equal length of substitution
    "inframe_delins" => 'inframe_delins',

    # for refSeq the same with call seq.
    "no-change" => 'no-change',

    # for span and splice
    "unknown-likely-deleterious" => 'unknown-likely-deleterious',

    # for ncRNA, utr, intron or intergenic
    "unknown" => 'unknown',

    # for no-call variant
    "unknown-no-call" => 'unknown-no-call',

    # the followings are replaced by 'unknown-likely-deleterious' in Voyager
    "SO:0001575" => 'splice_donor_variant',
    "SO:0001574" => 'splice_acceptor_variant',

    # the followings are replaced by 'unknown' in Voyager
    "SO:0001623" => '5_prime_UTR_variant',
    "SO:0001624" => '3_prime_UTR_variant',
    "SO:0001619" => 'nc_transcript_variant',
    "SO:0001627" => 'intron_variant'
);

%Name2SO = reverse(%SO2name);

=head2 new

=over

=item About : Creat a new annotation entry

=item Usage :

    my $beda = $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', short => 1 );

=item Args	- (all database files should be tabix indexed)

=over

=item Essential Args:

=over

=item I<db> <in.bed.gz>

=over

=item annotation database file. 

=back

=item I<tr> <in.trans.fas.gz>

=over

=item transcript sequence fasta file

=back

=item See L</DATABASE FORMAT> for more infomation. 

=back

=item Optional Args :

=over

=item Common options :

=over 

=item I<threads> [number]

=over

=item support paralleled annotation, give the number of threads be used for annotation.

=back

=item I<cytoBand> [cytoBand.bed.gz]

=over

=item add cytoBand information

=back

=item I<pfam> [pfam.tsv.gz]

=over

=item add pfam information

=back

=item I<prediction> [ensembl_prediction_db.tsv.gz]

=over

=item add sift, polyphen2 prediction and scores

=back

=item I<phyloP> [phyloP_scores.tsv.gz]

=over

=item add phyloP scores of all 3 datasets.

=back

=item I<dbSNP> [snp137.bed.gz]

=over

=item add rsID and dbSNP frequency information

=back

=item I<tgp> [tgp_phaseI_v3.bed.gz]

=over

=item add 1000 genomes allele frequency information

=back

=item I<cg54> [CG54.bed.gz]

=over

=item add 54 whole genomes allele frequency information from CompleteGenomics.

=back

=item I<wellderly> [wellderly.bed.gz]

=over

=item add wellderly's allele frequency information from CompleteGenomics.

=back

=item I<esp6500> [ESP6500.bed.gz]

=over

=item add ESP6500 allele frequency information from NHLBI

=back

=item I<short> [boolean]

=over

=item use 1 char peptide code instead of 3 chars format

=back

=item I<batch> [boolean]

=over

=item use batch mode annotation, default in daemon mode as an annotation engine.

=back

=item I<genes> [ "genes.list" | $rh_geneslist ]

=over

=item annotate transcripts for I<genes>. e.g. {"ABC" => 1, "DEF" => 1} or "genes.list" 

=back

=item I<trans> [ "trans.list" | $rh_translist ]

=over

=item annotate transcripts in I<trans>. e.g. {"NM_0012.1" => 1, "NM_0034.2" => 1} or "trans.list" 

=back

=item I<onlyPr> [boolean]

=over

=item limit to primary transcript for I<genes>, boolean option, default no limit.

=back

=item I<mmap> [boolean]

=over

=item allow annotating all other non "BEST" multiple-mapping records, boolen option, default not allowed.
      e.g. NM_0123.1 have 3 mapping location or alternate splicing, default only the "BEST" one will be 
      annotated. See L</DATABASE FORMAT> for the "BEST" definition.

=back

=back

=item Batch mode options :

=over

=item I<region> [region_string]

=over

=item limit to only annotate transcript in I<region>. e.g. "chr20:1234567-1234568", prior to I<regbed>.

=back
 
=item I<regbed> [BED format file]

=over

=item similar to I<region>, with allowing multiple regions. e.g. "in.region.bed". 

=back

=back

=item Notes

=over

=item Batch mode is designed for annotation of a complete list of variants 
      on same chromosome, read all information of the chr into memory, 
      and annotate all variants together in the order of chr coordinates.
      This mode can avoid frequent IO, brought by tabix searching, but need
      huge memory cost.

=back

=back

=back

=item Returns

=over

=item Annotation Engine object entry, please see L</load_anno> for more information.

=back

=back

=cut
sub new {
    my ( $class, @args ) = @_;
    my $self : shared;
    $self = shared_clone( {@args} );
    bless $self, ref($class) || $class;

    if (   ( !exists $self->{db} )
        or ( !-e $self->{db} )
        or ( !exists $self->{tr} )
        or ( !-e $self->{tr} ) )
    {
        $self->throw(
"Error filename or too few args, need db, tr file to be available at least."
        );
    }

    my $anno_all_opt = 1;

    my %open_args;
    if ( exists $self->{region} ) {
        $anno_all_opt = 0;
        $open_args{region} = $self->{region};
    }
    if ( !exists $self->{region} and exists $self->{regbed} ) {
        $anno_all_opt = 0;
        $open_args{regbed} = $self->{regbed};
    }
    if ( exists $self->{genes} ) {
        $anno_all_opt = 0;
        if ( ref( $self->{genes} ) eq 'HASH' ) {
            $open_args{genes} = $self->{genes};
        }
        else {
            open( GENE, $self->{genes} ) or $self->throw("$self->{genes} : $!");
            my %genes = map { s/\s+//g; $_ => 1 } <GENE>;
            close GENE;
            $open_args{genes} = \%genes;
        }
    }
    if ( exists $self->{trans} ) {
        $anno_all_opt = 0;
        if ( ref( $self->{trans} ) eq 'HASH' ) {
            $open_args{trans} = $self->{trans};
        }
        else {
            open( TRAN, $self->{trans} ) or $self->throw("$self->{trans} : $!");
            my %trans = map { s/\s+//g; $_ => 1 } <TRAN>;
            close TRAN;
            $open_args{trans} = \%trans;
        }
        $open_args{clean_trans} =
          { map { s/\-\d+$//; $_ => 1 } keys %{ $self->{trans} } };
    }
    if ( exists $self->{onlyPr} ) {
        $anno_all_opt = 0;
        $open_args{onlyPr} = $self->{onlyPr};
    }
    if ( exists $self->{mmap} ) {
        $open_args{mmap} = $self->{mmap};
    }
    if ( exists $self->{short} ) {
        %Code2Pep = %C1;
    }

    if ( exists $self->{cytoBand} ) {
        confess "Error: [$self->{cytoBand}] $!"
          if ( !-e $self->{cytoBand} or !-r $self->{cytoBand} );
        require GetCytoband;
        my $cytoBand_h = GetCytoBand->new( db => $self->{cytoBand} );
        $self->{cytoBand_h} = shared_clone($cytoBand_h);
    }

    if ( exists $self->{pfam} ) {
        confess "Error: [$self->{pfam}] $!"
          if ( !-e $self->{pfam} or !-r $self->{pfam} );
        require GetPfam;
        my $pfam_h = GetPfam->new( db => $self->{pfam} );
        $self->{pfam_h} = shared_clone($pfam_h);
    }

    if ( exists $self->{prediction} ) {
        confess "Error: [$self->{prediction}] $!"
          if ( !-e $self->{prediction} or !-r $self->{prediction} );
        require GetPrediction;
        my $prediction_h = GetPrediction->new( db => $self->{prediction} );
        $self->{prediction_h} = shared_clone($prediction_h);
    }

    if ( exists $self->{phyloP} ) {
        confess "Error: [$self->{phyloP}] $!"
          if ( !-e $self->{phyloP} or !-r $self->{phyloP} );
        require GetPhyloP;
        my $phyloP_h = GetPhyloP->new( db => $self->{phyloP} );
        $self->{phyloP_h} = shared_clone($phyloP_h);
    }

    if ( exists $self->{dbSNP} ) {
        confess "Error: [$self->{dbSNP}] $!"
          if ( !-e $self->{dbSNP} or !-r $self->{dbSNP} );
        require GetDBSNP;
        my $dbSNP_h = GetDBSNP->new( db => $self->{dbSNP} );
        $self->{dbSNP_h} = shared_clone($dbSNP_h);
    }

    if ( exists $self->{tgp} ) {
        confess "Error: [$self->{tgp}] $!"
          if ( !-e $self->{tgp} or !-r $self->{tgp} );
        require GetTgp;
        my $tgp_h = GetTgp->new( db => $self->{tgp} );
        $self->{tgp_h} = shared_clone($tgp_h);
    }

    if ( exists $self->{cg54} ) {
        confess "Error: [$self->{cg54}] $!"
          if ( !-e $self->{cg54} or !-r $self->{cg54} );
        require GetCGpub;
        my $cg54_h = GetCGpub->new( db => $self->{cg54} );
        $self->{cg54_h} = shared_clone($cg54_h);
    }

    if ( exists $self->{wellderly} ) {
        confess "Error: [$self->{wellderly}] $!"
          if ( !-e $self->{wellderly} or !-r $self->{wellderly} );
        require GetCGpub;
        my $wellderly_h = GetCGpub->new( db => $self->{wellderly} );
        $self->{wellderly_h} = shared_clone($wellderly_h);
    }

    if ( exists $self->{esp6500} ) {
        confess "Error: [$self->{esp6500}] $!"
          if ( !-e $self->{esp6500} or !-r $self->{esp6500} );
        require GetPfam;
        my $esp6500_h = GetPfam->new( db => $self->{esp6500} );
        $self->{esp6500_h} = shared_clone($esp6500_h);
    }

    $self->{trInfo} = shared_clone( $self->readtr() );
    $self->{annodb} = shared_clone( $self->load_anno(%open_args) )
      if ( exists $self->{batch} );
    return $self;
}

=head2 readtr

    About   : Read transcript information, and even sequences if in batch mode.
    Usage   : my $rtrSeqs = $beda->readtr( genes => $rh_NewGenes, trans => $rh_NewTrans );
    Args    : Optional args "genes" and "trans" only accept hash ref values.
	      if no args specified, it will load information based on the
	      configuration of BedAnno entry.
    Returns : A hash ref of trSeqs:
	      {
		$tr_acc => {
		    len      => $tr_len,
		    gene     => $gene_sym,

		    # optional tags:
		    prot     => $prot_acc,
		    plen     => $prot_len,
		    csta     => $cds_start_on_trSeq, # 1 based
		    csto     => $cds_end_on_trSeq,   # 1 based
		    seq      => $tr_sequence,
		    
		    X	     => 1,		     # inseqStop
		    U	     => 1,		     # selenocysteine
		    A	     => 1,		     # polyATail
		    altstart => $startCodons,	     # alt start codons
		},
		...
	      }

=cut
sub readtr {
    my $self = shift;
    my %opts = @_;
    if ( exists $opts{genes} and ref( $opts{genes} ) ne 'HASH' ) {
        $self->throw("Options arg 'genes' only accept hash ref as value.");
    }
    if ( exists $opts{trans} and ref( $opts{trans} ) ne 'HASH' ) {
        $self->throw("Options arg 'trans' only accept hash ref as value.");
    }

    open( FAS, "zcat -f $self->{tr} |" ) or confess "$self->{tr}: $!";
    local $/ = ">";
    my %seqs = ();
    while (<FAS>) {
        s/[\s>]+$//g;
        next if (/^\s*$/);
        my $hd = $1 if (s/^(\S+[^\n]*)\n//);
        confess "Error: trSeq parse error!" if ( !defined $hd );
        my @headers = split( /\s+/, $hd );
        confess "Error: trSeq header parse error!" if ( 7 > @headers );
        s/\s+//g;
        next if ( $headers[6] =~ /FAIL/ );    # skip failed transcript

        if ( !exists $opts{trans} and !exists $opts{genes} ) {
            next
              unless (
                ( !exists $self->{clean_trans} and !exists $self->{genes} )
                or (    exists $self->{clean_trans}
                    and exists $self->{clean_trans}->{ $headers[0] } )
                or (    exists $self->{genes}
                    and exists $self->{genes}->{ $headers[2] } )
              );
        }
        else {
            next
              unless (
                ( exists $opts{trans} and exists $opts{trans}{ $headers[0] } )
                or
                ( exists $opts{genes} and exists $opts{genes}{ $headers[2] } )
              );
        }
        $seqs{ $headers[0] }{seq} = $_
          if ( exists $self->{batch} );    # hash sequence when batch
        $seqs{ $headers[0] }{len}  = $headers[1];    # tx length
        $seqs{ $headers[0] }{gene} = $headers[2];    # gene symbol
        $seqs{ $headers[0] }{prot} = $headers[3]
          if ( $headers[3] ne "." );                 # prot acc.ver
        $seqs{ $headers[0] }{plen} = $headers[4]
          if ( $headers[4] ne "." );                 # prot length
        $seqs{ $headers[0] }{U} = 1
          if ( $headers[6] =~ /selenocysteine/ );    # selenocysteine
        $seqs{ $headers[0] }{X} = 1
          if ( $headers[6] =~ /inseqStop/ );         # inseqStop
        $seqs{ $headers[0] }{A} = 1
          if ( $headers[6] =~ /polyATail/ );         # polyATail

        if ( $headers[5] ne "." ) {                  # cds start, cds end in tx
            @{ $seqs{ $headers[0] } }{qw{csta csto}} =
              split( /,/, $headers[5] );
        }
        $seqs{ $headers[0] }{altstart} =
          { map { $_ => 1 } split( /;/, $headers[7] ) }
          if ( defined $headers[7] );
    }
    close FAS;
    return \%seqs;
}

=head2 load_anno

    About   : load all needed annotation infomation into memory for multi-process annotation
    Usage   : my $ranndb = $beda->load_anno( region => "chr20:1234567-1234568", trans => \%trans );
    Args    : Using %args to override class's properties: region, regbed, genes, trans
	      if no args, use the the class's properties as default.
    Returns : a localized merged anno db, The returned annotation database is a hash ref.
		{
		    $chr => [
			{
			    sta   => $start, (0 based)
			    sto   => $stop,  (1 based)
			    annos => {
				$anno_string => $offset, ...
			    }

			    detail => {
				$tid => {
				    gsym => $gsym,    (gene symbol)
				    gid	 => $gid,     (Entrez gene id)
				    gpSO => $gpSO,    (GeneParts SO)
				    blka => $blka,    (block attribute)
				    exin => $exin,    (exon intron number)
				    nsta => $nsta,    (n./r. left  of whole block)
				    nsto => $nsto,    (n./r. right of whole block)
				    csta => $csta,    (c.    left  of whole block)
				    csto => $csto,    (c.    right of whole block)
				    wlen => $wlen,    (length of whole block)
				    pr   => $pr,      (primary tag)
				    strand => $strd,  (strand)
				    offset => $offset,(offset of current block to whole block)
				    mismatch => $mismatch (non-equal length block descripter)
				}, ...
			    }
			}, ... 
		    ], ...
		}
	      Note: when variation hit one of the annotation entry, the anno_string will be parsed.
	      and the "detail" tag will be added then.

=cut
sub load_anno {
    my ( $self, %args ) = @_;
    my $cmd        = "$self->{db}";
    my $tabix_args = qq['$self->{db}'];
    if ( -e $self->{db} && $self->{db} =~ /\.gz/i ) {
        if ( exists $args{region} and defined $args{region} ) {
            $tabix_args .= qq[ '$args{region}'];
            $cmd = "tabix $tabix_args |";
        }
        elsif ( exists $args{regbed}
            and defined $args{regbed}
            and -e $args{regbed} )
        {
            $tabix_args = '-B ' . $tabix_args . qq[ '$args{regbed}'];
            $cmd        = "tabix $tabix_args |";
        }
        else {
            $cmd = "zcat -f '$self->{db}' |";
        }
    }
    open( ANNO, $cmd ) or $self->throw("$cmd: $!");

    # trans filter is always be ahead of genes
    my $prTag   = ( exists $args{onlyPr} ) ? 1 : 0;
    my $mmapTag = ( exists $args{mmap} )   ? 1 : 0;
    my $geneTag = ( exists $args{genes} )  ? 1 : 0;
    my $tranTag = ( exists $args{trans} )  ? 1 : 0;
    my $geneList = $args{genes} if ($geneTag);
    my $tranList = $args{trans} if ($tranTag);
    my %pureTran = map { s/\-\d+$//; $_ => 1 } keys %$tranList if ($tranTag);

    my $rannodb = {};
    while (<ANNO>) {
        s/\s+$//;
        my ( $chr, $start, $stop, $annostr ) = split(/\t/);
        my @annos = split( /; /, $annostr );

        #		 0-based	1-based
        my %ent = ( sta => $start, sto => $stop );
        foreach my $anno_ent (@annos) {
            my @cont = split( /\|/, $anno_ent );

            # no multiple mapping transcripts if !$mmapTag
            next
              if (  $cont[0] =~ /\-\d+$/
                and !$mmapTag
                and ( !$tranTag or !exists( $$tranList{ $cont[0] } ) ) );
            next
              if (  $prTag
                and ( $cont[13] ne 'Y' )
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

            my $ofst = $1
              if ( $anno_ent =~ s/\|(\d+)$// )
              or $self->throw("db format error: [$anno_ent]");
            $ent{annos}{$anno_ent} = $ofst;
        }

        next if ( !exists $ent{annos} );
        push( @{ $$rannodb{$chr} }, {%ent} );
    }
    close ANNO;

    $rannodb = region_merge($rannodb)
      if ( $prTag or $mmapTag or $geneTag or $tranTag );

    return $rannodb;
}

=head2 region_merge

    About   : merge consecutive same-entries regions
    Usage   : my $rannodb = region_merge($loaded_db);
    Args    : A hash ref of loaded_db.
    Returns : A hash ref of merged db.

=cut
sub region_merge {

    my $radb = shift;
    my %local_radb = %$radb;

    foreach my $chr (keys %local_radb) {
	my @annoents = @{$local_radb{$chr}};
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

	$local_radb{$chr} = [@annoents] if ($oricount > (scalar @annoents));
    }

    return \%local_radb;
}

# just parse annoents when mutation hits.
sub parse_annoent {
    my $annoent = shift;
    my %annoinfo = ();
    # tid         gsym  gid strand blka gpgo exin nsta nsto csta csto wlen mismatch    pr
    # NM_015658.3|NOC2L|26155|-|3U1E|205|EX19E|2817|2801|*508|*492|0|D,879582,879582,.|Y
    my @infos = split(/\|/, $annoent);
    my $tid = shift @infos;
    confess "Error format of anno ents [$annoent]" if (10 != @infos);
    my @tags  = qw(gsym gid strand blka gpSO exin nsta nsto csta csto wlen mismatch pr);
    @annoinfo{@tags} = @infos;
    return ($tid, \%annoinfo);
}

=head2 anno

    About   : Annotate single short variation by annotation db.
    Usage   : my $anno_ent = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );
	      or $anno_ent = $beda->anno( 'chr20', 1234568, 'AG', 'AGGG' );
    Args    : for CG's shell variants, need 5 args in UCSC coordinates (0-based start):
		chr id, chr start, chr end, reference, alternative.
	      for variants in VCF, need 4 args, which is lack of 
	      chr end, in 1-based coordinates.
    Returns : a hash ref of annotation informations, see varanno().

=cut
sub anno {
    my $self = shift;
    my $var = parse_var(@_);
    return $self->varanno($var);
}

=head2 varanno

    About   : generate all the needed annotation for a var entry
    Usage   : my $rAnnoRst = $beda->varanno($var);
    Args    : see parse_var()
    Returns : a hash ref:
                {
                    var => {

                        # the first part is from var parsing result.
                        # please see parse_var().
                        # import all the keys from original parsed var entry
			# and add the following keys by this method.

                        # information
                        cytoBand  => $cytoBand,
                        varTypeSO => $varTypeSO,
                        gHGVS     => $gHGVS,
                        refbuild  => $referenceBuild,

                        # Here's some optional parts which may be generated
                        # when extra resource is available:
                        dbsnp => {
                            $rsID => {
                                AN => $dbsnp_total_allele_count,
                                AF => $dbsnp_alt_allele_frequency,  # though ref
                            },
                            ...
                        },

                        tgp => {
                            AN => $tgp_total_allele_count,
                            AF => $tgp_alt_allele_frequency,
                        },

                        cg54 => {
                            AN => $cg54_total_allele_count,
                            AF => $cg54_alt_allele_frequency,
                        },

                        wellderly => {
                            AN => $wellderly_total_allele_count,
                            AF => $wellderly_alt_allele_frequency,
                        },

                        esp6500 => {
                            AN => $esp6500_total_allele_count,
                            AF => $esp6500_alt_allele_frequency,
                        },
                    },
                    trInfo => {
                        $tid => {
                            geneId        => $Entrez_Gene_ID,
                            geneSym       => $Gene_Symbol,
                            prot          => $Protein_Acc_Ver,
                            strd          => $strand,
                            rnaBegin      => $Begin_in_RNA_transcript,
                            rnaEnd        => $End_in_RNA_transcript,
                            protBegin     => $Begin_in_Protein,
                            protEnd       => $End_in_Protein,
                            c             => $cHGVS,
                            p             => $pHGVS,
                            cc            => $codon_change,
                            polar         => $polar_change,
                            r             => $imp_funcRegion,
                            func          => $imp_funcCode,
                            exin          => $exIntr_number,
                            genepart      => $GenePart,
                            genepartSO    => $GenePartSO,
                            genepartIndex => $GenePartIndex,
                            exonIndex   => $exonIndex,          # '.' for N/A
                            intronIndex => $intronIndex,        # '.' for N/A
                            funcSOname  => $FunctionImpact,
                            funcSO      => $FunctionImpactSO,
                            pfamId      => $PFAM_ID,
                            pfamName    => $PFAM_NAME,
                            phyloPpm    => $PhyloPscorePlacentalMammals,
                            phyloPpr    => $PhyloPscorePrimates,
                            phyloPve    => $PhyloPscoreVetebrates,
                            siftPref    => $SIFTpred,
                            siftScore   => $SIFTscore,
                            pp2divPred  => $Polyphen2HumDivPred,
                            pp2divScore => $Polyphen2HumDivScore,
                            pp2varPred  => $Polyphen2HumVarPred,
                            pp2varScore => $Polyphen2HumVarScore,
                        },
                        ...
                    }
                }

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

=head2 get_gHGVS

    About   : get genomics (chromosomal) HGVS string of variation
    Usage   : my $gHGVS = get_gHGVS($var);
    Args    : variation entry, after parse_var(), see parse_var().
    Returns : chromosomal HGVS string.

=cut
sub get_gHGVS {
    my $var = shift;
    my $gHGVS = 'g.';
    if ($var->{chr} =~ /M/) { # hit mito chromosome
	$gHGVS = 'm.';
    }

    $_ = $$var{guess};
    if ($_ eq 'snv') {
	$gHGVS .= $$var{pos}.$$var{ref}.'>'.$$var{alt};
    }
    elsif ($_ eq 'ins') {
	$gHGVS .= $$var{pos}.'_'.($$var{pos}+1).'ins'.(substr($$var{alt},1));
    }
    elsif ($_ eq 'del' or $_ eq 'delins') {
	# 1bp del/delins
	$gHGVS .= ($$var{pos}+1);
	if ($$var{reflen} > 2) {
	    $gHGVS .= '_'.($$var{pos} + $$var{reflen} - 1);
	}
	$gHGVS .= 'del'.(substr($$var{ref}, 1));
	$gHGVS .= 'ins'.(substr($$var{alt}, 1)) if (/delins/);
    }
    elsif ($_ eq 'rep') {
	$gHGVS .= ($$var{pos} + 1);
	if ($$var{ref_cn} == 1 and $$var{alt_cn} == 2) { # dup
	    if ($$var{replen} > 1) {
		$gHGVS .= '_' . ($$var{pos} + $$var{replen});
	    }
	    $gHGVS .= 'dup'. $$var{rep};
	}
	else {
	    $gHGVS .= $$var{rep}.'['.$$var{ref_cn}.'>'.$$var{alt_cn}.']';
	}
    }
    elsif ($_ eq 'ref') {
	$gHGVS .= '=';
    }
    else {
	confess "Can not recognize type $$var{guess}.";
    }

    return $gHGVS;
}


=head2 pairanno

    About   : get infomations from select_position(), (std grp) besides function,
	      codon-change, flank_region position, polar and HGVS strings:
		  for coding gene region gives c. HGVS and p. HGVS
		  for non-coding RNA gives n. HGVS, 
		  for non-available gives '.'.
	      function of mutation, follow the definitions mainly from dbSNP:
		  for intergenic or non-target region: using '.'
		  utr-3, utr-5, ncRNA, 
		  misstart, missense, nonsense, coding-synon, 
		  intron, splice-3, splice-5,
		  cds-indel, frameshift, stop-gain, stop-loss, init-loss, 
		  unknown for 'c.?'
		  abnormal-inseq-stop for non-end stop codon in refseq.
		  abnormal-intron for too short intron (<4bp) due to mapping deduced.
	      Sub-region and Exon/Intron number infomation
    Usage   : my $h_annos = $beda->pairanno($var, $sel);
    Args    : a var entry and a selected anno
    Returns : ( $tid, { c => c.HGVS, p => p.HGVS,  cc => condon-change,
			r => region, exin => exin, func => func,  polar => pol, 
			strd => [+-],
			flanks => { l => left_region, r => right_region }
	              } )

=cut
sub pairanno {
    my ($self, $var, $sel) = @_;

    # annotated by paired (sta, sto) cPos
    my ($tid, $cL, $cR) = @$sel;
    return ( $tid, {} ) if (!exists $$cL{cpos});

    my ( $ref, $alt, $pos, $rl, $al ) =
      ( $$var{ref}, $$var{alt}, $$var{pos}, $$var{reflen}, $$var{altlen} );
    if (exists $$var{'+'} and exists $$var{'+'}{p}) {
	# for diff strand position cases
	$ref = $$var{$$cL{strd}}{r};
	$alt = $$var{$$cL{strd}}{a};
	$pos = $$var{$$cL{strd}}{p};
	$rl  = $$var{$$cL{strd}}{rl};
	$al  = $$var{$$cL{strd}}{al};
    }

    if ($rl > 1 or $al > 1 and $rl != $al) { # get the real diff base in ref
	$ref = substr($ref, 1);
	$alt = substr($alt, 1);
	$rl --;
	$al --;
	$pos ++;
    }

    # get flank region string ready for fetchseq
    my $flanks = get_flanks($$var{chr}, $pos, $rl);

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
		strd   => $$cL{strd},
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
    # get cPos without c./n. 
    my $cP_L  = substr( $$cL{cpos}, 2 );
    my $cP_R  = substr( $$cR{cpos}, 2 );
    my $pre_L = substr( $$cL{cpos}, 0, 2 );
    my $pre_R = substr( $$cR{cpos}, 0, 2 );

    my $coding_max = int(length($$self{codonseq}{$query_tid})/3) if (exists $$self{codonseq}{$query_tid});

    $_ = $$var{guess};

    if ($_ eq 'snv') { # only one position
	$cHGVS = $$cL{cpos}.$t_ref.'>'.$t_alt; # c.123A>G
	$reg = $$cL{reg};
	$exin = $$cL{exin};
	($func, $cc, $pHGVS, $polar) = $self->parse_cPos($tid, $$cL{cpos}, $t_alt);
	if ($$cL{reg} =~ /^I/ and $$cL{bd} =~ /1/) {
	    $func = 'abnormal-intron';
	}
    }
    elsif ($_ eq 'ins') {

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
		    strd   => $$cL{strd},
		    bc     => $bc_cHGVS,
		    flanks => $flanks
		}
	    );
	}

	if ($$cL{reg} eq $$cR{reg}) { # non-edge insertion

	    $reg  = $$cL{reg};
	    $exin = $$cL{exin};

	    if ($reg =~ /^C/) { # coding region
		( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_L, $t_alt, 1);
	    }
	    elsif ($reg =~ /^I/) {
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
	    elsif ($reg =~ /^5/) {
		$func = 'utr-5';
	    }
	    elsif ($reg =~ /^3/) {
		$func = 'utr-3';
	    }
	    elsif ($reg =~ /^R/) {
		$func = 'ncRNA';
	    }
	    else { $self->throw("unexpected region string [$reg]"); }
	}
	else { # insertion on the edge of region
	    $reg = $$cL{reg}.'-'.$$cR{reg};
	    $exin = ( $$cL{exin} eq $$cR{exin} ) ? $$cL{exin} : $$cL{exin}.'-'.$$cR{exin};

	    my $indicator = substr($$cL{reg},0,1).substr($$cR{reg},0,1);

	    if ($indicator eq 'IC') { # 5' coding exon insertion
		( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_R, $t_alt, -1);
	    }
	    elsif ($indicator eq 'CI') { # 3' conding exon insertion
		( $pHGVS, $func ) = $self->get_aaInsInfo($query_tid, $cP_L, $t_alt, 1);
	    }
	    elsif ($indicator =~ /5/) {
		$func = 'utr-5';
	    }
	    elsif ($indicator =~ /3/) {
		$func = 'utr-3';
	    }
	    elsif ($indicator =~ /R/) {
		$func = 'ncRNA';
	    }
	    else { $self->throw("unexpected indicator [$indicator]"); }
	}
    }
    elsif ( $_ eq 'del' or $_ eq 'delins' ) {

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
	    $reg =
	      ( $$cL{reg} eq $$cR{reg} )
	      ? $$cL{reg}
	      : $$cL{reg} . '-' . $$cR{reg};
	    $exin =
	      ( $$cL{exin} eq $$cR{exin} )
	      ? $$cL{exin}
	      : $$cL{exin} . '-' . $$cR{exin};

	    if (($pre_L eq '--') and ($pre_R eq '++')) { # whole transcript is missing
		$pHGVS = 'p.0?';
		$func  = 'knockout'; # which may not be exists in variation calling result.
	    }
	    elsif ($pre_L eq '--') { # the 5' end of transcript is missing
		if ($$cR{reg} =~ /^5/) {
		    $func = 'utr-5';
		}
		elsif ($$cR{reg} =~ /^C/) {
		    $pHGVS = 'p.0?';
		    $func  = 'init-loss';
		}
		elsif ($$cR{reg} =~ /^3/) {
		    $pHGVS = 'p.0?';
		    $func  = 'knockout';
		}
		elsif ($$cR{reg} =~ /^R/) {
		    $func  = 'ncRNA';
		}
		elsif ($$cR{reg} =~ /^I/) {
		    if ($pre_R eq 'n.') {
			$func = 'ncRNA';
		    }
		    elsif ($cP_R =~ /^\-|^1\-\d/) {
			$func = 'utr-5';
		    }
		    elsif ($cP_R =~ /^\d+/) {
			$pHGVS = 'p.0?';
			$func  = 'init-loss';
		    }
		    elsif ($cP_R =~ /^\*/) {
			$pHGVS = 'p.0?';
			$func  = 'knockout';
		    }
		    else {
			$self->throw("unexpected cP_R [$cP_R]");
		    }
		}
		else {
		    $self->throw("unexpected region [$$cR{reg}]");
		}
	    }
	    else { # $pre_R eq '++' the 3'end of transcript is missing.
		if ($$cL{reg} =~ /^5/) {
		    $pHGVS = 'p.0?';
		    $func = 'knockout';
		}
		elsif ($$cL{reg} =~ /^C/) {
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
		elsif ($$cL{reg} =~ /^3/) {
		    $func  = 'utr-3';
		}
		elsif ($$cL{reg} =~ /^R/) {
		    $func  = 'ncRNA';
		}
		elsif ($$cL{reg} =~ /^I/) {
		    if ($pre_L eq 'n.') {
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
		else {
		    $self->throw("unexpected region [$$cL{reg}]");
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
		    strd   => $$cL{strd},
		    bc     => $bc_cHGVS,
		    flanks => $flanks
		}
	    );
	    
	}

	my $indicator = substr($$cL{reg},0,1).substr($$cR{reg},0,1);

	if ( $$cL{reg} eq $$cR{reg} ) {
	    $reg  = $$cL{reg};
	    $exin = $$cL{exin};

	    if ($indicator eq 'CC') {
		($pHGVS, $func) = $self->get_aaDelInfo($query_tid, $cP_L, $cP_R, $t_alt);
	    }
	    elsif ($indicator eq 'II') {
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
	    elsif ($indicator eq '55') {
		$func = 'utr-5';
	    }
	    elsif ($indicator eq '33') {
		$func = 'utr-3';
	    }
	    elsif ($indicator eq 'RR') {
		$func = 'ncRNA';
	    }
	    else {
		$self->throw("unexpected region [$$cL{reg}]");
	    }
	}
	else { # multiple region del/delins
	    $reg  = $$cL{reg} . '-' . $$cR{reg};
	    $exin = ( $$cL{exin} eq $$cR{exin} ) ? $$cL{exin} : $$cL{exin} . '-' . $$cR{exin};

	    if ($indicator eq 'CC') {
		($pHGVS, $func) = $self->get_aaDelInfo($query_tid, $cP_L, $cP_R, $t_alt);
	    }
	    elsif ($indicator eq '55') {
		$func = 'utr-5';
	    }
	    elsif ($indicator eq '33') {
		$func = 'utr-3';
	    }
	    elsif ($indicator eq 'RR') {
		$func = 'ncRNA';
	    }
	    elsif ($indicator eq 'II') {
		if ($$cL{bd} =~ /1/) {
		    $func = 'abnormal-intron';
		}
		elsif ($$cL{bd} =~ /0/ and $$cR{bd} =~ /0/) {
		    if ($pre_L eq 'n.') {
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
	    elsif ($indicator eq '5C') {
		$pHGVS = 'p.0?';
		$func  = 'init-loss';
	    }
	    elsif ($indicator eq 'C3') {
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
	    elsif ($indicator eq '53') {
		$pHGVS = 'p.0?';
		$func  = 'knockout';
	    }
	    else { # for other complex cases
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
    # when meet repeats, the input selection must come from
    # the standard group: (repeat element cds position), 
    # the backward compatible cds position is used to calculate
    # the bcHGVS, pHGVS and give the function prediction as well.
    # the standard cds position is only for the cHGVS.
    elsif ( $_ eq 'rep' ) {
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
		$bc_var{guess} = guess_type();
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
	    }
	}
    }
    elsif ($_ eq 'ref') {
        $cHGVS = ( $$var{chr} =~ /M/ ) ? 'm.=' : 'g.=';
        $reg =
          ( $$cL{reg} eq $$cR{reg} ) ? $$cL{reg} : $$cL{reg} . '-' . $$cR{reg};
        $exin =
          ( $$cL{exin} eq $$cR{exin} )
          ? $$cL{exin}
          : $$cL{exin} . '-' . $$cR{exin};
    }
    else { $self->throw("unexpected type of variation [$$var{guess}]."); }

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
	    strd   => $$cL{strd},
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
	my $add_aa = (exists $Code2Pep{$add_codon}) ? $Code2Pep{$add_codon} : "";
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
    Args    : query_tid is transcript id with out tail [-n], cpos is cds pos in hgvs c./n. format.
	      t_alt are strand specific alt char.
    Returns : the function code, codon-change, 'p.' format HGVS string, polar-change. see pairanno()

=cut
sub parse_cPos {
    my ($self, $query_tid, $cpos, $t_alt) = @_;
    my ($func, $cc, $pHGVS, $polar) = ('.') x 4;

    if ($cpos =~ /^c\.(\d+)$/) { # coding region
	my $cP = $1;
	my ($codon, $aa, $pol, $frame) = get_codon($$self{codonseq}, $query_tid, $cP);
	my $new_codon = $codon;
	if ($aa ne '') {
	    substr($new_codon, $frame, 1, $t_alt);
	    my $new_aa = (exists $Code2Pep{$new_codon}) ? $Code2Pep{$new_codon} : ".";
	    my $new_pol = (exists $Polar{$new_aa}) ? $Polar{$new_aa} : ".";
	    $func = $self->get_aafunc($aa, $new_aa, $query_tid, $cP);
	    $pHGVS = get_aaHGVS($aa, $new_aa, $func, $cP);

	    $pHGVS .= 'ext*?' if ($func eq 'stop-loss');

	    $cc = $codon.'=>'.$new_codon;
	    $polar = ($pol eq $new_pol) ? '.' : $pol.'=>'.$new_pol;
	    if ($aa eq '*' and $cP + 3 <= length($$self{codonseq}{$query_tid})) {
		$func = 'abnormal-inseq-stop';
	    }
	}
    }
    elsif ($cpos =~ /^[cr]\.[\-\*]?\d+([\-\+])(\d+)/) { # intron region
	if ($2 > 2) {
	    $func = 'intron';
	}
	elsif ($1 eq '+') {
	    $func = 'splice-5';
	}
	elsif ($1 eq '-') {
	    $func = 'splice-3';
	}
    }
    elsif ($cpos =~ /^c\.-\d+$/) { # 5' UTR
	$func = 'utr-5';
    }
    elsif ($cpos =~ /^c\.\*\d+$/) { # 3' UTR
	$func = 'utr-3';
    }
    elsif ($cpos =~ /^n\.\d+$/) { # ncRNA
	$func = 'ncRNA';
    }
    else { $func = 'unknown'; }

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
    if ($ins_aa =~ /^$/) { #frameshit or cds deletion without insertion 
	if ($pP_L == 1) {

	    $pHGVS = 'p.0?';
	    $func  = 'init-loss';
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
    elsif ($ins_aa =~ /\*/) { # stop gain after delins
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
    else {
	if ($todel_aa eq $ins_aa) { # frameshift or delins 

	    $pHGVS = ($fsOpt) ? 'p.' . $aa_R . $pP_R . 'fs*?' : 'p.=';
	    $func  = ($fsOpt) ? 'frameshift' : 'coding-synon';

	}
	elsif ($pP_L == 1 and substr($todel_aa, 0, 3) ne substr($ins_aa, 0, 3)) {

	    $pHGVS = 'p.0?';
	    $func  = 'init-loss';

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
	    $func  = 'init-loss';

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
    if ($func eq 'misstart') {
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
	return ((!exists $Code2Pep{$codon}) ? "" : $Code2Pep{$codon});
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
	$func = 'misstart';
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
	return ("", "", "", $frame);
    }
    my $codon = uc(substr($$rseq{$tid}, ($cpos - $frame - 1), 3));
    my $aa = (exists $Code2Pep{$codon}) ? $Code2Pep{$codon} : "";
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
    Returns : var with standard {sel}{std}, and possible backward compatible {sel}{bc}

=cut
sub select_position {
    my ($self, $var, $rcPos) = @_;

    $$var{sel} = {};
    if (!exists $$var{'+'} or (!exists $$var{'+'}{p} and $$var{guess} ne 'rep')) {
    # simple variation or same pos(strand) non-repeat variation which caused by complex delins or multiple samples
	my $anno_sels = [];

	$_ = $$var{guess};
	if ($_ eq 'snv') { # only select pos for the case: c.123A>G
	    if (exists $$rcPos{$$var{pos}}) {
		foreach my $tid (sort keys %{$$rcPos{$$var{pos}}}) {
		    push (@$anno_sels, [ $tid, $$rcPos{$$var{pos}}{$tid} ]);
		}
	    }
	}
	elsif ($_ eq 'ins') { # the pos and pos+1 are selected: c.123_124insTG
	    $anno_sels = $self->get_cover($$var{chr}, $$var{pos}, ($$var{pos}+1));
	}
	elsif ($_ eq 'del' or ($_ eq 'delins' and $$var{reflen} != $$var{altlen})) { # the pos+1 and pos+reflen-1 are selected
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
	elsif ($_ eq 'delins' and $$var{reflen} == $$var{altlen}) { # substitution case for CompleteGenomics.
	    $anno_sels = $self->get_cover(
		$$var{chr}, $$var{pos}, ($$var{pos} + $$var{reflen} - 1)
	    );
	}
	elsif ($_ eq 'ref') {
	    if ($$var{reflen} == 1 and exists $$rcPos{$$var{pos}}) {
		foreach my $tid (sort keys %{$$rcPos{$$var{pos}}}) {
		    push (@$anno_sels, [ $tid, $$rcPos{$$var{pos}}{$tid} ]);
		}
	    }
	    else {
		$anno_sels = $self->get_cover($$var{chr}, $$var{pos}, ($$var{pos} + $$var{reflen} - 1));
	    }
	}
	else { $self->throw("Error: unexpected guess for normal case: $$var{guess}"); }

	$$var{sel}{std} = $anno_sels;
    }
    elsif (exists $$var{'+'}{p}) { # strand different pos
	my ($f_annos, $r_annos) = ([], []);
	my @sel_std = ();
	$_ = $$var{guess}; # modified repeat caused ambiguous
	if ($_ eq 'ins') {
	    $f_annos = $self->get_cover( $$var{chr}, $$var{'+'}{p},
		( $$var{'+'}{p} + 1 ) );
	    $r_annos = $self->get_cover( $$var{chr}, $$var{'-'}{p},
		( $$var{'-'}{p} + 1 ) );
	}
	elsif ($_ eq 'del' or $_ eq 'delins') {
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
	else { $self->throw("Error: unexpected guess for different fr pos $$var{guess}"); }
	
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
    Usage   : $beda = BedAnno->new( db => 'in.bed.gz', codon => 'in.codon.fa'); 
	      @all_annos = ();
	      foreach chr {
		my @vars = (); 
		foreach pos, ref, @alt {
		  $var = parse_var($chr, $pos, $ref, $alt);
		  push (@vars, $var);
		}
		my $rChrAnnos = $beda->batch_anno( chr, \@vars); 
		push (@all_annos, @$rChrAnnos);
	      }
    Args    : chr, and an array ref of vars.
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
    elsif ($$var{guess} eq 'snv' or $$var{guess} eq 'ref') {
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

    $_ = $$rh{blka};
    if (/^C/ or /^5/) { # cds or utr-5
	if ($strandopt) {
	    $cpos = 'c.'.($$rh{csta} + $total_left_offset);
	}
	else {
	    $cpos = 'c.'.($$rh{csto} + $total_right_offset);
	}
    }
    elsif (/^I/) { # intron
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
		$cpos = 'n.'.$$rh{nsta}.$opt.($total_left_offset + 1);
	    }
	    else {
		my $opt = ($strandopt) ? '-' : '+';
		$cpos = 'n.'.$$rh{nsto}.$opt.($total_right_offset + 1);
	    }
	}
    }
    elsif (/^3/) { # utr-3
	if ($strandopt) {
	    my $p = $$rh{csta}; $p =~ s/^\*//;
	    $cpos = 'c.*'.($p + $total_left_offset);
	}
	else {
	    my $p = $$rh{csto}; $p =~ s/\*//;
	    $cpos = 'c.*'.($p + $total_right_offset);
	}
    }
    elsif (/^R/) { # ncRNA
	if ($strandopt) {
	    $cpos = 'n.'.($$rh{nsta} + $total_left_offset);
	}
	else {
	    $cpos = 'n.'.($$rh{nsto} + $total_right_offset);
	}
    }
    else { confess "unrecognized block attribution!"; }

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


=head2 get_cover_batch

    About   : get covered region in batch mode
    Usage   : my $cover_href = $beda->get_cover_batch( $chr, \@stasto );
    Args    : a chromosome id, and an array ref of [ [ $start, $stop ], ... ]
    Returns : a hash ref of:
		{
		    "$start-$stop" => [ 
					[ $tid, $left_cPos, $right_cPos ], ... 
				      ], ...
		}
	      Note: the pos-pair which do not hit any annotation blocks, will
		    not exist in the returned results.

=cut
sub get_cover_batch {
    my ($self, $chr, $stasto_aref) = @_;
    my $rAnnos = $$self{annodb};
    return {} if (!exists $$rAnnos{$chr});

    my @sorted_stasto = sort pairsort @$stasto_aref;
    return {} if (0 == @sorted_stasto);

    my %ret_cov = ();

    my $cur_blkId = 0;
    for (my $i = 0; $i < @{$$rAnnos{$chr}}; $i++) {
	my $cur_bedent = $$rAnnos{$chr}[$i];

	last if ($cur_blkId >= @sorted_stasto);
	# skip left anno blocks
	next if ($sorted_stasto[$cur_blkId][0] > $$cur_bedent{sto});

	# skip no hit pos-pair
	while ($sorted_stasto[$cur_blkId][1] <= $$cur_bedent{sta}) {
	    $cur_blkId ++;
	    return \%ret_cov if ($cur_blkId == @sorted_stasto);
	}

	# skip left anno blocks again for the new shifted pos-pair.
	next if ($sorted_stasto[$cur_blkId][0] > $$cur_bedent{sto});

	# hit the annotation blks
	my $pospair = join( "-", @{$sorted_stasto[$cur_blkId]} );
	$ret_cov{$pospair} = cal_covered( $$rAnnos{$chr}, $i, @{$sorted_stasto[$cur_blkId]} );
	$cur_blkId ++; 
	$i --; # keep annotation block to stay at current index
    }
    return \%ret_cov;
}

sub cal_covered {
    my ($rchrAnnos, $k, $start, $stop) = @_;

    my %min_left = ();
    my %min_right = ();
    my %left_most = ();
    my %right_most = ();

    while ($k < @$rchrAnnos and $$rchrAnnos[$k]{sta} < $stop) {
	my $bedent = $$rchrAnnos[$k];
	if (!exists $$bedent{detail}) {
	    foreach my $annoblk (keys %{$$bedent{annos}}) {
		my $offset1 = $$bedent{annos}{$annoblk};
		my ($tid, $ranno) = parse_annoent($annoblk);
		$$bedent{detail}{$tid} = $ranno;
		$$bedent{detail}{$tid}{offset} = $offset1;
	    }
	}
	
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

    my @ret_pairs = ();
    foreach my $tid (sort keys %left_most) {
	my ($left_cPos, $right_cPos);
	my ( $left_bed, $right_bed ) =
	  ( $$rchrAnnos[ $left_most{$tid} ], $$rchrAnnos[ $right_most{$tid} ] );

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

sub pairsort {
    $$a[0] <=> $$b[0] or $$a[1] <=> $$b[1]
}


=head2 parse_var

    About   : parse the variation directly by the ref and alt string
    Usage   : my $var = parse_var( $chr, $start, $end, $ref, $alt );
	   or my $var = parse_var( $chr, $pos, $ref, $alt );
    Returns : a hash ref of variation:
	      {
		chr	    => $chr,
		pos	    => $start,	    # 0-based start
		ref	    => $ref,
		alt	    => $alt,
		reflen	    => $ref_len,
		altlen	    => $alt_len,    # not exists if no-call
		guess	    => $varType,    # the output varType
		imp	    => $imp_varType,# the implicit varType
		sm	    => $sm,	    # single/multiple base indicator
					    # equal/non-equal length indicator

		# for hgvs naming convinient, reparse it
		# in forward and reverse strand separately,
		# If the result are the same, then only 
		# update the start, end, ref, alt to the
		# result. otherwise, the following '+', '-',
		# structure will be generated to reflect
		# the difference. they are all optional

		'+' => {

		  # This group assume the transcript is on forward strand,
		  # and give offsets and ref/alt string based on the rule
		  # with 'rep' annotation available
		  p  => $fpos,		    # forward strand pos, 0-based
		  r  => $fref,		    # forward strand ref string
		  a  => $falt,		    # forward strand alt string
		  rl => $freflen,	    # forward strand ref length
		  al => $faltlen,	    # forward strand alt length

		  # This group simplely trim off the leading same chars
		  # on forward strand, and then trim the same tail
		  bp  => $backward_fpos, 
		  br  => $backward_fref,
		  ba  => $backward_falt,
		  brl => $backward_freflen,
		  bal => $backward_faltlen,

		},

		'-' => {
		  same to '+', but for reverse strand
		},

	      # for repeat available variant
	      rep    => $repeat_element
	      replen => $repeat_element_length
	      ref_cn => $copy_number_in_ref,
	      alt_cn => $copy_number_in_alt,
	    }

=cut
sub parse_var {
    confess "Error: not enough args [",scalar(@_),"], need at least 4 args." if (4 > @_);
    my ($chr, $start, $end, $ref, $alt) = @_;
    my %var :shared;

    if ($end !~ /^\d+$/) { # from VCF v4.1 1-based start
	$alt = $ref;
	$ref = $end;
	$ref = normalise_seq($ref);
	$alt = normalise_seq($alt);
	my $rl = length($ref);
	if (substr($ref,0,1) eq substr($alt,0,1)) {
	    $ref = substr($ref,1);
	    $alt = substr($alt,1);
	}
	else {
	    $start -= 1; # change to 0-based start
	}
	$end = $start + $rl;
    }

    my $len_ref = $end - $start;            # chance to annotate long range
    my ($varType, $implicit_varType, $sm) = guess_type($len_ref, $ref, $alt);

    %var = (
        chr    => $chr,
        pos    => $start,
        ref    => $ref,
        alt    => $alt,
        reflen => $len_ref,
	guess  => $varType,
	imp    => $implicit_varType,
	sm     => $sm
    );
    $var{altlen} = length($alt) if ($alt ne '?');

    if ($guess eq 'no-call' or $implicit_varType ne 'delins') {
	return \%var;
    }
    
    return parse_complex( \%var );
}

sub normalise_seq {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq = "" if ($seq eq '.');
    $seq = uc($seq);
    if ($seq =~ /[^ACGTN]/ and $seq ne '?') {
	confess "Error: unrecognized pattern exists, no multiple alts please. [$seq]\n";
    }
    return $seq;
}

=head2 parse_complex
    
    About   : parse complex delins variants to recognize 
	      repeat and differ strand-pos var.
    Usage   : my $rguess = parse_complex( $var );
    Args    : delins variantion entry. have been uniform to CG's shell list format.
    Returns : see parse_var()

=cut
sub parse_complex {
    my $var = shift;
    my ($ref, $alt, $len_ref, $len_alt) = @{$var}{qw(ref alt reflen altlen)};

    my $get_rst = get_internal( $ref, $len_ref, $alt, $len_alt );
    if ($get_rst->{r} != $len_ref) {
	my ($imp_guess, $sm) = guess_type_by_length( $get_rst->{r}, $get_rst->{a} );
	$var->{imp} = $imp_guess;
	$var->{sm}  = $sm;
	if ( $get_rst->{'+'} != $get_rst->{'-'} ) {
	    for my $strd ( '+', '-' ) {
		$var->{$strd}->{bp} = $var->{pos} + $get_rst->{$strd};
		$var->{$strd}->{br} =
		  substr( $var->{ref}, $get_rst->{$strd}, $get_rst->{r} );
		$var->{$strd}->{ba} =
		  substr( $var->{alt}, $get_rst->{$strd}, $get_rst->{a} );
		$var->{$strd}->{brl} = $get_rst->{r};
		$var->{$strd}->{bal} = $get_rst->{a};
	    }
	}
	else {
	    $var->{pos} += $get_rst->{'+'};
	    $var->{reflen} = $get_rst->{r};
	    $var->{altlen} = $get_rst->{a};
	    $var->{ref} = substr( $var->{ref}, $get_rst->{'+'}, $get_rst->{r});
	    $var->{alt} = substr( $var->{alt}, $get_rst->{'+'}, $get_rst->{a});
	}
    }


    my %rst = ();
    @rst{ ("bcGuess", "bcRlen", "bcAlen") } = ($guess, $$get_rst{r}, $$get_rst{a});
    $rst{'+'}{bcOffst} = $$get_rst{'+'};
    $rst{'-'}{bcOffst} = $$get_rst{'-'};

    my $rc_ref = count_content($ref);
    my $rc_alt = count_content($alt);
    my @diff = map { $$rc_ref[$_] - $$rc_alt[$_] } (0 .. 5);

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
            $larger  = $ref;
            $llen    = $len_ref - 1;
            $smaller = $alt;
            $slen    = $len_alt - 1;
        }
        else {
            $larger  = $alt;
            $llen    = $len_alt - 1;
            $smaller = $ref;
            $slen    = $len_ref - 1;
        }

	my %has = ();
	for (my $rep = $llen; $rep > 0; $rep--) {
	    while ($larger =~ /([ACGTN]+)(?:\1){$rep}/g) {
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

# guess only by length for get_internal result
# the only delins without no-call left
# get_internal will recognize the real mutant
# region, discard the influnce by other 
# sample's result.
sub guess_type_by_length {
    my ($rl, $al) = @_;
    my ($imp_guess, $sm);
    if ($rl == 0) {
	$imp_guess = 'ins';
	$sm = 0;
    }
    elsif ($rl == 1) {
	if ($al == 1) {
	    $imp_guess = 'snv';
	}
	elsif ($al == 0) {
	    $imp_guess = 'del';
	}
	else {
	    $imp_guess = 'delins';
	}
	$sm = 1;
    }
    else {
	$imp_guess = ($al == 0) ? 'del' : 'delins';
	$sm = ($rl == $al) ? 3 : 2;
    }

    return ($imp_guess, $sm);
}

# input direct guess
sub guess_type {
    my ($len, $ref, $alt) = @_;
    # imp_varType: implicit variant type is for HGVS naming
    # varType    : to be output as varType name. can be as the key
    #		   to get the SO term by Name2SO hash.
    # sm	 : single or multiple bases tag
    #		    0 - for insertion, 0 base
    #		    1 - for single base variants
    #		    2 - for non-equal-length multiple bases variants
    #		    3 - for equal-length multiple bases delins
    my ($imp_varType, $varType, $sm);
    if ($len == 0) {
	$imp_varType = 'ins';
	$sm = 0;
    }
    elsif ($len == 1) {
	if ($ref eq $alt) {
	    $imp_varType = 'ref';
	}
	elsif ($alt eq '') {
	    $imp_varType = 'del';
	}
	elsif (1 == length($alt)) {
	    $imp_varType = 'snv';
	}
	else {
	    $imp_varType = 'delins';
	}
	$sm = 1;
    }
    elsif ($len > 1) {
	if ($ref eq $alt) {
	    $imp_varType = 'ref';
	}
	elsif ($alt eq '') {
	    $imp_varType = 'del';
	}
	else {
	    $imp_varType = 'delins';
	}

	if (length($ref) != length($alt)) {
	    $sm = 2; # non-equal-length subs
	}
	else {
	    $sm = 3; # equal-length delins
	}
    }
    $varType = ($alt eq '?') ? 'no-call' : $varType;
    return ($varType, $imp_varType, $sm);
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
    my ($new_ref_len, $new_alt_len);
    if ($shorter >= $loff + $roff) {
	$new_ref_len = $reflen - $loff - $roff;
	$new_alt_len = $altlen - $loff - $roff;
	return {
	    '+' => $loff,
	    '-' => $loff,
	    'r' => $new_ref_len,
	    'a' => $new_alt_len
	};
    }
    else {
	$new_ref_len = $reflen - $shorter;
	$new_alt_len = $altlen - $shorter;
	return {
	    '+' => $loff,
	    '-' => ($shorter - $roff),
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
    for (1 .. 5) {
        return 0 if ($$rc[$_] != $$rcs[$_] * $div);
    }
    return $div;
}

sub count_content {
    my $s = uc(shift);
    my $l = ($s =~ tr/ACGTN/12345/);
    my @count = ($l, 0, 0, 0, 0, 0);
    while ($s=~/(\d)/g) {
	$count[$1] ++;
    }
    return \@count;
}


=head2 fetchseq

    About   : get sequence from fasta db using samtools faidx
    Usage   : my $seq = fetchseq('db.fasta', $region_str);
	      my $rhash = fetchseq('db.fasta', \@regions);
    Args    : a list of region in format: (chr1:123-456, chr1:789-1000, chr2:234-567, or NM_01130.1:345 )
    Returns : a hash ref of { region => seq }

=cut
sub fetchseq {
    my ($fasta, $rRegs) = @_;
    if (!ref($rRegs)) {
	if ($rRegs =~ /[\S]\s+[\S]/) {
	    $rRegs = [split(/\s+/, $rRegs)];
	}
	else {
	    open (FASTA,"samtools faidx $fasta $rRegs |") or confess "samtools faidx: $!";
	    <FASTA>;
	    my @seq = <FASTA>;
	    my $seq = join("",@seq);
	    $seq =~ s/\s+//g;
	    close FASTA;
	    return $seq;
	}
    }
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

sub Code2Pep {
    my $code = shift;
    if (exists $Code2Pep{$code}) {
	return $Code2Pep{$code};
    }
    else {
	return '.';
    }
}

sub C3toC1 {
    my $c3 = shift;
    if (exists $C3toC1{$c3}) {
	return $C3toC1{$c3};
    }
    else {
	return '.';
    }
}

sub C1toC3 {
    my $c1 = shift;
    if (exists $C1toC3{$c1}) {
	return $C1toC3{$c1};
    }
    else {
	return '.';
    }
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

=head1 DATABASE FORMAT

    The Format

=head1 SEE ALSO

    HGVS     :  http://www.hgvs.org/mutnomen/recs.html
    Mutalyzer:  https://mutalyzer.nl

=head1 AUTHOR

liutao, E<lt>liutao@genomics.cnE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by liutao

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.14.2 or,
at your option, any later version of Perl 5 you may have available.


=cut
