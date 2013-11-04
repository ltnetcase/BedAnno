package BedAnno;

use strict;
use warnings;
use threads::shared; # if threads used before BedAnno, then this will be threaded
use Carp;
use Data::Dumper;

use Tabix;

our $VERSION = '0.35';

=head1 NAME

BedAnno - Perl module for annotating variation depend on the BED +1 format database.

=head2 VERSION v0.35

From version 0.32 BedAnno will change to support CG's variant shell list
and use ncbi annotation release 104 as the annotation database
with reformatted database format, and won't give any individual
annotation, so the individual_anno() is no longer available.
VCF4.1 format variant description (chr, pos, ref, alt) will also
be supported.

=head1 SYNOPSIS

  use BedAnno;
  my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz' );
  my $anno = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );

=head1 DESCRIPTION

This module need bgzipped BED+1 format database, combined with tabix index,
tabix is require for fast extract anno infomation randomly.
The fasta sequences are indexed by samtools faidx, sequence splicing also requires
samtools. If allele frequency information, prediction information or cytoBand 
information is needed, then extra dependencies will be required.

=head1 PARALLELED
    
This module can be used by multiple threads by:

    use threads;
    use threads::shared;
    use Thread::Queue;

to deal with large amount of variants.

=head1 Methods

=cut

our (
    %C3,      %C1,    %SO2Name, %func2SO,
    %Name2SO, %Polar, %C1toC3,  %AAnumber,
    $AAcount,
) : shared;

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
    
    UAA=>"X", UAG=>"X",

    UGA=>"U"
);

%C1toC3 = ();
foreach my $threebase (sort keys %C1) {
    $C1toC3{$C1{$threebase}} = $C3{$threebase};
}

$AAcount = scalar keys %C1toC3;
@AAnumber{sort keys %C1toC3} = (1 .. $AAcount);

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
    "SO:0000316" => 'CDS',
    "SO:0000204" => 'five_prime_UTR',
    "SO:0000205" => 'three_prime_UTR',
    "SO:0000655" => 'ncRNA',
    "SO:0000191" => 'interior_intron',
    "SO:0000448" => 'three_prime_UTR_intron',
    "SO:0000447" => 'five_prime_UTR_intron',
    "SO:0000163" => 'five_prime_cis_splice_site',
    "SO:0000164" => 'three_prime_cis_splice_site',
    "SO:0000167" => 'promoter',
    "SO:0000605" => 'intergenic_region',
    "span"       => 'span',

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
    "SO:0001627" => 'intron_variant',

    "annotation-fail"     => 'annotation-fail',
    "abnormal-intron"     => 'abnormal-intron',
    "abnormal-inseq-stop" => "abnormal-inseq-stop",
);

%Name2SO = reverse(%SO2Name);

%func2SO = (
    "abnormal-inseq-stop" => "abnormal-inseq-stop",
    "abnormal-intron"     => "abnormal-intron",
    "annotation-fail"     => 'annotation-fail',
    "cds-del"             => "SO:0001822",
    "cds-indel"           => "inframe_delins",
    "cds-ins"             => "SO:0001821",
    "cds-loss"		  => "unknown-likely-deleterious",
    "coding-synon"        => "SO:0001819",
    "init-loss"           => "SO:0001582",
    "no-change"           => "no-change",
    "splice"              => "unknown-likely-deleterious",
    "splice-3"            => "unknown-likely-deleterious",
    "splice-5"            => "unknown-likely-deleterious",
    "stop-gain"           => "SO:0001587",
    "stop-loss"           => "SO:0001578",
    "stop-retained"       => "SO:0001567",
    "unknown-no-call"     => 'unknown-no-call',
    "utr-3"               => "unknown",
    "utr-5"               => "unknown",
    altstart              => "SO:0001582",
    frameshift            => "SO:0001589",
    intron                => "unknown",
    knockout              => "SO:0001893",
    missense              => "SO:0001583",
    ncRNA                 => "unknown",
    nonsense              => "SO:0001587",
    promoter		  => "unknown",
    span		  => "unknown-likely-deleterious",
    unknown               => "unknown",
);

our %GenePartsOrder : shared;
@GenePartsOrder{
    (qw(CDS span five_prime_cis_splice_site
      three_prime_cis_splice_site ncRNA five_prime_UTR
      three_prime_UTR interior_intron five_prime_UTR_intron
      three_prime_UTR_intron abnormal-intron promoter
      annotation-fail intergenic_region), "")
} = ( 1 .. 15 );

our $REF_BUILD = 'GRCh37';

=head2 new

=over

=item About : Creat a new annotation entry

=item Usage :

    my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', batch => 1 );

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

    $self->throw("Error: please at least give 'db' and 'tr' path.")
	if (!exists $self->{db} or !exists $self->{tr});

    $self->set_db($self->{db});
    $self->set_tr($self->{tr});
    $self->set_refbuild($REF_BUILD);

    if ( exists $self->{cytoBand} ) {
	$self->set_cytoBand($self->{cytoBand});
    }

    if ( exists $self->{pfam} ) {
	$self->set_pfam($self->{pfam});
    }

    if ( exists $self->{prediction} ) {
	$self->set_prediction($self->{prediction});
    }

    if ( exists $self->{phyloP} ) {
	$self->set_phyloP($self->{phyloP});
    }

    if ( exists $self->{dbSNP} ) {
	$self->set_dbSNP($self->{dbSNP});
    }

    if (exists $self->{tgp}) {
	$self->set_tgp($self->{tgp});
    }

    if (exists $self->{esp6500}) {
	$self->set_esp6500($self->{esp6500});
    }

    if (exists $self->{cg54}) {
	$self->set_cg54($self->{cg54});
    }

    if (exists $self->{wellderly}) {
	$self->set_wellderly($self->{wellderly});
    }

    return $self;
}

=head2 set/get methods for properties

    List of Properties:
			get	    set
    db			o	    o
    tr			o	    o
    refbuild		o	    o
    tidb		o	    x
    annodb		o	    x
    trInfodb		o	    x
    cytoBand		o	    o
    cytoBand_h		o	    x
    pfam		o	    o
    pfam_h		o	    x
    prediction		o	    o
    prediction_h	o	    x
    phyloP		o	    o
    phyloP_h		o	    x
    dbSNP		o	    o
    dbSNP_h		o	    x
    tgp			o	    o
    tgp_h		o	    x
    esp6500		o	    o
    esp6500_h		o	    x
    cg54		o	    o
    cg54_h		o	    x
    wellderly		o	    o
    wellderly_h		o	    x


    e.g.    : $self->set_refbuild($refbuild);
	      my $refbuild = $self->get_refbuild();

=cut
sub set_refbuild {
    my $self = shift;
    my $custom_build_info = shift;
    $self->{refbuild} = $custom_build_info;
    return $self;
}

sub get_refbuild {
    my $self = shift;
    return $self->{refbuild};
}

sub set_db {
    my $self = shift;
    my $db = shift;
    if ( ! -e $db or ! -r $db ) {
	$self->throw("Error: cannot read $db.");
    }

    $self->{db} = $db;
    if ( ! -e $db.".tbi" or ! -r $db.".tbi" ) {
	$self->throw("Error: [$db] index (.tbi) file not found, please build it first.");
    }

    $self->{tidb} = shared_clone( Tabix->new( -data => $db ) );
    my %open_args;
    if ( exists $self->{region} ) {
        $open_args{region} = $self->{region};
    }
    if ( !exists $self->{region} and exists $self->{regbed} ) {
        $open_args{regbed} = $self->{regbed};
    }
    if ( exists $self->{genes} ) {
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
        $open_args{onlyPr} = $self->{onlyPr};
    }
    if ( exists $self->{mmap} ) {
        $open_args{mmap} = $self->{mmap};
    }

    $self->{annodb} = shared_clone( $self->load_anno(%open_args) )
      if ( exists $self->{batch} );
    return $self;
}

sub get_db {
    my $self = shift;
    return $self->{db};
}

sub get_tidb {
    my $self = shift;
    return $self->{tidb};
}

sub get_annodb {
    my $self = shift;
    if (exists $self->{annodb}) {
	return $self->{annodb};
    }
    else {
	return undef;
    }
}

sub set_tr {
    my $self = shift;
    my $tr = shift;
    if ( ! -e $tr or ! -r $tr ) {
	$self->throw("Error: cannot read $tr.");
    }
    $self->{tr} = $tr;
    $self->{trInfodb} = shared_clone( $self->readtr() );
    return $self;
}

sub get_tr {
    my $self = shift;
    return $self->{tr};
}

sub get_trInfodb {
    my $self = shift;
    return $self->{trInfodb};
}

sub set_cytoBand {
    my $self = shift;
    my $cytodb = shift;
    $self->{cytoBand} = $cytodb;
    require GetCytoBand if (!exists $self->{cytoBand_h});
    my $cytoBand_h = GetCytoBand->new( db => $cytodb );
    $self->{cytoBand_h} = shared_clone($cytoBand_h);
    return $self;
}

sub get_cytoBand {
    my $self = shift;
    return $self->{cytoBand} if (exists $self->{cytoBand});
    return undef;
}

sub get_cytoBand_h {
    my $self = shift;
    return $self->{cytoBand_h} if (exists $self->{cytoBand_h});
    return undef;
}

sub set_pfam {
    my $self = shift;
    my $pfamdb = shift;
    $self->{pfam} = $pfamdb;
    require GetPfam if (!exists $self->{pfam_h});
    my $pfam_h = GetPfam->new( db => $pfamdb );
    $self->{pfam_h} = shared_clone($pfam_h);
    return $self;
}

sub get_pfam {
    my $self = shift;
    return $self->{pfam} if (exists $self->{pfam});
    return undef;
}

sub get_pfam_h {
    my $self = shift;
    return $self->{pfam_h} if (exists $self->{pfam_h});
    return undef;
}

sub set_prediction {
    my $self = shift;
    my $predictiondb = shift;
    $self->{prediction} = $predictiondb;
    require GetPrediction if (!exists $self->{prediction_h});
    my $prediction_h = GetPrediction->new( db => $predictiondb );
    $self->{prediction_h} = shared_clone($prediction_h);
    return $self;
}

sub get_prediction {
    my $self = shift;
    return $self->{prediction} if (exists $self->{prediction});
    return undef;
}

sub get_prediction_h {
    my $self = shift;
    return $self->{prediction_h} if (exists $self->{prediction_h});
    return undef;
}

sub set_phyloP {
    my $self = shift;
    my $phyloPdb = shift;
    $self->{phyloP} = $phyloPdb;
    require GetPhyloP46wayScore if (!exists $self->{phyloP_h});
    my $phyloP_h = GetPhyloP46wayScore->new( db => $phyloPdb );
    $self->{phyloP_h} = shared_clone($phyloP_h);
    return $self;
}

sub get_phyloP {
    my $self = shift;
    return $self->{phyloP} if (exists $self->{phyloP});
    return undef;
}

sub get_phyloP_h {
    my $self = shift;
    return $self->{phyloP_h} if (exists $self->{phyloP_h});
    return undef;
}

sub set_dbSNP {
    my $self = shift;
    my $dbSNPdb = shift;
    $self->{dbSNP} = $dbSNPdb;
    require GetDBSNP if (!exists $self->{dbSNP_h});
    my $dbSNP_h = GetDBSNP->new( db => $dbSNPdb );
    $self->{dbSNP_h} = shared_clone($dbSNP_h);
    return $self;
}

sub get_dbSNP {
    my $self = shift;
    return $self->{dbSNP} if (exists $self->{dbSNP});
    return undef;
}

sub get_dbSNP_h {
    my $self = shift;
    return $self->{dbSNP_h} if (exists $self->{dbSNP_h});
    return undef;
}

sub set_tgp {
    my $self = shift;
    my $tgpdb = shift;
    $self->{tgp} = $tgpdb;
    require GetVcfAF if (!exists $self->{tgp_h} and !exists $self->{esp6500_h});
    my $tgp_h = GetVcfAF->new( db => $tgpdb );
    $self->{tgp_h} = shared_clone($tgp_h);
    return $self;
}

sub get_tgp {
    my $self = shift;
    return $self->{tgp} if (exists $self->{tgp});
    return undef;
}

sub get_tgp_h {
    my $self = shift;
    return $self->{tgp_h} if (exists $self->{tgp_h});
    return undef;
}

sub set_esp6500 {
    my $self = shift;
    my $esp6500db = shift;
    $self->{esp6500} = $esp6500db;
    require GetVcfAF if (!exists $self->{tgp_h} and !exists $self->{esp6500_h});
    my $esp6500_h = GetVcfAF->new( db => $esp6500db );
    $self->{esp6500_h} = shared_clone($esp6500_h);
    return $self;
}

sub get_esp6500 {
    my $self = shift;
    return $self->{esp6500} if (exists $self->{esp6500});
    return undef;
}

sub get_esp6500_h {
    my $self = shift;
    return $self->{esp6500_h} if (exists $self->{esp6500_h});
    return undef;
}

sub set_cg54 {
    my $self = shift;
    my $cg54db = shift;
    $self->{cg54} = $cg54db;
    require GetCGpub if (!exists $self->{cg54_h} and !exists $self->{wellderly_h});
    my $cg54_h = GetCGpub->new( db => $cg54db );
    $self->{cg54_h} = shared_clone($cg54_h);
    return $self;
}

sub get_cg54 {
    my $self = shift;
    return $self->{cg54} if (exists $self->{cg54});
    return undef;
}

sub get_cg54_h {
    my $self = shift;
    return $self->{cg54_h} if (exists $self->{cg54_h});
    return undef;
}

sub set_wellderly {
    my $self = shift;
    my $wellderlydb = shift;
    $self->{wellderly} = $wellderlydb;
    require GetCGpub if (!exists $self->{cg54_h} and !exists $self->{wellderly_h});
    my $wellderly_h = GetCGpub->new( db => $wellderlydb );
    $self->{wellderly_h} = shared_clone($wellderly_h);
    return $self;
}

sub get_wellderly {
    my $self = shift;
    return $self->{wellderly} if (exists $self->{wellderly});
    return undef;
}

sub get_wellderly_h {
    my $self = shift;
    return $self->{wellderly_h} if (exists $self->{wellderly_h});
    return undef;
}


sub TO_JSON {
    return { %{ shift() } };
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
                    
                    X        => 1,                   # inseqStop
                    U        => 1,                   # selenocysteine
                    A        => 1,                   # polyATail
                    altstart => {                    # altstart codons
                        $startCodons1 => 1,
                        $startCodons2 => 1,
                        ...
                    },
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
                            gid  => $gid,     (Entrez gene id)
                            gpSO => $gpSO,    (GeneParts SO)
                            blka => $blka,    (block attribute)
                            exin => $exin,    (exon intron number)
                            nsta => $nsta,    (n./r. left  of whole block)
                            nsto => $nsto,    (n./r. right of whole block)
                            csta => $csta,    (c.    left  of whole block)
                            csto => $csto,    (c.    right of whole block)
                            wlen => $wlen,    (length of whole block)
                            pr   => $pr,      (primary tag)
                            strd => $strd,    (strand)
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

    my @query_region = ();
    if ( exists $args{region} ) {
	my @regions = split(/\s+/, $args{region});
	foreach my $reg (@regions) {
	    next if ($reg eq "");
	    if ($reg =~ /^(\S+):(\d+)\-(\d+)$/) {
		my ($name, $beg, $end) = ($1, $2, $3);
		if ($beg == 0) {
		    $self->warn("Warning: region string should be 1 based [$reg], has been changed to 1 based");
		    $beg = 1;
		}
                push( @query_region, [ $name, ( $beg - 1 ), $end ] );
	    }
	    else {
		$self->throw("Error: unavailable region string [$reg].");
	    }
	}
    }
    elsif ( exists $args{regbed}
	and defined $args{regbed}
	and -e $args{regbed} )
    {
	open (BED, $args{regbed}) or $self->throw("Error: [$args{regbed}] $!");
	while (<BED>) {
	    chomp;
	    my @beditm = split(/\t/);
	    $self->throw("Error: bed format error.") if (3 > @beditm);
	    push (@query_region, [ @beditm[0, 1, 2] ]);
	}
	close BED;
    }

    my $read_all_opt = 0;
    my @all_querys   = ();
    if ( 0 == @query_region ) {
	$read_all_opt = 1;
        open( ANNO, "zcat -f $self->{db} |" )
          or $self->throw("Error: [$self->{db}] $!");
    }
    else {
        my @sorted_regions = sort {
            $a->[0] cmp $b->[0] or $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2]
        } @query_region;

	my ($cname, $cbeg, $cend) = @{$sorted_regions[0]};
	for (my $k = 1; $k < @sorted_regions; $k++) {
	    my ($dname, $dbeg, $dend) = @{$sorted_regions[$k]};
	    if ($dname ne $cname or $dbeg > $cend) {
                push( @all_querys,
                    $self->{tidb}->query( $cname, $cbeg, $cend ) );
		($cname, $cbeg, $cend) = @{$sorted_regions[$k]};
	    }
	    else {
		$cend = $dend;
	    }
	}
	my $q = $self->{tidb}->query( $cname, $cbeg, $cend );
	push( @all_querys, $q ) if (defined $q->{_});
    }

    # trans filter is always be ahead of genes
    my $prTag   = ( exists $args{onlyPr} ) ? 1 : 0;
    my $mmapTag = ( exists $args{mmap} )   ? 1 : 0;
    my $geneTag = ( exists $args{genes} )  ? 1 : 0;
    my $tranTag = ( exists $args{trans} )  ? 1 : 0;
    my $geneList = $args{genes} if ($geneTag);
    my $tranList = $args{trans} if ($tranTag);
    my %pureTran = map { s/\-\d+$//; $_ => 1 } keys %$tranList if ($tranTag);

    my $rannodb = {};
    while (1)
    {
	my $tb_ent;
	if ($read_all_opt) {
	    $tb_ent = <ANNO>;
	    last if (!defined $tb_ent or $tb_ent eq "");
	}
	else {
	    while (0 < @all_querys) {
		my $region_read = $self->{tidb}->read($all_querys[0]);
		if (!defined $region_read or $region_read eq "") {
		    shift (@all_querys);
		}
		else {
		    $tb_ent = $region_read;
		    last;
		}
	    }
	    last if (0 == @all_querys);
	}

        $tb_ent =~ s/\s+$//;
        my ( $chr, $start, $stop, $annostr ) = split(/\t/, $tb_ent);
	if ( !defined $annostr or $annostr eq "" ) {
	    $self->throw("Error: db format unmatched");
	}

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
    close ANNO if ($read_all_opt);

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

# involke parse_annoent to assign detail information to annodb.
sub assign_detail {
    my $self = shift;
    my $rannodb_k = shift;
    my %detail = ();
    foreach my $annoblk (sort keys %{$$rannodb_k{annos}}) {
        my $offset = $$rannodb_k{annos}{$annoblk};
        my ($tid, $ranno) = parse_annoent($annoblk);
        $detail{$tid} = $ranno;
        $detail{$tid}{offset} = $offset;
    }
    $rannodb_k->{detail} =
      ( is_shared($rannodb_k) ) ? shared_clone( {%detail} ) : {%detail};
    return $rannodb_k;
}

# just parse annoents when mutation hits.
sub parse_annoent {
    my $annoent  = shift;
    my %annoinfo = ();

# tid         gsym  gid strd blka gpgo exin nsta nsto csta csto wlen mismatch    pr
# NM_015658.3|NOC2L|26155|-|3U1E|205|EX19E|2817|2801|*508|*492|0|D,879582,879582,.|Y
    my @infos = split( /\|/, $annoent );
    confess "Error format of anno ents [$annoent]" if ( 14 != @infos );

    my $tid = shift @infos;
    my @tags =
      qw(gsym gid strd blka gpSO exin nsta nsto csta csto wlen mismatch pr);
    @annoinfo{@tags} = @infos;
    return ( $tid, \%annoinfo );
}


=head2 anno

    About   : Annotate single short variation by annotation db.
    Usage   : my $anno_ent = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );
              or $anno_ent = $beda->anno( 'chr20', 1234568, 'AG', 'AGGG' );
    Args    : for CG's shell variants, need 5 args in UCSC coordinates
              (0-based start), they are:
                chr id, chr start, chr end, reference, alternative.
              for variants in VCF, need 4 args, which is lack of 
                chr end, and "chr start" is in 1-based coordinates.
              for crawler: a input object with keys: 
		chr,begin,referenceSequence,variantSequence,[end].
                if end is specified, then use 0-based coordinates,
                otherwise 1-based (VCF) coordinates will be used.
    Returns : a hash ref of annotation informations, see varanno().

=cut
sub anno {
    my $self = shift;
    my $var = BedAnno::Var->new(@_);
    my ($annoEnt, $idx) = $self->varanno($var);
    return $annoEnt;
}


=head2 varanno

    About   : implicitly create a new BedAnno::Anno entry, and 
              assign all the needed annotation for to it.
    Usage   : ($rAnnoRst, $AEIndex) = $beda->varanno($var, $AEIndex);
    Args    : The BedAnno entry, BedAnno::Var entry and current dbidx, 
              current dbidx should be always 0, if used for non-batch mode.
    Returns : A BedAnno::Anno entry and current dbidx for nex query in batch.
            {
                var => {

                    # the first part is from var parsing result.
                    # please see "BedAnno::Var new()".
                    # import all the keys from original parsed var entry
                    # and add the following keys by this method.

                    # information
                    varTypeSO => $varTypeSO,
                    gHGVS     => $gHGVS,
                    refbuild  => $referenceBuild,

                    # Here's some optional parts which may be generated
                    # when extra resource is available:

                    cytoBand  => $cytoBand,

                    # For single position for now
                    phyloPpm    => $PhyloPscorePlacentalMammals,
                    phyloPpr    => $PhyloPscorePrimates,
                    phyloPve    => $PhyloPscoreVetebrates,

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
                        cdsBegin      => $Begin_in_CDS,            # cDot format
                        cdsEnd        => $End_in_CDS,              # cDot format
                        protBegin     => $Begin_in_Protein,
                        protEnd       => $End_in_Protein,
                        c             => $cHGVS,
                        p             => $pHGVS,
                        cc            => $codon_change,
                        polar         => $polar_change,
                        r             => $imp_funcRegion,
                        r_Begin       => $imp_beginfuncRegion,
                        r_End         => $imp_endfuncRegion,
                        func          => $imp_funcCode,
                        exin          => $exIntr_number,
                        ei_Begin      => $imp_Begin_exIntr_number,
                        ei_End        => $imp_End_exIntr_number,
                        genepart      => $GenePart,
                        genepartSO    => $GenePartSO,
                        genepartIndex => $GenePartIndex,
                        exonIndex     => $exonIndex,               # '.' for N/A
                        intronIndex   => $intronIndex,             # '.' for N/A
                        funcSOname    => $FunctionImpact,
                        funcSO        => $FunctionImpactSO,
                        trAlt         => $alt_string_on_transcript,
                        trRef         => $ref_string_on_transcript,
                        prAlt         => $protein_alt_sequence,
                        prRef         => $protein_ref_sequence,
			primaryTag    => $refstandard_primary_or_not,	# Y/N
                        preStart => {    # the position before the start of var
                            nDot => $rna_hgvs_pos,
                            cDot => $cds_hgvs_pos,
                            r    => $func_region,
                            exin => $exon_intron_number,
                        },
                        postEnd => {     # the position after the end of var
                            nDot => $rna_hgvs_pos,
                            cDot => $cds_hgvs_pos,
                            r    => $func_region,
                            exin => $exon_intron_number,
                        },
                        trRefComp => {

                            # some trRef components
                        },

                        # The following parts will be exists if extra resource
                        # is available.
                        pfamId      => $PFAM_ID,
                        pfamName    => $PFAM_NAME,
                        siftPred    => $SIFTpred,
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
    my ($self, $var, $AEIndex) = @_;
    $AEIndex ||= 0;

    if (exists $self->{cytoBand}) {
	$var->{cytoBand} = $self->{cytoBand_h}->getCB(@$var{qw(chr pos end)});
    }

    if (exists $self->{phyloP}) {
	if ($var->{sm} == 1) {
	    @$var{qw(phyloPpm phyloPpr phyloPve)} = 
		@{$self->{phyloP_h}->getPhyloP46wayScore($var->{chr}, ($var->{pos}+1))};
	}
    }

    if ( exists $self->{dbSNP} ) {
        if ( exists $var->{sep_snvs} ) {
            my @new_sqls  = ();
            my $cur_start = $var->{sep_snvs}->[0];
            for my $i ( 0 .. $#{ $var->{sep_snvs} } ) {
                if (   $i == $#{ $var->{sep_snvs} }
                    or $var->{sep_snvs}->[ $i + 1 ] - $cur_start > 1 )
                {
                    my $new_ref = substr(
                        $var->{ref},
                        ( $cur_start - $var->{pos} - 1 ),
                        ( $var->{sep_snvs}->[$i] - $cur_start + 1 )
                    );
                    my $new_alt = substr(
                        $var->{alt},
                        ( $cur_start - $var->{pos} - 1 ),
                        ( $var->{sep_snvs}->[$i] - $cur_start + 1 )
                    );
                    push(
                        @new_sqls,
                        [
                            $cur_start, $var->{sep_snvs}->[$i],
                            $new_ref,   $new_alt
                        ]
                    );
                    $cur_start = $var->{sep_snvs}->[ $i + 1 ]
                      if ( $i < $#{ $var->{sep_snvs} } );
                }
            }
            foreach my $rSE (@new_sqls) {
                my $rOneSql = $self->{dbSNP_h}->getRS( $var->{chr}, @$rSE );
                @{ $$var{dbsnp} }{ sort keys %$rOneSql } =
                  @$rOneSql{ sort keys %$rOneSql };
            }
        }
        else {
            $var->{dbsnp} =
              $self->{dbSNP_h}->getRS( @$var{qw(chr pos end ref alt)} );
        }
    }

    if (exists $self->{tgp}) {
	$var->{tgp} = $self->{tgp_h}->getAF(@$var{qw(chr pos end ref alt)});
    }

    if (exists $self->{esp6500}) {
	$var->{esp6500} = $self->{esp6500_h}->getAF(@$var{qw(chr pos end ref alt)});
    }

    if (exists $self->{cg54}) {
	$var->{cg54} = $self->{cg54_h}->getAF(@$var{qw(chr pos end ref alt)});
    }

    if (exists $self->{wellderly}) {
	$var->{wellderly} = $self->{wellderly_h}->getAF(@$var{qw(chr pos end ref alt)});
    }

    $var->{varTypeSO} = $Name2SO{$var->{guess}};
    $var->{refbuild}  = $self->{refbuild};
    $var->get_gHGVS();

    # Due to the bed format database, 
    # add flank left and right 1bp to query the database
    # to involve all the mismatches
    my $var_region =
      (     $var->{chr} . ':'
          . ( $var->{pos} - 1 ) . '-'
          . ( $var->{pos} + $var->{reflen} + 1 ) );

    my $localdb;
    if (!exists $self->{batch}) { # daemon mode
	my $locLoad = $self->load_anno( region => $var_region );
	$localdb = (exists $locLoad->{$var->{chr}}) ? $locLoad->{$var->{chr}} : [];
    }
    else {
	if ( !exists $self->{annodb}->{$var->{chr}} ) {
	    $self->warn("Warning: incomplete annotation database in batch mode, [$var->{chr}]");
	    $localdb = [];
	}
	else {
	    $localdb = $self->{annodb}->{$var->{chr}};
	}
    }

    my $annoEnt = BedAnno::Anno->new($var);
    $AEIndex = $annoEnt->getTrPosition($localdb, $AEIndex);
    $annoEnt = $self->getTrChange($annoEnt);
    $annoEnt = $self->finaliseAnno($annoEnt);

    return ($annoEnt, $AEIndex);
}

=head2 finaliseAnno

    About   : finalise the BedAnno::Anno entry by check all tag values,
              and uniform them for AE output usage, query transcript
              oringinated additional resources and add them into the data
              frame.
    Usage   : $beda->finaliseAnno($annEnt);
    Args    : BedAnno entry and a BedAnno::Anno entry
    Returns : A finalised BedAnno::Anno entry

=cut
sub finaliseAnno {
    my ($self, $annoEnt) = @_;
    
    if (exists $annoEnt->{trInfo}) {
	foreach my $tid (sort keys %{$annoEnt->{trInfo}}) {
	    my $trAnnoEnt = $annoEnt->{trInfo}->{$tid};
	    my $trdbEnt = $self->{trInfodb}->{$tid};
	    $trAnnoEnt->{prot} =
	      ( !exists $trdbEnt->{prot} or $trdbEnt->{prot} eq '.' )
	      ? ""
	      : $trdbEnt->{prot};

	    my ($strd) = ($trAnnoEnt->{strd} eq '+') ? 1 : 0;

	    # complete informations depend on already assigned values
	    my $genepartSO;
	    if ($trAnnoEnt->{rnaBegin} =~ /^\+/ or $trAnnoEnt->{rnaEnd} =~ /^\+/) {
		$genepartSO = 'span'; # 3'downstream span
	    }
	    elsif ( $trAnnoEnt->{r_Begin} eq $trAnnoEnt->{r_End} ) {
		$genepartSO =
		  ( $trAnnoEnt->{genepartSO} =~ /^\d+$/ )
		  ? sprintf( "SO:%07d", $trAnnoEnt->{genepartSO} )
		  : $trAnnoEnt->{genepartSO};

		$trAnnoEnt->{r} = $trAnnoEnt->{r_Begin};
		$trAnnoEnt->{exin} = $trAnnoEnt->{ei_Begin};

		if ($trAnnoEnt->{exin} =~ /^EX(\d+)/) {
		    $trAnnoEnt->{genepartIndex} = $1;
		    $trAnnoEnt->{exonIndex} = $1;
		}
		elsif ($trAnnoEnt->{exin} =~ /^IVS(\d+)/) {
		    $trAnnoEnt->{genepartIndex} = $1;
		    $trAnnoEnt->{intronIndex} = $1;
		}
		else {
		    $trAnnoEnt->{genepartIndex} = 0;
		}
	    }
	    elsif ( $trAnnoEnt->{ei_Begin} =~ /^IVS/ 
		and $trAnnoEnt->{ei_Begin} eq $trAnnoEnt->{ei_End} )
	    { # span splice or interior intron
		$genepartSO = 'span'; # need to discuss? whether to use interior
		$trAnnoEnt->{r} =
		  ($strd)
		  ? $trAnnoEnt->{r_Begin} . '-' . $trAnnoEnt->{r_End}
		  : $trAnnoEnt->{r_End} . '-' . $trAnnoEnt->{r_Begin};
		$trAnnoEnt->{exin}          = $trAnnoEnt->{ei_Begin};
		$trAnnoEnt->{genepartIndex} = $1 if ($trAnnoEnt->{ei_Begin} =~ /(\d+)/);
		$trAnnoEnt->{intronIndex}   = $trAnnoEnt->{genepartIndex};
	    }
	    else { # non-equal exon intron number or UTR/CDS span
		$genepartSO = 'span';
		$trAnnoEnt->{r} =
		  ($strd)
		  ? $trAnnoEnt->{r_Begin} . '-' . $trAnnoEnt->{r_End}
		  : $trAnnoEnt->{r_End} . '-' . $trAnnoEnt->{r_Begin};
		if ($trAnnoEnt->{ei_Begin} eq $trAnnoEnt->{ei_End}) {
		    $trAnnoEnt->{exin} = $trAnnoEnt->{ei_Begin};
		    $trAnnoEnt->{genepartIndex} = $1
		      if ( $trAnnoEnt->{ei_Begin} =~ /EX(\d+)/ );
		    $trAnnoEnt->{exonIndex} = $trAnnoEnt->{genepartIndex};
		}
		else { # for same index exon intron give this index
		    $trAnnoEnt->{exin} =
		      ($strd)
		      ? $trAnnoEnt->{ei_Begin} . '-' . $trAnnoEnt->{ei_End}
		      : $trAnnoEnt->{ei_End} . '-' . $trAnnoEnt->{ei_Begin};

		    if ( $trAnnoEnt->{exin} =~ /(\d+)\D+(\d+)/ and $1 eq $2 ) {
			$trAnnoEnt->{genepartIndex} = $1;
		    }
		}
	    }

	    confess "Error: unknown genepartSO [$genepartSO]."
	      if ( !exists $SO2Name{$genepartSO} );

	    $trAnnoEnt->{genepart} = $SO2Name{$genepartSO};
	    $trAnnoEnt->{genepartSO} = ( $genepartSO =~ /^SO:/ ) ? $genepartSO : "";
	    $trAnnoEnt->{genepartIndex} = '.'
	      if ( !exists $trAnnoEnt->{genepartIndex} );
	    $trAnnoEnt->{exonIndex} = '.' if ( !exists $trAnnoEnt->{exonIndex} );
	    $trAnnoEnt->{intronIndex} = '.'
	      if ( !exists $trAnnoEnt->{intronIndex} );

	    # change annotation-fail 's other value to empty?

	    # uniform rnaBegin end and protBegin end
	    $trAnnoEnt->{rnaBegin} =
	      ( $trAnnoEnt->{rnaBegin} =~ /^\+/ ) ? "" : $trAnnoEnt->{rnaBegin};
	    $trAnnoEnt->{rnaEnd} =
	      ( $trAnnoEnt->{rnaEnd} =~ /^\+/ ) ? "" : $trAnnoEnt->{rnaEnd};
	    $trAnnoEnt->{protBegin} =
	      ( !exists $trAnnoEnt->{protBegin} or $trAnnoEnt->{protBegin} eq '0' ) ? "" : $trAnnoEnt->{protBegin};
	    $trAnnoEnt->{protEnd} =
	      ( !exists $trAnnoEnt->{protEnd} or $trAnnoEnt->{protEnd} eq '0' ) ? "" : $trAnnoEnt->{protEnd};

	    $trAnnoEnt->{func} = 'unknown' if (!exists $trAnnoEnt->{func});
	    # uniform function group
	    confess "Error: unknown func code [$trAnnoEnt->{func}]."
	      if ( !exists $func2SO{ $trAnnoEnt->{func} } );
	    $trAnnoEnt->{funcSO}     = $func2SO{ $trAnnoEnt->{func} };
	    $trAnnoEnt->{funcSOname} = $SO2Name{ $trAnnoEnt->{funcSO} };
	    $trAnnoEnt->{funcSO} =
	      ( $trAnnoEnt->{funcSO} =~ /^SO:/ ) ? $trAnnoEnt->{funcSO} : "";

	    
	    # add additional resource
	    if (exists $trAnnoEnt->{prot} and $trAnnoEnt->{prot} ne "") {
		$trAnnoEnt->{protBegin} =
		  ( $trAnnoEnt->{protBegin} eq '0' ) ? "" : $trAnnoEnt->{protBegin};
		$trAnnoEnt->{protEnd} =
		  ( $trAnnoEnt->{protEnd} eq '0' ) ? "" : $trAnnoEnt->{protEnd};
		if ($trAnnoEnt->{protBegin} ne "" and $trAnnoEnt->{protEnd} ne "") {
		    if (exists $self->{pfam}) {
			my ($pb, $pe) = ($trAnnoEnt->{protBegin}, $trAnnoEnt->{protEnd});
			if ($pb > $pe) {
			    my $tmp = $pb;
			    $pb = $pe;
			    $pe = $tmp;
			}
			( $trAnnoEnt->{pfamId}, $trAnnoEnt->{pfamName} ) =
			  $self->{pfam_h}->getPfam( $trAnnoEnt->{prot}, $pb, $pe );
		    }
		}
		if ( exists $self->{prediction} ) {
		    if ( exists $trAnnoEnt->{p}
			and $trAnnoEnt->{p} =~ /^p\.[A-Z](\d+)([A-Z])$/ )
		    {
			my $rpred =
			  $self->{prediction_h}
			  ->getPredScore( $trAnnoEnt->{prot}, $1, $2 );
			if ( exists $rpred->{sift} ) {
			    $trAnnoEnt->{siftPred}  = $rpred->{sift}->[0];
			    $trAnnoEnt->{siftScore} = $rpred->{sift}->[1];
			}
			if ( exists $rpred->{polyphen2_humdiv} ) {
			    $trAnnoEnt->{pp2divPred} =
			      $rpred->{polyphen2_humdiv}->[0];
			    $trAnnoEnt->{pp2divScore} =
			      $rpred->{polyphen2_humdiv}->[1];
			}
			if ( exists $rpred->{polyphen2_humvar} ) {
			    $trAnnoEnt->{pp2varPred} =
			      $rpred->{polyphen2_humvar}->[0];
			    $trAnnoEnt->{pp2varScore} =
			      $rpred->{polyphen2_humvar}->[1];
			}

		    }
		    elsif ( exists $trAnnoEnt->{prRef}
			and exists $trAnnoEnt->{prAlt} )
		    {
			my $prReflen = length( $trAnnoEnt->{prRef} );
			my $prAltlen = length( $trAnnoEnt->{prAlt} );
			if (    $prReflen == 1
			    and $prAltlen == 1
			    and $trAnnoEnt->{protBegin} eq '1'
			    and $trAnnoEnt->{protEnd} eq '1' )
			{
			    my $init_pred =
			      $self->{prediction_h}
			      ->getPredScore( $trAnnoEnt->{prot}, 1,
				$trAnnoEnt->{prAlt} );

			    if ( exists $init_pred->{sift} ) {
				$trAnnoEnt->{siftPred} =
				  $init_pred->{sift}->[0];
				$trAnnoEnt->{siftScore} =
				  $init_pred->{sift}->[1];
			    }

			    if ( exists $init_pred->{polyphen2_humdiv} ) {
				$trAnnoEnt->{pp2divPred} =
				  $init_pred->{polyphen2_humdiv}->[0];
				$trAnnoEnt->{pp2divScore} =
				  $init_pred->{polyphen2_humdiv}->[1];
			    }

			    if ( exists $init_pred->{polyphen2_humvar} ) {
				$trAnnoEnt->{pp2varPred} =
				  $init_pred->{polyphen2_humvar}->[0];
				$trAnnoEnt->{pp2varScore} =
				  $init_pred->{polyphen2_humvar}->[1];
			    }

			}
		    }
		}
	    }
	}
    }

    return $annoEnt;
}

=head2 getTrChange

    About   : Calculate the transcript changes, based on TrPostition
    Usage   : $beda->getTrChange($annoEnt);
    Returns : assign the following tags in annoEnt
                trRef, prot, c, p, cc, polar, func
                prRef, prAlt

=cut
sub getTrChange {
    my ($self, $annoEnt) = @_;

    # do nothing if not hit any transcript
    if ( !exists $annoEnt->{trInfo}
        or 0 == ( scalar keys %{ $annoEnt->{trInfo} } ) )
    {
        return $annoEnt;
    }

    my $trInfodb = $self->{trInfodb};
    my @tids = sort keys %{$annoEnt->{trInfo}};

    my %trSeqs = ();
    if (exists $self->{batch}) {
	@trSeqs{@tids} = map {$_->{seq}} @$trInfodb{@tids};
    }
    else {
	# due to a bug in samtools faidx, when meet "|" 
	# in a fasta's header, and the id is the last id of input
	# then it will fail to generate correct result.
	# here we can only try to use single transcript,
	# It may become a little slow for annotation
	#my $rtrs = fetchseq($self->{tr}, \@tids);
	#%trSeqs = %$rtrs;
	
	foreach my $t (@tids) {
	    $trSeqs{$t} = fetchseq($self->{tr}, $t);
	}
    }

#    print STDERR Data::Dumper->Dump( [\@tids, \%trSeqs], ["tids", "trSeqs"] );

    foreach my $tid (@tids) {
	my $trdbEnt = $self->{trInfodb}->{$tid};
	my $trannoEnt = $annoEnt->{trInfo}->{$tid};

	my ( $unify_p, $unify_r, $unify_a, $unify_rl, $unify_al ) =
	  $annoEnt->{var}->getUnifiedVar( $trannoEnt->{strd} );

	# 1. check if annotation-fail
        if (   $trannoEnt->{genepartSO} eq 'annotation-fail'
            or $trannoEnt->{rnaBegin} eq '?'
            or $trannoEnt->{rnaEnd} eq '?' )
        {
            $trannoEnt->{func} = 'annotation-fail';
            next;
        }

	my $cdsOpt = (exists $trdbEnt->{prot} and $trdbEnt->{prot} ne ".") ? 1 : 0;
	my $strd   = ($trannoEnt->{strd} eq '+') ? 1 : 0;

	# fix the trRefComp when ended in 3'downstream

	# debug
#	print STDERR Data::Dumper->Dump( [$tid, $trannoEnt, $unify_r, $trSeqs{$tid}, $strd], ["tid", "trannoEnt", "unify_r", "trSeq", "strd"] );

	if (!exists $trSeqs{$tid} or $trSeqs{$tid} eq "") {
	    $self->throw("Error: your fasta database file may not complete. [$tid]");
	}
	my $trRef = getTrRef($trannoEnt, $unify_r, $trSeqs{$tid}, $strd);
	$trannoEnt->{trRef} = $trRef;
	my $trAlt = $trannoEnt->{trAlt};

	$trannoEnt->{prot} = $trdbEnt->{prot} if ($cdsOpt);
	my $f = ($cdsOpt) ? 'c.' : 'n.';
	my $chgvs_5 = ($cdsOpt) ? $trannoEnt->{cdsBegin} : $trannoEnt->{rnaBegin};
	my $chgvs_3 = ($cdsOpt) ? $trannoEnt->{cdsEnd}   : $trannoEnt->{rnaEnd};

	# [ aaPos, codon, aa, polar, frame ]
        my @cInfo5 =
          getCodonPos( $trdbEnt, $trSeqs{$tid}, $trannoEnt->{rnaBegin} );
        my @cInfo3 =
          getCodonPos( $trdbEnt, $trSeqs{$tid}, $trannoEnt->{rnaEnd} );

	# we just use position comparison to show transcript ref
	# stat, instead of sm in var, because there may exists
	# indel difference bewteen refgenome and refSeq
	my $cmpPos = $self->cmpPos( $chgvs_5, $chgvs_3 );


	# 2. check if no-call
	if ( $trAlt eq '?' ) {
	    $trannoEnt->{func} = 'unknown-no-call';
	    if ($cmpPos == 0) { # 1 bp
		$trannoEnt->{c} = $f . $chgvs_5 . $trRef . '>?';
		if ($cInfo5[0] > 0) {
                    my $aaOut = $cInfo5[2];
                    $trannoEnt->{p} =
                      'p.' . $aaOut . $cInfo5[0] . '?';
		    $trannoEnt->{cc} = $cInfo5[1].'=>?';
		    $trannoEnt->{polar} = $cInfo5[3].'=>?';
		    $trannoEnt->{protBegin} = $cInfo5[0];
		    $trannoEnt->{protEnd} = $cInfo5[0];
		}
	    }
	    elsif ($cmpPos == 1) { # multiple bp
		$trannoEnt->{c} = $f . $chgvs_5 . '_' . $chgvs_3 . 'del'
		    . $trRef . 'ins?';
		if ($cInfo5[0] > 0 and $cInfo3[0] > 0) {
		    if ($cInfo5[0] == $cInfo3[0]) {
			my $aaOut = $cInfo5[2];
			$trannoEnt->{p} =
			  'p.' . $aaOut . $cInfo5[0] . '?';
			$trannoEnt->{cc} = $cInfo5[1].'=>?';
			$trannoEnt->{polar} = $cInfo5[3].'=>?';
			$trannoEnt->{protBegin} = $cInfo5[0];
			$trannoEnt->{protEnd} = $cInfo5[0];
		    }
		    else {
			my $aa5 = $cInfo5[2];
			my $aa3 = $cInfo3[2];

                        $trannoEnt->{p} = 'p.'
                          . $aa5
                          . $cInfo5[0] . '_'
                          . $aa3
                          . $cInfo3[0]
                          . 'delins?';
			
			$trannoEnt->{protBegin} = $cInfo5[0];
			$trannoEnt->{protEnd} = $cInfo3[0];
		    }
		}
	    }
	    else { # ins : reverse the positions
		$trannoEnt->{c} = $f . $chgvs_3 . '_' . $chgvs_5 . 'ins?';
		if ($cInfo5[0] > 0 and $cInfo3[0] > 0) {
		    if ($cInfo5[0] == $cInfo3[0]) {
			my $aaOut = $cInfo5[2];
			$trannoEnt->{p} =
			  'p.' . $aaOut . $cInfo5[0] . 'delins?';
			$trannoEnt->{protBegin} = $cInfo5[0];
			$trannoEnt->{protEnd} = $cInfo5[0];
		    }
		    else {
			my $aa5 = $cInfo5[2];
			my $aa3 = $cInfo3[2];

                        $trannoEnt->{p} = 'p.'
                          . $aa3
                          . $cInfo3[0] . '_'
                          . $aa5
                          . $cInfo5[0]
                          . 'ins?';
			
			$trannoEnt->{protBegin} = $cInfo5[0];
			$trannoEnt->{protEnd} = $cInfo3[0];
		    }
		}
	    }

	    next;
	}

	# reparse the transcript originated var
	my $real_var = BedAnno::Var->new( $tid, 0, length($trRef), $trRef, $trAlt );
	if ($real_var->{guess} eq 'ref') {
	    $trannoEnt->{func} = 'no-change';
	    $trannoEnt->{c}    = $f . '=';
	    next;
	}
	
	my ($real_p, $real_r, $real_a, $real_rl, $real_al) = 
	    $real_var->getUnifiedVar('+');

	my ($trBegin, $trEnd);
	$trBegin = reCalTrPos_by_ofst( $trannoEnt, $real_p );
        $trEnd = reCalTrPos_by_ofst( $trannoEnt, ( $real_p + $real_rl - 1 ) );

        $chgvs_5 =
          ($cdsOpt)
          ? cPosMark( $trBegin, $trdbEnt->{csta}, $trdbEnt->{csto},
            $trdbEnt->{len} )
          : $trBegin;
        $chgvs_3 =
          ($cdsOpt)
          ? cPosMark( $trEnd, $trdbEnt->{csta}, $trdbEnt->{csto},
            $trdbEnt->{len} )
          : $trEnd;
        $cmpPos = $self->cmpPos( $chgvs_5, $chgvs_3 );

	# debug
#        print STDERR Data::Dumper->Dump(
#            [ $real_var, $trBegin,  $trEnd,  $chgvs_5,  $chgvs_3,  $cmpPos ],
#            [ "real_var", "trBegin", "trEnd", "chgvs_5", "chgvs_3", "cmpPos" ]
#        );

	# * check if a repeat case, and use cerntain chgvs string
	# * assign cHGVS string
	if ($real_var->{imp} eq 'rep') {
            my $trRep = $real_var->{rep};

            if (    $real_var->{ref_cn} == 1
                and $real_var->{alt_cn} == 2 )
            {
                if ( $real_var->{replen} == 1 ) {
                    $trannoEnt->{c} = $f . $chgvs_5 . 'dup' . $trRep;
                }
                else {
                    $trannoEnt->{c} =
                      $f . $chgvs_5 . '_' . $chgvs_3 . 'dup' . $trRep;
                }
            }
	    else {
		$trannoEnt->{c} =
		    $f
		  . $chgvs_5
		  . $trRep . '['
		  . $real_var->{ref_cn} . '>'
		  . $real_var->{alt_cn} . ']';
	    }
	}
	else {
	    if ( $cmpPos == 0 ) {    # 1 bp

		$trannoEnt->{c} = $f . $chgvs_5;
		if ($real_al == 1) {
		    $trannoEnt->{c} .= $real_r . '>' . $real_a; # 1bp alt
		}
		elsif ($real_al == 0) {
		    $trannoEnt->{c} .= 'del' . $real_r;
		}
		else {
		    $trannoEnt->{c} .= 'delins' . $real_a;
		}

	    }
	    elsif ( $cmpPos > 0 ) {    # multiple bp

		$trannoEnt->{c} =
		  $f . $chgvs_5 . '_' . $chgvs_3 . 'del';
		if ($real_al > 0) {
		    $trannoEnt->{c} .= 'ins' . $real_a;
		}
	    }
	    else {                      # ins : reverse the positions
		$trannoEnt->{c} =
		  $f . $chgvs_3 . '_' . $chgvs_5 . 'ins' . $real_a;
	    }
	}

        # 3. check if transcript-ablation
        if (( $trBegin =~ /^\-/ or $trBegin eq '1' ) 
		and ( $trEnd =~ /^\+/ or $trEnd eq "$trdbEnt->{len}" ) )
        {
            $trannoEnt->{func} = 'knockout';
	    # need prot begin end?

            next;
        }
	
	# * check if exon remapping introduced special case
	if ( $chgvs_5 =~ /^\+|\+d/ and $chgvs_3 =~ /^\+|\+d/ ) {
	    $trannoEnt->{func} = 'unknown';
	    next;
	}

	# 4. check if span different exon/intron/promoter/downstream
        if (
            (
                   $chgvs_3 =~ /^\+|\+d/
                or $trannoEnt->{ei_Begin} ne $trannoEnt->{ei_End}
            )
            and $cmpPos >= 0	    # not insertion at the edge
          )
        {
            $trannoEnt->{func} = 'span';

	    # need prot begin end?

            next;
        }
	
	# 5. check if all promoter
        if (    $trannoEnt->{r_Begin} eq 'PROM'
            and $trannoEnt->{r_Begin} eq $trannoEnt->{r_End} )
        {
            $trannoEnt->{func} = 'promoter';
	    
            next;
        }

	# 6. check if all intron
        if (    $trannoEnt->{ei_Begin} =~ /IVS/
            and $trannoEnt->{ei_Begin} eq $trannoEnt->{ei_End} )
        {
            if ( $trannoEnt->{genepartSO} eq 'abnormal-intron' ) {
                $trannoEnt->{func} = 'abnormal-intron';
            }
            elsif ( $trannoEnt->{r_Begin} =~ /^D/
                and $trannoEnt->{r_End} =~ /^A/ )
            {
                $trannoEnt->{func} = 'splice';
            }
            elsif ( $trannoEnt->{r_Begin} =~ /^D/ ) {
                $trannoEnt->{func} = 'splice-5';
            }
            elsif ( $trannoEnt->{r_End} =~ /^A/ ) {
                $trannoEnt->{func} = 'splice-3';
            }
            else {
                $trannoEnt->{func} = 'intron';
            }

            next;
        }


	# 7. check if all exon or ex/intron edge insertion
	#    count all this insertion to be exon region insertion
        if (
            (
                    $trannoEnt->{ei_Begin} eq $trannoEnt->{ei_End}
                and $trannoEnt->{ei_Begin} =~ /EX/
            )
            or
            ( $cmpPos < 0 and $trannoEnt->{ei_Begin} ne $trannoEnt->{ei_End} )
          )
	{ # any edge-insertion case will count on the exon change
	  # any coding-utr edge insertion case will count on the
	  # utr parts
	    if (!$cdsOpt) {
		if ($chgvs_5 =~ /^\+/) {
		    $trannoEnt->{func} = 'unknown';
		}
		elsif ($chgvs_3 =~ /^\-/) {
		    $trannoEnt->{func} = 'promoter';
		}
		elsif ($chgvs_5 eq '1' and $chgvs_3 eq $trdbEnt->{len}) {
		    $trannoEnt->{func} = 'knockout';
		}
		else {
		    $trannoEnt->{func} = 'ncRNA';
		}
		next;
	    }
	    else { # coding RNA
		if ($chgvs_3 =~ /\-u1$/) {
		    $trannoEnt->{func} = 'promoter';
		    next;
		}
		elsif ($chgvs_5 =~ /\+d1$/) { # 3' downstream edge insertion
		    $trannoEnt->{func} = 'unknown';
		}
		elsif ($chgvs_3 =~ /^\-/) { # besides utr-5 - cds edge
		    $trannoEnt->{func} = 'utr-5';
		    next;
		}
		elsif ($chgvs_5 =~ /^\*/) { # besides utr-3 - cds edge
		    $trannoEnt->{func} = 'utr-3';
		    next;
		}
		elsif ($chgvs_5 =~ /^\-/ and $chgvs_3 =~ /^\*/) {
                    my $utr5_len = $1 if ( $chgvs_5 =~ /^\-(\d+)$/ );
                    my $utr3_len = $1 if ( $chgvs_3 =~ /^\*(\d+)$/ );
                    if (    $utr5_len eq $trdbEnt->{csta}
                        and $utr3_len eq ( $trdbEnt->{len} - $trdbEnt->{csto} )
                      )
                    {
                        $trannoEnt->{func} = 'knockout';
                    }
                    else {
                        $trannoEnt->{func} = 'cds-loss';
                    }
                    $trannoEnt->{p} = 'p.0?';
                    next;
		}
		elsif ($chgvs_5 =~ /^\-/) {
		    $trannoEnt->{func} = 'init-loss';
		    $trannoEnt->{p}    = 'p.0?';
		    next;
		}
		elsif ($chgvs_3 eq '1-1') { # cds-begin - last 5'utr intron edge
		    $trannoEnt->{func} = 'utr-5';
		    next;
		}
                elsif ( $chgvs_5 =~ /^(\d+)\+1/
                    and $1 == ( $trannoEnt->{csto} - $trannoEnt->{csta} ) )
                {
                    $trannoEnt->{func} = 'utr-3';
                    next;
                }
		# here protein sequence var should be reparsed
		else {

		    $chgvs_5 = $1+1 if ($chgvs_5 =~ /^(\d+)\+1/);
		    $chgvs_3 = $1-1 if ($chgvs_3 =~ /^(\d+)\-1/);

		    # 5' end should be in coding region.
		    # leave this as debug information
                    $self->throw(
                        "Error: unavailable cDot HGVS. [$chgvs_5, $chgvs_3]")
                      if ( $chgvs_5 !~ /^\d+$/
                        or $chgvs_3 !~ /^\d+$/
                        or $chgvs_5 == 0
                        or $chgvs_3 == 0 );

		    @cInfo5 =
		      getCodon_by_cdsPos( $trdbEnt, $trSeqs{$tid}, $chgvs_5 );
		    @cInfo3 =
		      getCodon_by_cdsPos( $trdbEnt, $trSeqs{$tid}, $chgvs_3 );

		    my %translate_opts = ();
		    $translate_opts{mito} = 1 if ($tid =~ /^NM_MT-/);
		    $translate_opts{polyA} = 1 if (exists $trdbEnt->{A});
		    $translate_opts{nostop} = 1 if (exists $trdbEnt->{X});
		    my %altcodon_opts = %translate_opts;
                    delete $altcodon_opts{nostop}
                      if ( exists $altcodon_opts{nostop} );

                    if ( !exists $trdbEnt->{pseq} ) {
                        my $whole_cds = substr( $trSeqs{$tid}, $trdbEnt->{csta},
                            ( $trdbEnt->{csto} - $trdbEnt->{csta} ) );
                        my ( $pseq, $frame_next ) = translate(
                            $whole_cds, \%translate_opts
                        );
                        $trdbEnt->{pseq} = $pseq;
                    }

		    my ($prBegin, $prEnd, $prRef, $prAlt);

		    $prBegin = $cInfo5[0];

		    my $ready_to_add_5 = substr($cInfo5[1], 0, $cInfo5[4]);

		    my $diff_ra = $real_rl - $real_al;

		    # probably frame shift flag
		    my $frameshift_flag = ($diff_ra % 3 > 0) ? 1 : 0;

		    # end in cds or not.
                    my $end_in_cds_flag = (
                             $chgvs_3 =~ /^\*/
                          or $chgvs_3 > ( $trdbEnt->{csto} - 3 )
                    ) ? 0 : 1;

                    my $ready_to_add_3;
                    if ( !$end_in_cds_flag or $frameshift_flag ) {
                        $ready_to_add_3 = substr( $trSeqs{$tid}, $trEnd );
			# this variants's effect will be end at the terminal
			$prEnd = $trdbEnt->{plen} + 1; # terminal
                    }
                    else {    # chgvs_3 in cds region with no frameshift
                        $ready_to_add_3 = substr(
                            $cInfo3[1],
                            ( $cInfo3[4] + 1 ),
                            ( 2 - $cInfo3[4] )
                        );
			$prEnd = $cInfo3[0];
                    }

		    # Make protBegin and protEnd be coordinate to the 
		    # coordinates of nucl variants, extend the 
		    # frameshift's effect to the end of protein
		    $trannoEnt->{protBegin} = $prBegin;
		    $trannoEnt->{protEnd} = $prEnd;

		    # make the order of func assignment as follow:
		    # * total
		    # 1. init-loss
		    # 2. frameshift
		    # 3. coding-synon
		    # 
		    # * 1bp aa change
		    # 0. abnormal-inseq-stop
		    # 1. stop-retained
		    # 2. stop-loss
		    # 3. nonsense
		    # 4. missense
		    # 
		    # * other variants
		    # 0. abnormal-inseq-stop
		    # 1. stop-retained
		    # 2. frameshift
		    # 3. stop-loss
		    # 4. stop-gain
		    # 5. cds-ins
		    # 6. cds-del
		    # 7. cds-indel
		    
		    # extend the altered sequence to codon coordinates
		    my $codon_alt = $ready_to_add_5 . $real_a . $ready_to_add_3;

		    # init codon flag
                    my $hit_init_flag = ( $trBegin > $trdbEnt->{csta}
                          and $trBegin <= ( $trdbEnt->{csta} + 3 ) ) ? 1 : 0;
		    # stop codon flag
                    my $hit_stop_flag =
                    ( $trEnd > ( $trdbEnt->{csta} + 3 * $trdbEnt->{plen} ) )
                      ? 1
                      : 0;

		    # check if an alterstart
		    my %start_codons = ('ATG' => -1);
                    if ( exists $trdbEnt->{altstart} ) {
                        %start_codons =
                          map { $_ => -1 } keys %{ $trdbEnt->{altstart} };
                    }

		    my $init_synon = 0; # indicate start codon search result
		    my $init_frame_shift = 0; # indicate if init frameshift
		    my $current_startCodon = substr($trSeqs{$tid},$trdbEnt->{csta},3);
                    if ($hit_init_flag) {
			# check if altered sequence have new start codon
                        foreach my $startCodon ( sort keys %start_codons ) {
                            $start_codons{$startCodon} =
                              index( $codon_alt, $startCodon );
                            if ( $start_codons{$startCodon} > -1 ) {
                                $init_synon = 1;
                                $current_startCodon = $startCodon;

				# trim the leading non-coding region
                                $codon_alt  = substr( $codon_alt,
                                    $start_codons{$startCodon} );

				# init frameshift
                                if ( $start_codons{$startCodon} % 3 > 0 ) {
                                    $init_frame_shift = 1;
                                    my $original_var_frame = $diff_ra % 3;
                                    my $reparsed_frame =
                                      $start_codons{$startCodon} % 3;
                                    if ($original_var_frame == $reparsed_frame )
                                    { # correct to non-frameshift
                                        $frameshift_flag = 0;
                                    }
                                    elsif ( $original_var_frame == 0 ) {
					# change to a frameshift case
                                        $codon_alt .= substr(
                                            $trSeqs{$tid},
                                            (
                                                $trdbEnt->{csta} +
                                                  ( 3 * $cInfo3[0] )
                                            )
                                        );
                                        $frameshift_flag = 1;
					$prEnd = $trdbEnt->{plen} + 1;
					$trannoEnt->{protEnd} = $prEnd;
                                    }
                                }
                                last;
                            }
                        }

			if ($init_synon == 0) {
			    $trannoEnt->{func} = 'init-loss';
			    $trannoEnt->{p}    = 'p.0?';
			    if ($prBegin eq '1' and $prEnd eq '1') {
                                my $zero;
                                $trannoEnt->{prRef} = 'M';
                                my $altCodon = substr( $codon_alt, 0, 3 );
                                ( $trannoEnt->{prAlt}, $zero ) =
                                  translate( $altCodon, \%altcodon_opts );
			    }
			    next;
			}
			else {
			    $trannoEnt->{func} = 'altstart';
			    # can it continue to be assigned
			    # other function code and pHGVS
			    # or just finished as altstart?
			    #
			    # I choose to continue.
			}
                    }

		    # give out the first of single codon to be
		    my $codon_to_be = substr($codon_alt, 0, 3);

		    # debug
#		    print STDERR Data::Dumper->Dump( [$trdbEnt, $prBegin, $prEnd, $codon_alt], ["trdbEnt", "prBegin", "prEnd", "codon_alt"] );

                    $prRef = substr(
                        $trdbEnt->{pseq},
                        ( $prBegin - 1 ),
                        ( $prEnd - $prBegin + 1 )
                    );

		    # we can not allow inseq stop codon in altered sequence
		    # due to frameshift and ambiguity.
		    my $next_alt_frame;
                    ( $prAlt, $next_alt_frame ) =
                      translate( $codon_alt, \%altcodon_opts );


		    if ($hit_init_flag) {
			$prRef =~ s/^[A-Z]/M/; # change init pep to M
			# change alt-init to M
			$prAlt =~ s/^[A-Z]/M/ if ($init_synon);
		    }

		    $trannoEnt->{prRef} = $prRef;
		    $trannoEnt->{prAlt} = $prAlt;

		    # to indicate whether the alternate sequence 
		    # encode a stop codon or non-frameshift, or otherwise
		    # with a non stopped frameshift.
                    my $non_stop_flag = ( $prAlt !~ /\*$/
                          and ( $next_alt_frame or $frameshift_flag ) ) ? 1 : 0;

		    # To avoid large range perlre, before parse protein 
		    # var we should first deal with frameshift issues
		    # need to deal with non_stop_flag and frameshift and 
		    # imp ref case.
		    
		    # assign cc and polar to the same length, 
		    # single aa substitution, no matter which position.
		    if ($diff_ra == 0 and $cInfo5[0] == $cInfo3[0]) {
			# using trim due to terminal codon
			my ($aa_to_be, $polar_to_be) = getAAandPolar($codon_to_be);
			$trannoEnt->{cc} = $cInfo5[1].'=>'.$codon_to_be;
			$trannoEnt->{polar} = $cInfo5[3].'=>'.$polar_to_be;
		    }

		    my $prSimpleSame = 0;
                    for (
                        my $pp = 0 ;
                        $pp < length($prRef) and $pp < length($prAlt) ;
                        $pp++
                      )
                    {
                        if (
                            substr( $prRef, $pp, 1 ) eq
                            substr( $prAlt, $pp, 1 ) )
                        {
                            $prSimpleSame++;
                        }
                        else {
                            last;
                        }
                    }
		    my $no_parsed_pP = $prBegin - 1 + $prSimpleSame;
                    my $no_parsed_prStart =
                      substr( $trdbEnt->{pseq}, $no_parsed_pP, 1 );

		    # non stop frameshift with extra non-coded base
		    if ($non_stop_flag) {
                        $trannoEnt->{func} =
                          ( $no_parsed_pP >= $trdbEnt->{plen} )
                          ? 'stop-loss'
                          : 'frameshift';
                        $trannoEnt->{p} = 'p.' . $no_parsed_prStart
                          . ( $no_parsed_pP + 1 ) . 'fs*?';
			next;
		    }

		    # identical to reference (can be from any kind of vars)
		    if ($prRef eq $prAlt) {
                        if ($init_synon)
                        { # altstart
                            # do nothing
                        }
			elsif ($hit_stop_flag) {

			    # need this special assertion?
			    $trannoEnt->{func} = 'stop-retained';
			}
			else {
			    $trannoEnt->{func} = 'coding-synon';
			}
			$trannoEnt->{p}    = 'p.=';
			next;
		    }

		    # frameshift
		    if (!$end_in_cds_flag or $frameshift_flag) {

			$trannoEnt->{func} = 'frameshift';
			$trannoEnt->{p}    = 'p.'. $no_parsed_prStart 
			    . ( $no_parsed_pP + 1 ) . 'fs*';
			if ($prAlt =~ /\*$/) { # ext length estimated
                            $trannoEnt->{p} .=
                              ( length($prAlt) - $prSimpleSame );
			}
			else { # don't meet a stop codon
			    $trannoEnt->{p} .= '?';
			}
			next;
		    }

		    # parse the protein variants
		    # to recognize the repeat and adjust to correct position
                    my $prVar =
                      BedAnno::Var->new( $trdbEnt->{prot}, ( $prBegin - 1 ),
                        $prEnd, $prRef, $prAlt );
		    # 0-based
                    my ( $p_P, $p_r, $p_a, $prl, $pal ) =
                      $prVar->getUnifiedVar('+');

		    my $prStart = substr( $trdbEnt->{pseq}, $p_P, 1 );
                    my $prStop = substr( $trdbEnt->{pseq}, ($p_P + $prl - 1), 1 )
                      if ( $p_P > 0 or $prl > 0 );

		    # single substitution
                    if ( $prVar->{imp} eq 'snv' ) {
                        if ( $p_r eq 'X' ) {
			    $trannoEnt->{func} = 'abnormal-inseq-stop';
                        }
			elsif ( $p_r eq '.' ) { # N in refseq transcript
			    $trannoEnt->{func} = 'unknown';
			}
			elsif ( $p_r eq '*' ) {
			    $trannoEnt->{func} = 'stop-loss';
			}
                        elsif ( $p_a eq '*' ) {
			    $trannoEnt->{func} = 'nonsense';
                        }
			elsif ($init_synon) {
			    $trannoEnt->{func} = 'altstart';
			}
                        else {
                            $trannoEnt->{func} = 'missense';
                        }
                        $trannoEnt->{p} = 'p.' . $p_r . ($p_P + 1) . $p_a;
                        next;
                    }

		    # repeat
		    if ( $prVar->{imp} eq 'rep' ) {
			if ( $prVar->{ref_cn} < $prVar->{alt_cn} ) {
			    $trannoEnt->{func} = 'cds-ins';
			}
			else {
			    $trannoEnt->{func} = 'cds-del';
			}

			if ($prVar->{replen} == 1) {
			    $trannoEnt->{p} = 'p.' . $prVar->{rep} . ( $p_P + 1 );
			}
			else {
                            $trannoEnt->{p} = 'p.'
                              . $prStart
                              . ( $p_P + 1 ) . '_'
                              . substr( $prVar->{rep}, -1, 1 )
                              . ( $p_P + $prVar->{replen} );
			}
                        if (    $prVar->{ref_cn} == 1
                            and $prVar->{alt_cn} == 2 )
                        {
                            $trannoEnt->{p} .= 'dup';
                        }
                        else {
                            $trannoEnt->{p} .= '['
                              . $prVar->{ref_cn} . '>'
                              . $prVar->{alt_cn} . ']';
                        }

			next;
		    }

		    # insertion
		    if ( $prVar->{sm} == 0 ) {
			# insert into the 5'edge of start codon
			# this may caused by a insertion or delins
			# around init-codon to recreat a new init-codon
			# in the altered sequence
			if ($p_P == 0) {
			    # altstart don't fix this
			    # use pHGVS to indicate 
			    # altstart here
			    $trannoEnt->{p} = 'p.=';

			    # can it be utr-5?
			    $trannoEnt->{func} = 'utr-5';
			    next;
			}
			# insert between the last aa and stop codon
			# this may caused by a insertion or delins
			# aroud the stop codon
			elsif ($hit_stop_flag) {
			    $trannoEnt->{func} = 'stop-retained';
			}
			else {
			    $trannoEnt->{func} = 'cds-ins';
			}

                        $trannoEnt->{p} = 'p.'
                          . $prStop
                          . $p_P . '_'
                          . $prStart
                          . ( $p_P + 1 ) . 'ins'
                          . $p_a;
			next;
		    }

                    if ( $prVar->{sm} >= 1 ) {
                        if ( $p_r =~ /\*/ ) {
                            $trannoEnt->{func} = 'stop-loss';
                        }
                        elsif ($hit_stop_flag) {
                            $trannoEnt->{func} = 'stop-retained';
                        }
                        elsif ( $pal == 0 ) {
                            $trannoEnt->{func} = 'cds-del';
                        }
                        elsif ( $p_a =~ /\*/ ) {
                            $trannoEnt->{func} = 'stop-gain';
                        }
                        elsif ( $pal == $prl ) {
                            $trannoEnt->{func} = 'missense';
                        }
                        else {
                            $trannoEnt->{func} = 'cds-indel';
                        }

                        $trannoEnt->{p} = 'p.' . $prStart . ( $p_P + 1 );
                        if ( $prVar->{sm} == 1 ) {
                            $trannoEnt->{p} .= 'del';
                        }
                        else {
                            $trannoEnt->{p} .=
                              '_' . $prStop . ( $p_P + $prl ) . 'del';
                        }
                        $trannoEnt->{p} .= 'ins' . $p_a if ( $p_a ne "" );
                        next;
                    }
		}
	    }
        }

	# any other cases?
    }

    return $annoEnt;
}

=head2 cmpPos

    About   : judge the order or p1 and p2, because the insertion
              will have a reverted order of left & right position
    Usage   : my $cmpRst = BedAnno->cmpPos($p1, $p2);
    Args    : hgvs positio p1 and p2, with out 'c.' or 'n.' flag
    Return  : 0 for same, 1 for normal order, -1 for reversed order.

=cut
sub cmpPos {
    my $self = shift;
    my ($p1, $p2) = @_;
    return 0 if ($p1 eq $p2);
    my ($s1, $s2, $anc1, $anc2, $int_s1, $int_s2, $ofst1, $ofst2);

    my %order_s = (
	'-u' => 1,
	'-' => 2,
	''  => 3,
	'+' => 4,
	'*' => 5,
	'+d' => 6,
    );

    if ($p1 =~ /^([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$/) {
	$s1 = $1; $anc1 = $2; $int_s1 = $3; $ofst1 = $4;
    }
    if ($p2 =~ /^([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$/) {
	$s2 = $1; $anc2 = $2; $int_s2 = $3; $ofst2 = $4;
    }

    if ($order_s{$s1} < $order_s{$s2}) {
	return 1;
    }
    elsif ($order_s{$s1} > $order_s{$s2}) {
	return -1;
    }
    else { # $s1 eq $s2
	if ($s1 eq '-') {
	    if ($anc1 > $anc2) {
		return 1;
	    }
	    elsif ($anc1 < $anc2) {
		return -1;
	    }
	}
	else {
	    if ($anc1 < $anc2) {
		return 1;
	    }
	    elsif ($anc1 > $anc2) {
		return -1;
	    }
	}

	# anc1 eq anc2
	if ($order_s{$int_s1} < $order_s{$int_s2}) {
	    return 1;
	}
	elsif ($order_s{$int_s1} > $order_s{$int_s2}) {
	    return -1;
	}
	else { # int_s2 same
	    if ( $int_s1 =~ /\-/ ) {
		if ($ofst1 > $ofst2) {
		    return 1;
		}
		elsif ($ofst1 < $ofst2) {
		    return -1;
		}
	    }
	    else {
		if ($ofst1 < $ofst2) {
		    return 1;
		}
		elsif ($ofst1 > $ofst2) {
		    return -1;
		}
	    }
	}
    }
    return 0;
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
    local $/ = ">";
    if (!ref($rRegs)) { # a scalar string
	if ($rRegs =~ /[\S]\s+[\S]/) { # multiple regions
	    $rRegs = [split(/\s+/, $rRegs)];
	}
	else { # single region
	    open (FASTA,"samtools faidx $fasta $rRegs |") or confess "samtools faidx: $!";
	    <FASTA>; # shift first ">";
	    my $seq = <FASTA>;
	    close FASTA;
	    $seq =~ s/^[^\n]+//;
	    $seq =~ s/[\s>]+//g;
	    return uc($seq);
	}
    }
    my @localRegs = @$rRegs;
    my %seqs = ();
    while (0 < @localRegs) {
	my @set = ();
	if (1000 < @localRegs) { 
	    @set = splice (@localRegs, 0, 1000);
	}
	else {
	    @set = @localRegs;
	    @localRegs = ();
	}
	my $regions = join(" ", @set);
        open( FASTA, "samtools faidx $fasta $regions |" )
          or confess "samtools faidx: $!";
	my @all_seq = <FASTA>;
	close FASTA;
	shift @all_seq;
	foreach my $rec (@all_seq) {
	    chomp $rec;
	    my $hd = $1 if ($rec =~ s/^(\S+)[^\n]*//);
	    $rec =~ s/[\s>]+//g;
	    next if ($rec =~ /^$/);
	    if (!defined $hd) {
		confess "Error: exists empty header in fasta [$rec]";
	    }
	    $seqs{$hd} = uc($rec);
	}
    }
    return \%seqs;
}


=head2 getCodon_by_cdsPos

    About   : get codon position, codon string, aa string, and frame info
              for one single cds position
    Usage   : my ($pP, $codon, $aa, $polar, $frame) = getCodon_by_cdsPos($trdbEnt, $trSeq, $p);
    Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto
              for cds start/end position
              trSeq is transcript sequence
              p is cds position
    Returns : AA position  - 0 for not in cds region
              codon string - codon string, e.g. "ATA"
              aa char      - AA code, 1 bp mode, e.g. "E"
              polar        - Polar properties, e.g. "P+"
              frame        - current position's frame info,
                             -1 for not available.

=cut
sub getCodon_by_cdsPos {
    my ($trdbEnt, $trSeq, $p) = @_;
    if (   !exists $trdbEnt->{csta}
        or $trdbEnt->{csta} eq '.'
	or $p !~ /^\d+$/
        or $p <= 0
        or $p > ($trdbEnt->{csto} - $trdbEnt->{csta}) )
    {
        return ( 0, "", "", "", -1 );
    }

    $trSeq = substr( $trSeq, $trdbEnt->{csta},
        ( $trdbEnt->{csto} - $trdbEnt->{csta} ) );
    my $cds_p = $p; # > 0
    my $pP = int($cds_p / 3);
    if ($cds_p % 3 > 0) {
	$pP ++;
    }
    my $frame = 2 - ($pP * 3 - $cds_p);
    my $codon = uc(substr($trSeq, ($p - $frame - 1), 3));
    my $rtrans_opts = {};

    $rtrans_opts->{mito} = 1 if ($trdbEnt->{gene} =~ /^MT\-/); # chrM gene

    if ($pP <= $trdbEnt->{plen}) { # plen not involve terminal codon
	$rtrans_opts->{nostop} = 1;
    }
    elsif ($pP == $trdbEnt->{plen} + 1)
    { # terminal with polyA complement
	$codon .= 'A' x (3 - length($codon)) if (exists $trdbEnt->{A});
    }
    else {
	confess "Error: may be not correct protein length in transcript database.";
    }

    my ($aa, $polar);
    if ($pP == 1) {
	$aa = 'M';
	$polar = 'NP';
    }
    else {
	($aa, $polar) = getAAandPolar($codon, $rtrans_opts);
    }

    return ($pP, $codon, $aa, $polar, $frame);
}


=head2 getCodonPos

    About   : get codon position, codon string, aa string, and frame info
              for one single cds position
              See getCodon_by_cdsPos()

=cut
sub getCodonPos {
    my ($trdbEnt, $trSeq, $p) = @_;
    if (   !exists $trdbEnt->{csta}
        or $trdbEnt->{csta} eq '.'
	or $p !~ /^\d+$/
        or $p <= $trdbEnt->{csta}
        or $p > $trdbEnt->{csto} )
    {
        return ( 0, "", "", "", -1 );
    }

    my $cds_p = $p - $trdbEnt->{csta}; # > 0
    return getCodon_by_cdsPos($trdbEnt, $trSeq, $cds_p);
}

sub getAAandPolar {
    my $codon = shift;
    my $rtrans_opts = shift;
    my ($aa, $zero) = translate($codon, $rtrans_opts);
    my $polar = (exists $Polar{$aa}) ? $Polar{$aa} : '.';
    return ($aa, $polar);
}

# change the 1based nDot position into cDot format
sub cPosMark {
    #	1based	0-based 1based
    my ($trPos, $csta, $csto, $trlen) = @_;
    return $trPos if ($trPos eq "");

    my $cDot;
    if ($trPos =~ /^([\+\-]?)(\d+)([\+\-]?)(\d*)$/) {
        my ( $s1, $anchor, $s2, $ofst ) = ( $1, $2, $3, $4 );
        if ( $s1 eq '' and $s2 eq '' ) {
            $cDot = $anchor - $csta;
            $cDot -= 1 if ( $cDot <= 0 );
            if ( $anchor > $csto ) {
                $cDot = '*' . ( $anchor - $csto );
            }
        }
        elsif ( $s1 ne '' and $s2 ne '' ) {
            confess "Error: not a valid nDot string: [$trPos]";
        }
        elsif ( $s1 eq '-' ) {
            $cDot = ( ( $csta == 0 ) ? 1 : ( -$csta ) ) . '-u' . $anchor;
        }
        elsif ( $s1 eq '+' ) {
            $cDot = ( ( $csto == $trlen )
                ? ( $csto - $csta )
                : '*' . ( $trlen - $csto ) )
              . '+d'
              . $anchor;
        }
	elsif ( $s2 =~ /[\+\-]/ ) {
	    my $cDot_tmp = cPosMark($anchor, $csta, $csto, $trlen);
	    $cDot = $cDot_tmp.$s2.$ofst;
	}
	else {
	    confess "Error: unmatched nDot string: [$trPos]";
	}
    }
    return $cDot;
}

=head2 translate

    About   : Translate nucleotides to aa seqs
    Usage   : my ($aa_seq, $next_frame) = translate( $nucl, { mito => 1, polyA => 1 } );
    Args    : first args should be the nucleotide seq to be translated,
              the second args is a hash ref of optional args (all boolean):
              mito   : indicate it is for mDNA (mitochondrion)
              nostop : indicate there's no stop codon in aa seq,
                       translate 'UGA' to 'U', and other stop codon to 'X'
              polyA  : indicate extra A should be added to 3'end of codon,
                       to help encode a stop codon (usually used with mito).
    Returns : $aa_seq is the aa sequence, 
              $next_frame gives the next base's frame to the 3'end of sequence.

=cut
sub translate {
    my $nucl = shift;
    my $ropt = shift if ( @_ );
    my $mito = (defined $ropt and exists $ropt->{mito}) ? 1 : 0;
    my $nostop = (defined $ropt and exists $ropt->{nostop}) ? 1 : 0;
    my $polyA = (defined $ropt and exists $ropt->{polyA}) ? 1 : 0;

    my $lenN = length($nucl);
    my $frame_next = $lenN % 3;
    if ($frame_next != 0 and $polyA) {
	$nucl .= 'A' x (3 - $frame_next);
    }
    $nucl = uc($nucl);
    my $prot = "";
    for (my $i = 0; $i < length($nucl); $i+=3) {
	my $codon = substr($nucl,$i,3);
	my $aa = (exists $C1{$codon}) ? $C1{$codon} : '.';
	if ( $mito and $codon eq 'ATA' ) {
	    $aa = 'M';
	}
	if ($aa eq '*' and $mito and $codon eq 'TGA') {
	    $aa = 'W';
	}
	elsif ($aa eq '*' and $nostop) {
	    $aa = ($codon eq 'TGA') ? 'U' : 'X';
	}
	elsif ($aa eq '*') {
	    $prot .= $aa;
	    return ($prot, 0);
	}
	$prot .= $aa;
    }
    $prot =~ s/[UX]$/*/g if ($nostop);
    $prot =~ s/\.$// if ($frame_next != 0);
    return ($prot, $frame_next);
}

# new var generated from new() in BedAnnoVar the real transcript var
# will bring some new transcript start/end information
# here will recalculate it from the offset.
sub reCalTrPos_by_ofst {
    my ( $trannoEnt, $trRef_ofst ) = @_;
    return $trannoEnt->{rnaBegin} if ( $trRef_ofst == 0 );
    return $trannoEnt->{preStart}->{nDot}
      if ( $trRef_ofst == -1 and exists $trannoEnt->{preStart} );

    my @tag_sort = sort trRefSort keys %{ $trannoEnt->{trRefComp} };
    my $cumulate_len = 0;
    my ($cur_ex_start, $intOfst);
    if ($trannoEnt->{rnaBegin} =~ /^(\-?\d+)([\+\-]?\d*)$/) {
	$cur_ex_start = $1;
	$intOfst = $2;
    }
    if ($cur_ex_start =~ /^-/) {
	$intOfst = $cur_ex_start;
	$cur_ex_start = 1;
    }
    if ($intOfst eq '') {
	$intOfst = 0;
    }

    foreach my $exin (@tag_sort) {
	my $cur_blk_len;
	if ( $exin !~ /^EX/ ) {
            $cur_blk_len =
              ( $trannoEnt->{trRefComp}->{$exin}->[1] -
                  $trannoEnt->{trRefComp}->{$exin}->[0] );
	}
	else {
	    $cur_blk_len = $trannoEnt->{trRefComp}->{$exin};
	}

	my $cur_blk_ofst = ($trRef_ofst - $cumulate_len);

	if ( $cur_blk_ofst >= $cur_blk_len ) {
	    if ( $cur_blk_ofst == $cur_blk_len and $exin eq $tag_sort[-1] ) {
		return $trannoEnt->{postEnd}->{nDot};
	    }
	    else {
		$cumulate_len += $cur_blk_len;
		if ($exin =~ /^EX/) { # only cumulate ex pos
		    $cur_ex_start += $cur_blk_len;
		    $intOfst = 0;
		}
	    }
	}
	else { # hit current block
	    if ( $exin !~ /^EX/) {
		if ($intOfst < 0 and $cur_ex_start == 1) { # in promoter
		    confess "Error: non-zero cumulate_len in promoter [$cumulate_len]." if ($cumulate_len > 0);
		    return ($intOfst + $cur_blk_ofst);
		}
		elsif ($exin =~ /^Z/) {
		    return '+'.($cur_blk_ofst + 1);
		}
		else { # intron
		    my ($blk_begin, $blk_end);
		    if ($exin eq $tag_sort[0]) {
			$blk_begin = $trannoEnt->{rnaBegin};
		    }
		    else {
			$blk_begin = ($cur_ex_start - 1)."+1";
		    }
		    if ($exin eq $tag_sort[-1]) {
			$blk_end = $trannoEnt->{rnaEnd};
		    }
		    else {
			$blk_end   = $cur_ex_start."-1";
		    }

                    return getIntrPos( $blk_begin, $blk_end, $cur_blk_ofst,
                        $cur_blk_len );
		}
	    }
	    else {
		return ($cur_ex_start + $cur_blk_ofst);
	    }
	}
    }
    confess "Error: out of trRefComp range [$trRef_ofst].";
}

# calculate position in intron region
sub getIntrPos {
    my ($begin, $end, $ofst, $blk_len) = @_;
    
    my ($exPosBegin, $exPosEnd, $intPosBegin, $intPosEnd);
    if ($begin =~ /^(\d+)([\+\-]\d+)$/) {
	$exPosBegin = $1;
	$intPosBegin = $2;
    }
    if ($end =~ /^(\d+)([\+\-]\d+)$/) {
	$exPosEnd = $1;
	$intPosEnd = $2;
    }

    if ($exPosBegin eq $exPosEnd) {
        return
            $exPosBegin
          . ( ( $intPosBegin < 0 ) ? "" : '+' )
          . ( $intPosBegin + $ofst );
    }
    else {
	my $wlen = $intPosBegin - $intPosEnd + $blk_len - 2;
	my $lofst = $intPosBegin + $ofst - 1;
	if ($lofst < ($wlen / 2 - 1)) {
	    return ($exPosBegin . '+' . ($lofst + 1));
	}
	else {
	    return ($exPosEnd . '-' . ($wlen - $lofst));
	}
    }
}


=head2 getTrRef

    About   : generate concatenated transcript originated reference
    Usage   : my $trRef = getTrRef( $trannoEnt, $refgenome, $trSeq, $strd );
    Args    : trannoEnt - BedAnno::Anno->{trInfo}->{$tid}
              refgenome - Unified reference in BedAnno::Var
              trSeq     - whole transcript
              strd      - strand of transcript.
    Returns : transcript originated reference
    Notes   : Using sequence from transcript as the exon part,
              and using the sequence from reference genome
              as the intron part. and concatenate them.

=cut
sub getTrRef {
    my ( $trannoEnt, $refgenome, $trSeq, $strd ) = @_;
    my @tag_sort = sort trRefSort keys %{ $trannoEnt->{trRefComp} };
    my $trRef    = "";
    my $trStart  = getTrStart( $trannoEnt->{rnaBegin} );
    foreach my $exin (@tag_sort) {
        if ( $exin !~ /^EX/ ) {
            my $int_seq = substr(
                $refgenome,
                $trannoEnt->{trRefComp}->{$exin}->[0],
                (
                    $trannoEnt->{trRefComp}->{$exin}->[1] -
                      $trannoEnt->{trRefComp}->{$exin}->[0]
                )
            );
            $int_seq = BedAnno->rev_comp($int_seq) if ( !$strd );
            $trRef .= $int_seq;
        }
        else {
            $trRef .=
              substr( $trSeq, $trStart, $trannoEnt->{trRefComp}->{$exin} );
            $trStart += $trannoEnt->{trRefComp}->{$exin};
        }
    }
    return $trRef;
}

# get the next or current on-transcript position
# from any nDot format HGVS position string.
sub getTrStart {
    my $trStart = shift;
    if ($trStart =~ /^-/) {
	return 0;
    }
    elsif ($trStart =~ /^(\d+)\-/) {
	return ($1 - 1);
    }
    elsif ($trStart =~ /^(\d+)\+/) {
	return $1;
    }
    elsif ($trStart =~ /^(\d+)$/) {
	return ($1 - 1);
    }
    else {
	confess "Error: not available trStart [$trStart]";
    }
}

sub trRefSort {
    my ( $atag, $anum, $btag, $bnum );
    if ( $a =~ /^(\D+)(\d+)/ ) {
        $atag = $1;
        $anum = $2;
    }
    if ( $b =~ /^(\D+)(\d+)/ ) {
        $btag = $1;
        $bnum = $2;
    }
    confess "Error exon intron number [$a, $b]"
      if ( !defined $atag
        or !defined $anum
        or !defined $btag
        or !defined $bnum );
    $anum <=> $bnum or $atag cmp $btag;
}

=head2 batch_anno

    About   : The fastest way to annotate multiple snv and 1bp deletion variations,
              indel and other types also can be annotated, but no faster than annotated
              one by one.
    Usage   : $beda = BedAnno->new( db => 'in.bed.gz', tr => 'in.trans.fas', batch => 1);
              $rAnnoRst = $beda->batch_anno($rVars);
    Args    : an array ref of BedAnno::Var entries.
    Returns : an array ref of BedAnno::Anno entries, see varanno().

=cut
sub batch_anno {
    my ($self, $rVars) = @_;
    my @all_annoRst = ();
    my @sorted_vars =
      sort { $a->{chr} cmp $b->{chr} or $a->{pos} <=> $b->{pos} } @$rVars;

    my $cur_chr = "";
    my $anno_idx = 0;
    foreach my $var (@sorted_vars) {
	if ($var->{chr} ne $cur_chr) {
	    $cur_chr = $var->{chr};
	    $anno_idx = 0;
	}
	my $annoEnt;
	($annoEnt, $anno_idx) = $self->varanno($var, $anno_idx);
	push (@all_annoRst, $annoEnt);
    }
    return \@all_annoRst;
}

sub rev_comp {
    my $self = shift;
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


=head1 BedAnno::Var

    BedAnno::Var sub module

=cut
package BedAnno::Var;
use base qw(BedAnno);
use strict;
use warnings;
use Carp;
 
=head1 METHOD

=head2 new

    About   : Create a new object class entry, BedAnno::Var,
              parse the variation directly by the ref and alt string.
    Usage   : my $var = BedAnno::Var->new( $chr, $start, $end, $ref, $alt );
           or my $var = BedAnno::Var->new( $chr, $pos, $ref, $alt );
           or my $var = BedAnno::Var->new( $varInput );
    Args    : Input can be variable format
              1. 5 parameters format: CG shell list format: chr,start,end,ref,alt
              2. 4 parameters format: VCF format: chr,pos,ref,alt
              3. Crawler input object: A hash ref with nessesary keys: 
                 chr,begin,referenceSequence,variantSequence,  
                 optional key is "end", if end specified,
                 coordinates are treat as 0-based, otherwise, use 1-based (VCF)
    Returns : a new BedAnno::Var entry :
            {
                chr    => $chr,
                pos    => $start,          # 0-based start
                end    => $end,
                ref    => $ref,
                alt    => $alt,
                reflen => $ref_len,
                altlen => $alt_len,        # not exists if no-call
                guess  => $varType,        # the output varType
                imp    => $imp_varType,    # the implicit varType
                sm     => $sm,             # single/multiple base indicator
                                           # equal/non-equal length indicator

                # for hgvs naming convinient, reparse delins(guess),
                # in forward and reverse strand separately,
                # If the result are the same, then only
                # give the following optional
                # rescaled strand-same description group
                bp  => $bc_pos,       # backward compatible pos, 0-based
                br  => $bc_ref,       # backward compatible ref string
                ba  => $bc_alt,       # backward compatible alt string
                brl => $bc_reflen,    # backward compatible ref length
                bal => $bc_altlen,    # backward compatible alt length

                # otherwise, the following '+', '-',
                # structure will be generated to reflect
                # the difference. they are all optional

                '+' => {

                  # This group simplely trim off the leading same chars
                  # on forward strand, and then trim the same tail
                  bp  => $backward_fpos,
                  br  => $backward_fref,
                  ba  => $backward_falt,
                  brl => $backward_freflen,
                  bal => $backward_faltlen,

                },

                '-' => { 
                  # similar to '+', but for reverse strand 
                },

                # this group gives ref/alt string based on the rule
                # with 'rep' annotation available
                p      => $rep_left_pos,         # repeat pos, 0-based
                r      => $rep_ref,              # repeat ref string
                a      => $rep_alt,              # repeat alt string
                rl     => $rep_reflen,           # repeat ref length
                al     => $rep_altlen,           # repeat alt length
                rep    => $repeat_element,
                replen => $repeat_element_length,
                ref_cn => $copy_number_in_ref,
                alt_cn => $copy_number_in_alt,

                # for equal length long substitution
                # record the separated snvs positions
                # all positions are 1 based.
                sep_snvs => [ $snv_pos1, $snv_pos2, ... ],
            }

=cut
sub new {
    my $class = shift;
    my ($chr, $start, $end, $ref, $alt);
    if (ref($_[0])) {
        if (   !exists $_[0]->{chr}
            or !exists $_[0]->{begin}
            or !exists $_[0]->{variantSequence}
            or !exists $_[0]->{referenceSequence} )
        {
            confess "Error: unavailable object. need keys: ",
              "chr, start, alt, ref specified.";
        }
	
	$chr = $_[0]->{chr};
	$start = $_[0]->{begin};
	$alt = $_[0]->{variantSequence};
	$ref = $_[0]->{referenceSequence};
	$end = $_[0]->{end} if ( exists $_[0]->{end} );

    }
    else {
	confess "Error: not enough args, need at least 4 args." if (4 > @_);
	($chr, $start, $end, $ref, $alt) = @_;
    }


    if (!defined $end or $end !~ /^\d+$/) { # from VCF v4.1 1-based start
	if (defined $end) {
	    $alt = $ref;
	    $ref = $end;
	}

	$ref = normalise_seq($ref);
	$alt = normalise_seq($alt);
	if (substr($ref,0,1) eq substr($alt,0,1)) {
	    $ref = substr($ref,1);
	    $alt = substr($alt,1);
	}
	else {
	    $start -= 1; # change to 0-based start
	}
	my $rl = length($ref);
	$end = $start + $rl;
    }

    my %var;

    my $len_ref = $end - $start;            # chance to annotate long range
    my ($varType, $implicit_varType, $sm) = guess_type($len_ref, $ref, $alt);

    $chr = normalise_chr($chr);
    %var = (
        chr    => $chr,
        pos    => $start,
        ref    => uc($ref),
	end    => $end,
        alt    => uc($alt),
        reflen => $len_ref,
	guess  => $varType,
	imp    => $implicit_varType,
	sm     => $sm,
    );
    $var{altlen} = length($alt) if ($alt ne '?');

    my $var;
    $var = {%var};
    bless $var, ref($class) || $class;

    if ($varType eq 'no-call' or $implicit_varType ne 'delins') {
	return $var;
    }
    
    return $var->parse_complex();
}

sub TO_JSON {
    return { %{ shift() } };
}

=head2 getUnifiedVar

    About   : uniform the pos and ref/alt pair selection,
              after new(), give info for HGVS naming.
    Usage   : my @unified_desc = $var->getUnifiedVar($strd);
    Args    : BedAnno::Var entry and current strand for annotation.
    Returns : an array of ( 
                $pos,  # 0-based start pos
                $ref,  # reference bases
                $alt,  # called bases
                $reflen, # reference len
                $altlen )# called len, undef if no-call

=cut
sub getUnifiedVar {
    my $var = shift;
    my $strd = shift;

    my $consPos = $$var{pos};
    my $consRef = $$var{ref};
    my $consAlt = $$var{alt};
    my $consRL	= $$var{reflen};

    if (!exists $$var{altlen}) {
	return ($consPos, $consRef, $consAlt, $consRL);
    }

    my $consAL	= $$var{altlen};
    
    if (exists $var->{p}) { # rep
	$consPos = $$var{p};
	$consRef = $$var{r};
	$consAlt = $$var{a};
	$consRL	 = $$var{rl};
	$consAL	 = $$var{al};
    }
    elsif (exists $var->{bp}) { # complex bc strand same
	$consPos = $$var{bp};
	$consRef = $$var{br};
	$consAlt = $$var{ba};
	$consRL	 = $$var{brl};
	$consAL	 = $$var{bal};
    }
    elsif (exists $var->{$strd}) { # complex bc strand diff
	$consPos = $$var{$strd}{bp};
	$consRef = $$var{$strd}{br};
	$consAlt = $$var{$strd}{ba};
	$consRL	 = $$var{$strd}{brl};
	$consAL	 = $$var{$strd}{bal};
    }

    return ($consPos, $consRef, $consAlt, $consRL, $consAL);
}


sub normalise_chr {
    my $chr = shift;
    $chr =~ s/^chr//i;
    if ($chr eq '23') {
	$chr = 'X';
    }
    elsif ($chr eq '24') {
	$chr = 'Y';
    }
    elsif ($chr eq '25') {
	$chr = 'MT';
    }
    elsif ($chr =~ /^[xym]$/) {
	$chr = uc($chr);
	$chr .= 'T' if ($chr eq 'M');
    }
    elsif ($chr =~ /^M_/i) {
	$chr = 'MT';
    }
    return $chr;
}

sub normalise_seq {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq = "" if ($seq eq '.');
    $seq = uc($seq);
    if ($seq =~ /[^ACGTN]/ and $seq ne '?') {
	confess "Error: unrecognized pattern exists,",
	    "no multiple alts please. [$seq]";
    }
    return $seq;
}

=head2 parse_complex
    
    About   : parse complex delins variants to recognize 
              repeat and differ strand-pos var.
    Usage   : my $var = $var->parse_complex();
    Args    : variantion entry, which have been uniform to 
              CG's shell list format, with its 'guess':delins.
    Returns : see BedAnno::Var->new()

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
		$var->{$strd}->{brl} = $get_rst->{r};
		$var->{$strd}->{bal} = $get_rst->{a};
		$var->{$strd}->{br} =
		  substr( $var->{ref}, $get_rst->{$strd}, $get_rst->{r} );
		$var->{$strd}->{ba} =
		  substr( $var->{alt}, $get_rst->{$strd}, $get_rst->{a} );
	    }
	}
	else {
	    $var->{bp}  = $var->{pos} + $get_rst->{'+'};
	    $var->{brl} = $get_rst->{r};
	    $var->{bal} = $get_rst->{a};
	    $var->{br} = substr( $var->{ref}, $get_rst->{'+'}, $get_rst->{r} );
	    $var->{ba} = substr( $var->{alt}, $get_rst->{'+'}, $get_rst->{a} );
	}
    }

    if ($var->{sm} == 3) { # equal length long subs
	my @separate_snvs = ();
	for (my $i = 0; $i < $var->{reflen}; $i++) {
	    next if (substr($var->{ref}, $i, 1) eq substr($var->{alt}, $i, 1));
	    push (@separate_snvs,  $var->{pos} + $i + 1); # 1-based
	}
	$var->{sep_snvs} = [@separate_snvs];
    }

    my $rc_ref = count_content($ref);
    my $rc_alt = count_content($alt);
    my @diff = map { $$rc_ref[$_] - $$rc_alt[$_] } (0 .. $AAcount);

    # check if the sign of all base diff are consistent.
    if (check_sign( \@diff )) { # possible short tandom repeat variation
	my @absdiff = map {abs} @diff;
	my ($larger, $smaller, $llen, $slen);
        if ( $len_ref > $len_alt ) {
            $larger  = $ref;
            $llen    = $len_ref;
            $smaller = $alt;
            $slen    = $len_alt;
        }
        else {
            $larger  = $alt;
            $llen    = $len_alt;
            $smaller = $ref;
            $slen    = $len_ref;
        }

	my %has = ();
	for (my $rep = $llen; $rep > 0; $rep--) {
	    while ($larger =~ /([A-Z]+)(?:\1){$rep}/g) {
		next if (exists $has{$1});

                my $rep_el = $1;
                my $lofs   = length($`);    # $` is the prematched string

                my $cn = check_div( $rep_el, \@absdiff );
		if ($cn and check_insrep($larger, $smaller, $rep_el, $cn)) {
		    my $lenrep = length($rep_el);

                    @$var{qw(p rep replen)} =
                      ( ( $var->{pos} + $lofs ), $rep_el, $lenrep );

		    $rep += 1; # add the first copy of element

		    my $l_cn = $rep;
		    my $s_cn = $rep - $cn;
		    my $l_wholerep = $rep_el x $l_cn;
		    my $s_wholerep = $rep_el x $s_cn;
		    my $l_replen   = $lenrep * $l_cn;
		    my $s_replen   = $lenrep * $s_cn;

		    if ($llen == $len_ref) { # ref is longer
                        @$var{qw(ref_cn alt_cn r a rl al)} = (
                            $l_cn, $s_cn, $l_wholerep, $s_wholerep, $l_replen,
                            $s_replen
                        );
		    }
		    else { # alt is longer
                        @$var{qw(alt_cn ref_cn a r al rl)} = (
                            $l_cn, $s_cn, $l_wholerep, $s_wholerep, $l_replen,
                            $s_replen
                        );
		    }

		    $var->{imp} = 'rep';
		    $var->{sm}  = ($var->{rl} == 1) ? 1 : 2;

		    return $var;
		}
		$has{$rep_el} = 1;
	    }
	}
    }
    return $var;
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


=head2 guess_type

    About   : guess the varType directly from the input information.
    Usage   : my ($guess, $implicit_varType, $sm) = guess_type($len_ref, $ref, $alt);
    Args    : 1. length of reference (derived from end - start)
              2. reference sequence ('.' or empty for ins )
              3. called sequence ( '?' for no-call, '.' or empty for del )
    Returns : $guess is varType in output (ref,snv,ins,del,delins,no-call)
              $implicit_varType (ref,snv,ins,del,delins,rep)
              $sm is single/multiple/equal/non-equal-len indicator
                0 - ins case
                1 - single base ref case
                2 - multiple base, different length case
                3 - multiple base, equal length case 

=cut
sub guess_type {
    my ($len_ref, $ref, $alt) = @_;
    # imp_varType: implicit variant type is for HGVS naming
    # varType    : to be output as varType name. can be as the key
    #		   to get the SO term by Name2SO hash.
    # sm	 : single or multiple bases tag
    #		   0 - for insertion, 0 base
    #		   1 - for single base variants
    #		   2 - for non-equal-length multiple bases variants
    #		   3 - for equal-length multiple bases delins
    my ($imp_varType, $varType, $sm);
    if ($len_ref == 0) {
	if ($ref eq $alt) {
	    $imp_varType = 'ref';
	}
	else {
	    $imp_varType = 'ins';
	}
	$sm = 0;
    }
    elsif ($len_ref == 1) {
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
    elsif ($len_ref > 1) {
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
	    $sm = 2; # non-equal-length delins
	}
	else {
	    $sm = 3; # equal-length subs
	}
    }
    $varType = ($alt eq '?' or $alt eq 'N') ? 'no-call' : $imp_varType;
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
    for (1 .. $AAcount) {
        return 0 if ($$rc[$_] != $$rcs[$_] * $div);
    }
    return $div;
}

sub count_content {
    my $s = uc(shift);
    my $l = length($s);
    my @count = ($l, (0) x $AAcount);
    while ($s=~/(.)/g) {
	confess "unknown base [$1]" if (!exists $AAnumber{$1});
	$count[$AAnumber{$1}] ++;
    }
    return \@count;
}

=head2 get_gHGVS

    About   : get genomic (chromosomal) HGVS string of variation
    Usage   : my $gHGVS = $var->get_gHGVS();
    Args    : variation entry, after BedAnno::Var->new().
    Returns : chromosomal HGVS string.

=cut
sub get_gHGVS {
    my $var = shift;
    my $gHGVS = 'g.';
    if ($var->{chr} =~ /^M/) { # hit mito chromosome
	$gHGVS = 'm.';
    }

    my $imp = $var->{imp};
    my $sm  = $var->{sm};
    my ($pos, $ref, $alt, $reflen, $altlen) = $var->getUnifiedVar('+');

    if ($imp eq 'snv') {
	$gHGVS .= $pos.$ref.'>'.$alt;
    }
    elsif ($imp eq 'ins') {
	$gHGVS .= $pos.'_'.($pos+1).'ins'.$alt;
    }
    elsif ($imp eq 'del' or $imp eq 'delins') {
	# 1bp del/delins
	$gHGVS .= ($pos+1);
	if ($sm > 1) {
	    $gHGVS .= '_'.($pos + $reflen);
	}
	$gHGVS .= 'del'.$ref;
	$gHGVS .= 'ins'.$alt if ( $imp =~ /delins/);
    }
    elsif ($imp eq 'rep') {
	$gHGVS .= ($pos + 1);
	if ($var->{ref_cn} == 1 and $var->{alt_cn} == 2) { # dup
	    if ($sm > 1) {
		$gHGVS .= '_' . ($pos + $var->{replen});
	    }
	    $gHGVS .= 'dup'. $var->{rep};
	}
	else {
	    $gHGVS .= $var->{rep}.'['.$var->{ref_cn}.'>'.$var->{alt_cn}.']';
	}
    }
    elsif ($imp eq 'ref') {
	$gHGVS .= '=';
    }
    else {
	confess "Can not recognize type $imp.";
    }
    $var->{gHGVS} = $gHGVS;

    return $var;
}

=head1 BedAnno::Anno

    BedAnno::Anno sub package

=cut

package BedAnno::Anno;
use base qw(BedAnno BedAnno::Var);
use threads::shared;
use strict;
use warnings;

use Data::Dumper;
use Carp;

=head1 METHOD

=head2 new

    About   : Create BedAnno::Anno object
    Usage   : my $annoEnt = BedAnno::Anno->new($var);
    Args    : BedAnno::Var entry
    Returns : BedAnno::Anno entry

=cut
sub new {
    my $class = shift;
    my $var = shift;
    my $self;
    $self->{var} = $var;
    bless $self, ref($class) || $class;
    return $self;
}

sub TO_JSON {
    return { %{ shift() } };
}

=head2 decide_major

    About   : In finalise step, decide a major transcript to report for a var.
    Usage   : my $majorTranscriptVarName = $annoEnt->decide_major();
    Returns : A string in the following format:
	      If the transcript has a pName: 
		    mrnaAcc(geneSymbol): cName (pName), 
		    e.g. NM_145651.2(SCGB1C1): c.13C>T (p.R5C)
	      If the transcript does not have a pName: 
		    mrnaAcc(geneSymbol): cName
	      If only intergenic
		    chr: gName (intergenic)
		    e.g. chrX: g.220025A>T (intergenic)
    Note    : use the primaryTag to find reference standard or primary transcript,
	      if only one have the primaryTag "Y", then use it,
	      otherwise, sort by GenePart: 
		  1.  CDS
		  2.  span
		  3.  five_prime_cis_splice_site
		  4.  three_prime_cis_splice_site
		  5.  ncRNA
		  6.  five_prime_UTR
		  7.  three_prime_UTR
		  8.  interior_intron
		  9.  five_prime_UTR_intron
		  10. three_prime_UTR_intron
		  11. abnormal-intron
		  12. promoter
		  13. annotation-fail
		  14. intergenic_region
	      and choose the first one, if more than one transcript have the 
	      same reference standard and same GenePart, then choose the first
	      one which prior in name sorting.

=cut
sub decide_major {
    my $annoEnt = shift;
    if (   !exists $annoEnt->{trInfo}
        or 0 == scalar keys %{ $annoEnt->{trInfo} }
        or exists $annoEnt->{trInfo}->{""} )
    {
        my $intergenic =
          $annoEnt->{var}->{chr} . ": " . $annoEnt->{var}->{gHGVS};
        $intergenic = "chr" . $intergenic if ( $intergenic =~ /^[\dXYM]/ );
        $intergenic .= " (intergenic)";    # for only intergenic
        return $intergenic;
    }
    else {
	my %prTrs = ();
	my %nonPrTrs = ();
	foreach my $tid (keys %{$annoEnt->{trInfo}}) {
	    if ( !exists $annoEnt->{trInfo}->{$tid}->{genepart} ) {
		confess "Error: cannot decide major transcript before finaliseAnno";
	    }
	    if ( $annoEnt->{trInfo}->{$tid}->{primaryTag} eq "Y" ) {
		$prTrs{$tid} = $annoEnt->{trInfo}->{$tid}->{genepart};
	    }
	    else {
		$nonPrTrs{$tid} = $annoEnt->{trInfo}->{$tid}->{genepart};
	    }
	}

	my $majorTr;
	if ( 0 < keys %prTrs ) { # hit primary transcripts
	    if ( 1 == keys %prTrs ) {
		my ($majorTid) = keys %prTrs;
		$majorTr = $majorTid;
	    }
	    else {
		$majorTr = get_first_tr(\%prTrs);
	    }
	}
	else { # hit only non-primary transcripts
	    $majorTr = get_first_tr(\%nonPrTrs);
	}

        my $rTrEnt = $annoEnt->{trInfo}->{$majorTr};
        my $trhit  = $majorTr . "(" . $rTrEnt->{geneSym} . "): " . $rTrEnt->{c};
        $trhit .= " (" . $rTrEnt->{p} . ")"
          if ( exists $rTrEnt->{p} and $rTrEnt->{p} ne "" );
	return $trhit;
    }
}

sub get_first_tr {
    my $rGeneParts_hash = shift;
    my @ordered_Trs = sort {
        $GenePartsOrder{ $rGeneParts_hash->{$a} }
          <=> $GenePartsOrder{ $rGeneParts_hash->{$b} }
          or $a cmp $b
    } keys %$rGeneParts_hash;

    return $ordered_Trs[0];
}

=head2 reformatAnno

    About   : reformat the BedAnno::Anno entry to fit the need of crawler
    Usage   : my $crawler_need = $anno->reformatAnno();
    Note    : ID mapping changes are list the following:
	      var group:
	      refbuild -> referenceBuild
	      chr -> chr
	      pos -> begin
	      end -> end
	      ref -> referenceSequence
	      alt -> variantSequence
	      guess -> varType
	      varTypeSO -> varTypeSO
	      cytoBand -> cytoband
	      dbsnp => {rsID1,rsID2 ... } -> dbsnpIds
	      cg54 => AF -> CG54_AF
	      cg54 => AN -> CG54_AN
	      tgp => AF -> 1000G_AF
	      tgp => AN -> 1000G_AN
	      wellderly => AF -> Wellderly_AF
	      wellderly => AN -> Wellderly_AN
	      esp6500 => AF -> ESP6500_AF
	      esp6500 => AN -> ESP6500_AN
	      phyloPpm -> PhyloPscorePlacentalMammals
	      phyloPpr -> PhyloPscorePrimates
	      phyloPve -> PhyloPscoreVetebrates
	      gHGVS -> gHGVS
	      TranscriptVarName => $anno->decide_major();

	      trInfo group:
	      geneId -> GeneID
	      geneSym -> GeneSym
	      prot -> ProteinAccession
	      strd -> TranscriptOrientation
	      rnaBegin -> TranscriptBegin
	      rnaEnd -> TranscriptEnd
	      protBegin -> ProteinBegin
	      protEnd -> ProteinEnd
	      pfamId -> PFAM_ID
	      pfamName -> PFAM_NAME
	      genepart -> GenePart
	      genepartSO -> GenePartSO
	      genepartIndex -> GenePartIndex
	      exonIndex -> ExonNumber
	      intronIndex -> IntronNumer
	      funcSOname -> FunctionImpact
	      funcSO -> FunctionImpactSO
	      c -> cHGVS
	      p -> pHGVS
	      cc -> CodonChange
	      polar -> AAPolarityRef + AAPolarityVar
	      siftPred -> SIFTpred
	      siftScore -> SIFTscore
	      pp2divPred -> Polyphen2HumDivPred
	      pp2divScore -> Polyphen2HumDivScore
	      pp2varPred -> Polyphen2VarPred
	      pp2varScore -> Polyphen2VarScore

=cut
sub reformatAnno {
    my $anno = shift;
    my $var = $anno->{var};
    my $crawler_need = {
        var => {
            referenceBuild    => $var->{refbuild},
            chr               => $var->{chr},
            begin             => $var->{pos},
            end               => $var->{end},
            referenceSequence => $var->{ref},
            variantSequence   => $var->{alt},
            varType           => $var->{guess},
            varTypeSO         => ($var->{varTypeSO} =~ /SO:/) ? $var->{varTypeSO} : "",
            gHGVS             => $var->{gHGVS},
	    TranscriptVarName => $anno->decide_major(),
            cytoband => ( ( exists $var->{cytoBand} ) ? $var->{cytoBand} : "" ),
            dbsnpIds => (
                ( exists $var->{dbsnp} )
                ? join( ";", sort keys %{ $var->{dbsnp} } )
                : ""
            ),
            CG54_AF    => ( ( exists $var->{cg54} ) ? $var->{cg54}->{AF} : "" ),
            CG54_AN    => ( ( exists $var->{cg54} ) ? $var->{cg54}->{AN} : "" ),
            "1000G_AF" => ( ( exists $var->{tgp} )  ? $var->{tgp}->{AF}  : "" ),
            "1000G_AN" => ( ( exists $var->{tgp} )  ? $var->{tgp}->{AN}  : "" ),
            Wellderly_AF =>
              ( ( exists $var->{wellderly} ) ? $var->{wellderly}->{AF} : "" ),
            Wellderly_AN =>
              ( ( exists $var->{wellderly} ) ? $var->{wellderly}->{AN} : "" ),
            ESP6500_AF =>
              ( ( exists $var->{esp6500} ) ? $var->{esp6500}->{AF} : "" ),
            ESP6500_AN =>
              ( ( exists $var->{esp6500} ) ? $var->{esp6500}->{AN} : "" ),
            PhyloPscorePlacentalMammals =>
              ( ( exists $var->{phyloPpm} ) ? $var->{phyloPpm} : "" ),
            PhyloPscorePrimates =>
              ( ( exists $var->{phyloPpr} ) ? $var->{phyloPpr} : "" ),
            PhyloPscoreVetebrates =>
              ( ( exists $var->{phyloPve} ) ? $var->{phyloPve} : "" ),
        },
    };

    if ( !exists $anno->{trInfo} ) {
	my @trTags = qw(GeneID GeneSym ProteinAccession 
	  TranscriptOrientation TranscriptBegin TranscriptEnd 
	  ProteinBegin ProteinEnd 
	  PFAM_ID PFAM_NAME 
	  GenePartIndex ExonNumber IntronNumer 
	  FunctionImpact FunctionImpactSO 
	  cHGVS pHGVS CodonChange 
	  AAPolarityRef AAPolarityVar 
	  SIFTpred SIFTscore 
	  Polyphen2HumDivPred Polyphen2HumDivScore 
	  Polyphen2VarPred Polyphen2VarScore);

	my %empty_tr;
	@empty_tr{@trTags} = ("") x (scalar @trTags);
	$crawler_need->{trInfo}->{""} = {%empty_tr};
	$crawler_need->{trInfo}->{""}->{GenePart} = "intergenic_region";
	$crawler_need->{trInfo}->{""}->{GenePartSO} = "SO:0000605";
    }
    else {
	foreach my $trAcc (sort keys %{$anno->{trInfo}}) {
	    my $rTr = $anno->{trInfo}->{$trAcc};
	    my %trInfo = ();
            @trInfo{qw(GeneID GeneSym TranscriptOrientation)} =
              @$rTr{qw(geneId geneSym strd)};
            $trInfo{ProteinAccession} =
              ( exists $rTr->{prot} ) ? $rTr->{prot} : "";
	    $trInfo{TranscriptBegin} = 
	      ( exists $rTr->{rnaBegin} ) ? $rTr->{rnaBegin} : "";
	    $trInfo{TranscriptEnd} = 
	      ( exists $rTr->{rnaEnd} ) ? $rTr->{rnaEnd} : "";
	    $trInfo{ProteinBegin} = 
	      ( exists $rTr->{protBegin} ) ? $rTr->{protBegin} : "";
	    $trInfo{ProteinEnd} = 
	      ( exists $rTr->{protEnd} ) ? $rTr->{protEnd} : "";
	    $trInfo{PFAM_ID} = 
	      ( exists $rTr->{pfamId} ) ? $rTr->{pfamId} : "";
	    $trInfo{PFAM_NAME} = 
	      ( exists $rTr->{pfamName} ) ? $rTr->{pfamName} : "";
	    $trInfo{GenePart} = 
	      ( exists $rTr->{genepart} ) ? $rTr->{genepart} : "";
	    $trInfo{GenePartSO} = 
	      ( exists $rTr->{genepartSO} and $rTr->{genepartSO} =~ /SO:/ ) ? $rTr->{genepartSO} : "";
	    $trInfo{GenePartIndex} = 
	      ( exists $rTr->{genepartIndex} ) ? $rTr->{genepartIndex} : "";
	    $trInfo{ExonNumber} = 
	      ( exists $rTr->{exonIndex} ) ? $rTr->{exonIndex} : "";
	    $trInfo{IntronNumber} = 
	      ( exists $rTr->{intronIndex} ) ? $rTr->{intronIndex} : "";
	    $trInfo{FunctionImpact} = 
	      ( exists $rTr->{funcSOname} ) ? $rTr->{funcSOname} : "";
	    $trInfo{FunctionImpactSO} = 
	      ( exists $rTr->{funcSO} and $rTr->{funcSO} =~ /SO:/ ) ? $rTr->{funcSO} : "";
	    $trInfo{cHGVS} = 
	      ( exists $rTr->{c} ) ? $rTr->{c} : "";
	    $trInfo{pHGVS} = 
	      ( exists $rTr->{p} ) ? $rTr->{p} : "";
	    $trInfo{CodonChange} = 
	      ( exists $rTr->{cc} ) ? $rTr->{cc} : "";
	    $trInfo{SIFTpred} = 
	      ( exists $rTr->{siftPred} ) ? $rTr->{siftPred} : "";
	    $trInfo{SIFTscore} = 
	      ( exists $rTr->{siftScore} ) ? $rTr->{siftScore} : "";
	    $trInfo{Polyphen2HumDivPred} = 
	      ( exists $rTr->{pp2divPred} ) ? $rTr->{pp2divPred} : "";
	    $trInfo{Polyphen2HumDivScore} = 
	      ( exists $rTr->{pp2divScore} ) ? $rTr->{pp2divScore} : "";
	    $trInfo{Polyphen2VarPred} = 
	      ( exists $rTr->{pp2varPred} ) ? $rTr->{pp2varPred} : "";
	    $trInfo{Polyphen2VarScore} = 
	      ( exists $rTr->{pp2varScore} ) ? $rTr->{pp2varScore} : "";
	    $trInfo{AAPolarityRef} = "";
	    $trInfo{AAPolarityVar} = "";
	    if (exists $rTr->{polar} ) {
		my @pol = split(/=>/, $rTr->{polar});
		if (2 == @pol) {
		    $trInfo{AAPolarityRef} = $pol[0];
		    $trInfo{AAPolarityVar} = $pol[1];
		}
	    }

	    $crawler_need->{trInfo}->{$trAcc} = {%trInfo};
	}
    }
    return $crawler_need;
}

=head2 getTrPosition

    About   : Assign the BedAnno::Anno obj's trInfo with affected regions
    Usage   : my $AEIndex = $annoEnt->getTrPosition($rannodb, $AEIndex);
    Args    : $rannodb is a BedAnno::Anno object created by varanno(),
              $AEIndex is the current index for annodb searching.
    Returns : A new AEIndex for next query.
    Notes   : $AEIndex is used for same chr batch mode.
              assign the following tag to $annoEnt
              {
                trInfo => {
                    $tid => {
                        geneId,   geneSym, strd, primaryTag,
                        trAlt => $stranded_alt_string_with_ext_at_mismatches,

                        preStart => {
                            nDot, cDot, exin, r
                        },

                        postEnd => {
                            nDot, cDot, exin, r
                        },

                        trRefComp => {
                            $exon_number   => $transcript_exon_length,
                            $intron_number => [
                                $start_in_non_stranded_reference,
                                $stop_in_non_stranded_reference
                              ],
                            ...
                        },
                    },
                    ...
                }
              }

=cut
sub getTrPosition {
    my ($annoEnt, $rannodb, $aeIndex) = @_;

    my $var = $annoEnt->{var};
    # gather the covered entries
    my $new_aeIndex;

    my %tidExblk = ();
    my %rpreLeft = ();
    
    for (my $k = $aeIndex; $k < @$rannodb; $k ++) {
	if ($$rannodb[$k]{sto} < $var->{pos}) { # not reach var
	    $aeIndex ++;
	    next;
	}
	elsif ( $$rannodb[$k]{sta} > $var->{end} ) { # past var
	    last;
	}
	else { # covered by var
	    $new_aeIndex = $k if ($$rannodb[$k]{sta} <= $var->{pos}
		and $var->{pos} <= $$rannodb[$k]{sto} );
	    if (!exists $$rannodb[$k]{detail}) { # parse anno db
		$$rannodb[$k] = BedAnno->assign_detail($$rannodb[$k]);
	    }

	    foreach my $tid (keys %{$$rannodb[$k]{detail}}) {
		my $rtidDetail = $$rannodb[$k]{detail}{$tid};

		my $strd = ($rtidDetail->{strd} eq '+') ? 1 : 0;

                my ( $unify_p, $unify_r, $unify_a, $unify_rl, $unify_al ) =
                  $var->getUnifiedVar( $rtidDetail->{strd} );

		# skip non hitted block
                next
                  if ( $unify_p > $$rannodb[$k]{sto}
                    or ( $unify_p + $unify_rl ) < $$rannodb[$k]{sta} );


                my $total_left_ofst =
                  $$rtidDetail{offset} + ( $unify_p - $$rannodb[$k]{sta} );
                my $total_right_ofst = $total_left_ofst + $unify_rl;


                if (   !exists $annoEnt->{trInfo}
                    or !exists $annoEnt->{trInfo}->{$tid}
                    or !exists $annoEnt->{trInfo}->{$tid}->{trAlt} )
                {
		    # $tid added into trInfo
                    $annoEnt->{trInfo}->{$tid}->{trAlt} =
                      (       $unify_a =~ /^[ACGTN]+$/
                          and $rtidDetail->{strd} eq '-' )
                      ? BedAnno->rev_comp($unify_a)
                      : $unify_a;

		    my $tmp_trInfo = $annoEnt->{trInfo}->{$tid};

                    $tmp_trInfo->{geneId}     = $rtidDetail->{gid};
                    $tmp_trInfo->{geneSym}    = $rtidDetail->{gsym};
                    $tmp_trInfo->{strd}       = $rtidDetail->{strd};
                    $tmp_trInfo->{primaryTag} = $rtidDetail->{pr};
                    # assign left rna positions
		    # if no assignment then skip the right position assignment
                    if (
                        $annoEnt->cal_hgvs_pos(
                            tid       => $tid,
                            tidDetail => $rtidDetail,
                            offset    => $total_left_ofst,
                            LR        => 1,                  # left mode
                        )
                      )
                    {
                        if ( $rtidDetail->{blka} =~ /^PROM/ ) {

                            #			              use P0 for sort
                            $tmp_trInfo->{trRefComp}->{P0}->[0]
                              = 0;
                        }
                        elsif ( $rtidDetail->{exin} =~ /^IVS/ ) {
                            $tmp_trInfo->{trRefComp}
                              ->{ $rtidDetail->{exin} }->[0] = 0;
                        }
			elsif ( !$strd and $tmp_trInfo->{rnaEnd} =~ /^\+(\d+)/ ) {
			    $tmp_trInfo->{trRefComp}->{ Z999 }->[0] = 0;
			    $tmp_trInfo->{trRefComp}->{ Z999 }->[1] = $1;
			}
			# let right rna positions to record Exon block

			if ($total_left_ofst > 0) {
                            $rpreLeft{$tid} = $annoEnt->cal_hgvs_pos(
                                tid       => $tid,
                                tidDetail => $rtidDetail,
                                offset    => $total_left_ofst,
                                LR        => 0,                  # preLeft mode
                                noassign  => 1,
                            );
			}

			if (ref($rpreLeft{$tid})) {
			    if ($strd) {
				$tmp_trInfo->{preStart} = $rpreLeft{$tid};
			    }
			    else {
				$tmp_trInfo->{postEnd} = $rpreLeft{$tid};
			    }
			}
                    }
                    else {
			$rpreLeft{$tid} = $annoEnt->cal_hgvs_pos(
			    tid       => $tid,
			    tidDetail => $rtidDetail,
			    offset    => $total_left_ofst,
			    LR        => 0,                  # preLeft mode
			    noassign  => 1,
			);
                        next;
                    }
                }

		# $tid have been added into trInfo in the left end assignment
		my $trinfoEnt = $annoEnt->{trInfo}->{$tid};
		# assign right rna positions
                if (
                    $annoEnt->cal_hgvs_pos(
                        tid       => $tid,
                        tidDetail => $rtidDetail,
                        offset    => $total_right_ofst,
                        LR        => 0,                   # right mode
                    )
                  )
                {
		    # assign the trRefComp
		    my $tmp_exin = $rtidDetail->{exin};
                    if ( $tmp_exin =~ /^EX/ ) {
                        if (   !exists $tidExblk{$tid}
                            or !exists $tidExblk{$tid}{$tmp_exin}
                            or !exists $tidExblk{$tid}{$tmp_exin}{
				$rtidDetail->{nsta} . ',' . $rtidDetail->{nsto}
				}
                          )
                        {    # non Insertion on refseq, and first occured exon
                            $tidExblk{$tid}{$tmp_exin}{
				    $rtidDetail->{nsta} . ','
                                  . $rtidDetail->{nsto} } = 1;

                            my ( $cmpStart, $cmpEnd );
                            if ($strd) {
                                $cmpStart =
                                  (
                                    0 < BedAnno->cmpPos(
                                        $trinfoEnt->{rnaBegin},
                                        $rtidDetail->{nsta}
                                    )
                                  )
                                  ? $rtidDetail->{nsta}
                                  : $trinfoEnt->{rnaBegin}
                                  ;    # select the larger one
				$cmpEnd = 
                                  (
                                    0 < BedAnno->cmpPos(
                                        $trinfoEnt->{rnaEnd},
                                        $rtidDetail->{nsto}
                                    )
                                  )
                                  ? $trinfoEnt->{rnaEnd}
                                  : $rtidDetail->{nsto}
                                  ;    # select the smaller one
			    }
			    else {
                                $cmpStart =
                                  (
                                    0 < BedAnno->cmpPos(
                                        $trinfoEnt->{rnaBegin},
                                        $rtidDetail->{nsto}
                                    )
                                  )
                                  ? $rtidDetail->{nsto}
                                  : $trinfoEnt->{rnaBegin}
                                  ;    # select the larger one
                                $cmpEnd =
                                  (
                                    0 < BedAnno->cmpPos(
                                        $trinfoEnt->{rnaEnd},
                                        $rtidDetail->{nsta}
                                    )
                                  )
                                  ? $trinfoEnt->{rnaEnd}
                                  : $rtidDetail->{nsta}
                                  ;    # select the smaller one
                            }

                            if ( $cmpStart !~ /^\d+$/ or $cmpEnd !~ /^\d+$/ ) {
                                confess "Error: unknown bug cmpStart ",
				    "[$cmpStart], cmpEnd [$cmpEnd]";
                            }

                            $trinfoEnt->{trRefComp}->{$tmp_exin} +=
                              $cmpEnd - $cmpStart + 1;
                            if (    $tmp_exin =~ /E$/
                                and $strd
                                and $trinfoEnt->{rnaEnd} =~ /^\+(\d+)/ )
                            {
                                $trinfoEnt->{trRefComp}->{Z999}->[0] =
                                  $unify_rl - $1;
                                $trinfoEnt->{trRefComp}->{Z999}->[1] =
                                  $unify_rl;
                            }
                        }

                    }
		    else { # for intron region and PROM
			if ( $rtidDetail->{blka} =~ /^PROM/ ) {
			    $tmp_exin = 'P0'; # for sort
			}
                        if (   !exists $trinfoEnt->{trRefComp}
                            or !exists $trinfoEnt->{trRefComp}
                            ->{$tmp_exin} )
                        {
                            $trinfoEnt->{trRefComp}
                              ->{$tmp_exin}->[0] =
                              $$rannodb[$k]{sta} - $unify_p;
                        }
                        my $less =
                          ( $$rannodb[$k]{sto} > ($unify_p + $unify_rl) )
                          ? ($unify_p + $unify_rl)
                          : $$rannodb[$k]{sto};
                        $trinfoEnt->{trRefComp}->{$tmp_exin}
                          ->[1] = $less - $unify_p;
			if ( $tmp_exin eq 'P0' and !$strd ) {
			    $trinfoEnt->{trRefComp}->{$tmp_exin}->[1] = $unify_rl;
			}
		    }
                }
		my $rpostRight = $annoEnt->cal_hgvs_pos(
		    tid       => $tid,
		    tidDetail => $rtidDetail,
		    offset    => $total_right_ofst,
		    LR        => 1,                   # post right mode
		    noassign  => 1,
		);
		if (ref($rpostRight)) {
		    if ($strd) {
			$trinfoEnt->{postEnd} = $rpostRight;
		    }
		    else {
			$trinfoEnt->{preStart} = $rpostRight;
		    }
		}
	    }
        }
    }
    $new_aeIndex ||= $aeIndex;
    return $new_aeIndex;
}




=head2 cal_hgvs_pos

    About   : calculate nDot, cDot HGVS position, depend on given offset,
              assign trAlt string and nDot HGVS and cDot HGVS positions.
    Usage   : $annoEnt->cal_hgvs_pos(
                    offset => $offset, 
                    tid    => $tid,
                    LR     => $lr,
                    tidDetail => $rh_tidDetail,
              );
    Args    : ofst is total offset to the start(left) of currunt annoblk,
              tid is the transcript id for the detail entry
              tidDetail is the currunt annoblk detail
              LR indicate this offset is left or right pos,
                1 for left and assign sta group,
                0 for right and assign end group.
              "noassign" to indicate whether to assign those information
              to annoEnt, it's used return the useful information only,
              without change the annoEnt. 
    Returns : if noassign is used, then return a hash ref, which contains
                { nDot, cDot, exin, r } if successful.
              otherwise, 0 for no assigned status, 1 for successful assigned.
             
    Notes   : For position mirror on transcript, there are 2 other cases 
              than normal case:
              1. annotation fail, which can not be annotated in the region
                 of it except for its promoter region.
              2. block with length changing mismatches, or long substitution
                 mismatch, which contain the following three cases:
                 
                 a. insertion (I) on refSeq

                        +-------+---+--------+  refgenome
                         \       \ /        /
                          +-------+--------+    refSeq

                 b. deletion (D) on refSeq 

                          +-------+--------+    refgenome
                         /       / \        \
                        +-------+---+--------+  refSeq

                 c. delins (DI) on refSeq (equal/non-equal length)

                        +-------+---+--------+  refgenome
                        |       |  /        /
                        +-------+-+--------+    refSeq

                 Insertion will have an reversed start/end position on refSeq,
                 due to the 1-based position description system.
                 
                 Any position located in a non-zero length refgenome mismatch
                 block have to extend to total region of mismatched block,
                 and alternate "trAlt" value in annotation entry for this tid.

              This method assign the following tag to $annoEnt
                {
                    trInfo => {
                        $tid => {
                            rnaBegin, rnaEnd,  cdsBegin, cdsEnd,
                            ei_Begin, ei_End,  r_Begin,  r_End,
                            genepartSO, trAlt,
                          },
                        ...
                    }
                }
=cut
sub cal_hgvs_pos {
    my $annoEnt = shift;
    my %cal_args = @_;
    if (exists $cal_args{noassign} and $cal_args{noassign} == 0) {
	delete $cal_args{noassign};
    }
    
    my ( $ofst, $tid, $rtidDetail, $lr ) =
      @cal_args{qw(offset tid tidDetail LR)};
    my ( $nDot, $cDot ) = ( '', '' );
    my ( $lr_pair_nDot, $lr_pair_cDot ) = ( '', '' );
    my $strd = ($$rtidDetail{strd} eq '+') ? 1 : 0;
    my $trAlt = $annoEnt->{trInfo}->{$tid}->{trAlt} if (!exists $cal_args{noassign});

    my $lofst = $ofst;
    if ($lr) { # only one time assignment for one var, offset for left 
	my $rofst = $$rtidDetail{wlen} - $ofst - 1;
        if ($lofst < 0) {
	    if ($strd) {

		# 5'promoter is always available
		$nDot = $$rtidDetail{nsta} + $lofst;
		if ($$rtidDetail{csta} =~ /^(\S+)\-u(\d+)/) {
                    $cDot = $1 . '-u' . ( $2 - $lofst );
		}
	    }
	    else { # outside transcript region into 3'downstream
		$nDot = '+'.(0-$lofst);
		if ($$rtidDetail{csta} =~ /^\*?\d+$/) {
                    $cDot = $$rtidDetail{csta} . '+d' . ( 0 - $lofst );
		}
	    }
        }
	elsif ($lofst > $$rtidDetail{wlen}) {
	    if (!exists $cal_args{noassign}) {
		confess "Error: exceed the whole length [$lofst>$$rtidDetail{wlen}]";
	    }
	    else {
		return 0;
	    }
	}
	else { # 0 <= $lofst <= $$rtidDetail{wlen}
        
	    # 1. check if a mismatch block (only exon have this)
	    if ( $rtidDetail->{mismatch} ne "" ) {
		# for each kind of mismatch block,
		# hit the left edge of var, the refSeq ref will
		# contain start from the start of the mismatch
		# block, case I in database already has
		# a reversed coordinate between left and right.
		# transcript-alt string shoule be assigned to 
		# this tid
		$nDot = $rtidDetail->{nsta};
		$cDot = $rtidDetail->{csta};
		if (!exists $cal_args{noassign} and $lofst>0) {
		    # mismatch info record the transcript-stranded
		    # reference sequence
                    my ( $mType, $mStart, $mEnd, $strand_ref ) =
                      split( /,/, $rtidDetail->{mismatch} );
		    
		    if ($strd) { # cut left
			$trAlt = substr($strand_ref, 0, $lofst).$trAlt if ($trAlt ne '?');
		    }
		    else { # cut right
			$trAlt .= substr($strand_ref, -$lofst, $lofst) if ($trAlt ne '?');
		    }
		}
	    }
	    # 2. check if hit a normal block's right edge
	    elsif ( $lofst == $$rtidDetail{wlen} ) {
		# give back the chance for left edge parsing
		if (!exists $cal_args{noassign}) {
		    delete $annoEnt->{trInfo}->{$tid};
		}
		return 0;
	    }
	    # 3. check if annotation-fail
	    elsif ($rtidDetail->{gpSO} eq 'annotation-fail') {
		( $nDot, $cDot ) = ( '?', '?' );
	    }
	    # 4. check if a promoter
	    elsif ($rtidDetail->{blka} eq 'PROM') {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $lofst )
                  : ( $rtidDetail->{nsta} - $lofst );
		if ($rtidDetail->{csta} =~ /^(\S+)\-u(\d+)$/) {
                    $cDot =
                        ($strd)
                      ? ( $1 . '-u' . ( $2 - $lofst ) )
                      : ( $1 . '-u' . ( $2 + $lofst ) );
		}
	    }
	    # 5. check if an exon
	    elsif ($rtidDetail->{exin} =~ /EX/) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $lofst )
                  : ( $rtidDetail->{nsta} - $lofst );
		if ($rtidDetail->{csta} =~ /^(\*?)(\-?\d+)$/) {
		    $cDot = ($strd) ? $1.($2 + $lofst) : $1.($2 - $lofst);
		}
	    }
	    # 6. intron
	    elsif ($rtidDetail->{exin} =~ /IVS/) {
		my $half_length = $rtidDetail->{wlen} / 2;
		if ($lofst > $half_length) { # drop into latter part
                    if (   $rtidDetail->{nsto} =~ /^\d+\d+$/
                        or $rtidDetail->{csto} =~ /^\d+\d+$/ )
                    {
                        confess "Error format of database, ",
                          "intron right pos description error!\n";
                    }

		    if ($rtidDetail->{nsto} =~ /^(\d+\+?)(\-?\d+)$/) {
                        $nDot =
                            ($strd)
                          ? ( $1 . ( $2 - $rofst ) )
                          : ( $1 . ( $2 + $rofst ) );
		    }
		    if ($rtidDetail->{csto} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/) {
                        $cDot =
                            ($strd)
                          ? ( $1 . ( $2 - $rofst ) )
                          : ( $1 . ( $2 + $rofst ) );
		    }
		}
		else {
                    if (   $rtidDetail->{nsta} =~ /^\d+\d+$/
                        or $rtidDetail->{csta} =~ /^\d+\d+$/ )
                    {
                        confess "Error format of database, ",
                          "intron left pos description error!\n";
                    }

		    if ($rtidDetail->{nsta} =~ /^(\d+\+?)(\-?\d+)$/) {
                        $nDot =
                            ($strd)
                          ? ( $1 . ( $2 + $lofst ) )
                          : ( $1 . ( $2 - $lofst ) );
		    }
		    if ($rtidDetail->{csta} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/) {
                        $cDot =
                            ($strd)
                          ? ( $1 . ( $2 + $lofst ) )
                          : ( $1 . ( $2 - $lofst ) );
		    }
		}
	    }
	    else {
		confess "Error: what's this? $rtidDetail->{exin}\n";
	    }
	}

	if (!exists $cal_args{noassign}) {

	    # the Begin End assignment is always from 5'->3' except for insertion
	    if ($strd) {
		$annoEnt->{trInfo}->{$tid}->{rnaBegin} = $nDot;
		$annoEnt->{trInfo}->{$tid}->{cdsBegin} = $cDot;
		$annoEnt->{trInfo}->{$tid}->{ei_Begin} = $rtidDetail->{exin}; 
		$annoEnt->{trInfo}->{$tid}->{r_Begin}  = $rtidDetail->{blka};
	    }
	    else {
		$annoEnt->{trInfo}->{$tid}->{rnaEnd} = $nDot;
		$annoEnt->{trInfo}->{$tid}->{cdsEnd} = $cDot;
		$annoEnt->{trInfo}->{$tid}->{ei_End} = $rtidDetail->{exin};
		$annoEnt->{trInfo}->{$tid}->{r_End}  = $rtidDetail->{blka};
	    }
	    $annoEnt->{trInfo}->{$tid}->{genepartSO} = $rtidDetail->{gpSO};
	    $annoEnt->{trInfo}->{$tid}->{trAlt} = $trAlt;
	}
    }
    else { # can be multiple times assignment if large var, offset for right
	my $real_ofst = $lofst - 1; # end_ofst change
	my $rofst = $$rtidDetail{wlen} - $ofst;
	if ($lofst > $$rtidDetail{wlen}) {
	    my $ex_ofst = ($lofst - $$rtidDetail{wlen});
	    if ($strd) {
		$nDot = '+'.$ex_ofst;
                if ( $$rtidDetail{csto} =~ /^\*?\d+$/ ) {
                    $cDot = $$rtidDetail{csto} . '+d' . $ex_ofst;
                }
	    }
	    else { # outside 5'promoter region ?
                $nDot =
                    ( $$rtidDetail{nsto} =~ /^\-\d+$/ )
                  ? ( $$rtidDetail{nsto} - $ex_ofst )
                  : ( -$ex_ofst );
                if ( $$rtidDetail{csto} =~ /^(\S+)\-u(\d+)/ ) {
                    $cDot = $1 . '-u' . ( $2 - $ex_ofst );
                }
	    }
	}
	elsif ($lofst < 0) { # not proper called
	    if (!exists $cal_args{noassign}) {
		confess "Error: right preceed the block [$lofst < 0]";
	    }
	    else {
		return 0;
	    }
	}
	else { # 0 <= $lofst <= $$rtidDetail{wlen}

	    # 1. check if a mismatch block (exon)
	    if ($rtidDetail->{mismatch} ne "") {
		$nDot = $rtidDetail->{nsto};
		$cDot = $rtidDetail->{csto};
		if (!exists $cal_args{noassign} and $lofst < $rtidDetail->{wlen}) {
                    my ( $mType, $mStart, $mEnd, $strand_ref ) =
                      split( /,/, $rtidDetail->{mismatch} );
		    
		    my $cut_len = $rtidDetail->{wlen} - $lofst;
		    if ($strd) { # cut right
			$trAlt .= substr($strand_ref, -$cut_len, $cut_len) if ($trAlt ne '?');
		    }
		    else { # cut left
			$trAlt = substr($strand_ref, 0, $cut_len).$trAlt if ($trAlt ne '?');
		    }
		}
	    }
	    # 2. check if hit a normal block's left edge
	    elsif ($lofst == 0) {
		return 0;
	    }
	    # 3. check if annotation-fail
	    elsif ($rtidDetail->{gpSO} eq 'annotation-fail') {
		( $nDot, $cDot ) = ( '?', '?' );
	    }
	    # 4. check if a promoter
	    elsif ($rtidDetail->{blka} =~ /^PROM/) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $real_ofst )
                  : ( $rtidDetail->{nsta} - $real_ofst );
		if ($rtidDetail->{csta} =~ /^(\S+)\-u(\d+)$/) {
                    $cDot =
                        ($strd)
                      ? ( $1 . '-u' . ( $2 - $real_ofst ) )
                      : ( $1 . '-u' . ( $2 + $real_ofst ) );
		}
	    }
	    # 5. check if an exon
	    elsif ($rtidDetail->{exin} =~ /EX/) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $real_ofst )
                  : ( $rtidDetail->{nsta} - $real_ofst );
		if ($rtidDetail->{csta} =~ /^(\*?)(\-?\d+)$/) {
		    $cDot = ($strd) ? $1.($2 + $real_ofst) : $1.($2 - $real_ofst);
		}
	    }
	    # 6. intron
	    elsif ($rtidDetail->{exin} =~ /IVS/) {
		my $half_length = $rtidDetail->{wlen} / 2;
		if ($real_ofst > $half_length) { # drop into latter part
                    if (   $rtidDetail->{nsto} =~ /^\d+\d+$/
                        or $rtidDetail->{csto} =~ /^\d+\d+$/ )
                    {
                        confess "Error format of database, ",
                          "intron right pos description error!\n";
                    }

		    if ($rtidDetail->{nsto} =~ /^(\d+\+?)(\-?\d+)$/) {
                        $nDot =
                            ($strd)
                          ? ( $1 . ( $2 - $rofst ) )
                          : ( $1 . ( $2 + $rofst ) );
		    }
		    if ($rtidDetail->{csto} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/) {
                        $cDot =
                            ($strd)
                          ? ( $1 . ( $2 - $rofst ) )
                          : ( $1 . ( $2 + $rofst ) );
		    }
		}
		else {
                    if (   $rtidDetail->{nsta} =~ /^\d+\d+$/
                        or $rtidDetail->{csta} =~ /^\d+\d+$/ )
                    {
                        confess "Error format of database, ",
                          "intron left pos description error!\n";
                    }

		    if ($rtidDetail->{nsta} =~ /^(\d+\+?)(\-?\d+)$/) {
                        $nDot =
                            ($strd)
                          ? ( $1 . ( $2 + $real_ofst ) )
                          : ( $1 . ( $2 - $real_ofst ) );
		    }
		    if ($rtidDetail->{csta} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/) {
                        $cDot =
                            ($strd)
                          ? ( $1 . ( $2 + $real_ofst ) )
                          : ( $1 . ( $2 - $real_ofst ) );
		    }
		}
	    }
	    else {
		confess "Error: what's this? $rtidDetail->{exin}\n";
	    }
	}

	if (!exists $cal_args{noassign}) {
	    # the Begin End assignment is always from 5'->3' except for insertion
	    if ($strd) {
		$annoEnt->{trInfo}->{$tid}->{rnaEnd} = $nDot;
		$annoEnt->{trInfo}->{$tid}->{cdsEnd} = $cDot;
		$annoEnt->{trInfo}->{$tid}->{ei_End} = $rtidDetail->{exin};
		$annoEnt->{trInfo}->{$tid}->{r_End}  = $rtidDetail->{blka};
	    }
	    else {
		$annoEnt->{trInfo}->{$tid}->{rnaBegin} = $nDot;
		$annoEnt->{trInfo}->{$tid}->{cdsBegin} = $cDot;
		$annoEnt->{trInfo}->{$tid}->{ei_Begin} = $rtidDetail->{exin}; 
		$annoEnt->{trInfo}->{$tid}->{r_Begin}  = $rtidDetail->{blka};
	    }
	    $annoEnt->{trInfo}->{$tid}->{trAlt} = $trAlt;
	    
	    if (    exists $annoEnt->{trInfo}->{$tid}->{genepartSO}
		and $annoEnt->{trInfo}->{$tid}->{genepartSO} ne "span"
		and $rtidDetail->{gpSO} ne $annoEnt->{trInfo}->{$tid}->{genepartSO}
	      )
	    {
		$annoEnt->{trInfo}->{$tid}->{genepartSO} = "span";
	    }
	    elsif (!exists $annoEnt->{trInfo}->{$tid}->{genepartSO}) {
		carp "Warning: no genepartSO specified in left?";
		$annoEnt->{trInfo}->{$tid}->{genepartSO} = $rtidDetail->{gpSO};
	    }
	}
    }

    if (exists $cal_args{noassign}) {
	return {
	    nDot => $nDot,
	    cDot => $cDot,
	    exin => $rtidDetail->{exin},
	    r    => $rtidDetail->{blka},
	};
    }
    return 1;
}


1;
__END__

=head1 DATABASE FORMAT

The Format of annotation database is listed as following:

=head2 BED+1 FORMAT DATABASE

   Departed block with tag for annotation, Tag entries are separated by "; ",
   and Infos items in entry are separated by "|".

   Each entry contains the follwing infomation:

   1. Acc.Ver
   2. GeneID
   3. Gene Symbol
   4. Strand
                 5'=====|>>>|[=============]|>>>>>>|[==========]|>>>>>>|[=============]|>>>>|==3'
   5. BlockAttr  : PROM 5U2E D5U1 I5U1 A5U1 5U1 C1  DC1 IC1 AC1 C2E 3U1 D3U1 I3U1 A3U1 3U2E
   6. GenePartsSO: 167  204  163  447  164  204 316 163 191 164 316 205 163  448  164  448
   7. ExIn Num   :    . |EX1|      IVS1     |  EX2 |    IVS2    |  EX3 |    IVS3       |EX4E|
   8. nHGVS start for block before departing
   9. nHGVS end for block before departing
   10.cHGVS start for block before departing
   11.cHGVS end for block before departing
   12.Length for block before departing
   13.MismatchBlock :  $type,$gstart,$gstop,$gseq
                       ($gseq is in the strand of refseq, '.' for deletion)
   14.Primary Tag : Please see "PRIMARY TAG ASSIGNMENT"
   15.Offset to leftmost of non departing block.

=head2 TRANSCRIPT FASTA DATABASE

   One-line sequence fasta file
   ============================
   Header format is: ( separate by " ", with "." for unavailable value )

       >rnaAcc.ver rnaLen gene protAcc.ver protLen cdsSta,cdsEnd tags

=head2 PRIMARY TAG ASSIGNMENT

    Our database need to sort the refSeq record with same Acc, 
    or with same gene by the following rules:

    1. whether the transcript is on the primary assembly.
    2. whether the CDS is 3-codons (badcds will be put to the tail)
    3. concatenated CDS is longer
    4. concatenated Exon is longer
    5. union Exon with flank region is longer
    6. Chromosome ID number is smaller
    7. Position number on forward strand-chromosome is smaller

    if LRG_RefSeqGene file is used, then the primary tag "Y" will
    only assign to the reference standard transcripts's first 
    mapping entry. But for non-RefSeqGene gene/transcript,
    assign primary tag "Y" to the first record for same genes' all
    transcripts, which is the same with no LRG_RefSeqGene file case.

    For multiple-mapping of a same transcript, add postfix "-N" 
    (1..n-1) to the other records in the order of sort.

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
