package BedAnno;

use strict;
BEGIN {
   use Exporter   ();
   use vars       qw($VERSION @ISA @EXPORT @EXPORT_OK);
   # set the version for version checking
   $VERSION     = 2.00;
   @ISA         = qw(Exporter);
   @EXPORT      = qw();
   # your exported package globals go here,
   # as well as any optionally exported functions
   @EXPORT_OK   = qw($AAcount %AAnumber);
}
use warnings;
use Carp;
use Data::Dumper;
use IO::Uncompress::Gunzip qw($GunzipError);
use Time::HiRes qw(gettimeofday tv_interval);

use Tabix;

=head1 NAME

BedAnno - Perl module for annotating variation depend on the BED format database.

=head2 VERSION v2.00

From version 0.32 BedAnno will change to support CG's variant shell list
and use ncbi annotation release 104 as the annotation database
with reformatted database format, and won't give any individual
annotation, so the individual_anno() is no longer available.
VCF4.1 format variant description (chr, pos, ref, alt) will also
be supported.

=head1 SYNOPSIS

  use BedAnno;
  my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas' );
  my $anno = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );

=head1 DESCRIPTION

By using this module, we can get variants from whole-genome or exome-capture 
NGS genotyping result annotated. The information contains various possible 
HGVS string together with a most recent strandard HGVS mutation name. 
It can not annotate ambiguous variants (transition, transvertion, or unknown 
break point large deletion and duplication).

This module need bgzipped BED format database, combined with tabix index.
tabix perl module is required to be installed. If allele frequency information, 
prediction information or cytoBand information etc. are needed, then extra 
resource dependencies will be required.

=cut

our (
    %C3,                %C1,                 %SO2Name,
    %func2SO,           %Name2SO,            %Polar,
    %C1toC3,            $AAcount,            %AAnumber,    %canonicalSS,
);

our $CURRENT_MT  = 'NC_012920.1';
my $load_opt_vcfaf = 0;

%C3 = (
    AAA => "Lys", AAC => "Asn", AAG => "Lys", AAT => "Asn",
    ACA => "Thr", ACC => "Thr", ACG => "Thr", ACT => "Thr",
    AGA => "Arg", AGC => "Ser", AGG => "Arg", AGT => "Ser",
    ATA => "Ile", ATC => "Ile", ATG => "Met", ATT => "Ile",
    CAA => "Gln", CAC => "His", CAG => "Gln", CAT => "His",
    CCA => "Pro", CCC => "Pro", CCG => "Pro", CCT => "Pro",
    CGA => "Arg", CGC => "Arg", CGG => "Arg", CGT => "Arg",
    CTA => "Leu", CTC => "Leu", CTG => "Leu", CTT => "Leu",
    GAA => "Glu", GAC => "Asp", GAG => "Glu", GAT => "Asp",
    GCA => "Ala", GCC => "Ala", GCG => "Ala", GCT => "Ala",
    GGA => "Gly", GGC => "Gly", GGG => "Gly", GGT => "Gly",
    GTA => "Val", GTC => "Val", GTG => "Val", GTT => "Val",
    TAA => "*",   TAC => "Tyr", TAG => "*",   TAT => "Tyr",
    TCA => "Ser", TCC => "Ser", TCG => "Ser", TCT => "Ser",
    TGA => "*",   TGC => "Cys", TGG => "Trp", TGT => "Cys",
    TTA => "Leu", TTC => "Phe", TTG => "Leu", TTT => "Phe",

    TCN => "Ser", CCN => "Pro", ACN => "Thr", GTN => "Val",
    CTN => "Leu", GCN => "Ala", CGN => "Arg", GGN => "Gly",

    # inseq stop codon
    UAA => "X", UAG => "X",

    # selenocysteine
    UGA => "Sec"
);

%C1 = (
    AAA => "K", AAC => "N", AAG => "K", AAT => "N",
    ACA => "T", ACC => "T", ACG => "T", ACT => "T",
    AGA => "R", AGC => "S", AGG => "R", AGT => "S",
    ATA => "I", ATC => "I", ATG => "M", ATT => "I",
    CAA => "Q", CAC => "H", CAG => "Q", CAT => "H",
    CCA => "P", CCC => "P", CCG => "P", CCT => "P",
    CGA => "R", CGC => "R", CGG => "R", CGT => "R",
    CTA => "L", CTC => "L", CTG => "L", CTT => "L",
    GAA => "E", GAC => "D", GAG => "E", GAT => "D",
    GCA => "A", GCC => "A", GCG => "A", GCT => "A",
    GGA => "G", GGC => "G", GGG => "G", GGT => "G",
    GTA => "V", GTC => "V", GTG => "V", GTT => "V",
    TAA => "*", TAC => "Y", TAG => "*", TAT => "Y",
    TCA => "S", TCC => "S", TCG => "S", TCT => "S",
    TGA => "*", TGC => "C", TGG => "W", TGT => "C",
    TTA => "L", TTC => "F", TTG => "L", TTT => "F",

    TCN => "S", CCN => "P", ACN => "T", GTN => "V",
    CTN => "L", GCN => "A", CGN => "R", GGN => "G",

    UAA => "X", UAG => "X",

    UGA => "U"
);

%C1toC3 = ();
foreach my $threebase ( sort keys %C1 ) {
    $C1toC3{ $C1{$threebase} } = $C3{$threebase};
}

$AAcount                       = scalar keys %C1toC3;
@AAnumber{ sort keys %C1toC3 } = ( 1 .. $AAcount );
$AAnumber{'?'}                 = $AAcount + 1;
$AAnumber{'.'}                 = $AAcount + 2;

%Polar = (
    Ala => "NP", Arg => "P+", Asn => "P0", Asp => "P-",
    Cys => "P0", Gln => "P0", Glu => "P-", Gly => "P0",
    His => "P+", Ile => "NP", Leu => "NP", Lys => "P+",
    Met => "NP", Phe => "NP", Pro => "NP", Sec => "NP",
    Ser => "P0", Thr => "P0", Trp => "NP", Tyr => "P0",
    Val => "NP",

    'X' => '.', '*' => '.'
);

@Polar{ @C1{ ( sort keys %C1 ) } } = @Polar{ @C3{ ( sort keys %C1 ) } };

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

    # newly added Functions in 1.01 for complementary of splicing variants
    "SO:0001572" => 'exon_loss_variant',
    "SO:0001787" => 'splice_donor_5th_base_variant',
    "SO:0001630" => 'splice_region_variant',           # exon: 3bp + intron: 8bp
    "SO:0001995" => 'extended_intronic_splice_region_variant',    # intron: 10bp

    # the followings are replaced by 'unknown' in Voyager
    "SO:0001623" => '5_prime_UTR_variant',
    "SO:0001624" => '3_prime_UTR_variant',
    "SO:0001619" => 'nc_transcript_variant',
    "SO:0001627" => 'intron_variant',

    "annotation-fail"     => 'annotation-fail',
    "abnormal-fs-site"    => 'abnormal-fs-site',
    "abnormal-intron"     => 'abnormal-intron',
    "abnormal-inseq-stop" => "abnormal-inseq-stop",
);


%Name2SO = reverse(%SO2Name);

%func2SO = (
    "abnormal-fs-site"    => 'abnormal-fs-site',
    "abnormal-inseq-stop" => "abnormal-inseq-stop",
    "abnormal-intron"     => "abnormal-intron",
    "annotation-fail"     => 'annotation-fail',
    "cds-del"             => "SO:0001822",
    "cds-indel"           => "inframe_delins",
    "cds-ins"             => "SO:0001821",
    "cds-loss"            => "unknown-likely-deleterious",
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
    promoter              => "unknown",
    span                  => "unknown-likely-deleterious",
    unknown               => "unknown",

    # newly added in 1.01 for splicing variants complementary
    # may exists in alt_func keys, alt_funcSO, alt_funcSOname
    "exon-loss"     => "SO:0001572",
    "splice-5-5th"  => "SO:0001787",
    "splice-region" => "SO:0001630",
    "splice-ext"    => "SO:0001995",
);

%canonicalSS = (
    D => "GT",
    A => "AG",
);

our %GenePartsOrder;
@GenePartsOrder{
    (
        qw(CDS span five_prime_cis_splice_site
          three_prime_cis_splice_site ncRNA five_prime_UTR
          three_prime_UTR interior_intron five_prime_UTR_intron
          three_prime_UTR_intron abnormal-intron promoter
          annotation-fail intergenic_region), ""
    )
} = ( 1 .. 15 );

our $REF_BUILD = 'GRCh37';

=head1 Methods

=head2 new

=over

=item About : Creat a new annotation entry

=item Usage :

    my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', batch => 1 );

=item Args    - (all database files should be tabix indexed)

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

=item I<rmsk> [rmsk.bed.gz]

=over

=item add repeat tag information

=back

=item I<gwas> [gwasCatalog_snp137.bed.gz]

=over

=item add gwas information depend on only position

=back

=item I<pfam> [pfam.tsv.gz]

=over

=item add pfam information

=back

=item I<prediction> [ensembl_prediction_db.tsv.gz]

=over

=item add sift, polyphen2 prediction and scores

=back

=item I<condel> [condel config path]

=over

=item compute condel scores and prediction based on sift and polyphen2 HumVar prediction

=back

=item I<cosmic> [Cosmic_v67_241013.bed.gz]

=over

=item add cosmic annotation information

=back

=item I<phyloP> [phyloP_scores.tsv.gz]

=over

=item add phyloP scores of all 3 datasets.

=back

=item I<dbSNP> [snp137.bed.gz]

=over

=item add rsID and dbSNP frequency information

=back

=item I<tgp> [tgp_phase3_small_vars.vcf.gz]

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

=item I<exac> [ExAC.r0.2.sites.vep.vcf.gz]

=over

=item add ExAC allele frequency information from 61,486 unrelated individuals.

=back

=item I<gnomAD> [gnomad.exomes.r2.0.1.sites.vcf.gz]

=over

=item add gnomAD allele frequency information of 123,136 exomes and 15,496 genomes from unrelated individuals sequenced

=back

=item I<customdb_XX> [custom db in the same format with esp6500's]

=over

=item add customdb XX frequency infomation, XX is the ID. (multiple db will require multiple customdb option)


=back

=item I<quiet>

=over

=item Suppress warning messege to output.

=back

=item I<batch> [boolean]

=over

=item use batch mode annotation, default in daemon mode as an annotation engine.

=back

=item I<genome> [ "refgenome.fa.rz" ]

=over

=item reference genome fasta, razipped and samtools faidxed for use.

=back

=item I<genes> [ "genes.list" | $rh_geneslist ]

=over

=item annotate transcripts for I<genes>. e.g. {"ABC" => 1, "DEF" => 1} or "genes.list" 

=back

=item I<trans> [ "trans.list" | $rh_translist ]

=over

=item annotate transcripts in I<trans>. e.g. {"NM_0012.1" => 1, "NM_0034.2" => 1} or "trans.list" 

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
    my $self;
    $self = {@args};
    bless $self, ref($class) || $class;

    $self->throw("Error: please at least give 'db' and 'tr' path.")
      if ( !exists $self->{db} or !exists $self->{tr} );

    my $debugOpt = ( exists $self->{debug} ) ? 1 : 0;
    my $t0;
    if ($debugOpt) {
        $t0 = [gettimeofday];
    }

    if ( exists $self->{genes} ) {
        if ( !ref( $self->{genes} ) ) {
            open( GENE, $self->{genes} ) or $self->throw("$self->{genes} : $!");
            my %genes = map { s/\s+//g; $_ => 1 } <GENE>;
            close GENE;
            $self->{genes} = \%genes;
        }
    }
    if ( exists $self->{trans} ) {
        if ( !ref( $self->{trans} ) ) {
            open( TRAN, $self->{trans} ) or $self->throw("$self->{trans} : $!");
            my %trans = map { s/\s+//g; $_ => 1 } <TRAN>;
            close TRAN;
            $self->{trans} = \%trans;
        }
        my $rclean_trans =
          { map { s/\-\d+$//; $_ => 1 } keys %{ $self->{trans} } };
        $self->{clean_trans} = $rclean_trans;
    }

    $self->set_db( $self->{db} );

    my $t1;
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->new [db load] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $self->set_tr( $self->{tr} );
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->new [tr load] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $self->set_refbuild($REF_BUILD);

    if ( exists $self->{genome} ) {
        $self->set_genome( $self->{genome} );
    }

    if ( exists $self->{cytoBand} ) {
        $self->set_cytoBand( $self->{cytoBand} );
    }

    if ( exists $self->{pfam} ) {
        $self->set_pfam( $self->{pfam} );
    }

    if ( exists $self->{prediction} ) {
        $self->set_prediction( $self->{prediction} );
    }

    if ( exists $self->{condel} ) {
        $self->set_condel( $self->{condel} );
    }

    if ( exists $self->{phyloP} ) {
        $self->set_phyloP( $self->{phyloP} );
    }

    if ( exists $self->{cosmic} ) {
        $self->set_cosmic( $self->{cosmic} );
    }

    if ( exists $self->{dbSNP} ) {
        $self->set_dbSNP( $self->{dbSNP} );
    }

    if ( exists $self->{tgp} ) {
        $self->set_tgp( $self->{tgp} );
    }

    if ( exists $self->{esp6500} ) {
        $self->set_esp6500( $self->{esp6500} );
    }

    if ( exists $self->{exac} ) {
        $self->set_exac( $self->{exac} );
    }

    if ( exists $self->{gnomAD} ) {
        $self->set_gnomAD( $self->{gnomAD} );
    }

    foreach my $dbk ( sort keys %$self ) {
        if ( $dbk =~ /^customdb_(\S+)/ ) {
            my $dbID = $1;
            $self->set_customdb( $self->{$dbk}, $dbID );
        }
    }

    if ( exists $self->{cg54} ) {
        $self->set_cg54( $self->{cg54} );
    }

    if ( exists $self->{wellderly} ) {
        $self->set_wellderly( $self->{wellderly} );
    }

    if ( exists $self->{rmsk} ) {
        $self->set_rmsk( $self->{rmsk} );
    }

    if ( exists $self->{gwas} ) {
        $self->set_gwas( $self->{gwas} );
    }

    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->new [others load] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    return $self;
}

=head2 set/get methods for properties

    List of Properties:
                        get          set
    db                  o            o
    tr                  o            o
    refbuild            o            o
    tidb                o            x
    annodb              o            x
    trInfodb            o            x
    cytoBand            o            o
    cytoBand_h          o            x
    rmsk                o            o
    rmsk_h              o            x
    gwas                o            o
    gwas_h              o            x
    pfam                o            o
    pfam_h              o            x
    prediction          o            o
    prediction_h        o            x
    phyloP              o            o
    phyloP_h            o            x
    cosmic              o            o
    cosmic_h            o            x
    dbSNP               o            o
    dbSNP_h             o            x
    tgp                 o            o
    tgp_h               o            x
    esp6500             o            o
    esp6500_h           o            x
    exac                o            o
    gnomAD              o            o
    exac_h              o            x
    cg54                o            o
    cg54_h              o            x
    wellderly           o            o
    wellderly_h         o            x

    e.g.    : $beda->set_refbuild($refbuild);
              my $refbuild = $beda->get_refbuild();

=cut

sub set_refbuild {
    my $self              = shift;
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
    my $db   = shift;
    if ( !-e $db or !-r $db ) {
        $self->throw("Error: cannot read $db.");
    }

    $self->{db} = $db;
    if ( !-e $db . ".tbi" or !-r $db . ".tbi" ) {
        $self->throw(
            "Error: [$db] index (.tbi) file not found, please build it first.");
    }

    $self->{tidb} = Tabix->new( -data => $db );

    return $self if ( !exists $self->{batch} );

    my %open_args;
    if ( exists $self->{region} ) {
        $open_args{region} = $self->{region};
    }
    if ( !exists $self->{region} and exists $self->{regbed} ) {
        $open_args{regbed} = $self->{regbed};
    }
    if ( exists $self->{genes} ) {
        $open_args{genes} = $self->{genes};
    }
    if ( exists $self->{trans} and exists $self->{clean_trans} ) {
        $open_args{trans}       = $self->{trans};
        $open_args{clean_trans} = $self->{clean_trans};
    }
    if ( exists $self->{mmap} ) {
        $open_args{mmap} = $self->{mmap};
    }

    $self->{annodb} = $self->load_anno(%open_args);

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
    if ( exists $self->{annodb} ) {
        return $self->{annodb};
    }
    else {
        return undef;
    }
}

sub set_tr {
    my $self = shift;
    my $tr   = shift;
    if ( !-e $tr or !-r $tr ) {
        $self->throw("Error: cannot read $tr.");
    }
    $self->{tr} = $tr;
    my %load_opts = ();
    $load_opts{genes} = $self->{genes} if ( exists $self->{genes} );
    $load_opts{trans} = $self->{trans} if ( exists $self->{trans} );
    $self->{trInfodb} = $self->readtr(%load_opts);
    return $self;
}

sub set_genome {
    my $self      = shift;
    my $genome_rz = shift;
    if ( !-e $genome_rz or !-r $genome_rz ) {
        $self->throw("Error: cannot read $genome_rz.");
    }
    require Faidx if ( !exists $self->{genome_h} );
    $self->{genome}   = $genome_rz;
    $self->{genome_h} = Faidx->new($genome_rz);
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
    my $self   = shift;
    my $cytodb = shift;
    $self->{cytoBand} = $cytodb;
    require GetCytoBand if ( !exists $self->{cytoBand_h} );
    my $cytoBand_h = GetCytoBand->new( db => $cytodb );
    $self->{cytoBand_h} = $cytoBand_h;
    return $self;
}

sub get_cytoBand {
    my $self = shift;
    return $self->{cytoBand} if ( exists $self->{cytoBand} );
    return undef;
}

sub get_cytoBand_h {
    my $self = shift;
    return $self->{cytoBand_h} if ( exists $self->{cytoBand_h} );
    return undef;
}

sub set_rmsk {
    my $self   = shift;
    my $rmskdb = shift;
    $self->{rmsk} = $rmskdb if ( defined $rmskdb );
    require GetRepeatTag if ( !exists $self->{rmsk_h} );
    my $rmsk_h = GetRepeatTag->new( db => $self->{rmsk} );
    $self->{rmsk_h} = $rmsk_h;
    return $self;
}

sub get_rmsk {
    my $self = shift;
    return $self->{rmsk} if ( exists $self->{rmsk} );
    return undef;
}

sub get_rmsk_h {
    my $self = shift;
    return $self->{rmsk_h} if ( exists $self->{rmsk_h} );
    return undef;
}

sub set_gwas {
    my $self   = shift;
    my $gwasdb = shift;
    $self->{gwas} = $gwasdb if ( defined $gwasdb );
    require GetGWAS if ( !exists $self->{gwas_h} );
    my $gwas_h = GetGWAS->new( db => $self->{gwas} );
    $self->{gwas_h} = $gwas_h;
    return $self;
}

sub get_gwas {
    my $self = shift;
    return $self->{gwas} if ( exists $self->{gwas} );
    return undef;
}

sub get_gwas_h {
    my $self = shift;
    return $self->{gwas_h} if ( exists $self->{gwas_h} );
    return undef;
}

sub set_pfam {
    my $self   = shift;
    my $pfamdb = shift;
    $self->{pfam} = $pfamdb;
    require GetPfam if ( !exists $self->{pfam_h} );
    my $pfam_h = GetPfam->new( db => $pfamdb );
    $self->{pfam_h} = $pfam_h;
    return $self;
}

sub get_pfam {
    my $self = shift;
    return $self->{pfam} if ( exists $self->{pfam} );
    return undef;
}

sub get_pfam_h {
    my $self = shift;
    return $self->{pfam_h} if ( exists $self->{pfam_h} );
    return undef;
}

sub set_prediction {
    my $self         = shift;
    my $predictiondb = shift;
    $self->{prediction} = $predictiondb;
    require GetPrediction if ( !exists $self->{prediction_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $prediction_h = GetPrediction->new( db => $predictiondb, %common_opts );
    $self->{prediction_h} = $prediction_h;
    return $self;
}

sub get_prediction {
    my $self = shift;
    return $self->{prediction} if ( exists $self->{prediction} );
    return undef;
}

sub get_prediction_h {
    my $self = shift;
    return $self->{prediction_h} if ( exists $self->{prediction_h} );
    return undef;
}

sub set_condel {
    my $self         = shift;
    my $condelConfig = shift;
    $self->{condel} = $condelConfig;
    require CondelPred if ( !exists $self->{condel_h} );
    my $condel_h = CondelPred->new($condelConfig);
    $self->{condel_h} = $condel_h;
    return $self;
}

sub get_condel {
    my $self = shift;
    return $self->{condel} if ( exists $self->{condel} );
    return undef;
}

sub get_condel_h {
    my $self = shift;
    return $self->{condel_h} if ( exists $self->{condel_h} );
    return undef;
}

sub set_phyloP {
    my $self     = shift;
    my $phyloPdb = shift;
    $self->{phyloP} = $phyloPdb;
    require GetPhyloP46wayScore if ( !exists $self->{phyloP_h} );
    my $phyloP_h = GetPhyloP46wayScore->new( db => $phyloPdb );
    $self->{phyloP_h} = $phyloP_h;
    return $self;
}

sub get_phyloP {
    my $self = shift;
    return $self->{phyloP} if ( exists $self->{phyloP} );
    return undef;
}

sub get_phyloP_h {
    my $self = shift;
    return $self->{phyloP_h} if ( exists $self->{phyloP_h} );
    return undef;
}

sub set_dbSNP {
    my $self    = shift;
    my $dbSNPdb = shift;
    $self->{dbSNP} = $dbSNPdb;
    require GetDBSNP if ( !exists $self->{dbSNP_h} );
    my $dbSNP_h = GetDBSNP->new( db => $dbSNPdb );
    $self->{dbSNP_h} = $dbSNP_h;
    return $self;
}

sub get_dbSNP {
    my $self = shift;
    return $self->{dbSNP} if ( exists $self->{dbSNP} );
    return undef;
}

sub get_dbSNP_h {
    my $self = shift;
    return $self->{dbSNP_h} if ( exists $self->{dbSNP_h} );
    return undef;
}

sub get_cosmic {
    my $self = shift;
    return $self->{cosmic} if ( exists $self->{cosmic} );
    return undef;
}

sub set_cosmic {
    my $self      = shift;
    my $cosmic_db = shift;
    $self->{cosmic} = $cosmic_db;
    require GetCOSMIC if ( !exists $self->{cosmic_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $cosmic_h = GetCOSMIC->new( db => $cosmic_db, %common_opts );
    $self->{cosmic_h} = $cosmic_h;
    return $self;
}

sub get_cosmic_h {
    my $self = shift;
    return $self->{cosmic_h} if ( exists $self->{cosmic_h} );
    return undef;
}

sub set_tgp {
    my $self  = shift;
    my $tgpdb = shift;
    $self->{tgp} = $tgpdb;
    require GetTGP if ( !exists $self->{tgp_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $tgp_h = GetTGP->new( db => $tgpdb, %common_opts );
    $self->{tgp_h} = $tgp_h;
    return $self;
}

sub get_tgp {
    my $self = shift;
    return $self->{tgp} if ( exists $self->{tgp} );
    return undef;
}

sub get_tgp_h {
    my $self = shift;
    return $self->{tgp_h} if ( exists $self->{tgp_h} );
    return undef;
}

sub set_esp6500 {
    my $self      = shift;
    my $esp6500db = shift;
    $self->{esp6500} = $esp6500db;
    require GetVcfAF if ( !$load_opt_vcfaf );
    $load_opt_vcfaf = 1;
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $esp6500_h = GetVcfAF->new( db => $esp6500db, %common_opts );
    $self->{esp6500_h} = $esp6500_h;
    return $self;
}

sub get_esp6500 {
    my $self = shift;
    return $self->{esp6500} if ( exists $self->{esp6500} );
    return undef;
}

sub set_exac {
    my $self    = shift;
    my $exac_db = shift;
    $self->{exac} = $exac_db;
    require GetExAC if ( !exists $self->{exac_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $exac_h = GetExAC->new( db => $exac_db, %common_opts );
    $self->{exac_h} = $exac_h;
    return $self;
}

sub get_exac {
    my $self = shift;
    return $self->{exac} if ( exists $self->{exac} );
    return undef;
}

sub get_exac_h {
    my $self = shift;
    return $self->{exac_h} if ( exists $self->{exac_h} );
    return undef;
}

sub set_gnomAD {
    my $self    = shift;
    my $gnomAD_db = shift;
    $self->{gnomAD} = $gnomAD_db;
    require GetGAD if ( !exists $self->{gnomAD_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $gnomAD_h = GetGAD->new( db => $gnomAD_db, %common_opts );
    $self->{gnomAD_h} = $gnomAD_h;
    return $self;
}

sub get_gnomAD {
    my $self = shift;
    return $self->{gnomAD} if ( exists $self->{gnomAD} );
    return undef;
}

sub get_gnomAD_h {
    my $self = shift;
    return $self->{gnomAD_h} if ( exists $self->{gnomAD_h} );
    return undef;
}

sub set_customdb {
    my $self = shift;
    my ( $cusdb, $dbID ) = @_;
    require GetVcfAF if ( !$load_opt_vcfaf );
    $load_opt_vcfaf = 1;
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $cusdb_h = GetVcfAF->new( db => $cusdb, %common_opts );
    $self->{ "cusdb_" . $dbID . "_h" } = $cusdb_h;
    return $self;
}

sub get_esp6500_h {
    my $self = shift;
    return $self->{esp6500_h} if ( exists $self->{esp6500_h} );
    return undef;
}

sub set_cg54 {
    my $self   = shift;
    my $cg54db = shift;
    $self->{cg54} = $cg54db;
    require GetCGpub
      if ( !exists $self->{cg54_h} and !exists $self->{wellderly_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $cg54_h = GetCGpub->new( db => $cg54db, %common_opts );
    $self->{cg54_h} = $cg54_h;
    return $self;
}

sub get_cg54 {
    my $self = shift;
    return $self->{cg54} if ( exists $self->{cg54} );
    return undef;
}

sub get_cg54_h {
    my $self = shift;
    return $self->{cg54_h} if ( exists $self->{cg54_h} );
    return undef;
}

sub set_wellderly {
    my $self        = shift;
    my $wellderlydb = shift;
    $self->{wellderly} = $wellderlydb;
    require GetCGpub
      if ( !exists $self->{cg54_h} and !exists $self->{wellderly_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    my $wellderly_h = GetCGpub->new( db => $wellderlydb, %common_opts );
    $self->{wellderly_h} = $wellderly_h;
    return $self;
}

sub get_wellderly {
    my $self = shift;
    return $self->{wellderly} if ( exists $self->{wellderly} );
    return undef;
}

sub get_wellderly_h {
    my $self = shift;
    return $self->{wellderly_h} if ( exists $self->{wellderly_h} );
    return undef;
}

sub TO_JSON {
    return { %{ shift() } };
}

sub DESTROY {
    my $self = shift;
    if ( exists $self->{tidb} and defined $self->{tidb} ) {
        $self->{tidb}->DESTROY() if ( $self->{tidb}->can('DESTROY') );
        delete $self->{tidb};
    }

    foreach my $dbhk ( sort keys %$self ) {
        if ( $dbhk =~ /\S+_h$/ and defined $self->{$dbhk} ) {
            $self->{$dbhk}->DESTROY() if ( $self->{$dbhk}->can('DESTROY') );
            delete $self->{$dbhk};
        }
    }

    return;
}

=head2 write_using

    About   : write the current using database information to files
    Usage   : $beda->write_using( $file, $type );
    Args    : file gives the filename to output, and type is one of the following:
		g:  gene symbol list
		t:  transcript acc.ver list
		c:  the complete annotation region, in bed format,
		    which can be used as the variation calling region.
		b:  standard bed format of only exon region
		a:  BED format, exon region with '+1' annotation, 
		    oneline per one exon for one transcript, this
		    format allow redundancy exists, and not sorted by
		    chromosomal coordinates, but by the transcript acc.ver
		    and the exon number, this file is mainly for 
		    transcript originated statistics for only exon.

=cut

sub write_using {
    my ( $self, $file, $type ) = @_;
    open( F, ">", $file ) or $self->throw("$file: $!");

    my $annodb_load;
    if ( exists $self->{annodb} ) {
        $annodb_load = $self->{annodb};
    }
    else {
        my %open_args;
        if ( exists $self->{region} ) {
            $open_args{region} = $self->{region};
        }
        if ( !exists $self->{region} and exists $self->{regbed} ) {
            $open_args{regbed} = $self->{regbed};
        }
        if ( exists $self->{genes} ) {
            $open_args{genes} = $self->{genes};
        }
        if ( exists $self->{trans} and exists $self->{clean_trans} ) {
            $open_args{trans}       = $self->{trans};
            $open_args{clean_trans} = $self->{clean_trans};
        }
        if ( exists $self->{mmap} ) {
            $open_args{mmap} = $self->{mmap};
        }
        $annodb_load = $self->load_anno(%open_args);
    }

    my ( %genes, %trans, @beds, %anno, @complete ) = ();
    foreach my $chr ( sort keys %$annodb_load ) {
        foreach my $bedent ( @{ $annodb_load->{$chr} } ) {
            if ( $type eq 'c' ) {
                push( @complete, [ $chr, $$bedent{sta}, $$bedent{sto} ] );
                next;
            }
            my $exOpt = 0;
            foreach my $annoblk ( keys %{ $$bedent{annos} } ) {

# tid         gsym  gid strd blka gpgo exin nsta nsto csta csto wlen mismatch    pr
# NM_015658.3|NOC2L|26155|-|3U1E|205|EX19E|2817|2801|*508|*492|0|D,879582,879582,.|Y
                my @info = split( /\|/, $annoblk );

                if ( $type eq 'g' ) {
                    $genes{ $info[1] } = 1;
                }
                elsif ( $type eq 't' ) {
                    $trans{ $info[0] } = 1;
                }
                elsif ( $type eq 'b' ) {
                    $exOpt = 1 if ( $info[6] =~ /^EX/ );
                }
                elsif ( $type eq 'a' ) {
                    if ( $info[6] =~ /^EX/ ) {
                        if (   !exists( $anno{ $info[0] } )
                            or !exists( $anno{ $info[0] }{ $info[6] } ) )
                        {
                            $anno{ $info[0] }{ $info[6] }{chr} = $chr;
                            $anno{ $info[0] }{ $info[6] }{sta} = $$bedent{sta};
                            $anno{ $info[0] }{ $info[6] }{sto} = $$bedent{sto};
                            $anno{ $info[0] }{ $info[6] }{blk} = $annoblk;
                        }
                        elsif (
                            $$bedent{sto} > $anno{ $info[0] }{ $info[6] }{sto} )
                        {
                            $anno{ $info[0] }{ $info[6] }{sto} = $$bedent{sto};
                        }
                        else {
                            next;
                        }
                    }
                }
                else { $self->throw("type: $type not regconized."); }
            }
            push( @beds, [ $chr, $$bedent{sta}, $$bedent{sto} ] )
              if ( $type eq 'b' and $exOpt );
        }
    }

    if ( $type eq 'g' ) {
        say F join( "\n", ( sort keys %genes ) );
    }
    elsif ( $type eq 't' ) {
        say F join( "\n", ( sort keys %trans ) );
    }
    elsif ( $type eq 'b' ) {    # merge bed regions due to pre-sorted db.
        if ( 0 == @beds ) {
            $self->warn("no region in db.");
            return;
        }
        my ( $pre_chr, $pre_sta, $pre_sto ) = @{ $beds[0] };
        for ( my $i = 1 ; $i < @beds ; $i++ ) {
            my ( $cur_chr, $cur_sta, $cur_sto ) = @{ $beds[$i] };
            if ( ( $cur_chr ne $pre_chr ) or ( $cur_sta > $pre_sto ) ) {
                say F join( "\t", $pre_chr, $pre_sta, $pre_sto );
                ( $pre_chr, $pre_sta, $pre_sto ) =
                  ( $cur_chr, $cur_sta, $cur_sto );
            }
            elsif ( $cur_sta == $pre_sto ) {
                $pre_sto = $cur_sto;
            }
            else {
                $self->throw("Error: bad db, non-departed beds, or no-sort!");
            }
        }
        say F join( "\t", $pre_chr, $pre_sta, $pre_sto );
    }
    elsif ( $type eq 'a' ) {
        foreach my $nm ( sort keys %anno ) {
            foreach my $ex ( sort exsort keys %{ $anno{$nm} } ) {
                say F join( "\t",
                    $anno{$nm}{$ex}{chr}, $anno{$nm}{$ex}{sta},
                    $anno{$nm}{$ex}{sto}, $anno{$nm}{$ex}{blk} );
            }
        }
    }
    elsif ( $type eq 'c' ) {
        if ( 0 == @complete ) {
            $self->warn("no region in db.");
            return;
        }
        my ( $pre_chr, $pre_sta, $pre_sto ) = @{ $complete[0] };
        for ( my $i = 1 ; $i < @complete ; $i++ ) {
            my ( $cur_chr, $cur_sta, $cur_sto ) = @{ $complete[$i] };
            if ( ( $cur_chr ne $pre_chr ) or ( $cur_sta > $pre_sto ) ) {
                say F join( "\t", $pre_chr, $pre_sta, $pre_sto );
                ( $pre_chr, $pre_sta, $pre_sto ) =
                  ( $cur_chr, $cur_sta, $cur_sto );
            }
            elsif ( $cur_sta == $pre_sto ) {
                $pre_sto = $cur_sto;
            }
            else {
                $self->throw(
                    "Error: bad db, non-departed complete, or no-sort!");
            }
        }
        say F join( "\t", $pre_chr, $pre_sta, $pre_sto );
    }
    else { $self->throw("type: $type not regconized."); }

    close F;
    return;
}

sub exsort {
    my ( $sym, $anum, $bnum );
    if ( $a =~ /^EX([\+\-\*]?)(\d+)[EP]?$/ ) {
        $sym  = $1;
        $anum = $2;
    }
    if ( $b =~ /^EX[\+\-\*]?(\d+)[EP]?$/ ) {
        $bnum = $1;
    }
    confess "ExIn number format error. [$a, $b]"
      if ( !defined $anum or !defined $bnum );
    if ( !defined $sym or $sym !~ /\-/ ) {
        $anum <=> $bnum;
    }
    else {
        $bnum <=> $anum;
    }
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
		    # start is 1 based, to coordinate with region string format
		}
	      Note: the pos-pair which do not hit any annotation blocks, will
		    not exist in the returned results.

=cut

sub get_cover_batch {
    my ( $self, $chr, $stasto_aref ) = @_;

    $chr =~ s/^chr//i;
    if ( $chr =~ /^M/i ) {
        $chr = 'MT';
    }
    my $rAnnos;
    if ( !exists $self->{annodb} ) {
        my %open_args;
        if ( exists $self->{region} ) {
            $open_args{region} = $self->{region};
        }
        if ( !exists $self->{region} and exists $self->{regbed} ) {
            $open_args{regbed} = $self->{regbed};
        }
        if ( exists $self->{genes} ) {
            $open_args{genes} = $self->{genes};
        }
        if ( exists $self->{trans} and exists $self->{clean_trans} ) {
            $open_args{trans}       = $self->{trans};
            $open_args{clean_trans} = $self->{clean_trans};
        }
        if ( exists $self->{mmap} ) {
            $open_args{mmap} = $self->{mmap};
        }
        my $rwhole = $self->load_anno(%open_args);
        $rAnnos = $rwhole->{$chr} if ( exists $rwhole->{$chr} );
    }
    else {
        $rAnnos = $self->{annodb}->{$chr} if ( exists $self->{annodb}->{$chr} );
    }

    my @sorted_stasto = sort pairsort @$stasto_aref;
    return {} if ( 0 == @sorted_stasto );
    if ( !defined $rAnnos ) {
        $self->warn("Warning: no available annotation items in curdb for $chr")
          if ( !exists $self->{quiet} );
        return {};
    }

    my %ret_cov = ();

    my $cur_blkId = 0;
    foreach my $rgn (@sorted_stasto) {

        my $pseudo_var =
          BedAnno::Var->new( $chr, ( $rgn->[0] - 1 ), $rgn->[1], "=", "?" );
        my $anno_var = BedAnno::Anno->new($pseudo_var);

        $cur_blkId = $anno_var->getTrPosition( $rAnnos, $cur_blkId );

        my $pospair = join( "-", @$rgn );
        $ret_cov{$pospair} = $self->get_hitted_blk($anno_var);

    }
    return \%ret_cov;
}

sub pairsort {
    $$a[0] <=> $$b[0] or $$a[1] <=> $$b[1];
}

sub get_hitted_blk {
    my $self        = shift;
    my $anno_var    = shift;
    my @hitted_blks = ();

    # hit the annotation blks
    if ( exists $anno_var->{trInfo} ) {
        foreach my $tid ( sort keys %{ $anno_var->{trInfo} } ) {

            if (    $anno_var->{trInfo}->{$tid}->{rnaBegin} ne '?'
                and $anno_var->{trInfo}->{$tid}->{rnaEnd} ne '?' )
            {
                my ( $begin_hash, $end_hash );
                $begin_hash = {
                    gsym => $anno_var->{trInfo}->{$tid}->{geneSym},
                    gid  => $anno_var->{trInfo}->{$tid}->{geneId},
                    reg  => $anno_var->{trInfo}->{$tid}->{r_Begin},
                    exin => $anno_var->{trInfo}->{$tid}->{ei_Begin},
                    strd => $anno_var->{trInfo}->{$tid}->{strd},
                    cpos => 'n.' . $anno_var->{trInfo}->{$tid}->{rnaBegin}
                };
                $end_hash = {
                    gsym => $anno_var->{trInfo}->{$tid}->{geneSym},
                    gid  => $anno_var->{trInfo}->{$tid}->{geneId},
                    reg  => $anno_var->{trInfo}->{$tid}->{r_End},
                    exin => $anno_var->{trInfo}->{$tid}->{ei_End},
                    strd => $anno_var->{trInfo}->{$tid}->{strd},
                    cpos => 'n.' . $anno_var->{trInfo}->{$tid}->{rnaEnd}
                };

                if ( exists $anno_var->{trInfo}->{$tid}->{cdsBegin}
                    and $anno_var->{trInfo}->{$tid}->{cdsBegin} ne '' )
                {
                    $begin_hash->{cpos} =
                      'c.' . $anno_var->{trInfo}->{$tid}->{cdsBegin};
                    $end_hash->{cpos} =
                      'c.' . $anno_var->{trInfo}->{$tid}->{cdsEnd};
                }

                if ( $begin_hash->{strd} eq '+' ) {
                    push( @hitted_blks, [ $tid, $begin_hash, $end_hash ] );
                }
                else {
                    push( @hitted_blks, [ $tid, $end_hash, $begin_hash ] );
                }
            }
        }
    }
    return \@hitted_blks;
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
                    csta     => $cds_start_on_trSeq, # 0 based
                    csto     => $cds_end_on_trSeq,   # 1 based
                    seq      => $tr_sequence,
                    
		    nfs      => { $fslead => $fsbase, ... },
		    cfs      => { $cfslead => $fsbase, ... },
                    X        => 1,                   # inseqStop
                    U        => 1,                   # selenocysteine
                    A        => 1,                   # polyATail
                    altstart => {                    # altstart codons
                        $startCodons1 => 1,
                        $startCodons2 => 1,
                        ...
                    },

                    # the following two keys will be assigned
                    # to mRNA when needed
                    cseq     => $codonSequence,
                    pseq     => $proteinSequence,
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

    my $fas_h;
    if ($self->{tr} =~ /\.gz$/) {
        $fas_h = new IO::Uncompress::Gunzip $self->{tr}, AUTOCLOSE => 1, MultiStream => 1
            or confess "Error: [$self->{tr}] $GunzipError\n";
    }
    else {
        open($fas_h, $self->{tr}) or confess "Error: [$self->{tr}] $!";
    }
    local $/ = ">";
    my %seqs = ();
    while (<$fas_h>) {
        s/[\s>]+$//g;
        next if (/^\s*$/);
        my $hd = $1 if (s/^(\S+[^\n]*)\n//);
        confess "Error: trSeq parse error!" if ( !defined $hd );
        my @headers = split( /\s+/, $hd, 10 );
        confess "Error: trSeq header parse error!" if ( 7 > @headers );
        s/\s+//g;

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
        $seqs{ $headers[0] }{seq}  = $_;
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

        my @fss = ();
        if ( $headers[6] =~ /altstart/ ) {
            $seqs{ $headers[0] }{altstart} =
              { map { $_ => 1 } split( /;/, $headers[7] ) };
        }
        elsif ( defined $headers[7] ) {              # frameshift at 8th item
            @fss = split( /\|/, $headers[7] );
        }

        if ( defined $headers[8] ) {                 # frameshift at 9th item
            @fss = split( /\|/, $headers[8] );
        }

        if ( 0 < scalar @fss ) {
            confess "Error: [$headers[0]] no cds but with fs!"
              if ( $headers[5] eq "." );
            $seqs{ $headers[0] }{nfs} =
              { map { my @fsl = split(/;/); $fsl[0] => $fsl[1] } @fss };

            $seqs{ $headers[0] }{cfs} = {
                map {
                    my @cfsl = split(/;/);
                    ( $cfsl[0] - $seqs{ $headers[0] }{csta} ) => $cfsl[1]
                } @fss
            };
        }
    }
    close($fas_h);
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
        my @regions = split( /\s+/, $args{region} );
        foreach my $reg (@regions) {
            next if ( $reg eq "" );
            if ( $reg =~ /^(\S+):(\-?\d+)\-(\d+)$/ ) {
                my ( $name, $beg, $end ) = ( $1, $2, $3 );
                if ( $beg <= 0 ) {
                    $self->warn(
"Warning: region string should be 1 based [$reg], has been changed to 1 based"
                    ) if ( exists $self->{debug} );
                    $beg = 1;
                }
                $name =~ s/^chr//i;
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
        open( BED, $args{regbed} ) or $self->throw("Error: [$args{regbed}] $!");
        while (<BED>) {
            chomp;
            my @beditm = split(/\t/);
            $self->throw("Error: bed format error.") if ( 3 > @beditm );
            $beditm[0] =~ s/^chr//i;
            if ( $beditm[0] =~ /^M/i ) {
                $beditm[0] = "MT";
            }
            push( @query_region, [ @beditm[ 0, 1, 2 ] ] );
        }
        close BED;
    }

    my $read_all_opt = 0;
    my @all_querys   = ();
    my $annodb_h;
    if ( 0 == @query_region ) {
        $read_all_opt = 1;
        $annodb_h = new IO::Uncompress::Gunzip $self->{db}, AUTOCLOSE => 1, MultiStream => 1
            or confess "Error: [$self->{db}] $GunzipError\n";
    }
    else {
        my @sorted_regions = sort {
                 $a->[0] cmp $b->[0]
              or $a->[1] <=> $b->[1]
              or $a->[2] <=> $b->[2]
        } @query_region;

        my ( $cname, $cbeg, $cend ) = @{ $sorted_regions[0] };
        for ( my $k = 1 ; $k < @sorted_regions ; $k++ ) {
            my ( $dname, $dbeg, $dend ) = @{ $sorted_regions[$k] };
            if ( $dname ne $cname or $dbeg > $cend ) {

                my $query_ent = $self->{tidb}->query( $cname, $cbeg, $cend );
                if ( defined $query_ent->{_} ) {
                    push( @all_querys, $query_ent );
                }
                ( $cname, $cbeg, $cend ) = @{ $sorted_regions[$k] };
            }
            else {
                $cend = $dend;
            }
        }
        my $q = $self->{tidb}->query( $cname, $cbeg, $cend );
        if ( defined $q->{_} ) {
            push( @all_querys, $q );
        }
    }

    # trans filter is always be ahead of genes
    my $mmapTag = ( exists $args{mmap} )  ? 1 : 0;
    my $geneTag = ( exists $args{genes} ) ? 1 : 0;
    my $tranTag = ( exists $args{trans} ) ? 1 : 0;
    my $geneList = $args{genes} if ($geneTag);
    my $tranList = $args{trans} if ($tranTag);
    my %pureTran = map { s/\-\d+$//; $_ => 1 } keys %$tranList if ($tranTag);

    my $rannodb = {};
    while (1) {
        my $tb_ent;
        if ($read_all_opt) {
            $tb_ent = <$annodb_h>;
            last if ( !defined $tb_ent or $tb_ent eq "" );
        }
        else {
            while ( 0 < @all_querys ) {
                my $region_read = $self->{tidb}->read( $all_querys[0] );
                if ( !defined $region_read or $region_read eq "" ) {
                    shift(@all_querys);
                }
                else {
                    $tb_ent = $region_read;
                    last;
                }
            }
            last if ( 0 == @all_querys );
        }

        $tb_ent =~ s/\s+$//;
        my ( $chr, $start, $stop, $annostr ) = split( /\t/, $tb_ent );
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
    if ($read_all_opt) {
        close($annodb_h);
    }

    $rannodb = region_merge($rannodb)
      if ( $mmapTag or $geneTag or $tranTag );

    return $rannodb;
}

=head2 region_merge

    About   : merge consecutive same-entries regions
    Usage   : my $rannodb = region_merge($loaded_db);
    Args    : A hash ref of loaded_db.
    Returns : A hash ref of merged db.

=cut

sub region_merge {

    my $radb       = shift;
    my %local_radb = %$radb;

    foreach my $chr ( keys %local_radb ) {
        my @annoents = @{ $local_radb{$chr} };
        my $oricount = scalar(@annoents);
        next if ( $oricount <= 1 );

        # merge from bottom to top
        my $curid = $oricount - 1;
      ANNOENT: while ( $curid > 0 ) {
            my $latter = $annoents[$curid];
            my $former = $annoents[ $curid - 1 ];
            $curid--;

            next if ( $$latter{sta} != $$former{sto} );
            my @former_ann = keys %{ $$former{annos} };
            my $formern    = scalar @former_ann;
            next if ( $formern != ( scalar keys %{ $$latter{annos} } ) );
            foreach my $ann (@former_ann) {
                next ANNOENT if ( !exists $$latter{annos}{$ann} );
            }

            # splice the latter one and correct the former one
            $$former{sto} = $$latter{sto};
            splice( @annoents, ( $curid + 1 ), 1 );
        }

        $local_radb{$chr} = [@annoents] if ( $oricount > ( scalar @annoents ) );
    }

    return \%local_radb;
}

# involke parse_annoent to assign detail information to annodb.
sub assign_detail {
    my $self      = shift;
    my $rannodb_k = shift;
    my %detail    = ();
    foreach my $annoblk ( sort keys %{ $$rannodb_k{annos} } ) {
        my $offset = $$rannodb_k{annos}{$annoblk};
        my ( $tid, $ranno ) = parse_annoent($annoblk);
        $detail{$tid} = $ranno;
        $detail{$tid}{offset} = $offset;
    }
    $rannodb_k->{detail} = {%detail};
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
    $annoinfo{gpSO} = 'abnormal-intron'
      if ( $annoinfo{gpSO} eq 'abnormal_intron' );
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

    my $debugOpt = ( exists $self->{debug} ) ? 1 : 0;
    my $t0;
    if ($debugOpt) {
        $t0 = [gettimeofday];
    }
    my $var = BedAnno::Var->new(@_);

    my $t1;
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->anno [new BedAnno::Var] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }
    my ( $annoEnt, $idx ) = $self->varanno($var);

    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->anno [varanno] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }
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

                    varName => $var_mutation_name,

                    # information
                    varTypeSO => $varTypeSO,
                    gHGVS     => $gHGVS,
                    refbuild  => $referenceBuild,

                    # Here's some optional parts which may be generated
                    # when extra resource is available:

                    cytoBand  => $cytoBand,
                    reptag    => $repeatTag,
                    gwas      => $ref_gwas_ret,

                    # For single position for now
                    phyloPpm    => $PhyloPscorePlacentalMammals,
                    phyloPpr    => $PhyloPscorePrimates,
                    phyloPve    => $PhyloPscoreVetebrates,

                    reptag	=> $RepeatMaskerTag,

                    dbsnp => {
                        $rsID => {
                            AN => $dbsnp_total_allele_count,
                            AF => $dbsnp_alt_allele_frequency,  # though ref
                        },
                        ...
                    },

                    cosmic => $ref_cosmic_return,

                    tgp => $ref_tgp_return,
                    exac => $ref_exac_return,

                    esp6500 => {
                        AN => $esp6500_total_allele_count,
                        AF => $esp6500_alt_allele_frequency,
                    },

                    cusdb_XX => {
                    AN => $custom_db_allele_count,
                    AF => $custom_db_allele_frequency,
                    },
                    ...

                },
                trInfo => {
                    $tid => {
                        trVarName     => $transcriptVariantName,
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
                        p3	          => $threeletter_pHGVS,
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
                        componentIndex => $componentIndex,
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

                        # for some of splice variants, there may exists 
                        # the following function information
                        alt_func        => $alternative_func_code,
                        alt_funcSO      => $alternative_variant_SO_id,
                        alt_funcSOname  => $alt_variant_SO_name,

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
                        condelPred  => $Condelpred,
                        condelScore => $Condelscore,

                    },
                    ...
                }
            }

=cut

sub varanno {
    my ( $self, $var, $AEIndex ) = @_;
    $AEIndex ||= 0;

    my $debugOpt = ( exists $self->{debug} ) ? 1 : 0;
    my $t0;
    if ($debugOpt) {
        $t0 = [gettimeofday];
    }

    my ($unified_start, $unified_ref, $unified_alt, $unified_rl, $no_use);
    if ($var->{guess} eq 'ref') {
        ($unified_start, $unified_ref, $unified_alt, $unified_rl, $no_use) = @$var{qw(pos ref alt reflen altlen)};
    }
    else {
        ($unified_start, $unified_ref, $unified_alt, $unified_rl, $no_use) = $var->getUnifiedVar('-', 1);
    }
    my $unified_end = $unified_start + $unified_rl;

    if ( exists $self->{cytoBand} ) {
        $var->{cytoBand} = $self->{cytoBand_h}->getCB( $var->{chr}, $unified_start, $unified_end );
    }

    if ( exists $self->{rmsk} ) {
        $var->{reptag} = $self->{rmsk_h}->getRepTag( $var->{chr}, $unified_start, $unified_end );
    }

    if ( exists $self->{gwas} ) {
        $var->{gwas} = $self->{gwas_h}->getGWAS( $var->{chr}, $unified_start, $unified_end );
    }

    if ( exists $self->{phyloP} ) {
        if ( $var->{sm} == 1 ) {
            @$var{qw(phyloPpm phyloPpr phyloPve)} =
              @{ $self->{phyloP_h}
                  ->getPhyloP46wayScore( $var->{chr}, ( $var->{pos} + 1 ) ) };
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
              $self->{dbSNP_h}->getRS( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
        }
    }

    if ( exists $self->{tgp} ) {
        $var->{tgp} = $self->{tgp_h}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    if ( exists $self->{esp6500} ) {
        $var->{esp6500} =
          $self->{esp6500_h}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    if ( exists $self->{exac} ) {
        $var->{exac} = $self->{exac_h}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    if ( exists $self->{gnomAD} ) {
        $var->{gnomAD} = $self->{gnomAD_h}->getGAD( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    foreach my $dbhk ( sort keys %$self ) {
        if ( $dbhk =~ /^cusdb_(\S+)_h/ and defined $self->{$dbhk} ) {
            my $dbID = $1;
            $var->{"cusdb_$dbID"} =
              $self->{$dbhk}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
        }
    }

    if ( exists $self->{cg54} ) {
        $var->{cg54} = $self->{cg54_h}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    if ( exists $self->{wellderly} ) {
        $var->{wellderly} =
          $self->{wellderly_h}->getAF( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    if ( exists $self->{cosmic} ) {
        $var->{cosmic} =
          $self->{cosmic_h}->getCOSMIC( $var->{chr}, $unified_start, $unified_end, $unified_ref, $unified_alt );
    }

    my $t1;
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [var extradb sql] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $var->{varTypeSO} = $Name2SO{ $var->{guess} };
    $var->{refbuild}  = $self->{refbuild};
    $var->get_gHGVS();
    if ( $var->{gHGVS} =~ /\]$/ ) {
        $var->{standard_gHGVS} = gen_standard_gphgvs( $var->{gHGVS} );
        $var->{alt_gHGVS}      = gen_alt_ghgvs( $var->{gHGVS} );
    }
    if ( $var->{gHGVS} =~ /del[ACGTN]+ins/ ) {
        $var->{standard_gHGVS} = $var->{gHGVS};
        $var->{standard_gHGVS} =~ s/del[ACGTN]+ins/delins/;
    }

    # Due to the bed format database,
    # add flank left and right 1bp to query the database
    # to involve all the mismatches

    my $var_region =
      (     $var->{chr} . ':'
          . ( ( $var->{pos} > 0 ) ? ( $var->{pos} - 1 ) : 0 ) . '-'
          . ( $var->{pos} + $var->{reflen} + 1 ) );

    my $localdb;
    if ( !exists $self->{batch} ) {    # daemon mode
        my %locOpts = ( region => $var_region );
        if ( exists $self->{trans} ) {
            $locOpts{trans} = $self->{trans};
        }
        if ( exists $self->{genes} ) {
            $locOpts{genes} = $self->{genes};
        }
        my $locLoad = $self->load_anno(%locOpts);
        $localdb =
          ( exists $locLoad->{ $var->{chr} } ) ? $locLoad->{ $var->{chr} } : [];
    }
    else {
        if ( !exists $self->{annodb}->{ $var->{chr} } ) {
            $self->warn(
"Warning: incomplete annotation database in batch mode, [$var->{chr}]"
            ) if ($debugOpt);
            $localdb = [];
        }
        else {
            $localdb = $self->{annodb}->{ $var->{chr} };
        }
    }

    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [localdb load] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    my $annoEnt = BedAnno::Anno->new($var);

    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [new BedAnno::Anno] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $AEIndex = $annoEnt->getTrPosition( $localdb, $AEIndex );
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [getTrPosition] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $annoEnt = $self->getTrChange($annoEnt);
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [getTrChange] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    $annoEnt = $self->finaliseAnno($annoEnt);
    if ($debugOpt) {
        $t1 = [gettimeofday];
        print STDERR "BedAnno->varanno [finaliseAnno] ... "
          . tv_interval( $t0, $t1 ) . "\n";
        $t0 = $t1;
    }

    return ( $annoEnt, $AEIndex );
}

=head2 P1toP3

    About : Change 1 letter format of pHGVS string to 3 letter format
    Usage : my $p3 = P1toP3($p1);

=cut

sub P1toP3 {
    my $p1 = shift;
    if ( $p1 =~ /^p\.([A-Z\*])(\d+)([A-Z\*])((fs\*.+)?)$/ ) {

        # missense, frameshift
        return 'p.' . $C1toC3{$1} . $2 . $C1toC3{$3} . $4;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)([A-Z\*])((ext\*.+)?)$/ ) {

        # missense, frameshift
        return 'p.' . $C1toC3{$1} . $2 . $C1toC3{$3} . $4;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)(del|dup|\[.*\])$/ ) {

        # 1 bp deletion, duplication, repeats
        return 'p.' . $C1toC3{$1} . $2 . $3;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)delins([A-Z\*]+)$/ ) {

        # 1 bp delins
        my $former  = 'p.' . $C1toC3{$1} . $2 . 'delins';
        my @singles = split( //, $3 );
        my $latter  = join( "", ( map { $C1toC3{$_} } @singles ) );
        return $former . $latter;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)(del|dup|\[.*\])$/ ) {

        # long deletion, duplication, repeats
        return 'p.' . $C1toC3{$1} . $2 . '_' . $C1toC3{$3} . $4 . $5;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)((del)?ins)([A-Z\*]+)$/ ) {

        # insertion, long delins
        my $former  = 'p.' . $C1toC3{$1} . $2 . '_' . $C1toC3{$3} . $4 . $5;
        my @singles = split( //, $7 );
        my $latter  = join( "", ( map { $C1toC3{$_} } @singles ) );
        return $former . $latter;
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)((delins)?)\?$/ ) {

        # 1 base no-call
        return 'p.' . $C1toC3{$1} . $2 . $3 . '?';
    }
    elsif ( $p1 =~ /^p\.([A-Z\*])(\d+)_([A-Z\*])(\d+)(delins\?)$/ ) {

        # long no-call
        return 'p.' . $C1toC3{$1} . $2 . '_' . $C1toC3{$3} . $4 . $5;
    }
    elsif ( $p1 =~ /^p\.\(([A-Z\*])(\d+)=\)$/ ) {
        return 'p.(' . $C1toC3{$1} . $2 . '=)';
    }
    else {
        return $p1;
    }
}

sub _intronSO {
    my $r = shift;
    if ( $r =~ /^I5U/ ) {
        return 'SO:0000447';
    }
    elsif ( $r =~ /^I3U/ ) {
        return 'SO:0000448';
    }
    elsif ( $r =~ /^I/ ) {
        return 'SO:0000191';
    }
    else {
        confess "Not intron region tag [$r]";
    }
}

sub _checkDistEdge {
    my ( $tranno, $dbk ) = @_;
    my $edge_pos;
    if ( $dbk->{nsta} =~ /^(\d+)[\+\-]/ ) {
        $edge_pos = $1;
    }

    if ( defined $edge_pos ) {
        if (
            (
                $tranno->{rnaBegin} =~ /^\d+$/
                and abs( $tranno->{rnaBegin} - $edge_pos ) < 3
            )
            or ( $tranno->{rnaEnd} =~ /^\d+$/
                and abs( $tranno->{rnaEnd} - $edge_pos ) < 3 )
          )
        {
            return 1;
        }
    }
    return 0;
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
    my ( $self, $annoEnt ) = @_;

    if ( exists $annoEnt->{trInfo} ) {

        # use this extended local db to check if exon variant locates
        # in the 3bp range of exon-intron edge.
        my $av = $annoEnt->{var};
        my $ext_region =
            $av->{chr} . ':'
          . ( ( $av->{pos} - 3 > 0 ) ? ( $av->{pos} - 3 ) : 0 ) . '-'
          . ( $av->{pos} + $av->{reflen} + 3 );
        my $ext_localdb = $self->load_anno( region => $ext_region );
        $ext_localdb = $ext_localdb->{ $av->{chr} };   # no need to check exists

        foreach my $tid ( sort keys %{ $annoEnt->{trInfo} } ) {

            my $trAnnoEnt = $annoEnt->{trInfo}->{$tid};

            if ( $tid =~ /^N[MR]_MT-/ ) {
                $trAnnoEnt->{geneSym} =~ s/^MT-//;
            }

            my $qtid = $tid;
            $qtid =~ s/\-\d+$//;
            my $trdbEnt = $self->{trInfodb}->{$qtid};
            $trAnnoEnt->{prot} =
              ( !exists $trdbEnt->{prot} or $trdbEnt->{prot} eq '.' )
              ? ""
              : $trdbEnt->{prot};

            my ($strd) = ( $trAnnoEnt->{strd} eq '+' ) ? 1 : 0;

            # complete informations depend on already assigned values
            my $genepartSO;
            if (   $trAnnoEnt->{rnaBegin} eq '?'
                or $trAnnoEnt->{rnaEnd} eq '?'
                or $trAnnoEnt->{genepartSO} eq "annotation-fail" )
            {
                # do nothing?
                $genepartSO                  = $trAnnoEnt->{genepartSO};
                $trAnnoEnt->{r}              = '?';
                $trAnnoEnt->{exin}           = '?';
                $trAnnoEnt->{componentIndex} = '';
                $trAnnoEnt->{exonIndex}      = $1
                  if ( $trAnnoEnt->{ei_Begin} =~ /EX(\d+)/ );
                $trAnnoEnt->{intronIndex} = $1
                  if ( $trAnnoEnt->{ei_Begin} =~ /IVS(\d+)/ );
            }
            elsif ( $trAnnoEnt->{r_Begin} eq $trAnnoEnt->{r_End} ) {
                $genepartSO =
                  ( $trAnnoEnt->{genepartSO} =~ /^\d+$/ )
                  ? sprintf( "SO:%07d", $trAnnoEnt->{genepartSO} )
                  : $trAnnoEnt->{genepartSO};

                $trAnnoEnt->{r} = $trAnnoEnt->{r_Begin}
                  if ( !exists $trAnnoEnt->{r} );
                $trAnnoEnt->{exin} = $trAnnoEnt->{ei_Begin}
                  if ( !exists $trAnnoEnt->{exin} );

                if ( $trAnnoEnt->{exin} =~ /^EX(\d+)/ ) {
                    $trAnnoEnt->{componentIndex} = $1;
                    $trAnnoEnt->{exonIndex}      = $1;
                }
                elsif ( $trAnnoEnt->{exin} =~ /^IVS(\d+)/ ) {
                    $trAnnoEnt->{componentIndex} = $1;
                    $trAnnoEnt->{intronIndex}    = $1;
                }
                else {    # PROM
                    $trAnnoEnt->{componentIndex} = 0;
                }
            }
            elsif ( $trAnnoEnt->{ei_Begin} =~ /^IVS/
                and $trAnnoEnt->{ei_Begin} eq $trAnnoEnt->{ei_End} )
            {             # span splice or interior intron
                if (    $trAnnoEnt->{r_Begin} =~ /^D/
                    and $trAnnoEnt->{r_End} =~ /^A/ )
                {
                    $genepartSO = 'span';
                }
                elsif ( $trAnnoEnt->{r_Begin} =~ /^D/ ) {
                    $genepartSO = "SO:0000163";    # 5'
                }
                elsif ( $trAnnoEnt->{r_End} =~ /^A/ ) {
                    $genepartSO = "SO:0000164";    # 3'
                }
                elsif ( $trAnnoEnt->{r_Begin} =~ /^A/ )
                {                                  # insertion on intron edge
                    $genepartSO = _intronSO( $trAnnoEnt->{r_End} );
                }
                elsif ( $trAnnoEnt->{r_End} =~ /^D/ ) {
                    $genepartSO = _intronSO( $trAnnoEnt->{r_Begin} );
                }
                else {
                    $self->confess("Unknown Error in span case intron!");
                }

                $trAnnoEnt->{r} =
                  $trAnnoEnt->{r_Begin} . '-' . $trAnnoEnt->{r_End};
                $trAnnoEnt->{exin}           = $trAnnoEnt->{ei_Begin};
                $trAnnoEnt->{componentIndex} = $1
                  if ( $trAnnoEnt->{ei_Begin} =~ /(\d+)/ );
                $trAnnoEnt->{intronIndex} = $trAnnoEnt->{componentIndex};
            }
            else {
                $genepartSO                  = 'span';
                $trAnnoEnt->{componentIndex} = '';
                $trAnnoEnt->{componentIndex} = $1
                  if ( $trAnnoEnt->{ei_Begin} =~ /(\d+)/ )
                  ;    # index the begin part
                $trAnnoEnt->{r} =
                  $trAnnoEnt->{r_Begin} . '-' . $trAnnoEnt->{r_End};
                if ( $trAnnoEnt->{ei_Begin} eq $trAnnoEnt->{ei_End} ) {
                    $trAnnoEnt->{exin}      = $trAnnoEnt->{ei_Begin};
                    $trAnnoEnt->{exonIndex} = $trAnnoEnt->{componentIndex}
                      if ( $trAnnoEnt->{componentIndex} ne '' );
                }
                else {
                    $trAnnoEnt->{exin} =
                      $trAnnoEnt->{ei_Begin} . '-' . $trAnnoEnt->{ei_End};
                    $trAnnoEnt->{exonIndex} = $1
                      if ( $trAnnoEnt->{ei_Begin} =~ /EX(\d+)/ );
                    $trAnnoEnt->{intronIndex} = $1
                      if ( $trAnnoEnt->{ei_Begin} =~ /IVS(\d+)/ );
                    if (    $trAnnoEnt->{r_Begin} eq 'PROM'
                        and $trAnnoEnt->{r_End} =~ /^5U/ )
                    {
                        $trAnnoEnt->{componentIndex} = 0;
                    }
                    else {    # non-equal exon intron number or UTR/CDS span
                        if ( $trAnnoEnt->{trRef} eq '' ) { # edge insertion case
                            $trAnnoEnt->{exonIndex} = $1
                              if ( $trAnnoEnt->{exin} =~ /EX(\d+)/ );
                            $trAnnoEnt->{intronIndex} = $1
                              if ( $trAnnoEnt->{exin} =~ /IVS(\d+)/ );
                        }
                    }
                }
            }

            confess "Error: unknown genepartSO [$genepartSO]."
              if ( !exists $SO2Name{$genepartSO} );

            $trAnnoEnt->{exin} = '.' if ( !exists $trAnnoEnt->{exin} );
            $trAnnoEnt->{r}    = '.' if ( !exists $trAnnoEnt->{r} );
            $trAnnoEnt->{genepart} = $SO2Name{$genepartSO};
            $trAnnoEnt->{genepartSO} =
              ( $genepartSO =~ /^SO:/ ) ? $genepartSO : "";
            $trAnnoEnt->{componentIndex} = ''
              if ( !exists $trAnnoEnt->{componentIndex} );
            $trAnnoEnt->{exonIndex} = '.'
              if ( !exists $trAnnoEnt->{exonIndex} );
            $trAnnoEnt->{intronIndex} = '.'
              if ( !exists $trAnnoEnt->{intronIndex} );

            # change annotation-fail 's other value to empty?

            # uniform protBegin end
            $trAnnoEnt->{protBegin} = (
                !exists $trAnnoEnt->{protBegin}
                  or $trAnnoEnt->{protBegin} eq '0'
            ) ? "" : $trAnnoEnt->{protBegin};
            $trAnnoEnt->{protEnd} =
              ( !exists $trAnnoEnt->{protEnd} or $trAnnoEnt->{protEnd} eq '0' )
              ? ""
              : $trAnnoEnt->{protEnd};

            $trAnnoEnt->{func} = 'unknown' if ( !exists $trAnnoEnt->{func} );

            # uniform function group
            confess "Error: unknown func code [$trAnnoEnt->{func}]."
              if ( !exists $func2SO{ $trAnnoEnt->{func} } );
            $trAnnoEnt->{funcSO}     = $func2SO{ $trAnnoEnt->{func} };
            $trAnnoEnt->{funcSOname} = $SO2Name{ $trAnnoEnt->{funcSO} };
            $trAnnoEnt->{funcSO} =
              ( $trAnnoEnt->{funcSO} =~ /^SO:/ ) ? $trAnnoEnt->{funcSO} : "";

            if (
                exists $trAnnoEnt->{c}
                and (  $trAnnoEnt->{c} =~ /\d+\+5[ACGT]>[ACGT\?]$/
                    or $trAnnoEnt->{c} =~ /\d+\+5del[ACGT]?$/ )
              )
            {    # only single position
                $trAnnoEnt->{alt_func} = 'splice-5-5th';
            }
            elsif ( $trAnnoEnt->{func} !~ /splice/
                and $trAnnoEnt->{genepart} ne "annotation-fail"
                and !exists $trAnnoEnt->{alt_func} )
            {
                if (
                    $trAnnoEnt->{func} eq 'span'
                    and (
                           $trAnnoEnt->{ei_Begin} =~ /^IVS/
                        or $trAnnoEnt->{ei_End} =~ /^IVS/
                        or (    $trAnnoEnt->{r_Begin} eq 'PROM'
                            and $trAnnoEnt->{ei_End} !~ /^EX1E?$/ )
                        or (    $trAnnoEnt->{r_End} eq '3D'
                            and $trAnnoEnt->{ei_Begin} !~ /^EX(\d+)E$/ )
                    )
                  )
                {
                    $trAnnoEnt->{alt_func} =
                      'splice-region';    # or just ommit this?
                }
                else {
                    if ( $trAnnoEnt->{rnaBegin} =~ /\d+[\+\-](\d+)/ ) {
                        if ( $1 <= 8 ) {
                            $trAnnoEnt->{alt_func} = 'splice-region';
                        }
                        elsif ( $1 <= 10 ) {
                            $trAnnoEnt->{alt_func} = 'splice-ext';
                        }
                    }

                    if ( $trAnnoEnt->{rnaEnd} =~ /\d+[\+\-](\d+)/ ) {
                        if ( $1 <= 8 ) {
                            $trAnnoEnt->{alt_func} = 'splice-region';
                        }
                        elsif ( $1 <= 10 ) {
                            $trAnnoEnt->{alt_func} = 'splice-ext'
                              if ( !exists $trAnnoEnt->{alt_func} );
                        }
                    }

                    if ( !exists $trAnnoEnt->{alt_func} )    # on exon region
                    {
                        for ( my $k = 0 ; $k < @$ext_localdb ; $k++ ) {
                            $$ext_localdb[$k] =
                              $self->assign_detail( $$ext_localdb[$k] )
                              if ( !exists $$ext_localdb[$k]{detail} );

                            next if ( !exists $$ext_localdb[$k]{detail}{$tid} );
                            my $rtd = $$ext_localdb[$k]{detail}{$tid};
                            if ( $rtd->{blka} =~ /^[AD]/ ) {   # hit splice site
                                    # calculate distance from exon-intron edge
                                if ( _checkDistEdge( $trAnnoEnt, $rtd ) ) {
                                    $trAnnoEnt->{alt_func} = 'splice-region';
                                    last;
                                }
                            }
                        }
                    }
                }
            }

            if ( exists $trAnnoEnt->{alt_func} ) {
                $trAnnoEnt->{alt_funcSO} = $func2SO{ $trAnnoEnt->{alt_func} };
                $trAnnoEnt->{alt_funcSOname} =
                  $SO2Name{ $trAnnoEnt->{alt_funcSO} };
                $trAnnoEnt->{alt_funcSO} =
                  ( $trAnnoEnt->{alt_funcSO} =~ /^SO:/ )
                  ? $trAnnoEnt->{alt_funcSO}
                  : "";
            }

            # add additional resource
            if ( exists $trAnnoEnt->{prot} and $trAnnoEnt->{prot} ne "" ) {
                $trAnnoEnt->{protBegin} =
                  ( $trAnnoEnt->{protBegin} eq '0' )
                  ? ""
                  : $trAnnoEnt->{protBegin};
                $trAnnoEnt->{protEnd} =
                  ( $trAnnoEnt->{protEnd} eq '0' ) ? "" : $trAnnoEnt->{protEnd};
                if (    $trAnnoEnt->{protBegin} ne ""
                    and $trAnnoEnt->{protEnd} ne "" )
                {
                    if ( exists $self->{pfam} ) {
                        my ( $pb, $pe ) =
                          ( $trAnnoEnt->{protBegin}, $trAnnoEnt->{protEnd} );
                        if ( $pb > $pe ) {
                            my $tmp = $pb;
                            $pb = $pe;
                            $pe = $tmp;
                        }
                        ( $trAnnoEnt->{pfamId}, $trAnnoEnt->{pfamName} ) =
                          $self->{pfam_h}
                          ->getPfam( $trAnnoEnt->{prot}, $pb, $pe );
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

                    if ( exists $self->{condel_h}
                        and defined $self->{condel_h} )
                    {
                        my $rcondel_info = $self->{condel_h}->pred($trAnnoEnt);
                        if ( exists $rcondel_info->{condelPred}
                            and $rcondel_info->{condelPred} ne
                            "not_computable_was" )
                        {
                            $trAnnoEnt->{condelPred} =
                              $rcondel_info->{condelPred};
                            $trAnnoEnt->{condelScore} =
                              $rcondel_info->{condelScore};
                        }
                    }
                }
            }

            $trAnnoEnt->{standard_pHGVS} =
              gen_standard_gphgvs( $trAnnoEnt->{p} )
              if (  exists $trAnnoEnt->{p}
                and defined $trAnnoEnt->{p}
                and $trAnnoEnt->{p} =~ /\]$/ );

            $trAnnoEnt->{p3} = P1toP3( $trAnnoEnt->{p} )
              if ( exists $trAnnoEnt->{p} );
            $trAnnoEnt->{standard_p3} = P1toP3( $trAnnoEnt->{standard_pHGVS} )
              if ( exists $trAnnoEnt->{standard_pHGVS} );
            $trAnnoEnt->{alt_p3} = P1toP3( $trAnnoEnt->{alt_pHGVS} )
              if ( exists $trAnnoEnt->{alt_pHGVS} );

            # add standard_cHGVS for delXXXinsXXX format cHGVS, for the need of standardization
            if ( exists $trAnnoEnt->{c} and $trAnnoEnt->{c} =~ /del[ACGTN]+ins/) {
                $trAnnoEnt->{standard_cHGVS} = $trAnnoEnt->{c};
                $trAnnoEnt->{standard_cHGVS} =~ s/del[ACGTN]+ins/delins/;
            }

            # 1.23 add standard_cHGVS for gene flanking numbering cHGVS.
            if ( exists $trAnnoEnt->{c} and $trAnnoEnt->{c} =~ /-u|\+d/ ) {
                $trAnnoEnt->{standard_cHGVS} = $trAnnoEnt->{c};
            }
            if ( exists $trAnnoEnt->{standard_cHGVS} ) {
                $trAnnoEnt->{standard_cHGVS} =~ s/(?<=\*|-)(\d+)(-u|\+d)(\d+)/$1+$3/eg;
                $trAnnoEnt->{standard_cHGVS} =~ s/\d+-u/-/g;  # for no-5'utr case
                $trAnnoEnt->{standard_cHGVS} =~ s/\d+\+d/*/g; # for no-3'utr case
            }
            if ( exists $trAnnoEnt->{alt_cHGVS} and $trAnnoEnt->{alt_cHGVS} =~ /-u|\+d/ ) {
                $trAnnoEnt->{alt_cHGVS} =~ s/(?<=\*|-)(\d+)(-u|\+d)(\d+)/$1+$3/eg;
                $trAnnoEnt->{alt_cHGVS} =~ s/\d+-u/-/g;  # for no-5'utr case
                $trAnnoEnt->{alt_cHGVS} =~ s/\d+\+d/*/g; # for no-3'utr case
            }

            $trAnnoEnt->{trVarName} = _getTrVarName( $tid, $trAnnoEnt );
        }
    }

    $annoEnt->{var}->{varName} = decide_major($annoEnt);

    return $annoEnt;
}

# The standard gHGVS and pHGVS can be generated directly from
# original HGVS string, but cHGVS need more specification
sub gen_standard_gphgvs {
    my $rep_hgvs = shift;
    if ( $rep_hgvs =~ /^([gm]\.)(\d+)([ACGTN]+)\[(\d+)>(\d+)\]$/ ) {
        my ( $sym, $anchor, $rep, $ref_cn, $alt_cn ) = ( $1, $2, $3, $4, $5 );
        my $replen = length($rep);
        if ( $alt_cn - $ref_cn == 1 ) {    # convert to dup
            if ( $replen == 1 ) {
                return $sym . ( $anchor + $ref_cn - 1 ) . 'dup' . $rep;
            }
            else {
                return
                    $sym
                  . ( $anchor + ( $ref_cn - 1 ) * $replen ) . '_'
                  . ( $anchor + $ref_cn * $replen - 1 ) . 'dup'
                  . $rep;
            }
        }
        else {
            return $sym . $anchor . $rep . '[' . $alt_cn . ']';
        }
    }
    elsif ( $rep_hgvs =~ /^p\.(\D)(\d+)((_\D)(\d+))?\[(\d+)>(\d+)\]$/ ) {
        my ( $anchor1_char, $anchor1 ) = ( $1, $2 );
        my ( $anchor2_char, $anchor2 ) = ( $4, $5 );
        my ( $ref_cn,       $alt_cn )  = ( $6, $7 );
        $anchor2_char ||= "";
        $anchor2      ||= "";

        if ( $alt_cn - $ref_cn == 1 ) {
            if ( !defined $3 or $3 eq '' ) {
                return
                    'p.'
                  . $anchor1_char
                  . ( $anchor1 + $ref_cn - 1 ) . 'dup';
            }
            else {
                my $replen = $anchor2 - $anchor1 + 1;
                return
                    'p.'
                  . $anchor1_char
                  . ( $anchor1 + ( $ref_cn - 1 ) * $replen )
                  . $anchor2_char
                  . ( $anchor1 + $ref_cn * $replen - 1 ) . 'dup';
            }
        }
        else {
            return
                "p."
              . $anchor1_char
              . $anchor1
              . $anchor2_char
              . $anchor2 . '['
              . $alt_cn . ']';
        }
    }
    else {
        return $rep_hgvs;
    }
}

# only ghgvs can generate alt HGVS directly from repeat format string
sub gen_alt_ghgvs {
    my $rep_hgvs = shift;

    if ( $rep_hgvs =~ /^([gm]\.)(\d+)([ACGTN]+)\[(\d+)>(\d+)\]$/ ) {
        my ( $sym, $anchor, $rep, $ref_cn, $alt_cn ) = ( $1, $2, $3, $4, $5 );
        my $replen = length($rep);

        if ( $ref_cn > $alt_cn ) {    # deletion
            if ( $ref_cn - $alt_cn == 1 and $replen == 1 ) {
                return
                    $sym
                  . ( $anchor + $ref_cn * $replen - 1 ) . 'del'
                  . $rep;
            }
            else {
                return
                    $sym
                  . ( $anchor + $alt_cn * $replen ) . '_'
                  . ( $anchor + $ref_cn * $replen - 1 ) . 'del'
                  . ( $rep x ( $ref_cn - $alt_cn ) );
            }
        }
        else {    # insertion
            return
                $sym
              . ( $anchor + $ref_cn * $replen - 1 ) . '_'
              . ( $anchor + $ref_cn * $replen ) . 'ins'
              . ( $rep x ( $alt_cn - $ref_cn ) );
        }
    }

    return $rep_hgvs;
}

=head2 decide_major

    About   : In finalise step, decide a major transcript to report for a var.
    Usage   : my $majorTranscriptVarName = decide_major($annoEnt);
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
        my $gHGVS =
          ( exists $annoEnt->{var}->{standard_gHGVS} )
          ? $annoEnt->{var}->{standard_gHGVS}
          : $annoEnt->{var}->{gHGVS};
        my $intergenic = $annoEnt->{var}->{chr} . ": " . $gHGVS;
        $intergenic = "chr" . $intergenic if ( $intergenic =~ /^[\dXYM]/ );
        $intergenic .= " (intergenic)";    # for only intergenic

        if ( $intergenic =~ /^chrMT/ ) {
            $intergenic =~ s/^chrMT/$CURRENT_MT/;
        }
        return $intergenic;
    }
    else {
        my %prTrs    = ();
        my %nonPrTrs = ();
        foreach my $tid ( keys %{ $annoEnt->{trInfo} } ) {
            if ( !exists $annoEnt->{trInfo}->{$tid}->{genepart} ) {
                confess
                  "Error: cannot decide major transcript before finaliseAnno";
            }
            if ( $annoEnt->{trInfo}->{$tid}->{primaryTag} eq "Y" ) {
                $prTrs{$tid} = $annoEnt->{trInfo}->{$tid}->{genepart};
            }
            else {
                $nonPrTrs{$tid} = $annoEnt->{trInfo}->{$tid}->{genepart};
            }
        }

        my $majorTr;
        if ( 0 < keys %prTrs ) {    # hit primary transcripts
            if ( 1 == keys %prTrs ) {
                my ($majorTid) = keys %prTrs;
                $majorTr = $majorTid;
            }
            else {
                $majorTr = _get_first_tr( \%prTrs );
            }
        }
        else {                      # hit only non-primary transcripts
            $majorTr = _get_first_tr( \%nonPrTrs );
        }

        my $rTrEnt = $annoEnt->{trInfo}->{$majorTr};
        return $annoEnt->{trInfo}->{$majorTr}->{trVarName};
    }
}

sub _get_first_tr {
    my $rGeneParts_hash = shift;
    my @ordered_Trs     = sort {
        $GenePartsOrder{ $rGeneParts_hash->{$a} }
          <=> $GenePartsOrder{ $rGeneParts_hash->{$b} }
          or $a cmp $b
    } keys %$rGeneParts_hash;

    return $ordered_Trs[0];
}

sub _getTrVarName {
    my ( $trID, $rTrEnt ) = @_;
    my $trhit;

    if ( $trID =~ /N[MR]_MT-/ ) {    # mitochondrial mutation
        $trhit = $rTrEnt->{geneSym};
    }
    else {
        $trhit = $trID . "(" . $rTrEnt->{geneSym} . ")";
    }

    if ( $rTrEnt->{func} eq 'annotation-fail' ) {
        $trhit .= ": annotation-fail";
    }
    else {
        my $cHGVS =
          ( exists $rTrEnt->{standard_cHGVS} )
          ? $rTrEnt->{standard_cHGVS}
          : $rTrEnt->{c};
        $trhit .= ": " . $cHGVS;
    }

    if ( exists $rTrEnt->{p} and $rTrEnt->{p} ne "" ) {
        my $pHGVS =
          ( exists $rTrEnt->{standard_pHGVS} )
          ? $rTrEnt->{standard_pHGVS}
          : $rTrEnt->{p};
        $trhit .= " (" . $pHGVS . ")";
    }
    return $trhit;
}

sub _getCodingSeq {
    my $trdbEnt = shift;
    return "" if ( !exists $trdbEnt->{csta} );
    my $codingSeq = substr( $trdbEnt->{seq}, $trdbEnt->{csta},
        ( $trdbEnt->{csto} - $trdbEnt->{csta} ) );
    return $codingSeq if ( !exists $trdbEnt->{cfs} );
    foreach my $fsld ( sort { $b <=> $a } keys %{ $trdbEnt->{cfs} } ) {
        if (   $fsld <= 0
            or ( $fsld + $trdbEnt->{cfs}->{$fsld} ) <= 0
            or $fsld >= ( length($codingSeq) - 1 ) )
        {
            confess "Error: [$trdbEnt->{gene}] frameshift position error.";
        }
        if ( $trdbEnt->{cfs}->{$fsld} < 0 ) {
            substr(
                $codingSeq,
                $fsld, 0,
                substr(
                    $codingSeq,
                    ( $fsld + $trdbEnt->{cfs}->{$fsld} ),
                    abs( $trdbEnt->{cfs}->{$fsld} )
                )
            );
        }
        else {
            substr( $codingSeq, $fsld, $trdbEnt->{cfs}->{$fsld}, "" );
        }
    }
    return $codingSeq;
}

# let unavailable pr pos to be 0
sub _getCoveredProd {
    my ( $trdbEnt, $chgvs_5, $chgvs_3 ) = @_;
    if (
        $chgvs_3 =~ /^\-/       # 5'UTR
        or $chgvs_5 =~ /^\*/    # 3'UTR
      )
    {
        return ( 0, 0 );
    }

    my ( $prb,        $pre )     = ();
    my ( $prb_anchor, $prb_sig ) = ();
    my ( $pre_anchor, $pre_sig ) = ();

    if ( $chgvs_5 =~ /^(\d+)([\+\-])/ ) {
        $prb_anchor = $1;
        $prb_sig    = $2;
    }

    if ( $chgvs_3 =~ /^(\d+)([\+\-])/ ) {
        $pre_anchor = $1;
        $pre_sig    = $2;
    }

    if (    defined $prb_anchor
        and defined $prb_sig
        and defined $pre_anchor
        and defined $pre_sig )
    {
        if (
            $prb_anchor == $pre_anchor    # assume no 1 bp exon
            or (    $prb_anchor == ( $pre_anchor - 1 )
                and $prb_sig eq '+'
                and $pre_sig eq '-' )
          )
        {
            return ( 0, 0 );              # single intron
        }
    }

    if ( $chgvs_5 =~ /^\-/ ) {
        $prb = 1;
    }
    elsif ( $chgvs_5 =~ /^(\d+)/ ) {
        my $tmp_cdsp = $1;
        ( $prb, $tmp_cdsp ) = _calPosFrame($tmp_cdsp);
    }
    else {
        $prb = 0;
    }

    if ( $chgvs_3 =~ /^\*/ ) {
        $pre = $trdbEnt->{plen} + 1;
    }
    elsif ( $chgvs_3 =~ /^(\d+)/ ) {
        my $tmp_cdsp2 = $1;
        ( $pre, $tmp_cdsp2 ) = _calPosFrame($tmp_cdsp2);
    }
    else {
        $pre = 0;
    }

    return ( $prb, $pre );
}

=head2 getTrChange

    About   : Calculate the transcript changes, based on TrPostition
    Usage   : $beda->getTrChange($annoEnt);
    Returns : assign the following tags in annoEnt
                trRef, prot, c, p, cc, polar, func
                prRef, prAlt

=cut

sub getTrChange {
    my ( $self, $annoEnt ) = @_;

    # do nothing if not hit any transcript
    if ( !exists $annoEnt->{trInfo}
        or 0 == ( scalar keys %{ $annoEnt->{trInfo} } ) )
    {
        return $annoEnt;
    }

    my %getseq_cache = ();
    foreach my $tid ( sort keys %{ $annoEnt->{trInfo} } ) {
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

        # the database hash of tr don't hash the annotation failed transcript
        my $qtid = $tid;
        $qtid =~ s/\-\d+$//;    # trim the multiple mapping indicator in trAcc
        if ( !exists $self->{trInfodb}->{$qtid}
            or $self->{trInfodb}->{$qtid}->{seq} eq "" )
        {
            $self->throw(
                "Error: your fasta database file may not complete. [$qtid]");
        }
        my $trdbEnt = $self->{trInfodb}->{$qtid};
        my $trSeq   = $trdbEnt->{seq};

        $trdbEnt->{cseq} = _getCodingSeq($trdbEnt)
          if ( !exists $trdbEnt->{cseq} and exists $trdbEnt->{csta} );

        my $cdsOpt =
          (       exists $trdbEnt->{prot}
              and $trdbEnt->{prot} ne "."
              and $trannoEnt->{cdsBegin} ne '' ) ? 1 : 0;
        my $strd = ( $trannoEnt->{strd} eq '+' ) ? 1 : 0;

        # fix the trRefComp when ended in 3'downstream

        # debug
        #print STDERR Data::Dumper->Dump(
        #    [$tid, $trannoEnt, $unify_r, $trSeq, $strd],
        #    ["tid", "trannoEnt", "unify_r", "trSeq", "strd"] );

        my $trRef = getTrRef( $trannoEnt, $unify_r, $trSeq, $strd );
        $trannoEnt->{trRef} = $trRef;

        my $trAlt = $trannoEnt->{trAlt};

        $trannoEnt->{prot} = $trdbEnt->{prot}
          if ( exists $trdbEnt->{prot} and $trdbEnt->{prot} ne "." );
        my $f = ($cdsOpt) ? 'c.' : 'n.';
        my $chgvs_5 =
          ($cdsOpt) ? $trannoEnt->{cdsBegin} : $trannoEnt->{rnaBegin};
        my $chgvs_3 = ($cdsOpt) ? $trannoEnt->{cdsEnd} : $trannoEnt->{rnaEnd};
        my $cds_len = 0;
        if ($cdsOpt) {
            ( $trannoEnt->{protBegin}, $trannoEnt->{protEnd} ) =
              _getCoveredProd( $trdbEnt, $chgvs_5, $chgvs_3 );
            $cds_len = $trdbEnt->{csto} - $trdbEnt->{csta};
        }

        # [ aaPos, codon, aa, polar, frame, [framealt] ]
        my ( $rcInfo5, $rcInfo3 ) =
          _getPairedCodon( $trdbEnt, $trannoEnt->{rnaBegin},
            $trannoEnt->{rnaEnd} );

        # we just use position comparison to show transcript ref
        # stat, instead of sm in var, because there may exists
        # indel difference bewteen refgenome and refSeq
        my $cmpPos = $self->cmpPos( $chgvs_5, $chgvs_3 );

        # 2. check if no-call
        if ( $trAlt eq '?' ) {
            $trannoEnt->{func} = 'unknown-no-call';
            if ( $cmpPos == 0 ) {    # 1 bp
                $trannoEnt->{c} = $f . $chgvs_5 . $trRef . '>?';
            }
            elsif ( $cmpPos == 1 ) {    # multiple bp
                $trannoEnt->{c} =
                    $f
                  . $chgvs_5 . '_'
                  . $chgvs_3 . 'del'
                  . ( ( $trRef =~ /^[ACGTN]+$/ ) ? $trRef : "" ) . 'ins?';
            }
            else {                      # ins : reverse the positions
                $trannoEnt->{c} = $f . $chgvs_3 . '_' . $chgvs_5 . 'ins?';
            }

            if ( $rcInfo5->[0] == $rcInfo3->[0] and $rcInfo5->[0] > 0 )
            {                           # single codon
                my $aaOut = $rcInfo5->[2];
                $trannoEnt->{p}     = 'p.' . $aaOut . $rcInfo5->[0] . '?';
                $trannoEnt->{cc}    = $rcInfo5->[1] . '=>?';
                $trannoEnt->{polar} = $rcInfo5->[3] . '=>?';
            }
            elsif ( $rcInfo3->[0] > 0 and $rcInfo5->[0] > 0 ) { # multiple codon

                my $aa5 = $rcInfo5->[2] . $rcInfo5->[0];
                my $aa3 = $rcInfo3->[2] . $rcInfo3->[0];

                my $diopt = ( $rcInfo5->[0] < $rcInfo3->[0] ) ? 1 : 0;

                $trannoEnt->{p} = 'p.'
                  . ( $diopt ? $aa5  : $aa3 ) . '_'
                  . ( $diopt ? $aa3  : $aa5 )
                  . ( $diopt ? 'del' : '' ) . 'ins?';
            }

            next;
        }

        if ( $trRef =~ /=/ ) {    # reference sequence too long
            if ( $trannoEnt->{r_Begin} ne $trannoEnt->{r_End} ) {
                $trannoEnt->{func} = 'span';
            }
            else {
                $trannoEnt->{func} = 'unknown';
            }
            $trannoEnt->{c} = $f . $chgvs_5 . '_' . $chgvs_3 . 'del';
            $trannoEnt->{c} .= 'ins' . $trAlt if ( $trAlt ne "" );
            next;
        }
        if ( $trRef eq $trAlt ) {
            $trannoEnt->{func} = 'no-change';
            $trannoEnt->{c} = $f;
            if ($cmpPos == 0) {
                $trannoEnt->{c} .= $chgvs_5 . $trRef . '=';
            }
            elsif ($cmpPos > 0) {
                $trannoEnt->{c} .= $chgvs_5 . "_" . $chgvs_3 . "=";
            }
            else {
                # hgvs 2.1511 don't refer to this kind of no-change, which
                # come from del/ins mismatch from refgenome and refseq
                # with 0 length refSeq.
                $trannoEnt->{c} .= "=";
            }
            next;
        }

        my ( $trBegin, $trEnd, $real_var, $rUnified ) =
          $self->trWalker( $tid, $trannoEnt );
        my ( $real_p, $real_r, $real_a, $real_rl, $real_al ) = @$rUnified;

        $chgvs_5 =
          ($cdsOpt)
          ? _cPosMark( $trBegin, $trdbEnt->{csta}, $trdbEnt->{csto},
            $trdbEnt->{len} )
          : $trBegin;
        $chgvs_3 =
          ($cdsOpt)
          ? _cPosMark( $trEnd, $trdbEnt->{csta}, $trdbEnt->{csto},
            $trdbEnt->{len} )
          : $trEnd;
        $cmpPos = $self->cmpPos( $chgvs_5, $chgvs_3 );

        # check if a rep get across the utr/cds edge and can be curated
        # to a variant in cds region, add a new key uncurated_cHGVS
        # to restore the uncurated version of cHGVS
        if ( $real_var->{imp} eq 'rep' and $cdsOpt
                and ( $chgvs_5 =~ /^\d+$/ xor $chgvs_3 =~ /^\d+$/ )
                and ( $chgvs_5 =~ /^\-\d+$/ xor $chgvs_3 =~ /^\*\d+$/ )
        )
        {
            my $lesser_len = ($real_rl < $real_al) ? $real_rl : $real_al;
            if (
                ( $chgvs_5 =~ /^\-(\d+)$/ and $1 <= $lesser_len )
                or ( $chgvs_5 =~ /^(\d+)$/ and $cdsOpt
                    and ( $cds_len - $1 ) <= $lesser_len )
              )
            {
                my $offset = $1;
                my $opt_53 = ( $chgvs_5 =~ /^\-/ ) ? 1 : 0;
                if ( !$opt_53 ) {
                    $offset = $cds_len - $offset - 1;
                }
                my $change_cn = 0;
                $change_cn++
                  while ( $offset > $change_cn * $real_var->{replen} );

                # restore uncurated version of cHGVS
                $trannoEnt->{uncurated_cHGVS} =
                  $f . $chgvs_5 . '_' . $chgvs_3 . 'del';
                $trannoEnt->{uncurated_cHGVS} .=
                  ( $real_r =~ /^[ACGTN]+$/ and length($real_r) < 50 )
                  ? $real_r
                  : "";

                if ( $real_al > 0 ) {
                    $trannoEnt->{uncurated_cHGVS} .= 'ins' . $real_a;
                }

                my $trPosChanged = $change_cn * $real_var->{replen};
                $trBegin += $trPosChanged;
                $real_var = BedAnno::Var->new(
                    $tid, 0,
                    $real_rl - $trPosChanged,
                    substr( $real_r, $trPosChanged ),
                    substr( $real_a, $trPosChanged )
                );

                ( $real_p, $real_r, $real_a, $real_rl, $real_al ) =
                  $real_var->getUnifiedVar('+');
                $chgvs_5 =
                  ($cdsOpt)
                  ? _cPosMark(
                    $trBegin,         $trdbEnt->{csta},
                    $trdbEnt->{csto}, $trdbEnt->{len}
                  )
                  : $trBegin;
                $chgvs_3 =
                  ($cdsOpt)
                  ? _cPosMark(
                    $trEnd,           $trdbEnt->{csta},
                    $trdbEnt->{csto}, $trdbEnt->{len}
                  )
                  : $trEnd;
                $cmpPos = $self->cmpPos( $chgvs_5, $chgvs_3 );
            }
        }

 # debug
 #        print STDERR Data::Dumper->Dump(
 #            [ $real_var, $trBegin,  $trEnd,  $chgvs_5,  $chgvs_3,  $cmpPos ],
 #            [ "real_var", "trBegin", "trEnd", "chgvs_5", "chgvs_3", "cmpPos" ]
 #        );

        # * check if a repeat case, and use cerntain chgvs string
        # * assign cHGVS string
        if (
            $real_var->{imp} eq 'rep'

            # not get across the edge of cds/non-cds region
            and !( $chgvs_5 !~ /^\d+$/ xor $chgvs_3 !~ /^\d+$/ )

            # not get across the edge of intron/exon region or promoter / 3'd
            and !(
                $chgvs_5 !~ /\d+[\+\-][ud]?\d+/ xor $chgvs_3 !~
                /\d+[\+\-][ud]?\d+/
            )

            # we currently assume repeat variant won't get across more than
            # two region, so ommit other case test
          )
        {
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

            # 0.4.9
            # add a new key: alt_cHGVS, using for query database which
            # do not have current standard cHGVS string, especially
            # for repeat case
            #
            # 0.5.0
            # generate standard_cHGVS for this kind of variants
            # prior to alt_cHGVS, change alt_cHGVS to only ins/del
            # format notation

            # give this opt to indicate whether to give standard_cHGVS
            # a duplication format notation
            my $dupopt =
              (       $real_var->{ref_cn} > 1
                  and $real_var->{alt_cn} - $real_var->{ref_cn} == 1 )
              ? 1
              : 0;

            if (   $real_var->{ref_cn} > $real_var->{alt_cn}
                or $dupopt )
            {    # deletion and duplication
                if ( !$dupopt ) {
                    $trannoEnt->{standard_cHGVS} =
                      $f . $chgvs_5 . $trRep . '[' . $real_var->{alt_cn} . ']';
                }

                my $changed_cn =
                  abs( $real_var->{ref_cn} - $real_var->{alt_cn} );
                my $changed_type = ($dupopt) ? 'dup' : 'del';
                my $changed_cont = $trRep x $changed_cn;
                my $changed_len  = length($changed_cont);
                if ( $changed_cn == 1 and $real_var->{replen} == 1 ) {
                    my $single_change = $f . $chgvs_3 . $changed_type . $trRep;
                    if ($dupopt) {
                        $trannoEnt->{standard_cHGVS} = $single_change;
                    }
                    else {
                        $trannoEnt->{alt_cHGVS} = $single_change;
                    }
                }
                else {
                    my $renew_offset =
                        ($dupopt)
                      ? ( $real_var->{replen} * ( $real_var->{ref_cn} - 1 ) )
                      : ( $real_var->{replen} * $real_var->{alt_cn} );

                    my $renew_cHGVS;
                    if ( $chgvs_5 =~ /^\d+$/ ) {    # cds / ncRNA exon
                        $renew_cHGVS =
                            $f
                          . ( $chgvs_5 + $renew_offset ) . '_'
                          . $chgvs_3
                          . $changed_type
                          . $changed_cont;
                    }
                    elsif ( $chgvs_5 =~ /^([\*\+]?)(\-?\d+)$/ ) {

                        # utr region / ncRNA promoter 3'd region
                        my $sig       = $1;
                        my $start_pos = $2;
                        if ( $sig eq '' ) {    # 5 utr
                            $renew_cHGVS =
                                $f
                              . ( $start_pos + $renew_offset ) . '_'
                              . $chgvs_3
                              . $changed_type
                              . $changed_cont;
                        }
                        else {                 # 3 utr
                            $renew_cHGVS =
                                $f
                              . $sig
                              . ( $start_pos + $renew_offset ) . '_'
                              . $chgvs_3
                              . $changed_type
                              . $changed_cont;
                        }
                    }
                    elsif ( $chgvs_5 =~ /^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$/ ) {

                        # intron or promoter or 3'downstream
                        my $anchor5 = $1;
                        my $sig5    = $2;
                        my $offset5 = $3;
                        my $u5_opt  = 0;
                        if ( $offset5 =~ /u/ ) {
                            $u5_opt = 1;
                            $offset5 =~ s/u//;
                        }
                        if ( $chgvs_3 =~ /^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$/ ) {
                            my $anchor3 = $1;
                            my $sig3    = $2;
                            my $offset3 = $3;
                            my $u3_opt  = 0;
                            if ( $offset3 =~ /u/ ) {
                                $u3_opt = 1;
                                $offset3 =~ s/u//;
                            }
                            if ( $anchor3 eq $anchor5 and $sig5 eq $sig3 ) {
                                my $anc_5 = $offset5 + $renew_offset;
                                $anc_5 =~ s/^-/-u/ if ($u5_opt);
                                $renew_cHGVS =
                                    $f
                                  . $anchor5
                                  . $sig5
                                  . $anc_5 . '_'
                                  . $chgvs_3
                                  . $changed_type
                                  . $changed_cont;
                            }
                            elsif ( $anchor3 ne $anchor5
                                and $sig5 eq '+'
                                and $sig3 eq '' )
                            {    # assume same intron
                                my $half_offset =
                                  ( $real_rl - ( $offset5 + $offset3 ) ) / 2;

                                if ( $renew_offset < $half_offset ) {
                                    $renew_cHGVS =
                                        $f
                                      . $anchor5
                                      . $sig5
                                      . ( $offset5 + $renew_offset ) . '_'
                                      . $chgvs_3
                                      . $changed_type
                                      . $changed_cont;
                                }
                                else {
                                    $renew_cHGVS =
                                        $f
                                      . $anchor3
                                      . $sig3
                                      . ( $offset3 - $changed_len + 1 )
                                      . '_'
                                      . $chgvs_3
                                      . $changed_type
                                      . $changed_cont;
                                }
                            }
                            else {
                                $self->warn(
                                    "Warning: repeat get across more than",
                                    " 3 different region: $tid:$trannoEnt->{c}"
                                ) if ( !exists $self->{quiet} );
                            }
                        }
                        else {
                            $self->warn(
                                "Warning: not both intron while parsing",
                                " alt_cHGVS : $tid:$trannoEnt->{c}"
                            ) if ( !exists $self->{quiet} );
                        }
                    }
                    else {
                        $self->warn(
                            "Warning: Unknown chgvs5 while ",
                            "parsing alt_cHGVS : $tid: $chgvs_5"
                        ) if ( !exists $self->{quiet} );
                    }

                    if ( defined $renew_cHGVS ) {
                        if ($dupopt) {
                            $trannoEnt->{standard_cHGVS} = $renew_cHGVS;
                        }
                        else {
                            $trannoEnt->{alt_cHGVS} = $renew_cHGVS;
                        }
                    }
                }
            }

            if (    $real_var->{ref_cn} < $real_var->{alt_cn}
                and $real_var->{alt_cn} > 2 )
            {    # non-duplication insertion

                my $ins_cn   = $real_var->{alt_cn} - $real_var->{ref_cn};
                my $ins_cont = $trRep x $ins_cn;

                $trannoEnt->{standard_cHGVS} =
                  $f . $chgvs_5 . $trRep . '[' . $real_var->{alt_cn} . ']'
                  if ( !$dupopt );

                if ( $chgvs_3 =~ /^\d+$/ ) {    # cds / ncRNA exon
                    if (
                        (
                                $cdsOpt
                            and $chgvs_3 + 1 <= $cds_len
                        )
                        or ( !$cdsOpt and $chgvs_3 + 1 <= $trdbEnt->{len} )
                      )
                    {
                        $trannoEnt->{alt_cHGVS} =
                            $f
                          . $chgvs_3 . '_'
                          . ( $chgvs_3 + 1 ) . 'ins'
                          . $ins_cont;
                    }
                    else {
                        my $outlatter;
                        if ($cdsOpt) {
                            $outlatter = "*1";
                        }
                        else {
                            $outlatter = "+1";
                        }
                        $trannoEnt->{alt_cHGVS} =
                          $f . $chgvs_3 . '_' . $outlatter . 'ins' . $ins_cont;
                    }
                }
                elsif ( $chgvs_3 =~ /^([\*\+]?)(\-?\d+)$/ ) { # utr or ncRNA ext
                    my $sig     = $1;
                    my $ins_pos = $2;
                    if ( $ins_pos eq '-1' ) {
                        $trannoEnt->{alt_cHGVS} = $f . '-1_1ins' . $ins_cont;
                    }
                    else {
                        $trannoEnt->{alt_cHGVS} =
                            $f
                          . $chgvs_3 . '_'
                          . $sig
                          . ( $ins_pos + 1 ) . 'ins'
                          . $ins_cont;
                    }
                }
                elsif ( $chgvs_3 =~ /^(\-?[^\-\+]+)(\+?d?)(\-?u?\d+)$/ )
                {    # intron
                    my $anchor = $1;
                    my $sig    = $2;
                    my $ofst   = $3;
                    my $u3opt  = 0;
                    if ( $ofst =~ /u/ ) {
                        $u3opt = 1;
                        $ofst =~ s/u//;
                    }
                    if ( $ofst eq '-1' ) {
                        $trannoEnt->{alt_cHGVS} =
                          $f . $chgvs_3 . '_' . $anchor . 'ins' . $ins_cont;
                    }
                    else {
                        my $anc3 = $ofst + 1;
                        $anc3 =~ s/^-/-u/ if ($u3opt);
                        $trannoEnt->{alt_cHGVS} =
                            $f
                          . $chgvs_3 . '_'
                          . $anchor
                          . $sig
                          . $anc3 . 'ins'
                          . $ins_cont;
                    }
                }
                else {
                    $self->warn(
                        "Warning: Unknown chgvs5 while ",
                        "parsing alt_cHGVS : $tid: $chgvs_3"
                    ) if ( !exists $self->{quiet} );
                }
            }
        }
        else {
            if ( $cmpPos == 0 ) {    # 1 bp

                $trannoEnt->{c} = $f . $chgvs_5;
                if ( $real_al == 1 ) {
                    $trannoEnt->{c} .= $real_r . '>' . $real_a;    # 1bp alt
                }
                elsif ( $real_al == 0 ) {
                    $trannoEnt->{c} .=
                      'del' . ( ( $real_r =~ /^[ACGTN]+$/ ) ? $real_r : "" );
                }
                else {
                    $trannoEnt->{c} .= 'del'
                      . ( ( $real_r =~ /^[ACGTN]+$/ ) ? $real_r : "" ) . 'ins'
                      . $real_a;
                }

            }
            elsif ( $cmpPos > 0 ) {    # multiple bp

                $trannoEnt->{c} = $f . $chgvs_5 . '_' . $chgvs_3 . 'del';
                $trannoEnt->{c} .=
                  ( $real_r =~ /^[ACGTN]+$/ and length($real_r) < 50 )
                  ? $real_r
                  : "";

                if ( $real_al > 0 ) {
                    $trannoEnt->{c} .= 'ins' . $real_a;
                }
            }
            else {                     # ins : reverse the positions
                $trannoEnt->{c} =
                  $f . $chgvs_3 . '_' . $chgvs_5 . 'ins' . $real_a;
            }
        }

        # 3. check if transcript-ablation
        if (    ( $trBegin =~ /^\-/ or $trBegin eq '1' )
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

        if (
                exists $trannoEnt->{preStart}
            and exists $trannoEnt->{postEnd}
            and $trannoEnt->{preStart}->{exin} ne $trannoEnt->{postEnd}->{exin}
            and $trannoEnt->{preStart}->{exin} =~ /^IVS|^\./
            and (  $trannoEnt->{postEnd}->{exin} =~ /^IVS/
                or $trannoEnt->{postEnd}->{r} =~ /^3D/ )
          )
        {
            $trannoEnt->{alt_func} = 'exon-loss';    # highest priority
                 # due to the possibility of only delete one exon
                 # without any intron region, which can be further
                 # analysis the protein change, so we continue
                 # without skip the following steps
        }

        # 4. check if span different exon/intron/promoter/downstream
        if (
            (
                   $chgvs_3 =~ /^\+|\+d/
                or $trannoEnt->{ei_Begin} ne $trannoEnt->{ei_End}
            )
            and $cmpPos >= 0    # not insertion at the edge
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
            else {
                if (    exists $self->{genome_h}
                    and defined $self->{genome_h}
                    and $trannoEnt->{r_Begin} =~ /^[AD]/
                    and $trannoEnt->{r_End} eq $trannoEnt->{r_Begin}
                    and length($trRef) == length($trAlt) )
                {    # only in splice site substitution, to check if conanical
                    my $cis_tag = substr( $trannoEnt->{r_Begin}, 0, 1 );
                    my $ext_len = 2 - length($trRef);
                    $self->throw("Error length for trRef") if ( $ext_len < 0 );
                    my ( $Ladded, $Radded ) = ( "", "" );
                    my ( $gchr, $gstart, $gend ) = (
                        $annoEnt->{var}->{chr},
                        $annoEnt->{var}->{pos},
                        $annoEnt->{var}->{end}
                    );

                    $gchr = "chr" . $gchr if ( $gchr !~ /^chr/ );

                    if ( $ext_len == 1 ) {    # only can be 1
                        if (
                            (
                                $trannoEnt->{strd} eq '+'
                                and (  $trannoEnt->{rnaBegin} =~ /\+2$/
                                    or $trannoEnt->{rnaBegin} =~ /\-1$/ )
                            )
                            or (
                                $trannoEnt->{strd} eq '-'
                                and (  $trannoEnt->{rnaEnd} =~ /\+1$/
                                    or $trannoEnt->{rnaEnd} =~ /\-2$/ )
                            )
                          )
                        {
                            my $extpos = $gstart;    # extend 1 bp left
                            my $rgn_tmp = $gchr . ":" . $extpos . "-" . $extpos;
                            if ( exists $getseq_cache{$rgn_tmp} ) {
                                $Ladded = $getseq_cache{$rgn_tmp};
                            }
                            else {
                                my $toAdd1 =
                                  $self->{genome_h}->getseq($rgn_tmp);
                                if ( defined $toAdd1 ) {
                                    $Ladded = uc($toAdd1);
                                    $getseq_cache{$rgn_tmp} = $Ladded;
                                }
                                else {
                                    $self->warn( "Warning: Can not get string from genome for $rgn_tmp");
                                }
                            }
                        }
                        else {
                            my $extpos2 = $gend + 1;    # extend 1 bp right
                            my $rgn_tmp2 =
                              $gchr . ":" . $extpos2 . "-" . $extpos2;
                            if ( exists $getseq_cache{$rgn_tmp2} ) {
                                $Radded = $getseq_cache{$rgn_tmp2};
                            }
                            else {
                                my $toAdd2 =
                                  $self->{genome_h}->getseq($rgn_tmp2);
                                if ( defined $toAdd2 ) {
                                    $Radded = uc($toAdd2);
                                    $getseq_cache{$rgn_tmp2} = $Radded;
                                }
                                else {
                                    $self->warn( "Warning: Can not get string from genome for $rgn_tmp2");
                                }
                            }
                        }

                        if ( $trannoEnt->{strd} eq '-' ) { # change to tr strand
                            my $tmp = $Ladded;
                            $Ladded = $self->rev_comp($Radded);
                            $Radded = $self->rev_comp($tmp);
                        }
                    }

                    my ( $toCheckRef, $toCheckAlt ) = (
                        $Ladded . uc($trRef) . $Radded,
                        $Ladded . uc($trAlt) . $Radded
                    );

                    if (    $toCheckRef ne $canonicalSS{$cis_tag}
                        and $toCheckAlt eq $canonicalSS{$cis_tag} )
                    {
                        $trannoEnt->{func} = 'no-change';
                        # here may need to keep the original format of hgvs. 
                        # $trannoEnt->{c}    = 'c.=';
                        next;
                    }
                }

                if (    $trannoEnt->{r_Begin} =~ /^D/
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
        {    # any edge-insertion case will count on the exon change
                # any coding-utr edge insertion case will count on the
                # utr parts
            if ( !$cdsOpt ) {
                if ( $chgvs_5 =~ /^\+/ ) {
                    $trannoEnt->{func} = 'unknown';
                }
                elsif ( $chgvs_3 =~ /^\-/ ) {
                    $trannoEnt->{func} = 'promoter';
                }
                elsif ( $chgvs_5 eq '1' and $chgvs_3 eq $trdbEnt->{len} ) {
                    $trannoEnt->{func} = 'knockout';
                }
                else {
                    $trannoEnt->{func} = 'ncRNA';
                }
                next;
            }
            else {    # coding RNA
                if ( $chgvs_3 =~ /\-u1$/ ) {
                    $trannoEnt->{func} = 'promoter';
                    next;
                }
                elsif ( $chgvs_5 =~ /\+d1$/ ) {   # 3' downstream edge insertion
                    $trannoEnt->{func} = 'unknown';
                }
                elsif ( $chgvs_3 =~ /^\-/ ) {     # besides utr-5 - cds edge
                    $trannoEnt->{func} = 'utr-5';
                    next;
                }
                elsif ( $chgvs_5 =~ /^\*/ ) {     # besides utr-3 - cds edge
                    $trannoEnt->{func} = 'utr-3';
                    next;
                }
                elsif ( $chgvs_5 =~ /^\-/ and $chgvs_3 =~ /^\*/ ) {
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
                elsif ( $chgvs_3 eq '1-1' )
                {    # cds-begin - last 5'utr intron edge
                    $trannoEnt->{func} = 'utr-5';
                    next;
                }
                elsif ( $chgvs_5 =~ /^(\d+)\+1/
                    and $1 == $cds_len )
                {
                    $trannoEnt->{func} = 'utr-3';
                    next;
                }

                # here protein sequence var should be reparsed
                elsif ( $chgvs_5 =~ /^\-(\d+)/ ) {
                    my $u5_len = $1;
                    if ( $real_var->{imp} eq 'rep' ) {
                        if ( $u5_len > ( $real_rl - $real_al ) ) {
                            $trannoEnt->{func} = 'utr-5';
                        }
                        else {
                            $trannoEnt->{func} = 'unknown';
                        }
                    }
                    else {
                        $trannoEnt->{func} = 'init-loss';
                        $trannoEnt->{p}    = 'p.0?';
                    }
                    next;
                }
                else {
                    if ( $chgvs_3 =~ /^\*(\d+)/ ) {

                        # variant get across the 3' end of cds
                        my $u3_len = $1;
                        if ( $real_var->{imp} eq 'rep'
                            and ( $u3_len > ( $real_rl - $real_al ) ) )
                        {
                            $trannoEnt->{func} = 'utr-3';

                            # hgvs 2.1511 don't refer to this kind of
                            # coding synon, just keep it.
                            $trannoEnt->{p}    = 'p.(=)';
                            next;
                        }
                    }

                    $chgvs_5 = $1 + 1
                      if ( $cmpPos < 0 and $chgvs_5 =~ /^(\d+)\+1$/ );
                    $chgvs_3 = $1 - 1
                      if ( $cmpPos < 0 and $chgvs_3 =~ /^(\d+)\-1$/ );
                    $trBegin = $1 + 1
                      if ( $cmpPos < 0 and $trBegin =~ /^(\d+)\+1$/ );
                    $trEnd = $1 - 1
                      if ( $cmpPos < 0 and $trEnd =~ /^(\d+)\-1$/ );

                      if ( $chgvs_5 !~ /^\d+$/ ) {
                          # special case of transcript - refgenome differing
                          if ($chgvs_5 =~ /\d([\+\-])(\d+)$/ and my ($c5s, $c5i) = ($1, $2) and $chgvs_3 =~ /\d([\+\-])(\d+)$/ and my ($c3s, $c3i) = ($1, $2)) {
                              if ($c5s eq $c3s and $c5i < 3 and $c3i < 3) {
                                  $trannoEnt->{func} = 'splice';
                                  $trannoEnt->{func} .= ($c5s eq '+') ? '-5' : '-3';
                              }
                              else {
                                  $trannoEnt->{func} = 'intron';
                              }
                              next;
                          }
                          else {
                              $self->throw("Error: unavailable 5' cHGVS. $chgvs_5 annoEnt:".Dumper($trannoEnt))
                                if ( $chgvs_5 == 0 );
                          }
                      }
                    # 5' end should be in coding region.
                    # leave this as debug information
                    $self->throw("Error: unavailable 5' cHGVS. $chgvs_5 annoEnt:".Dumper($trannoEnt))
                      if ( $chgvs_5 == 0 );

                    ( $rcInfo5, $rcInfo3 ) =
                      _getPairedCodon( $trdbEnt, $trBegin, $trEnd );

                    if ( $real_a eq "" )
                    {    # for deletion case try to predict the pr change
                        if ( 5 < @$rcInfo5 ) {
                            splice( @$rcInfo5, 4, 1 );
                        }
                        if ( 5 < @$rcInfo3 ) {
                            splice( @$rcInfo3, 4, 1 );
                        }

                        if (    $rcInfo5->[0] > 0
                            and $rcInfo3->[0] > 0
                            and $rcInfo5->[4] == -1
                            and $rcInfo3->[4] == -1 )
                        {
                            my $hit_whole_deleted_fs = 0;
                            foreach my $fsld ( keys %{ $trdbEnt->{nfs} } ) {
                                if ( $fsld == ( $trBegin - 1 )
                                    and ( $fsld + $trdbEnt->{nfs}->{$fsld} ) ==
                                    $trEnd )
                                {
                                    $hit_whole_deleted_fs = 1;
                                    last;
                                }
                            }
                            if ($hit_whole_deleted_fs) {
                                $trannoEnt->{protBegin} = $rcInfo5->[0];
                                $trannoEnt->{protEnd}   = $rcInfo3->[0];
                                $trannoEnt->{func}      = 'no-change';
                                next;
                            }
                        }
                    }

                    if (    $rcInfo5->[0] > 0
                        and $rcInfo5->[4] == -1
                        and $rcInfo3->[0] > 0
                        and $rcInfo3->[4] == -1 )
                    {
                        $trannoEnt->{protBegin} = $rcInfo5->[0];
                        $trannoEnt->{protEnd}   = $rcInfo3->[0];
                        $trannoEnt->{func}      = 'abnormal-fs-site';
                        next;
                    }

                    my $real_rl_hidden_ofst = 0;
                    if ( exists $trdbEnt->{nfs} ) {
                        foreach my $fsld ( keys %{ $trdbEnt->{nfs} } ) {
                            my $fsEd = $fsld + $trdbEnt->{nfs}->{$fsld};
                            if (
                                (
                                        $trdbEnt->{nfs}->{$fsld} > 0
                                    and $fsld >= $trBegin
                                    and $fsEd <= $trEnd
                                )
                                or (
                                    $trdbEnt->{nfs}->{$fsld} < 0
                                    and (
                                        (
                                                $fsEd >= $trBegin
                                            and $fsld < $trEnd
                                        )
                                        or
                                        ( $fsld >= $trEnd and $fsEd < $trBegin )
                                    )
                                )
                              )
                            {    # contain fs site or in duplicated fs site
                                $real_rl_hidden_ofst +=
                                  $trdbEnt->{nfs}->{$fsld};
                            }
                            else {
                                if (    $trdbEnt->{nfs}->{$fsld} > 0
                                    and $fsld < $trBegin
                                    and $trBegin <= $fsEd )
                                {    # begin in deleted fs site
                                    $real_rl_hidden_ofst +=
                                      ( $fsEd - $trBegin );
                                }
                                if (    $trdbEnt->{nfs}->{$fsld} > 0
                                    and $fsld < $trEnd
                                    and $trEnd <= $fsEd )
                                {    # end in deleted fs site
                                    $real_rl_hidden_ofst += ( $trEnd - $fsld );
                                }
                                if (
                                    $trdbEnt->{nfs}->{$fsld} < 0
                                    and (
                                        (
                                                $fsld >= $trBegin
                                            and $fsEd < $trBegin
                                        )
                                        or
                                        ( $fsld >= $trEnd and $fsEd < $trEnd )
                                    )
                                  )
                                {    # start or end in duplicated fs site
                                    $real_rl_hidden_ofst +=
                                      $trdbEnt->{nfs}->{$fsld};
                                }
                            }
                        }
                    }

                    my %translate_opts = ();
                    $translate_opts{mito}  = 1 if ( $tid =~ /^NM_MT-/ );
                    $translate_opts{polyA} = 1 if ( exists $trdbEnt->{A} );
                    $translate_opts{nostop} = 1
                      if ( exists $trdbEnt->{X} or exists $trdbEnt->{U} );
                    my %altcodon_opts = %translate_opts;
                    delete $altcodon_opts{nostop}
                      if ( exists $altcodon_opts{nostop} );

                    if ( !exists $trdbEnt->{pseq} ) {
                        my ( $pseq, $frame_next ) =
                          translate( $trdbEnt->{cseq}, \%translate_opts );
                        $trdbEnt->{pseq} = $pseq;
                    }

                    my ( $prBegin, $prEnd, $prRef, $prAlt );

                    $prBegin = $rcInfo5->[0];

                    my $ready_to_add_5 =
                      substr( $rcInfo5->[1], 0, $rcInfo5->[4] );

                    my $diff_ra = $real_rl - $real_al - $real_rl_hidden_ofst;

                    # probably frame shift flag
                    my $frameshift_flag = ( $diff_ra % 3 > 0 ) ? 1 : 0;

                    my $start_in_cds_flag = ($chgvs_5 > ( $cds_len - 3 )) ? 0 : 1;

                    # end in cds or not.
                    my $end_in_cds_flag = (
                             $chgvs_3 =~ /^\*/
                          or $chgvs_3 > ( $cds_len - 3 )
                    ) ? 0 : 1;

                    my $ready_to_add_3;
                    if ( !$end_in_cds_flag or $frameshift_flag ) {
                        $ready_to_add_3 = _cdsubstr( $trdbEnt, $trEnd );

                        # this variants's effect will be end at the terminal
                        $prEnd = $trdbEnt->{plen} + 1;    # terminal
                        $ready_to_add_3 .=
                          "A" x ( 3 - ( length( $trdbEnt->{cseq} ) % 3 ) )
                          if ( exists $trdbEnt->{A} );
                    }
                    else {    # chgvs_3 in cds region with no frameshift
                        $ready_to_add_3 = substr(
                            $rcInfo3->[1],
                            ( $rcInfo3->[4] + 1 ),
                            ( 2 - $rcInfo3->[4] )
                        );
                        $prEnd = $rcInfo3->[0];
                    }

                    if (
                        ( $trdbEnt->{plen} + 1 ) != length( $trdbEnt->{pseq} ) )
                    {
                        $self->warn(
                            "[Critical Warning] : annotation for $tid may not be correct, \n",
                            "                     due to a mis-translated base in its \n",
                            "                     cds region, usually a frameshift on \n",
                            "                     TGA to GAX, please ignore any annotation \n",
                            "                     on latter part behind that frameshift.\n",
                            "                     protein annotation for $tid will be skipped\n",
                            "                     and the function of it will be annotation-fail\n",
                            "                     until the problem is fixed."
                        ) if ( !exists $self->{quiet} );
                        $trannoEnt->{func} = 'annotation-fail';
                        next;
                    }

                    # Make protBegin and protEnd be coordinate to the
                    # coordinates of nucl variants, extend the
                    # frameshift's effect to the end of protein
                    $trannoEnt->{protBegin} = $prBegin;
                    $trannoEnt->{protEnd}   = $prEnd;

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

                    # here only issues variant start/stop on exon region
                    # and 5' start position on the coding region
                    # init codon flag
                    my $hit_init_flag = ( $rcInfo5->[0] == 1 ) ? 1 : 0;

                    # stop codon flag
                    my $hit_stop_flag =
                      ( $trEnd > ( $trdbEnt->{csta} + 3 * $trdbEnt->{plen} ) )
                      ? 1
                      : 0;

                    # check if an alterstart
                    my %start_codons = ( 'ATG' => -1 );
                    if ( exists $trdbEnt->{altstart} ) {
                        %start_codons =
                          map { $_ => 1 } keys %{ $trdbEnt->{altstart} };
                    }

                    my $init_synon = 0;    # indicate start codon search result
                    if ($hit_init_flag) {

                        # check if altered sequence have new start codon
                        foreach my $startCodon ( sort keys %start_codons ) {
                            if ( substr( $codon_alt, 0, 3 ) eq $startCodon ) {
                                $init_synon = 1;
                                last;
                            }
                        }

                        if ( $init_synon == 0 ) {
                            if ( $trAlt =~ /N/ ) {
                                my $all_hit = 1;
                                foreach my $base ( 'A', 'C', 'G', 'T' ) {
                                    my $subalt_codon = $codon_alt;
                                    $subalt_codon =~ s/N/$base/;
                                    if ( !exists $start_codons{$subalt_codon} )
                                    {
                                        $all_hit = 0;
                                        last;
                                    }
                                }

                                $trannoEnt->{cc} =
                                  $rcInfo5->[1] . "=>" . $codon_alt
                                  if ( 3 == length($codon_alt) );
                                if ( $all_hit == 1 ) {
                                    $trannoEnt->{func}  = "altstart";
                                    $trannoEnt->{prAlt} = "M";
                                    $init_synon         = 1;
                                }
                                else {
                                    $trannoEnt->{func} = 'unknown-no-call';
                                }
                            }
                            else {
                                $trannoEnt->{func} = 'init-loss';
                                $trannoEnt->{p}    = 'p.0?';
                                if ( $prBegin eq '1' and $prEnd eq '1' ) {
                                    my $zero;
                                    $trannoEnt->{prRef} = 'M';
                                    my $altCodon = substr( $codon_alt, 0, 3 );
                                    ( $trannoEnt->{prAlt}, $zero ) =
                                      translate( $altCodon, \%altcodon_opts );
                                }
                                next;
                            }
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
                    my $codon_to_be = substr( $codon_alt, 0, 3 );

                    if (    $prBegin - 1 == $trdbEnt->{plen}
                        and $prEnd == $prBegin )
                    {
                        $prRef = '*';
                    }
                    elsif ( $prBegin - 1 <= $trdbEnt->{plen} ) {
                        $prRef = substr(
                            $trdbEnt->{pseq},
                            ( $prBegin - 1 ),
                            ( $prEnd - $prBegin + 1 )
                        );
                    }
                    else {
                        $self->throw(
                            "Error: Cannot get prRef from $tid, [$trdbEnt->{plen}, $prBegin, $prEnd]\n",
                            "      Var: [$annoEnt->{var}->{chr}, $annoEnt->{var}->{pos}, $annoEnt->{var}->{end},",
                            "$annoEnt->{var}->{ref}, $annoEnt->{var}->{alt}]"
                        );
                    }

                    # we can not allow inseq stop codon in altered sequence
                    # due to frameshift and ambiguity.

                    my $next_alt_frame;
                    ( $prAlt, $next_alt_frame ) =
                      translate( $codon_alt, \%altcodon_opts );

                    if ($hit_init_flag) {
                        $prRef =~ s/^[A-Z]/M/;    # change init pep to M
                                                  # change alt-init to M
                        $prAlt =~ s/^[A-Z]/M/ if ($init_synon);
                    }

                    $trannoEnt->{prRef} = $prRef;
                    $trannoEnt->{prAlt} ||= $prAlt;

                    $prAlt = $trannoEnt->{prAlt};    # uniform prAlt

                    # to indicate whether the alternate sequence
                    # encode a stop codon or non-frameshift, or otherwise
                    # with a non stopped frameshift.
                    my $non_stop_flag = (
                        $prAlt !~ /\*$/ and ( $prRef eq "*"
                            or $next_alt_frame
                            or $frameshift_flag )
                    ) ? 1 : 0;

                    # To avoid large range perlre, before parse protein
                    # var we should first deal with frameshift issues
                    # need to deal with non_stop_flag and frameshift and
                    # imp ref case.

                    # assign cc and polar to the same length,
                    # single aa substitution, no matter which position.
                    if ( $diff_ra == 0 and $rcInfo5->[0] == $rcInfo3->[0] ) {

                        # using trim due to terminal codon
                        my ( $aa_to_be, $polar_to_be ) =
                          _getAAandPolar( $codon_to_be, \%altcodon_opts );

                        $trannoEnt->{cc} = $rcInfo5->[1] . '=>' . $codon_to_be;
                        if ( $prRef ne $prAlt ) {
                            $trannoEnt->{polar} =
                              $rcInfo5->[3] . '=>' . $polar_to_be;
                        }
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
                    my $no_parsed_prAlt_Start =
                      substr( $prAlt, $prSimpleSame, 1 );

                    # non stop frameshift with extra non-coded base
                    if ($non_stop_flag) {
                        my $local_fs_ext;
                        if ( $no_parsed_pP >= $trdbEnt->{plen} ) {
                            $trannoEnt->{func} = 'stop-loss';
                            $local_fs_ext = 'ext*?';
                        }
                        else {
                            $trannoEnt->{func} = 'frameshift';
                            $local_fs_ext = 'fs*?';
                        }
                        $trannoEnt->{p} = 'p.'
                          . $no_parsed_prStart
                          . ( $no_parsed_pP + 1 )
                          . $no_parsed_prAlt_Start . $local_fs_ext;
                        $trannoEnt->{protBegin} += $prSimpleSame;
                        next;
                    }

                    # identical to reference (can be from any kind of vars)
                    if ( $prRef eq $prAlt ) {
                        if ($init_synon) {    # altstart
                                              # do nothing
                        }
                        elsif ($hit_stop_flag) {

                            # need this special assertion?
                            $trannoEnt->{func} = 'stop-retained';
                        }
                        else {
                            $trannoEnt->{func} = 'coding-synon';
                        }

                        if (1 == length($prRef)) {
                            $trannoEnt->{p} = 'p.(' . $prRef . $no_parsed_pP . '=)';
                        }
                        else {
                            $trannoEnt->{p} = 'p.(' . $prBegin . '_' . $no_parsed_pP  . '=)';
                        }

                        next;
                    }

                    # frameshift
                    if ( !$end_in_cds_flag or $frameshift_flag ) {

                        $trannoEnt->{func} = ($start_in_cds_flag and $frameshift_flag) ? 'frameshift' : 'stop-loss';
                        $trannoEnt->{p}    = 'p.'
                          . $no_parsed_prStart
                          . ( $no_parsed_pP + 1 )
                          . $no_parsed_prAlt_Start;
                        my $ext_length = ( length($prAlt) - $prSimpleSame );
                        if ( $ext_length > 0 or $prAlt !~ /\*$/ ) {
                            $trannoEnt->{p} .= ($start_in_cds_flag and $frameshift_flag) ? 'fs*' : 'ext*';
                            if ( $prAlt =~ /\*$/ ) {    # ext length estimated
                                $trannoEnt->{p} .= ($start_in_cds_flag and $frameshift_flag) ? $ext_length : ($ext_length - 1);
                            }
                            else {    # don't meet a stop codon
                                $trannoEnt->{p} .= '?';
                            }
                        }
                        else {
                            $trannoEnt->{func} = 'stop-retained';
                            $trannoEnt->{p} = 'p.(*' . ( $no_parsed_pP + 1 ) . '=)';
                        }
                        $trannoEnt->{protBegin} += $prSimpleSame;
                        next;
                    }

                    my $ins_stop_tag = 0;
                    if ( $prRef !~ /\*/ and $prAlt ne '*' and $prAlt =~ /\*$/ )
                    {                 # stop gain
                        $prEnd = $trdbEnt->{plen};   # extent to the end of prot
                        $prRef = substr(
                            $trdbEnt->{pseq},
                            ( $prBegin - 1 ),
                            ( $prEnd - $prBegin + 1 )
                        );
                        $prAlt =~ s/\*$//;
                        $ins_stop_tag = 1;
                    }

                    if ( $prAlt eq '?' ) {
                        $trannoEnt->{func} = 'unknown-no-call';
                        $trannoEnt->{p} = 'p.';
                        if ($prBegin == $prEnd) {
                            $trannoEnt->{p} .= $prRef . $prBegin . '?';
                        }
                        elsif ($prBegin < $prEnd) {
                            $trannoEnt->{p} .= substr($prRef,0,1) . $prBegin 
                                            . '_' . substr($prRef,-1,1) . $prEnd . 'delins?';
                        }
                        else {
                            $trannoEnt->{p} .= substr($trdbEnt->{pseq}, $prEnd - 1, 1) . $prEnd
                                            . '_' . substr($trdbEnt->{pseq}, $prBegin - 1, 1) . $prBegin . 'ins?';
                        }
                        next;
                    }

                    # parse the protein variants
                    # to recognize the repeat and adjust to correct position
                    my $prVar =
                      $self->prWalker( $trdbEnt->{pseq}, $prBegin, $prEnd,
                        $prRef, $prAlt );

                    # 0-based
                    my ( $p_P, $p_r, $p_a, $prl, $pal ) =
                      $prVar->getUnifiedVar('+');

                    my $prStart = substr( $trdbEnt->{pseq}, $p_P, 1 );
                    my $prStop =
                      substr( $trdbEnt->{pseq}, ( $p_P + $prl - 1 ), 1 )
                      if ( $p_P > 0 or $prl > 0 );

                    # single substitution
                    if ( $prVar->{imp} eq 'snv' ) {
                        if ( $p_r eq 'X' ) {
                            $trannoEnt->{func} = 'abnormal-inseq-stop';
                        }
                        elsif ( $p_r eq '.' ) {    # N in refseq transcript
                            $trannoEnt->{func} = 'unknown';
                        }
                        elsif ( $p_a eq '?' ) {    # substitution with N
                            $trannoEnt->{func} = 'unknown-no-call';
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
                        $trannoEnt->{p} = 'p.' . $p_r . ( $p_P + 1 ) . $p_a;
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

                        if ( $prVar->{replen} == 1 ) {
                            $trannoEnt->{p} =
                              'p.' . $prVar->{rep} . ( $p_P + 1 );
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

                        # add a new key alt_pHGVS for querying
                        # previous database, do exactly like alt_cHGVS
                        # but will be much more simple.

                        my $phgvs_5 = $p_P + 1;
                        my $phgvs_3 =
                          $p_P + ( $prVar->{ref_cn} * $prVar->{replen} );
                        my $repPrSta = substr( $prVar->{rep}, 0,  1 );
                        my $repPrEnd = substr( $prVar->{rep}, -1, 1 );
                        if ( $prVar->{ref_cn} > $prVar->{alt_cn} )
                        {    # deletion or duplication in previous definition
                            if ( ( $prVar->{ref_cn} - $prVar->{alt_cn} ) == 1
                                and $prVar->{replen} == 1 )
                            {
                                $trannoEnt->{alt_pHGVS} =
                                  'p.' . $prVar->{rep} . $phgvs_3 . 'del';
                            }
                            else {
                                $trannoEnt->{alt_pHGVS} = 'p.'
                                  . $repPrSta
                                  . ( $phgvs_5 +
                                      $prVar->{alt_cn} * $prVar->{replen} )
                                  . '_'
                                  . $repPrEnd
                                  . $phgvs_3 . 'del';
                            }
                        }
                        else {    # insertion
                            my $ins_cn   = $prVar->{alt_cn} - $prVar->{ref_cn};
                            my $ins_cont = $prVar->{rep} x $ins_cn;
                            my $prPost =
                              substr( $trdbEnt->{pseq}, $phgvs_3, 1 );
                            $trannoEnt->{alt_pHGVS} =
                                'p.'
                              . $repPrEnd
                              . $phgvs_3 . '_'
                              . $prPost
                              . ( $phgvs_3 + 1 ) . 'ins'
                              . $ins_cont;
                        }

                        next;
                    }

                    # when inserted bases contains stop codon but not framshift
                    if ( $p_a =~ /^\*$/ ) {    # prAlt is only a stop codon
                        $trannoEnt->{p} = 'p.' . $prStart . ( $p_P + 1 ) . '*';
                        $trannoEnt->{func} = 'stop-gain';
                        next;
                    }

                    # insertion
                    if ( $prVar->{sm} == 0 ) {

                        # insert into the 5'edge of start codon
                        # this may caused by a insertion or delins
                        # around init-codon to recreat a new init-codon
                        # in the altered sequence
                        if ( $p_a eq '?' ) {    # substitution with N
                            $trannoEnt->{func} = 'unknown-no-call';
                        }
                        elsif ( $p_P == 0 ) {

                            # altstart don't fix this
                            # use pHGVS to indicate
                            # altstart here
                            #
                            # hgvs 2.1511 don't refer to this kind of 
                            # coding synon, so just keep it.
                            $trannoEnt->{p} = 'p.(=)';

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
                            $trannoEnt->{func} =
                              ($ins_stop_tag) ? 'stop-gain' : 'cds-ins';
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
                        if ( $p_a eq '?' ) {    # substitution with N
                            $trannoEnt->{func} = 'unknown-no-call';
                        }
                        elsif ( $p_r =~ /\*/ ) {
                            $trannoEnt->{func} = 'stop-loss';
                        }
                        elsif ($hit_stop_flag) {
                            $trannoEnt->{func} = 'stop-retained';
                        }
                        elsif ($ins_stop_tag) {
                            $trannoEnt->{func} = 'stop-gain';
                        }
                        elsif ( $pal == 0 ) {
                            $trannoEnt->{func} = 'cds-del';
                        }
                        elsif ( $pal == $prl ) {
                            $trannoEnt->{func} = 'missense';
                        }
                        else {
                            $trannoEnt->{func} = 'cds-indel';
                        }

                        $trannoEnt->{p} = 'p.' . $prStart . ( $p_P + 1 );
                        if ($ins_stop_tag) {
                            if ($p_a eq "") {
                                $trannoEnt->{p} .= '*'
                            }
                            else {
                                my $Ptmp = substr($trdbEnt->{pseq}, $p_P + $pal, 1);
                                my $PendTmp = $p_P + $pal + 1;
                                $trannoEnt->{p} .= '_' . $Ptmp . $PendTmp . 'delins' . $p_a . '*';
                            }
                            next;
                        } elsif ( $prVar->{sm} == 1 ) {
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

=head2 prWalker

    About   : act as trWalker, but on protein sequence
    Usage   : $prVar = $beda->prWalker( $prSeq,  $prBegin, $prEnd, $prRef, $prAlt );
    Args    : prSeq   - whole protein sequence
	      prBegin - protein variant Begin position (1 based)
	      prEnd   - protein variant End position (1 based)
	      prRef   - protein variant reference sequence
	      prAlt   - protein variant alt sequence

=cut

sub prWalker {
    my ( $self, $prSeq, $prBegin, $prEnd, $prRef, $prAlt ) = @_;
    my $prVar =
      BedAnno::Var->new( "nouse", ( $prBegin - 1 ), $prEnd, $prRef, $prAlt );
    my ( $pr_P, $pr_r, $pr_a, $pr_rl, $pr_al ) = $prVar->getUnifiedVar('+');

    if (   ( defined $pr_al and $pr_rl == $pr_al )
        or ( 0 != index( $pr_r, $pr_a ) and 0 != index( $pr_a, $pr_r ) ) )
    {
        return $prVar;
    }

    $prBegin = $pr_P + 1;
    $prEnd   = $pr_P + $pr_rl;

    my ( $ref_sta, $ref_sto, $renew_prRef, $renew_prAlt ) =
      walker( $prBegin, $prEnd, $prSeq, $pr_r, $pr_a, $pr_rl, $pr_al );

    if ( $prBegin != $ref_sta or $prEnd != $ref_sto ) {
        $prVar = BedAnno::Var->new( "nouse", ( $ref_sta - 1 ),
            $ref_sto, $renew_prRef, $renew_prAlt );
    }

    return $prVar;
}

=head2 trWalker

    About   : walk around the variant position to find possible
              repeat start/end, and return the recalculated
              trBegin and trEnd, together with the real
              transcript originated variants and unified property
              Current implementation won't walk around in the following
              cases:
              1. no-call
              2. annotation-fail
              3. snv or mnp
              4. non-exon region
              5. span case or any edge case
              6. delins without repeat.
    Usage   : ($trBegin, $trEnd, $real_var) = 
              $beda->trWalker($tid, $rtrinfo);
    Args    : trAcc and hash ref of trInfo annotation for trAcc,
          

=cut

sub trWalker {
    my ( $self, $tid, $rtrinfo ) = @_;

    # reparse the transcript originated var
    my $real_var = BedAnno::Var->new( $tid, 0, length( $rtrinfo->{trRef} ),
        $rtrinfo->{trRef}, $rtrinfo->{trAlt} );
    my @Unified = $real_var->getUnifiedVar('+');
    my ( $real_p, $real_r, $real_a, $real_rl, $real_al ) = @Unified;

    my $trBegin = reCalTrPos_by_ofst( $rtrinfo, $real_p );
    my $trEnd = reCalTrPos_by_ofst( $rtrinfo, ( $real_p + $real_rl - 1 ) );

    if (
        # skip snv and mnp
        ( $real_rl == $real_al )

        # skip span case or any edge case
        or ( $rtrinfo->{ei_Begin} ne $rtrinfo->{ei_End} )

        # skip non exon begin/end position
        or ( $trBegin !~ /^\d+$/ or $trEnd !~ /^\d+$/ )

        # skip delins without repeat
        or (    0 != index( $real_r, $real_a )
            and 0 != index( $real_a, $real_r ) )
      )
    {
        # no correction
        return ( $trBegin, $trEnd, $real_var, \@Unified );
    }

    my $qtid = $tid;
    $qtid =~ s/\-\d+$//;
    my $trSeq = $self->{trInfodb}->{$qtid}->{seq};
    my ( $ref_sta, $ref_sto, $trRef, $trAlt ) =
      walker( $trBegin, $trEnd, $trSeq, $real_r, $real_a, $real_rl, $real_al );

    if (   $ref_sta ne $trBegin
        or $ref_sto ne $trEnd )
    {
        $trBegin = $ref_sta;
        $trEnd   = $ref_sto;

        $real_var =
          BedAnno::Var->new( $tid, 0, length($trRef), $trRef, $trAlt );
        @Unified = $real_var->getUnifiedVar('+');
    }

    return ( $trBegin, $trEnd, $real_var, \@Unified );
}

sub walker {
    my ( $ref_sta, $ref_sto, $whole_seq, $ref, $alt, $reflen, $altlen ) = @_;
    my $seqlen = length($whole_seq);

    my $ori_walker;
    my $track_opt;

    if ( $reflen < $altlen ) {
        $ref_sta    = $ref_sto + 1;
        $ori_walker = substr( $alt, $reflen );
        $track_opt  = 0;                         # walk on alt
    }
    else {
        $ref_sta += $altlen;
        $ori_walker = substr( $ref, $altlen );
        $track_opt  = 1;                         # walk on ref
    }

    # walk to 3' most
    my @cur_walker = split( //, $ori_walker );
    my $walk_forward_step = 0;
    for ( my $p = $ref_sto ; $p < $seqlen ; $p++ ) {
        if ( substr( $whole_seq, $p, 1 ) eq $cur_walker[0] ) {
            push( @cur_walker, shift(@cur_walker) );
            if ( $p == $seqlen - 1 ) {
                $walk_forward_step = $seqlen - $ref_sto;
                last;
            }
        }
        else {
            $walk_forward_step = $p - $ref_sto;
            last;
        }
    }

    $ref_sto += $walk_forward_step;    # the final end on reference track

    # use the 3' most element as the final matched difference
    my $match_target = join( "", @cur_walker );

    my $forward_footprint =
      substr( $whole_seq, $ref_sta - 1, $walk_forward_step );

    # walk to 5' most
    @cur_walker = split( //, $ori_walker );
    my $walk_back_step = 0;
    for ( my $q = $ref_sta - 1 ; $q > 0 ; $q-- ) {
        if ( substr( $whole_seq, $q - 1, 1 ) eq $cur_walker[-1] ) {
            unshift( @cur_walker, pop(@cur_walker) );
            if ( $q == 1 ) {
                $walk_back_step = $ref_sta - 1;
                last;
            }
        }
        else {
            $walk_back_step = $ref_sta - $q - 1;
            last;
        }
    }

    my $back_most = $ref_sta - $walk_back_step;
    my $backward_footprint =
      substr( $whole_seq, $ref_sta - $walk_back_step - 1, $walk_back_step );

    my $cur_short_track = $backward_footprint . $forward_footprint;
    my $cur_long_track  = $cur_short_track . $match_target;
    my $target_sta      = index( $cur_long_track, $match_target );
    $ref_sta = $back_most + $target_sta;

    my ( $renew_ref, $renew_alt );
    if ($track_opt) {
        $renew_ref = substr( $cur_long_track,  $target_sta );
        $renew_alt = substr( $cur_short_track, $target_sta );
    }
    else {
        $renew_ref = substr( $cur_short_track, $target_sta );
        $renew_alt = substr( $cur_long_track,  $target_sta );
    }

    return ( $ref_sta, $ref_sto, $renew_ref, $renew_alt );
}

=head2 cmpPos

    About   : judge the order for p1 and p2, because the insertion
              will have a reverted order of left & right position
    Usage   : my $cmpRst = BedAnno->cmpPos($p1, $p2);
    Args    : hgvs positio p1 and p2, with out 'c.' or 'n.' flag
    Return  : 0 for same, 1 for normal order, -1 for reversed order.

=cut

sub cmpPos {
    my $self = shift;
    my ( $p1, $p2 ) = @_;
    return 0 if ( $p1 eq $p2 );
    my ( $s1, $s2, $anc1, $anc2, $int_s1, $int_s2, $ofst1, $ofst2 );

    my %order_s = (
        '-u' => 1,
        '-'  => 2,
        ''   => 3,
        '+'  => 4,
        '*'  => 5,
        '+d' => 6,
    );

    if ( $p1 =~ /^([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$/ ) {
        $s1     = $1;
        $anc1   = $2;
        $int_s1 = $3;
        $ofst1  = $4;
    }

    if ( $p2 =~ /^([\-\+\*]?)(\d+)([\+\-]?[ud]?)(\d*)$/ ) {
        $s2     = $1;
        $anc2   = $2;
        $int_s2 = $3;
        $ofst2  = $4;
    }

    if ( $order_s{$s1} < $order_s{$s2} ) {
        return 1;
    }
    elsif ( $order_s{$s1} > $order_s{$s2} ) {
        return -1;
    }
    else {    # $s1 eq $s2
        if ( $s1 eq '-' ) {
            if ( $anc1 > $anc2 ) {
                return 1;
            }
            elsif ( $anc1 < $anc2 ) {
                return -1;
            }
        }
        else {
            if ( $anc1 < $anc2 ) {
                return 1;
            }
            elsif ( $anc1 > $anc2 ) {
                return -1;
            }
        }

        # anc1 eq anc2
        if ( $order_s{$int_s1} < $order_s{$int_s2} ) {
            return 1;
        }
        elsif ( $order_s{$int_s1} > $order_s{$int_s2} ) {
            return -1;
        }
        else {    # int_s2 same
            if ( $int_s1 =~ /\-/ ) {
                if ( $ofst1 > $ofst2 ) {
                    return 1;
                }
                elsif ( $ofst1 < $ofst2 ) {
                    return -1;
                }
            }
            else {
                if ( $ofst1 < $ofst2 ) {
                    return 1;
                }
                elsif ( $ofst1 > $ofst2 ) {
                    return -1;
                }
            }
        }
    }
    return 0;
}

sub _getAAandPolar {
    my $codon       = shift;
    my $rtrans_opts = shift;
    my ( $aa, $zero ) = translate( $codon, $rtrans_opts );
    my $polar =
      ( exists $Polar{$aa} ) ? $Polar{$aa} : ( ( $codon =~ /N/ ) ? '?' : '.' );
    return ( $aa, $polar );
}

# change the 1based nDot position into cDot format
sub _cPosMark {

    #	1based	0-based 1based
    my ( $trPos, $csta, $csto, $trlen ) = @_;
    return $trPos if ( $trPos eq "" );

    my $cDot;
    if ( $trPos =~ /^([\+\-]?)(\d+)([\+\-]?)(\d*)$/ ) {
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
            $cDot = (
                  ( $csto == $trlen )
                ? ( $csto - $csta )
                : '*' . ( $trlen - $csto )
              )
              . '+d'
              . $anchor;
        }
        elsif ( $s2 =~ /[\+\-]/ ) {
            my $cDot_tmp = _cPosMark( $anchor, $csta, $csto, $trlen );
            $cDot = $cDot_tmp . $s2 . $ofst;
        }
        else {
            confess "Error: unmatched nDot string: [$trPos]";
        }
    }
    return $cDot;
}

sub _getCPosFrame {
    my ( $trdbEnt, $p ) = @_;
    return (-1) if ( !exists $trdbEnt->{csta} or $p !~ /^\d+$/ );
    my $cds_p = $p - $trdbEnt->{csta};
    return _getCPosFrame_by_cdsPos( $trdbEnt, $cds_p );
}

=head2 _getCPosFrame_by_cdsPos

    About   : calculate codon position and frame info by
              involving frameshift case in consideration.
    Usage   : my @cPosFrame = _getCPosFrame_by_cdsPos( $trdbEnt, $cds_p );
    Args    : cds_p is relative position to the first bp in start codon on trans
    Returns : an array in the format of:
              ( $stat, $ra_posrefs1, $ra_posrefs2 )
              ra_posrefs* is an array ref for the following array
              [ $codonPos, $frame ]

              for non-cds region case stat is -1, without ra_posrefs.
              for normal case stat is 1, 
                  ra_posrefs1 is the information of input position,
                  without ra_posrefs2 specified.
              for deleted frame case stat is 0,
                  ra_posrefs1 is the information of position on 5' 
                  side next to the deleted frame, and ra_posrefs2
                  is on 3' side.
              for duplicated frame case stat is 2,
                  ra_posrefs1 is the 5' most hit of input position,
                  ra_posrefs2 is the corresponding position on the 
                  3' duplicated region.

=cut

# give codon number and frame with involving frameshift case
sub _getCPosFrame_by_cdsPos {
    my ( $trdbEnt, $p ) = @_;
    if (   !exists $trdbEnt->{csta}
        or $trdbEnt->{csta} eq '.'
        or $p !~ /^\d+$/
        or $p <= 0
        or $p > ( $trdbEnt->{csto} - $trdbEnt->{csta} ) )
    {
        return (-1);
    }
    if ( exists $trdbEnt->{cfs} ) {
        my $total_fs = 0;
        foreach my $fsld ( sort { $a <=> $b } keys %{ $trdbEnt->{cfs} } ) {
            if ( $fsld < $p and ( $fsld + $trdbEnt->{cfs}->{$fsld} ) < $p )
            {    # normal
                $total_fs += $trdbEnt->{cfs}->{$fsld};
            }
            elsif ( $fsld < $p and $fsld + $trdbEnt->{cfs}->{$fsld} >= $p )
            {    # in deleted frames
                my $codon_pre = $fsld - $total_fs;
                my ( $pP_pre, $frame_pre ) = _calPosFrame($codon_pre);
                my ( $pP_lat, $frame_lat );
                if ( $frame_pre < 2 ) {
                    $pP_lat    = $pP_pre;
                    $frame_lat = $frame_pre + 1;
                }
                else {
                    $pP_lat    = $pP_pre + 1;
                    $frame_lat = 0;
                }
                return ( 0, [ $pP_pre, $frame_pre ], [ $pP_lat, $frame_lat ] );
            }
            elsif ( $fsld >= $p and $fsld + $trdbEnt->{cfs}->{$fsld} < $p )
            {    # in dup frames
                my $codon_1 = $p - $total_fs;
                my $codon_2 = $p - $total_fs + abs( $trdbEnt->{cfs}->{$fsld} );
                my ( $pP_1, $frame_1 ) = _calPosFrame($codon_1);
                my ( $pP_2, $frame_2 ) = _calPosFrame($codon_2);
                return ( 2, [ $pP_1, $frame_1 ], [ $pP_2, $frame_2 ] );
            }
            else {    # after the current position
                last;
            }
        }
        $p -= $total_fs;
    }
    return ( 1, [ ( _calPosFrame($p) ) ] );
}

# calculate codon number and frame directly be cds position
sub _calPosFrame {
    my $cds_p = shift;
    my $pP    = int( $cds_p / 3 );
    if ( $cds_p % 3 > 0 ) {
        $pP++;
    }
    my $frame = 2 - ( $pP * 3 - $cds_p );
    return ( $pP, $frame );
}

# generate codon together with AA and Polar by codon number
sub _genCodonInfo {
    my ( $trdbEnt, $pP ) = @_;

    $trdbEnt->{cseq} = _getCodingSeq($trdbEnt) if ( !exists $trdbEnt->{cseq} );
    if ( $pP <= 0 or ( $pP - 1 ) > $trdbEnt->{plen} ) {
        confess "Error: [$trdbEnt->{gene}] no codon info for $pP position.";
    }
    my $codon = substr( $trdbEnt->{cseq}, ( $pP - 1 ) * 3, 3 );
    my $rtrans_opts = {};

    $rtrans_opts->{mito} = 1 if ( $trdbEnt->{gene} =~ /^MT\-/ );    # chrM gene

    if ( $pP <= $trdbEnt->{plen} ) {    # plen not involve terminal codon
        $rtrans_opts->{nostop} = 1;
    }
    elsif ( $pP == $trdbEnt->{plen} + 1 ) {    # terminal with polyA complement
        $codon .= 'A' x ( 3 - length($codon) ) if ( exists $trdbEnt->{A} );
    }

    return ( $codon, _getAAandPolar( $codon, $rtrans_opts ) );
}

=head2 _cdsubstr

    About   : substr from transcript seqeunce involving frameshift changing.
    Usage   : my $codonStr = _cdsubstr( $trdbEnt, $start, $length, $replace );
    Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto.
              start is 0 based start position of sub string on transcript seq.
              length is the length of target region on transcript seq.
              replace is the alternative.
    Returns : A substring cut from the transcript sequence, with frameshift involved.
              Behave like function "substr", but for replace mode, 
              it returns the whole transcript seq after replacement.

=cut

sub _cdsubstr {
    my ( $trdbEnt, $start, $len, $replace ) = @_;
    if ( !defined $len ) {
        $len = $trdbEnt->{len} - $start;
    }
    if ( !exists $trdbEnt->{csta} or !exists $trdbEnt->{nfs} ) {
        if ( defined $replace ) {
            my $trSeq = $trdbEnt->{seq};
            substr( $trSeq, $start, $len, $replace );
            return $trSeq;
        }
        else {
            return substr( $trdbEnt->{seq}, $start, $len );
        }
    }
    else {
        if ( $start < 0 ) {
            $start = $trdbEnt->{len} + $start;
        }
        if ( $len < 0 ) {
            $len = $trdbEnt->{len} + $len - $start;
        }
        my @start_ofst = _calfsOfst( $trdbEnt, $start );
        my @stop_ofst = _calfsOfst( $trdbEnt, ( $start + $len - 1 ) );

        $trdbEnt->{cseq} = _getCodingSeq($trdbEnt)
          if ( !exists $trdbEnt->{cseq} );
        my $new_wholeseq =
            substr( $trdbEnt->{seq}, 0, $trdbEnt->{csta} )
          . $trdbEnt->{cseq}
          . substr( $trdbEnt->{seq}, $trdbEnt->{csto} );

        my $returned_seq;
        my $max_len = -1;
        foreach my $sta_ofst (@start_ofst) {
            my $new_start = $start - $sta_ofst;
            foreach my $sto_ofst (@stop_ofst) {
                my $new_stop = $start + $len - 1 - $sto_ofst;
                next if ( $new_start > $new_stop );
                $new_stop += 1;
                my $new_len = $new_stop - $new_start;
                if ( $new_len > $max_len ) {
                    if ( defined $replace ) {
                        $returned_seq = $new_wholeseq;
                        substr( $returned_seq, $new_start, $new_len, $replace );
                        $max_len = $new_len;
                    }
                    else {
                        $returned_seq =
                          substr( $new_wholeseq, $new_start, $new_len );
                        $max_len = $new_len;
                    }
                }
            }
        }

        # only return the longest substring
        return $returned_seq;
    }
}

sub _calfsOfst {
    my ( $trdbEnt, $p ) = @_;    # 0 based
    if ( !exists $trdbEnt->{csta} or !exists $trdbEnt->{nfs} ) {
        return 0;
    }

    my $ofst = 0;
    foreach my $fsld ( sort { $a <=> $b } keys %{ $trdbEnt->{nfs} } ) {
        my $fsEd = $fsld + $trdbEnt->{nfs}->{$fsld};
        if ( $fsld <= $p and $fsEd < $p ) {    # before p
            $ofst += $trdbEnt->{nfs}->{$fsld};
        }
        elsif ( $fsld <= $p and $fsEd >= $p ) {
            $ofst += ( $fsld + $trdbEnt->{nfs}->{$fsld} - $p );
            return ($ofst);
        }
        elsif ( $fsld > $p and $fsEd <= $p ) {
            return ( $ofst, ( $ofst + $trdbEnt->{nfs}->{$fsld} ) );
        }
        else {
            last;
        }
    }
    return ($ofst);
}

=head2 _getPairedCodon

    About   : get codon position, codon string, aa string, and frame info
              for a pair of transcript position
    Usage   : my ($rcinfo5, $rcinfo3) = _getPairedCodon( $trdbEnt, $p5, $p3 );
    Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto
              for cds start/end position, besides many other feature tags.
              p5 and p3 is a pair of positions which give a region or single
              position or maybe an insertion anchor.
    Returns : codon info array ref for p5 and p3. The array ref:
              [
                  AA position  - 0 for not in cds region.
                  codon string - codon string, e.g. "ATA".
                  aa char      - AA code, 1 bp mode, e.g. "E".
                  polar        - Polar properties, e.g. "P+".
                  frame        - current position's frame info,
                                 -1 for not available or in fs-site.
		  [frame-alt]  - frame info around fs-site.
              ]

=cut

sub _getPairedCodon {
    my ( $trdbEnt, $p5, $p3 ) = @_;
    my @c5 = _getCPosFrame( $trdbEnt, $p5 );
    my @c3 = _getCPosFrame( $trdbEnt, $p3 );
    if ( $c5[0] == 1 and $c3[0] == 1 ) {    # all normal
        if ( exists $trdbEnt->{nfs} ) {
            my $contain_fs = 0;
            foreach my $fsld ( keys %{ $trdbEnt->{nfs} } ) {
                if ( $p5 <= $fsld and $p3 >= $fsld + $trdbEnt->{nfs}->{$fsld} )
                {
                    $contain_fs = 1;
                    last;
                }
            }
            if ( $contain_fs == 1 ) {
                return (
                    [
                        $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                        -1, $c5[1]->[1]
                    ],
                    [
                        $c3[1]->[0], ( _genCodonInfo( $trdbEnt, $c3[1]->[0] ) ),
                        -1, $c3[1]->[1]
                    ]
                );
            }
        }
        return (
            [
                $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                $c5[1]->[1]
            ],
            [
                $c3[1]->[0], ( _genCodonInfo( $trdbEnt, $c3[1]->[0] ) ),
                $c3[1]->[1]
            ]
        );
    }
    elsif ( $c5[0] == 0 and $c3[0] == 0 ) {    # all in deleted frame
        if ( $c5[1]->[0] == $c3[1]->[0] ) {    # same frame
            return (
                [ $c5[2]->[0], ( _genCodonInfo( $trdbEnt, $c5[2]->[0] ) ), -1 ],
                [ $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ), -1 ]
            );
        }
        else {                                 # different frames
            return (
                [
                    $c5[2]->[0], ( _genCodonInfo( $trdbEnt, $c5[2]->[0] ) ),
                    -1, $c5[2]->[1]
                ],
                [
                    $c3[1]->[0], ( _genCodonInfo( $trdbEnt, $c3[1]->[0] ) ),
                    -1, $c3[1]->[1]
                ],
            );
        }
    }
    elsif ( $c5[0] == 2 and $c3[0] == 2 ) {    # all in duplicated frame
        if ( _chkCodonPosRel( \@c5, \@c3 ) ) {    # connected
            return (
                [
                    $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                    -1, $c5[1]->[1]
                ],
                [
                    $c3[2]->[0], ( _genCodonInfo( $trdbEnt, $c3[2]->[0] ) ),
                    -1, $c3[2]->[1]
                ]
            );
        }
        else {    # non-connected, this case don't exists in current db,
                  # we will only take the first two returned value
                  # and give up the latter two.
            return (
                [
                    $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                    -1, $c5[1]->[1]
                ],
                [
                    $c3[1]->[0], ( _genCodonInfo( $trdbEnt, $c3[1]->[0] ) ),
                    -1, $c3[1]->[1]
                ],

                [
                    $c5[2]->[0], ( _genCodonInfo( $trdbEnt, $c5[2]->[0] ) ),
                    -1, $c5[2]->[1]
                ],
                [
                    $c3[2]->[0], ( _genCodonInfo( $trdbEnt, $c3[2]->[0] ) ),
                    -1, $c3[2]->[1]
                ]
            );
        }
    }
    else {
        my ( $r5_ret, $r3_ret );
        if ( $c5[0] == 1 ) {
            $r5_ret = [
                $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                $c5[1]->[1]
            ];
        }
        elsif ( $c5[0] < 0 ) {
            $r5_ret = [ 0, "", "", "", -1 ];
        }
        elsif ( $c5[0] == 0 ) {
            $r5_ret = [
                $c5[2]->[0], ( _genCodonInfo( $trdbEnt, $c5[2]->[0] ) ),
                -1, $c5[2]->[1]
            ];
        }
        elsif ( $c5[0] == 2 ) {
            $r5_ret = [
                $c5[1]->[0], ( _genCodonInfo( $trdbEnt, $c5[1]->[0] ) ),
                -1, $c5[1]->[1]
            ];
        }

        if ( $c3[0] == 1 ) {
            $r3_ret = [
                $c3[-1]->[0], ( _genCodonInfo( $trdbEnt, $c3[-1]->[0] ) ),
                $c3[-1]->[1]
            ];
        }
        elsif ( $c3[0] < 0 ) {
            $r3_ret = [ 0, "", "", "", -1 ];
        }
        elsif ( $c3[0] == 0 ) {
            $r3_ret = [
                $c3[1]->[0], ( _genCodonInfo( $trdbEnt, $c3[1]->[0] ) ),
                $c3[1]->[1]
            ];
        }
        elsif ( $c3[0] == 2 ) {
            $r3_ret = [
                $c3[2]->[0], ( _genCodonInfo( $trdbEnt, $c3[2]->[0] ) ),
                $c3[2]->[1]
            ];
        }

        return ( $r5_ret, $r3_ret );
    }
}

# check if region in duplicated frames are connected
sub _chkCodonPosRel {
    my ( $rsta, $rsto ) = @_;
    if (    $rsta->[2]->[0] == $rsto->[1]->[0]
        and $rsta->[2]->[1] <= $rsto->[1]->[1] + 1 )
    {
        # same codon connected frame
        return 1;
    }
    elsif ( $rsta->[2]->[0] == $rsto->[1]->[0] + 1
        and $rsta->[2]->[1] == 0
        and $rsto->[1]->[1] == 2 )
    {
        # neighbor codon consecutive frame
        return 1;
    }
    elsif ( $rsta->[2]->[0] < $rsto->[1]->[0] ) {

        # overlapped cases
        return 1;
    }
    else {
        return 0;
    }
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
    my $nucl   = shift;
    my $ropt   = shift if (@_);
    my $mito   = ( defined $ropt and exists $ropt->{mito} ) ? 1 : 0;
    my $nostop = ( defined $ropt and exists $ropt->{nostop} ) ? 1 : 0;
    my $polyA  = ( defined $ropt and exists $ropt->{polyA} ) ? 1 : 0;

    my $lenN       = length($nucl);
    my $frame_next = $lenN % 3;
    if ( $frame_next != 0 and $polyA ) {
        $nucl .= 'A' x ( 3 - $frame_next );
    }
    $nucl = uc($nucl);
    my $prot = "";
    for ( my $i = 0 ; $i < length($nucl) ; $i += 3 ) {
        my $codon = substr( $nucl, $i, 3 );
        my $aa =
          ( exists $C1{$codon} )
          ? $C1{$codon}
          : ( ( $codon =~ /N/ ) ? '?' : '.' );
        if ( $mito and $codon eq 'ATA' ) {
            $aa = 'M';
        }
        if ( $aa eq '*' and $mito and $codon eq 'TGA' ) {
            $aa = 'W';
        }
        elsif ( $aa eq '*' and $nostop ) {
            $aa = ( $codon eq 'TGA' ) ? 'U' : 'X';
        }
        elsif ( $aa eq '*' ) {
            $prot .= $aa;
            return ( $prot, 0 );
        }
        $prot .= $aa;
    }
    $prot =~ s/[UX]$/*/g if ($nostop);
    $prot =~ s/\.$// if ( $frame_next != 0 );
    return ( $prot, $frame_next );
}

# new var generated from new() in BedAnnoVar the real transcript var
# will bring some new transcript start/end information
# here will recalculate it from the offset.
sub reCalTrPos_by_ofst {
    my ( $trannoEnt, $trRef_ofst ) = @_;
    return $trannoEnt->{rnaBegin} if ( $trRef_ofst == 0 );
    return $trannoEnt->{preStart}->{nDot}
      if ( $trRef_ofst == -1 and exists $trannoEnt->{preStart} );

    my @tag_sort     = sort trRefSort keys %{ $trannoEnt->{trRefComp} };
    my $cumulate_len = 0;
    my ( $cur_ex_start, $intOfst );
    if ( $trannoEnt->{rnaBegin} =~ /^(\-?\d+)([\+\-]?\d*)$/ ) {
        $cur_ex_start = $1;
        $intOfst      = $2;
    }
    if ( $intOfst =~ /^\+/ ) {
        $cur_ex_start += 1;
    }

    if ( $cur_ex_start =~ /^-/ ) {
        $intOfst      = $cur_ex_start;
        $cur_ex_start = 1;
    }
    if ( $intOfst eq '' ) {
        $intOfst = 0;
    }

    foreach my $exin (@tag_sort) {
        my $cur_blk_len;
        if ( $exin !~ /^EX/ ) {
            $cur_blk_len =
              ( 2 > @{ $trannoEnt->{trRefComp}->{$exin} } )
              ? 0
              : ( $trannoEnt->{trRefComp}->{$exin}->[1] -
                  $trannoEnt->{trRefComp}->{$exin}->[0] );
        }
        else {
            $cur_blk_len = $trannoEnt->{trRefComp}->{$exin};
        }

        my $cur_blk_ofst = ( $trRef_ofst - $cumulate_len );

        if ( $cur_blk_ofst >= $cur_blk_len ) {
            if ( $cur_blk_ofst == $cur_blk_len and $exin eq $tag_sort[-1] ) {
                return $trannoEnt->{postEnd}->{nDot};
            }
            else {
                $cumulate_len += $cur_blk_len;
                if ( $exin =~ /^EX/ ) {    # only cumulate ex pos
                    $cur_ex_start += $cur_blk_len;
                    $intOfst = 0;
                }
            }
        }
        else {                             # hit current block
            if ( $exin !~ /^EX/ ) {
                if ( $intOfst < 0 and $cur_ex_start == 1 ) {    # in promoter
                    confess
"Error: non-zero cumulate_len in promoter [$cumulate_len]."
                      if ( $cumulate_len > 0 );
                    return ( $intOfst + $cur_blk_ofst );
                }
                elsif ( $exin =~ /^Z/ ) {
                    return '+' . ( $cur_blk_ofst + 1 );
                }
                else {                                          # intron
                    my ( $blk_begin, $blk_end );
                    if ( $exin eq $tag_sort[0] ) {
                        $blk_begin = $trannoEnt->{rnaBegin};
                    }
                    else {
                        $blk_begin = ( $cur_ex_start - 1 ) . "+1";
                    }
                    if ( $exin eq $tag_sort[-1] ) {
                        $blk_end = $trannoEnt->{rnaEnd};
                    }
                    else {
                        $blk_end = $cur_ex_start . "-1";
                    }

                    return getIntrPos( $blk_begin, $blk_end, $cur_blk_ofst,
                        $cur_blk_len );
                }
            }
            else {
                return ( $cur_ex_start + $cur_blk_ofst );
            }
        }
    }
    confess "Error: out of trRefComp range [$trRef_ofst].";
}

# calculate position in intron region
sub getIntrPos {
    my ( $begin, $end, $ofst, $blk_len ) = @_;

    my ( $exPosBegin, $exPosEnd, $intPosBegin, $intPosEnd );
    if ( $begin =~ /^(\d+)([\+\-]\d+)$/ ) {
        $exPosBegin  = $1;
        $intPosBegin = $2;
    }
    if ( $end =~ /^(\d+)([\+\-]\d+)$/ ) {
        $exPosEnd  = $1;
        $intPosEnd = $2;
    }

    if ( $exPosBegin eq $exPosEnd ) {
        return
            $exPosBegin
          . ( ( $intPosBegin < 0 ) ? "" : '+' )
          . ( $intPosBegin + $ofst );
    }
    else {
        my $wlen  = $intPosBegin - $intPosEnd + $blk_len - 2;
        my $lofst = $intPosBegin + $ofst - 1;
        if ( $lofst < ( $wlen / 2 - 1 ) ) {
            return ( $exPosBegin . '+' . ( $lofst + 1 ) );
        }
        else {
            return ( $exPosEnd . '-' . ( $wlen - $lofst ) );
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
            return '=' if ( $refgenome eq '=' );
            my $int_seq =
              ( 2 == scalar @{ $trannoEnt->{trRefComp}->{$exin} } )
              ? substr(
                $refgenome,
                $trannoEnt->{trRefComp}->{$exin}->[0],
                (
                    $trannoEnt->{trRefComp}->{$exin}->[1] -
                      $trannoEnt->{trRefComp}->{$exin}->[0]
                )
              )
              : "";
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
    if ( $trStart =~ /^-/ ) {
        return 0;
    }
    elsif ( $trStart =~ /^(\d+)\-/ ) {
        return ( $1 - 1 );
    }
    elsif ( $trStart =~ /^(\d+)\+/ ) {
        return $1;
    }
    elsif ( $trStart =~ /^(\d+)$/ ) {
        return ( $1 - 1 );
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
    my ( $self, $rVars ) = @_;
    my @all_annoRst = ();
    my @sorted_vars =
      sort { $a->{chr} cmp $b->{chr} or $a->{pos} <=> $b->{pos} } @$rVars;

    my $cur_chr  = "";
    my $anno_idx = 0;
    foreach my $var (@sorted_vars) {
        if ( $var->{chr} ne $cur_chr ) {
            $cur_chr  = $var->{chr};
            $anno_idx = 0;
        }
        my $annoEnt;
        ( $annoEnt, $anno_idx ) = $self->varanno( $var, $anno_idx );
        push( @all_annoRst, $annoEnt );
    }
    return \@all_annoRst;
}

sub rev_comp {
    my $self = shift;
    my $Seq  = shift;
    $Seq = reverse($Seq);
    $Seq =~ tr/ATCG/TAGC/;    # only deal with 'A', 'T', 'G', 'C'
    return $Seq;
}

sub throw {
    my ( $self, @msg ) = @_;
    confess @msg, "\n";
}

sub warn {
    my ( $self, @msg ) = @_;
    carp @msg, "\n";
}

1;

=head1 DATABASE FORMAT

The Format of annotation database is listed as following:

=head2 BED FORMAT DATABASE

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

       >rnaAcc.ver rnaLen gene protAcc.ver protLen cdsSta,cdsEnd tags [altStartCodons] [frameshiftSite]
   
    "tags" are a string of multiple properties separated by "|", which are:

	<alignStat>[|altstart][|selenocysteine][|inseqStop][|polyATail]

    when "tags" contain "altstart", alternative start codons will be list in 
    "altStartCodons", separated by ";"

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

liutao E<lt>liut@geneplus.org.cnE<gt>

=head1 COPYRIGHT AND LICENSE

Please check LICENSE file for detail

=cut

