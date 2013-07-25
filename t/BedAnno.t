# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl BedAnno.t'

our $data;

BEGIN {
    unless ( grep /blib/, @INC ) {
        chdir 't' if -d 't';
        unshift @INC, '../lib' if -d '../lib';
        $data = '../data';
    }
}
$data ||= "data";

use Test::Most;
BEGIN { use_ok('BedAnno') }

test_tabix("$data/test_db.bed.gz");
unlink("$data/test_db.bed.gz.tbi");
my $beda = BedAnno->new(
    db     => "$data/hg19_refseq_anno.bed.gz",
    codon  => "$data/hg19_coding_seq.fas.gz",
    regbed => "$data/test.in.bed",
    trans  => "$data/trans.list"
);

#explain "The database are:", $$beda{annodb};

my $snv_anno = {
    c      => "c.82G>T",
    p      => "p.Glu28*",
    cc     => "GAA=>TAA",
    r      => "C1",
    exin   => "EX1",
    func   => "nonsense",
    polar  => "P-=>.",
    bc     => '.',
    flanks => {
        l    => "chr1:169454921-169454922",
        r    => "chr1:169454924-169454925",
        strd => '-'
    }
};

my $del1_anno = {
    c      => 'c.757delA',
    p      => 'p.Ile253Serfs*?',
    cc     => '.',
    r      => 'C1E',
    exin   => 'EX1E',
    func   => 'frameshift',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:248525637-248525638",
        r    => "chr1:248525640-248525641",
        strd => '+'
    }
};

my $splice_anno = {
    c      => 'c.1408C[7>6]',
    p      => '.',
    cc     => '.',
    r      => 'IC3',
    exin   => 'IVS3',
    func   => 'abnormal-intron',
    polar  => '.',
    bc     => 'c.1413+1delC',
    flanks => {
        l    => "chr8:24811063-24811064",
        r    => "chr8:24811066-24811067",
        strd => '-'
    }
};

my $splice1_anno = {
    c      => 'c.617-3_617-2dupTA',
    p      => '.',
    cc     => '.',
    r      => 'IC6',
    exin   => 'IVS7',
    func   => 'splice-3',
    polar  => '.',
    bc     => 'c.617-2_617-1insTA',
    flanks => {
        l    => "chr21:43803307-43803308",
        r    => "chr21:43803309-43803310",
        strd => '-'
    }
};

my $ins_anno = {
    c      => 'c.338_339insTGTTACGCAGGAGAC',
    p      => 'p.Thr113_Ala114insValThrGlnGluThr',
    cc     => '.',
    r      => 'C2',
    exin   => 'EX2',
    func   => 'cds-indel',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr3:195518111-195518112",
        r    => "chr3:195518113-195518114",
        strd => '-'
    }
};

my $intr_anno = {
    c      => 'c.83-12787_83-12786insTGTTACGCAGGAGAC',
    p      => '.',
    cc     => '.',
    r      => 'IC1',
    exin   => 'IVS1',
    func   => 'intron',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr3:195518111-195518112",
        r    => "chr3:195518113-195518114",
        strd => '-'
    }
};

my $delins_anno = {
    c      => 'c.1192_1195delCAACinsC',
    p      => 'p.Gln398_Arg399delinsArg',
    cc     => '.',
    r      => 'C5',
    exin   => 'EX7',
    func   => 'cds-indel',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr22:38119753-38119754",
        r    => "chr22:38119759-38119760",
        strd => '+'
    }
};

my $rep_anno = {
    c      => 'c.21GGC[6>4]',
    p      => 'p.Ala11_Ala13delinsAla',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX1',
    func   => 'cds-indel',
    polar  => '.',
    bc     => 'c.33_38delGGCGGC',
    flanks => {
        l    => "chr3:49395672-49395673",
        r    => "chr3:49395680-49395681",
        strd => '-'
    }
};

my $n_anno = {
    c      => 'c.2376+47_2376+48delCTinsCTTGGNCT',
    p      => '.',
    cc     => '.',
    r      => 'IC4',
    exin   => 'IVS5',
    func   => 'intron',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => 'chr8:61713129-61713130',
        r    => 'chr8:61713133-61713134',
        strd => '+'
    }
};

my $snv_5utr_anno = {
    c      => 'c.-1T>G',
    p      => '.',
    cc     => '.',
    r      => '5U1',
    exin   => 'EX2',
    func   => 'utr-5',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:861319-861320",
        r    => "chr1:861322-861323",
        strd => '+'
    }
};

my $snv_intron_3_anno = {
    c      => 'c.73-1G>T',
    p      => '.',
    cc     => '.',
    r      => 'IC1',
    exin   => 'IVS2',
    func   => 'splice-3',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:865532-865533",
        r    => "chr1:865535-865536",
        strd => '+'
    }
};

my $snv_intron_3_new_anno = {
    c      => 'c.73-2A>T',
    p      => '.',
    cc     => '.',
    r      => 'IC1',
    exin   => 'IVS2',
    func   => 'splice-3',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:865531-865532",
        r    => "chr1:865534-865535",
        strd => '+'
    }
};

my $snv_cdna_scoden_anno = {
    c      => 'c.1A>C',
    p      => 'p.0?',
    cc     => 'ATG=>CTG',
    r      => 'C1',
    exin   => 'EX2',
    func   => 'altstart',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861323-861324",
        strd => '+'
    }
};

my $snv_cdna_ecoden_anno = {
    c      => 'c.2046A>T',
    p      => 'p.*682Cysext*?',
    cc     => 'TGA=>TGT',
    r      => 'C13E',
    exin   => 'EX14E',
    func   => 'stop-loss',
    polar  => '.=>P0',
    bc     => '.',
    flanks => {
        l    => "chr1:879531-879532",
        r    => "chr1:879534-879535",
        strd => '+'
    }
};

my $del_cdna_scoden_anno = {
    c      => 'c.2_6delTGTCC',
    p      => 'p.0?',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX2',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:861321-861322",
        r    => "chr1:861328-861329",
        strd => '+'
    }
};
my $del_intron_anno = {
    c      => 'c.948+1_949-1delGCTG',
    p      => '.',
    cc     => '.',
    r      => 'IC5',
    exin   => 'IVS5',
    func   => 'splice',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr6_cox_hap2:2892856-2892857",
        r    => "chr6_cox_hap2:2892862-2892863",
        strd => '+'
    }
};

my $del_intron1_anno = {
    c      => 'c.946_949-1delGCTGCTG',
    p      => '.',
    cc     => '.',
    r      => 'C5-IC5',
    exin   => 'EX5-IVS5',
    func   => 'splice',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr6_cox_hap2:2892853-2892854",
        r    => "chr6_cox_hap2:2892862-2892863",
        strd => '+'
    }
};

my $del_intron2_anno = {
    c      => 'c.948_949delTGCTGG',
    p      => 'p.Ala316_Ala317delinsAlaSerfs*?',
    cc     => '.',
    r      => 'C5-C6',
    exin   => 'EX5-EX6',
    func   => 'frameshift',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr6_cox_hap2:2892855-2892856",
        r    => "chr6_cox_hap2:2892863-2892864",
        strd => '+'
    }
};

my $del_extron_anno = {
    c      => 'c.1_3delATG',
    p      => 'p.0?',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX1',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:948952-948953",
        r    => "chr1:948957-948958",
        strd => '+'
    }
};
my $del_extron2_anno = {
    c      => 'c.1525-2_1536+2delAGGTTCCCAAAGAGGT',
    p      => '.',
    cc     => '.',
    r      => 'IC11-IC12',
    exin   => 'IVS11-IVS12',
    func   => 'splice',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr15:42694318-42694319",
        r    => "chr15:42694336-42694337",
        strd => '+'
    }
};

my $del_5utr_coding_anno = {
    c      => 'c.-2_3delCCATG',
    p      => 'p.0?',
    cc     => '.',
    r      => '5U1E-C1',
    exin   => 'EX1',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:948950-948951",
        r    => "chr1:948957-948958",
        strd => '+'
    }
};

my $del_3utr_coding_anno = {
    c      => 'c.498_*2delAGG',
    p      => 'p.166del*',
    cc     => '.',
    r      => 'C2E-3U1E',
    exin   => 'EX2E',
    func   => 'stop-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:949856-949857",
        r    => "chr1:949861-949862",
        strd => '+'
    }
};

my $del_in_out_coding_anno = {
    c      => 'c.1_4delATGG',
    p      => 'p.0?',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX1',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr19:32083185-32083186",
        r    => "chr19:32083191-32083192",
        strd => '+'
    }
};

my $ins_intron_splice_anno = {
    c      => 'c.255-2_255-1insT',
    p      => '.',
    cc     => '.',
    r      => 'IC2',
    exin   => 'IVS3',
    func   => 'splice-3',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:866416-866417",
        r    => "chr1:866418-866419",
        strd => '+'
    }
};

my $ins_utr_coding_anno = {
    c      => 'c.-1_1insAG',
    p      => '.',
    cc     => '.',
    r      => '5U1-C1',
    exin   => 'EX2',
    func   => 'utr-5',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861322-861323",
        strd => '+'
    }
};

my $ins_intron_exon_anno = {
    c      => 'c.255-1_255insT',
    p      => 'p.Arg85fs*?',
    cc     => '.',
    r      => 'IC2-C3',
    exin   => 'IVS3-EX4',
    func   => 'frameshift',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:866417-866418",
        r    => "chr1:866419-866420",
        strd => '+'
    }
};
my $ins_intron_exon2_anno = {
    c      => 'c.255-1_255insTAA',
    p      => 'p.Arg85delinsSerLys',
    cc     => '.',
    r      => 'IC2-C3',
    exin   => 'IVS3-EX4',
    func   => 'cds-indel',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:866417-866418",
        r    => "chr1:866419-866420",
        strd => '+'
    }
};

my $ins_intron_exon3_anno = {
    c      => 'c.2377-1_2377insTAA',
    p      => 'p.Gln792_Gln792ins*',
    cc     => '.',
    r      => 'IC4-C5',
    exin   => 'IVS5-EX6',
    func   => 'stop-gain',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr8:61714085-61714086",
        r    => "chr8:61714087-61714088",
        strd => '+'
    }
};

my $delins_intron_splice_anno = {
    c      => 'c.255-1_256delGAGinsG',
    p      => '.',
    cc     => '.',
    r      => 'IC2-C3',
    exin   => 'IVS3-EX4',
    func   => 'splice',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:866416-866417",
        r    => "chr1:866421-866422",
        strd => '+'
    }
};

my $delins_utr_coding_anno = {
    c      => 'c.-1_2delTATinsAGCTA',
    p      => 'p.0?',
    cc     => '.',
    r      => '5U1-C1',
    exin   => 'EX2',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:861319-861320",
        r    => "chr1:861324-861325",
        strd => '+'
    }
};

my $delins_exon_intron_anno = {
    c      => 'c.306-1_306delGCinsCCA',
    p      => '.',
    cc     => '.',
    r      => 'IC3-C4',
    exin   => 'IVS4-EX5',
    func   => 'splice',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr1:871149-871150",
        r    => "chr1:871153-871154",
        strd => '+'
    }
};

my $delins_in_out_anno = {
    c      => '--1_c.3delCATGinsGAGA',
    p      => 'p.0?',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX1',
    func   => 'init-loss',
    polar  => '.',
    bc     => '.',
    flanks => {
        l    => "chr19:32083184-32083185",
        r    => "chr19:32083190-32083191",
        strd => '+'
    }
};

my $rep_insertion_dup_anno = {
    c      => 'c.38dupC',
    p      => 'p.Ala13fs*?',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX1',
    func   => 'frameshift',
    polar  => '.',
    bc     => 'c.38_39insC',
    flanks => {
        l    => "chr3:49395672-49395673",
        r    => "chr3:49395674-49395675",
        strd => '-'
    }
};

my $rep_insertion_dup1_anno = {
    c      => 'c.-1T[1>4]',
    p      => '.',
    cc     => '.',
    r      => '5U1-C1',
    exin   => 'EX2',
    func   => 'utr-5',
    polar  => '.',
    bc     => 'c.-1_1insTTT',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861322-861323",
        strd => '+'
    }
};

my $rep_insertion_dup_utr_coding_anno = {
    c      => 'c.-1TA[1>4]',
    p      => 'p.=',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX2',
    func   => 'coding-synon',
    polar  => '.',
    bc     => 'c.1_2insTATATA',
    flanks => {
        l    => "chr1:861321-861322",
        r    => "chr1:861323-861324",
        strd => '+'
    }
};

my $rep_insertion_dup_extron_intron_anno = {
    c      => 'c.1064GGT[1>4]',
    p      => '.',
    cc     => '.',
    r      => 'IC9',
    exin   => 'IVS10',
    func   => 'splice-5',
    polar  => '.',
    bc     => 'c.1064+2_1064+3insGGTGGTGGT',
    flanks => {
        l    => "chr1:877869-877870",
        r    => "chr1:877871-877872",
        strd => '+'
    }
};

my $rep_insertion_dup_in_out_anno = {
    c      => 'c.446AG[1>3]',
    p      => '.',
    cc     => '.',
    r      => '.',
    exin   => '.',
    func   => '.',
    polar  => '.',
    bc     => 'c.447_++1insAGAG',
    flanks => {
        l    => "chr19:32083944-32083945",
        r    => "chr19:32083946-32083947",
        strd => '+'
    }
};

my $rep_insertion_cnv_coding_anno = {
    c      => 'c.58A[1>4]',
    p      => 'p.Arg19_Ile20insLys',
    cc     => '.',
    r      => 'C1',
    exin   => 'EX2',
    func   => 'cds-indel',
    polar  => '.',
    bc     => 'c.58_59insAAA',
    flanks => {
        l    => "chr1:861378-861379",
        r    => "chr1:861380-861381",
        strd => '+'
    }
};

my $rep_insertion_cnv_intron_anno = {
    c      => 'c.305+140T[1>4]',
    p      => '.',
    cc     => '.',
    r      => 'IC3',
    exin   => 'IVS4',
    func   => 'intron',
    polar  => '.',
    bc     => 'c.305+140_305+141insTTT',
    flanks => {
        l    => "chr1:866608-866609",
        r    => "chr1:866610-866611",
        strd => '+'
    }
};

my $rep_insertion_cnv_utr_anno = {
    c      => 'c.*11G[1>6]',
    p      => '.',
    cc     => '.',
    r      => '3U1E',
    exin   => 'EX14E',
    func   => 'utr-3',
    polar  => '.',
    bc     => 'c.*11_*12insGGGGG',
    flanks => {
        l    => "chr1:879543-879544",
        r    => "chr1:879545-879546",
        strd => '+'
    }
};

test_anno( $snv_anno,  "NM_006996.2",    "chr1", 169454923, "C",  "A" );
test_anno( $del1_anno, "NM_001004696.1", "chr1", 248525638, "CA", "C" );
test_anno( $ins_anno, "NM_018406.6", "chr3", 195518112, "T",
    "TGTCTCCTGCGTAACA" );
test_anno( $intr_anno, "NM_004532.5", "chr3", 195518112, "T",
    "TGTCTCCTGCGTAACA" );
test_anno( $splice_anno, "NM_006158.4", "chr8", 24811064, "AGGGGGGG",
    "AGGGGGG" );
test_anno( $splice1_anno, "NM_024022.2", "chr21", 43803308, "CTA", "CTATA" );
test_anno( $delins_anno, "NM_001039141.2", "chr22", 38119754, "TCAAC", "TC" );
test_anno( $rep_anno, "NM_000581.2", "chr3", 49395673, "GGCCGCCGCCGCCGCCGCC",
    "GGCCGCCGCCGCC" );
test_anno( $n_anno, "NM_017780.3", "chr8", 61713130, "ACT", "ACTTGGNCT" );
test_anno( $snv_5utr_anno,         "NM_152486.2", "chr1", 861320, "TT", "TG" );
test_anno( $snv_intron_3_anno,     "NM_152486.2", "chr1", 865533, "AG", "AT" );
test_anno( $snv_intron_3_new_anno, "NM_152486.2", "chr1", 865533, "AG", "TG" );
##test_anno($snv_intron_5_anno, "NM_152486.2", "chr1", 861394, "GT", "GA");
test_anno( $snv_cdna_scoden_anno, "NM_152486.2", "chr1", 861321, "TA", "TC" );
test_anno( $snv_cdna_ecoden_anno, "NM_152486.2", "chr1", 879532, "GA", "GT" );
test_anno( $del_cdna_scoden_anno, "NM_152486.2", "chr1", 861322, "ATGTCC",
    "A" );
test_anno( $del_intron_anno, "NM_000247.1-2", "chr6_cox_hap2", 2892857,
    "TGCTG", "T" );
test_anno( $del_intron1_anno, "NM_000247.1-2", "chr6_cox_hap2", 2892854,
    "TGCTGCTG", "T" );
test_anno( $del_intron2_anno, "NM_000247.1-2", "chr6_cox_hap2", 2892856,
    "CTGCTGG", "C" );
test_anno( $del_extron_anno, "NM_005101.3", "chr1", 948953, "CATG", "C" );
test_anno( $del_extron2_anno, "NM_000070.2", "chr15", 42694319,
    "AAGGTTCCCAAAGAGGT", "A" );
test_anno( $del_5utr_coding_anno, "NM_005101.3", "chr1", 948951, "GCCATG",
    "G" );
test_anno( $del_3utr_coding_anno, "NM_005101.3", "chr1", 949857, "AAGG", "A" );
test_anno( $del_in_out_coding_anno, "NM_001205273.1", "chr19", 32083186,
    "CATGG", "C" );
test_anno( $ins_intron_splice_anno, "NM_152486.2", "chr1", 866417, "A", "AT" );
test_anno( $ins_utr_coding_anno,    "NM_152486.2", "chr1", 861321, "T", "TAG" );
test_anno( $ins_intron_exon_anno,   "NM_152486.2", "chr1", 866418, "G", "GT" );
test_anno( $ins_intron_exon2_anno, "NM_152486.2", "chr1", 866418, "G", "GTAA" );
test_anno( $ins_intron_exon3_anno, "NM_017780.3", "chr8", 61714086, "G",
    "GTAA" );
test_anno( $delins_intron_splice_anno, "NM_152486.2", "chr1", 866417, "AGAG",
    "AG" );
test_anno( $delins_utr_coding_anno, "NM_152486.2", "chr1", 861320, "TTAT",
    "TAGCTA" );
test_anno( $delins_exon_intron_anno, "NM_152486.2", "chr1", 871150, "AGCC",
    "ACCAC" );
test_anno( $delins_in_out_anno, "NM_001205273.1", "chr19", 32083185, "ACATG",
    "AGAGA" );
test_anno( $rep_insertion_dup_anno, "NM_000581.2", "chr3", 49395673, "GG",
    "GGG" );
test_anno( $rep_insertion_dup1_anno, "NM_152486.2", "chr1", 861320, "TT",
    "TTTTT" );
test_anno( $rep_insertion_dup_utr_coding_anno,
    "NM_152486.2", "chr1", 861320, "TTA", "TTATATATA" );
test_anno( $rep_insertion_dup_extron_intron_anno,
    "NM_152486.2", "chr1", 877867, "CGGT", "CGGTGGTGGTGGT" );
test_anno( $rep_insertion_dup_in_out_anno,
    "NM_001205273.1", "chr19", 32083943, "TAG", "TAGAGAG" );
test_anno( $rep_insertion_cnv_coding_anno,
    "NM_152486.2", "chr1", 861378, "AA", "AAAAA" );
test_anno( $rep_insertion_cnv_intron_anno,
    "NM_152486.2", "chr1", 866608, "TT", "TTTTT" );
test_anno( $rep_insertion_cnv_utr_anno, "NM_152486.2", "chr1", 879543, "GG",
    "GGGGGGG" );

done_testing();
exit 0;

sub test_anno {
    my ( $expect, $tid, @fourargs ) = @_;
    my $ranno = $beda->anno(@fourargs);
    if (
        !is_deeply(
            $$ranno{info}{$tid}, $expect,
            "for [ $tid," . join( ",", @fourargs ) . " ]"
        )
      )
    {
        explain "The anno infomations are: ", $ranno;
    }
}

# test if the tabix is available
sub test_tabix {
    my ($file) = @_;
    my $cmd;
    $cmd = "tabix -p bed $file";
    system($cmd);
    if ( !is( $?, 0, "tabix .. $cmd" ) ) {
        exit 1;
    }
}

