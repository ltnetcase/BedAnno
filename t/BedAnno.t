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

unlink("$data/test_db.bed.gz.tbi");
test_tabix("$data/test_db.bed.gz");
my $beda = BedAnno->new(
    db     => "$data/test_db.bed.gz",
    codon  => "$data/test_coding_seq.fas.gz",
    regbed => "$data/test.in.bed",
    trans  => "$data/trans.list"
);

#explain "The database are:", $$beda{annodb};

my $ref_anno = {
    c      => "g.=",
    p      => ".",
    cc     => ".",
    r      => "C1",
    exin   => "EX1",
    func   => ".",
    polar  => ".",
    bc     => '.',
    strd   => '-',
    flanks => {
        l    => "chr1:169454921-169454922",
        r    => "chr1:169454924-169454925"
    }
};

my $snv_anno = {
    c      => "c.82G>T",
    p      => "p.Glu28*",
    cc     => "GAA=>TAA",
    r      => "C1",
    exin   => "EX1",
    func   => "nonsense",
    polar  => "P-=>.",
    bc     => '.',
    strd   => '-',
    flanks => {
        l    => "chr1:169454921-169454922",
        r    => "chr1:169454924-169454925"
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
    strd   => '+',
    flanks => {
        l    => "chr1:248525637-248525638",
        r    => "chr1:248525640-248525641"
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
    strd   => '-',
    flanks => {
        l    => "chr8:24811063-24811064",
        r    => "chr8:24811066-24811067"
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
    strd   => '-',
    flanks => {
        l    => "chr21:43803307-43803308",
        r    => "chr21:43803309-43803310"
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
    strd   => '-',
    flanks => {
        l    => "chr3:195518111-195518112",
        r    => "chr3:195518113-195518114",
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
    strd   => '-',
    flanks => {
        l    => "chr3:195518111-195518112",
        r    => "chr3:195518113-195518114"
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
    strd   => '+',
    flanks => {
        l    => "chr22:38119753-38119754",
        r    => "chr22:38119759-38119760"
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
    strd   => '-',
    flanks => {
        l    => "chr3:49395672-49395673",
        r    => "chr3:49395680-49395681"
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
    strd   => '+',
    flanks => {
        l    => 'chr8:61713129-61713130',
        r    => 'chr8:61713133-61713134'
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
    strd   => '+',
    flanks => {
        l    => "chr1:861319-861320",
        r    => "chr1:861322-861323"
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
    strd   => '+',
    flanks => {
        l    => "chr1:865532-865533",
        r    => "chr1:865535-865536"
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
    strd   => '+',
    flanks => {
        l    => "chr1:865531-865532",
        r    => "chr1:865534-865535"
    }
};

my $snv_cdna_scoden_anno = {
    c      => 'c.1A>C',
    p      => 'p.0?',
    cc     => 'ATG=>CTG',
    r      => 'C1',
    exin   => 'EX2',
    func   => 'misstart',
    polar  => '.',
    bc     => '.',
    strd   => '+',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861323-861324"
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
    strd   => '+',
    flanks => {
        l    => "chr1:879531-879532",
        r    => "chr1:879534-879535"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861321-861322",
        r    => "chr1:861328-861329"
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
    strd   => '+',
    flanks => {
        l    => "chr6_cox_hap2:2892856-2892857",
        r    => "chr6_cox_hap2:2892862-2892863"
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
    strd   => '+',
    flanks => {
        l    => "chr6_cox_hap2:2892853-2892854",
        r    => "chr6_cox_hap2:2892862-2892863"
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
    strd   => '+',
    flanks => {
        l    => "chr6_cox_hap2:2892855-2892856",
        r    => "chr6_cox_hap2:2892863-2892864"
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
    strd   => '+',
    flanks => {
        l    => "chr1:948952-948953",
        r    => "chr1:948957-948958"
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
    strd   => '+',
    flanks => {
        l    => "chr15:42694318-42694319",
        r    => "chr15:42694336-42694337"
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
    strd   => '+',
    flanks => {
        l    => "chr1:948950-948951",
        r    => "chr1:948957-948958"
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
    strd   => '+',
    flanks => {
        l    => "chr1:949856-949857",
        r    => "chr1:949861-949862"
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
    strd   => '+',
    flanks => {
        l    => "chr19:32083185-32083186",
        r    => "chr19:32083191-32083192"
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
    strd   => '+',
    flanks => {
        l    => "chr1:866416-866417",
        r    => "chr1:866418-866419"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861322-861323"
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
    strd   => '+',
    flanks => {
        l    => "chr1:866417-866418",
        r    => "chr1:866419-866420"
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
    strd   => '+',
    flanks => {
        l    => "chr1:866417-866418",
        r    => "chr1:866419-866420"
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
    strd   => '+',
    flanks => {
        l    => "chr8:61714085-61714086",
        r    => "chr8:61714087-61714088"
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
    strd   => '+',
    flanks => {
        l    => "chr1:866416-866417",
        r    => "chr1:866421-866422"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861319-861320",
        r    => "chr1:861324-861325"
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
    strd   => '+',
    flanks => {
        l    => "chr1:871149-871150",
        r    => "chr1:871153-871154"
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
    strd   => '+',
    flanks => {
        l    => "chr19:32083184-32083185",
        r    => "chr19:32083190-32083191"
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
    strd   => '-',
    flanks => {
        l    => "chr3:49395672-49395673",
        r    => "chr3:49395674-49395675"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861320-861321",
        r    => "chr1:861322-861323"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861321-861322",
        r    => "chr1:861323-861324"
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
    strd   => '+',
    flanks => {
        l    => "chr1:877869-877870",
        r    => "chr1:877871-877872"
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
    strd   => '+',
    flanks => {
        l    => "chr19:32083944-32083945",
        r    => "chr19:32083946-32083947"
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
    strd   => '+',
    flanks => {
        l    => "chr1:861378-861379",
        r    => "chr1:861380-861381"
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
    strd   => '+',
    flanks => {
        l    => "chr1:866608-866609",
        r    => "chr1:866610-866611"
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
    strd   => '+',
    flanks => {
        l    => "chr1:879543-879544",
        r    => "chr1:879545-879546"
    }
};

my $snv_parse={
   'alt' => 'A',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'snv',
   'pos' => 14410,
   'ref' => 'C',
   'reflen' => 1
};
my $insert_parse={
   '+' => {
     'ba' => 'CGAATAGCTA',
     'bal' => 10,
     'bp' => 14410,
     'br' => 'CTAG',
     'brl' => 4
   },
   '-' => {
     'ba' => 'CGAATAGCTA',
     'bal' => 10,
     'bp' => 14410,
     'br' => 'CTAG',
     'brl' => 4
   },
   'alt' => 'CGAATAGCTA',
   'altlen' => 10,
   'chr' => 'chr1',
   'guess' => 'delins',
   'pos' => 14410,
   'ref' => 'CTAG',
   'reflen' => 4
 }
;
my $del_parse={
   'alt' => 'C',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'del',
   'pos' => 14410,
   'ref' => 'CTAGATCG',
   'reflen' => 8
};
my $rep_parse={
   '+' => {
     'ba' => 'GTAGTAGTAG',
     'bal' => 10,
     'bp' => 14413,
     'br' => 'G',
     'brl' => 1
   },
   '-' => {
     'ba' => 'CTAGTAGTAG',
     'bal' => 10,
     'bp' => 14410,
     'br' => 'C',
     'brl' => 1
   },
   'alt' => 'CTAGTAGTAGTAG',   'alt_cn' => 4,
   'altlen' => 13,
   'chr' => 'chr1',
   'guess' => 'rep',
   'pos' => 14410,
   'ref' => 'CTAG',
   'ref_cn' => 1,
   'reflen' => 4,
   'rep' => 'TAG',
   'replen' => 3
};
my $dup_parse={
  '+' => {
     'ba' => 'ATA',
     'bal' => 3,
     'bp' => 14412,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'CTA',
     'bal' => 3,
     'bp' => 14410,
     'br' => 'C',
     'brl' => 1
   },
   'alt' => 'CTATA',
   'alt_cn' => 2,
   'altlen' => 5,
   'chr' => 'chr1',
   'guess' => 'rep',
   'pos' => 14410,
   'ref' => 'CTA',
   'ref_cn' => 1,
   'reflen' => 3,
   'rep' => 'TA',
   'replen' => 2
};

my $delins_parse={
   '+' => {
     'ba' => 'GT',
     'bal' => 2,
     'bp' => 14410,
     'br' => 'CTAGA',
     'brl' => 5
   },
   '-' => {
     'ba' => 'GT',
     'bal' => 2,
     'bp' => 14410,
     'br' => 'CTAGA',
     'brl' => 5
   },
   'alt' => 'GT',
   'altlen' => 2,
   'chr' => 'chr1',
   'guess' => 'delins',
   'pos' => 14410,
   'ref' => 'CTAGA',
   'reflen' => 5
 };
my $n_parse={
   'alt' => 'CNTAN',
   'altlen' => 5,
   'chr' => 'chr1',
   'guess' => 'ins',
   'pos' => 14410,
   'ref' => 'C',
   'reflen' => 1
 };
my $splice_parse_1={
   '+' => {
     'ba' => 'G',
     'bal' => 1,
     'bp' => 24811070,
     'br' => 'GG',
     'brl' => 2
   },
   '-' => {
     'ba' => 'A',
     'bal' => 1,
     'bp' => 24811064,
     'br' => 'AG',
     'brl' => 2
   },
   'alt' => 'AGGGGGG',
   'alt_cn' => 6,
   'altlen' => 7,
   'chr' => 'chr8',
   'guess' => 'rep',
   'pos' => 24811064,
   'ref' => 'AGGGGGGG',
   'ref_cn' => 7,
   'reflen' => 8,
   'rep' => 'G',
   'replen' => 1
};
my $splice_parse_2={
   '+' => {
     'ba' => 'GA',
     'bal' => 2,
     'bp' => 24811068,
     'br' => 'G',
     'brl' => 1
   },
   '-' => {
     'ba' => 'GA',
     'bal' => 2,
     'bp' => 24811068,
     'br' => 'G',
     'brl' => 1
   },
   'alt' => 'GA',
   'altlen' => 2,
   'chr' => 'chr8',
   'guess' => 'ins',
   'pos' => 24811068,
   'ref' => 'G',
   'reflen' => 1
};
my $splice1_parse_1={
'+' => {
     'ba' => 'ATA',
     'bal' => 3,
     'bp' => 43803310,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'CTA',
     'bal' => 3,
     'bp' => 43803308,
     'br' => 'C',
     'brl' => 1
   },
   'alt' => 'CTATA',
   'alt_cn' => 2,
   'altlen' => 5,
   'chr' => 'chr21',
   'guess' => 'rep',
   'pos' => 43803308,
   'ref' => 'CTA',
   'ref_cn' => 1,
   'reflen' => 3,
   'rep' => 'TA',
   'replen' => 2
};

my $splice1_parse_2={
   'alt' => 'T',
   'altlen' => 1,
   'chr' => 'chr21',
   'guess' => 'del',
   'pos' => 43803309,
   'ref' => 'TAG',
   'reflen' => 3
};

my $cds_parse_1={
   '+' => {
     'ba' => 'CC',
     'bal' => 2,
     'bp' => 49395684,
     'br' => 'C',
     'brl' => 1
   },
   '-' => {
     'ba' => 'GC',
     'bal' => 2,
     'bp' => 49395683,
     'br' => 'G',
     'brl' => 1
   },
   'alt' => 'GGCCGCCGCCGCC',
   'altlen' => 13,
   'chr' => 'chr3',
   'guess' => 'delins',
   'pos' => 49395673,
   'ref' => 'GGCCGCCGCCGC',
   'reflen' => 12
};
my $cds_parse_2={
   '+' => {
     'ba' => 'C',
     'bal' => 1,
     'bp' => 49395679,
     'br' => 'CGCCGC',
     'brl' => 6
   },
   '-' => {
     'ba' => 'C',
     'bal' => 1,
     'bp' => 49395678,
     'br' => 'CCGCCG',
     'brl' => 6
   },
   'alt' => 'GCC',
   'altlen' => 3,
   'chr' => 'chr3',
   'guess' => 'delins',
   'pos' => 49395677,
   'ref' => 'GCCGCCGC',
   'reflen' => 8
};

my $utr5_parse_1={
   '+' => {
     'ba' => 'G',
     'bal' => 1,
     'bp' => 861321,
     'br' => 'T',
     'brl' => 1
   },
   '-' => {
     'ba' => 'G',
     'bal' => 1,
    'bp' => 861321,
     'br' => 'T',
     'brl' => 1
   },
   'alt' => 'G',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'snv',
   'pos' => 861321,
   'ref' => 'T',
   'reflen' => 1
};
my $utr5_parse_2={ 
   'alt' => 'T',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'del',
   'pos' => 861320,
   'ref' => 'TTTA',
   'reflen' => 4
};

my $cdna_scoden_parse_1={
   '+' => {
     'ba' => 'C',
     'bal' => 1,
     'bp' => 861322,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'C',
     'bal' => 1,
     'bp' => 861322,
     'br' => 'A',
     'brl' => 1
   },
   'alt' => 'C',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'snv',
   'pos' => 861322,
   'ref' => 'A',
   'reflen' => 1
};
my $cdna_scoden_parse_2={
   'alt' => 'GC',
   'altlen' => 2,
   'chr' => 'chr1',
   'guess' => 'delins',
   'pos' => 861322,
   'ref' => 'A',
   'reflen' => 1
};
my $cdna_ecoden_parse_1={
   '+' => {
     'ba' => 'T',
     'bal' => 1,
     'bp' => 879533,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'T',
     'bal' => 1,
     'bp' => 879533,
     'br' => 'A',
     'brl' => 1
   },
   'alt' => 'T',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'snv',
   'pos' => 879533,
   'ref' => 'A',
   'reflen' => 1
};
my $cdna_ecoden_parse_2={
   '+' => {
     'ba' => 'AAAA',
     'bal' => 4,
     'bp' => 879533,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'GAAA',
     'bal' => 4,
     'bp' => 879532,
     'br' => 'G',
     'brl' => 1
   },
   'alt' => 'GAAAA',
   'alt_cn' => 4,
   'altlen' => 5,
  'chr' => 'chr1',
   'guess' => 'rep',
   'pos' => 879532,
   'ref' => 'GA',
   'ref_cn' => 1,
   'reflen' => 2,
   'rep' => 'A',
   'replen' => 1
};

my $intron_parse_1={
   'alt' => 'T',
   'altlen' => 1,
   'chr' => 'chr6_cox_hap2',
   'guess' => 'del',
   'pos' => 2892857,
   'ref' => 'TGCTG',
   'reflen' => 5
};
my $intron_parse_2={
   '+' => {
     'ba' => 'TA',
     'bal' => 2,
     'bp' => 2892857,
     'br' => 'T',
     'brl' => 1
   },
   '-' => {
     'ba' => 'TA',
     'bal' => 2,
     'bp' => 2892857,
     'br' => 'T',
     'brl' => 1
   },
   'alt' => 'TA',
   'altlen' => 2,
   'chr' => 'chr6_cox_hap2',
   'guess' => 'ins',
   'pos' => 2892857,
   'ref' => 'T',
   'reflen' => 1
};

my $utr5_del_parse_1={
   'alt' => 'G',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'del',
   'pos' => 948951,
   'ref' => 'GCCATG',
   'reflen' => 6
};
my $utr5_del_parse_2={
   'alt' => 'C',
   'altlen' => 1,
   'chr' => 'chr1',
   'guess' => 'del',
   'pos' => 948953,
   'ref' => 'CATG',
   'reflen' => 4
};

my $ins_intron_exon3_parse_1={
 'alt' => 'GTAA',
   'altlen' => 4,
   'chr' => 'chr8',
   'guess' => 'ins',
   'pos' => 61714086,
   'ref' => 'G',
   'reflen' => 1
};
my $ins_intron_exon3_parse_2={
   'alt' => 'TGGA',
   'altlen' => 4,
   'chr' => 'chr8',
   'guess' => 'ins',
   'pos' => 61714087,
   'ref' => 'T',
   'reflen' => 1
};

my $rep_insertion_dup_parse_1={
   '+' => {
     'ba' => 'GAGAG',
     'bal' => 5,
     'bp' => 32083945,
    'br' => 'G',
     'brl' => 1
   },
   '-' => {
     'ba' => 'TAGAG',
     'bal' => 5,
     'bp' => 32083943,
     'br' => 'T',
     'brl' => 1
   },
   'alt' => 'TAGAGAG',
   'alt_cn' => 3,
   'altlen' => 7,
   'chr' => 'chr19',
   'guess' => 'rep',
   'pos' => 32083943,
   'ref' => 'TAG',
   'ref_cn' => 1,
   'reflen' => 3,
   'rep' => 'AG',
   'replen' => 2
};
my $rep_insertion_dup_parse_2={
   '+' => {
     'ba' => 'AAAAA',
     'bal' => 5,
     'bp' => 32083944,
     'br' => 'A',
     'brl' => 1
   },
   '-' => {
     'ba' => 'TAAAA',
     'bal' => 5,
     'bp' => 32083943,
     'br' => 'T',
     'brl' => 1
   },
   'alt' => 'TAAAAA',
   'alt_cn' => 5,
   'altlen' => 6,
   'chr' => 'chr19',
   'guess' => 'rep',
   'pos' => 32083943,
   'ref' => 'TA',
   'ref_cn' => 1,
   'reflen' => 2,
   'rep' => 'A',
   'replen' => 1
};

my $snv_gHGVS="g.14410C>A";
my $insert_gHGVS="g.14411_14413delTAGinsGAATAGCTA";
my $del_gHGVS="g.14411_14417delTAGATCG";
my $rep_gHGVS="g.14411TAG[1>4]";
my $dup_HGVS="g.14411_14412dupTA";
my $delins_HGVS="g.14410_14414delCTAGAinsGT";
my $n_HGVS="g.14410_14411insNTAN";
my $splice_1_gHGVS="g.24811065G[7>6]";
my $splice_2_gHGVS="g.24811068_24811069insA";
my $splice1_1_gHGVS="g.43803309_43803310dupTA";
my $splice1_2_gHGVS="g.43803310_43803311delAG";
my $cds_1_gHGVS="g.49395674_49395684delGCCGCCGCCGCinsGCCGCCGCCGCC";
my $cds_2_gHGVS="g.49395678_49395684delCCGCCGCinsCC";
my $utr5_1_gHGVS="g.861321T>G";
my $utr5_2_gHGVS="g.861321_861323delTTA";
my $cdna_scoden_1_gHGVS="g.861322A>C";
my $cdna_scoden_2_gHGVS="g.861322delAinsGC";
my $cdna_ecoden_1_gHGVS="g.879533A>T";
my $cdna_ecoden_2_gHGVS="g.879533A[1>4]";
my $intron_1_gHGVS="g.2892858_2892861delGCTG";
my $intron_2_gHGVS="g.2892857_2892858insA";
my $utr5_del_1_gHGVS="g.948952_948956delCCATG";
my $utr5_del_2_gHGVS="g.948954_948956delATG";
my $ins_intron_exon3_1_gHGVS="g.61714086_61714087insTAA";
my $ins_intron_exon3_2_gHGVS="g.61714087_61714088insGGA"; 
my $rep_insertion_dup_1_gHGVS="g.32083944AG[1>3]";
my $rep_insertion_dup_2_gHGVS="g.32083944A[1>5]";

my $snv_varanno={
   'info' => {
     'NR_024540.1' => {
       'bc' => '.',
       'c' => 'n.1721G>T',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14408-14409',
         'r' => 'chr1:14411-14412'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
  },
   'var' => {
     'alt' => 'A',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'snv',
     'pos' => 14410,
     'ref' => 'C',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1721',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $insert_varanno={
   'info' => {
     'NR_024540.1' => {
       'bc' => '.',
       'c' => 'n.1718_1720delCTAinsTAGCTATTC',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14409-14410',
         'r' => 'chr1:14414-14415'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E',
     }
   },
   'var' => {
     '+' => {
       'ba' => 'CGAATAGCTA',
       'bal' => 10,
       'bp' => 14410,
       'br' => 'CTAG',
       'brl' => 4
     },
     '-' => {
       'ba' => 'CGAATAGCTA',
       'bal' => 10,
       'bp' => 14410,
       'br' => 'CTAG',
       'brl' => 4
     },
     'alt' => 'CGAATAGCTA',
     'altlen' => 10,
     'chr' => 'chr1',
     'guess' => 'delins',
     'pos' => 14410,
     'ref' => 'CTAG',
     'reflen' => 4,
     'sel' => {
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1718',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $del_varanno={
   'info' => {
     'NR_024540.1' => {
       'bc' => '.',
       'c' => 'n.1714_1720delCGATCTA',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14409-14410',
         'r' => 'chr1:14418-14419'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
   },
   'var' => {
     'alt' => 'C',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'del',
     'pos' => 14410,
     'ref' => 'CTAGATCG',
     'reflen' => 8,
     'sel' => {
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1714',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
   }
 };
my $rep_varanno={
   'info' => {
     'NR_024540.1' => {
       'bc' => 'n.1720_1721insCTACTACTA',
       'c' => 'n.1718CTA[1>4]',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14409-14410',
         'r' => 'chr1:14411-14412'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'GTAGTAGTAG',
       'bal' => 10,
       'bp' => 14413,
       'br' => 'G',
       'brl' => 1
     },
     '-' => {
       'ba' => 'CTAGTAGTAG',
       'bal' => 10,
       'bp' => 14410,
       'br' => 'C',
       'brl' => 1
     },
     'alt' => 'CTAGTAGTAGTAG',
     'alt_cn' => 4,
     'altlen' => 13,
     'chr' => 'chr1',
     'guess' => 'rep',
     'pos' => 14410,
     'ref' => 'CTAG',
     'ref_cn' => 1,
     'reflen' => 4,
     'rep' => 'TAG',
     'replen' => 3,
     'sel' => {
       'bc' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1721',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ],
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1718',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $dup_varanno={
   'info' => {
     'NR_024540.1' => {
       'bc' => 'n.1720_1721insTA',
       'c' => 'n.1719_1720dupTA',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14409-14410',
         'r' => 'chr1:14411-14412'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'ATA',
       'bal' => 3,
       'bp' => 14412,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'CTA',
       'bal' => 3,
       'bp' => 14410,
       'br' => 'C',
       'brl' => 1
     },
     'alt' => 'CTATA',
     'alt_cn' => 2,
     'altlen' => 5,
     'chr' => 'chr1',
     'guess' => 'rep',
     'pos' => 14410,
     'ref' => 'CTA',
     'ref_cn' => 1,
     'reflen' => 3,
     'rep' => 'TA',
     'replen' => 2,
     'sel' => {
       'bc' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1721',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ],
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1719',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $delins_varanno={
'info' => {
     'NR_024540.1' => {
       'bc' => '.',
       'c' => 'n.1717_1721delTCTAGinsAC',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14408-14409',
         'r' => 'chr1:14415-14416'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'GT',
       'bal' => 2,
       'bp' => 14410,
       'br' => 'CTAGA',
       'brl' => 5
     },
     '-' => {
       'ba' => 'GT',
       'bal' => 2,
       'bp' => 14410,
       'br' => 'CTAGA',
       'brl' => 5
     },
     'alt' => 'GT',
     'altlen' => 2,
     'chr' => 'chr1',
     'guess' => 'delins',
     'pos' => 14410,
     'ref' => 'CTAGA',
     'reflen' => 5,
     'sel' => {
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1721',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1717',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
          }
         ]
       ]
     }
   }
};
my $n_varanno={
'info' => {
     'NR_024540.1' => {
       'bc' => '.',
       'c' => 'n.1720_1721insNTAN',
       'cc' => '.',
       'exin' => 'EX11E',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr1:14409-14410',
         'r' => 'chr1:14411-14412'
       },
       'func' => 'ncRNA',
       'p' => '.',
       'polar' => '.',
       'r' => 'R11E'
     }
   },
   'var' => {
     'alt' => 'CNTAN',
     'altlen' => 5,
     'chr' => 'chr1',
     'guess' => 'ins',
     'pos' => 14410,
     'ref' => 'C',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NR_024540.1',
           {
             'bd' => 0,
             'cpos' => 'n.1721',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'n.1720',
             'exin' => 'EX11E',
             'gsym' => 'WASH7P',
             'reg' => 'R11E',
             'strd' => '-'
           }
         ]
       ]
     }
 }
};
my $splice_1_varanno={
'info' => {
     'NM_006158.4' => {
       'bc' => 'c.1413+1delC',
       'c' => 'c.1408C[7>6]',
       'cc' => '.',
       'exin' => 'IVS3',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr8:24811063-24811064',
         'r' => 'chr8:24811066-24811067'
       },
       'func' => 'abnormal-intron',
       'p' => '.',
       'polar' => '.',
       'r' => 'IC3'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'G',
       'bal' => 1,
       'bp' => 24811070,
       'br' => 'GG',
       'brl' => 2
     },
     '-' => {
       'ba' => 'A',
       'bal' => 1,
       'bp' => 24811064,
       'br' => 'AG',
       'brl' => 2
     },
     'alt' => 'AGGGGGG',
     'alt_cn' => 6,
     'altlen' => 7,
     'chr' => 'chr8',
     'guess' => 'rep',
     'pos' => 24811064,
     'ref' => 'AGGGGGGG',
     'ref_cn' => 7,
     'reflen' => 8,
     'rep' => 'G',
     'replen' => 1,
     'sel' => {
       'bc' => [
         [
           'NM_006158.4',
           {
             'bd' => 1,
             'cpos' => 'c.1413+1',
             'exin' => 'IVS3',
             'gsym' => 'NEFL',
             'reg' => 'IC3',
             'strd' => '-'
           }
         ]
       ],
       'std' => [
         [
           'NM_006158.4',
           {
             'bd' => 0,
             'cpos' => 'c.1408',
             'exin' => 'EX3',
             'gsym' => 'NEFL',
             'reg' => 'C3',
             'strd' => '-'
           }
         ]
       ]
     }
   }
 };
my $splice_2_varanno={
   'info' => {
     'NM_006158.4' => {
       'bc' => '.',
       'c' => 'c.1410_1411insT',
       'cc' => '.',
       'exin' => 'EX3',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr8:24811067-24811068',
         'r' => 'chr8:24811069-24811070'
       },
       'func' => 'frameshift',
       'p' => 'p.Pro471fs*?',
       'polar' => '.',
       'r' => 'C3'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'GA',
       'bal' => 2,
       'bp' => 24811068,
       'br' => 'G',
       'brl' => 1
     },
     '-' => {
       'ba' => 'GA',
       'bal' => 2,
       'bp' => 24811068,
       'br' => 'G',
       'brl' => 1
     },
     'alt' => 'GA',
     'altlen' => 2,
     'chr' => 'chr8',
     'guess' => 'ins',
     'pos' => 24811068,
     'ref' => 'G',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_006158.4',
           {
             'bd' => 0,
             'cpos' => 'c.1411',
             'exin' => 'EX3',
             'gsym' => 'NEFL',
             'reg' => 'C3',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'c.1410',
             'exin' => 'EX3',
             'gsym' => 'NEFL',
             'reg' => 'C3',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $splice1_1_varanno={
   'info' => {
     'NM_024022.2' => {
       'bc' => 'c.617-2_617-1insTA',
       'c' => 'c.617-3_617-2dupTA',
       'cc' => '.',
       'exin' => 'IVS7',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr21:43803307-43803308',
         'r' => 'chr21:43803309-43803310'
       },
       'func' => 'splice-3',
       'p' => '.',
       'polar' => '.',
       'r' => 'IC6'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'ATA',
       'bal' => 3,
       'bp' => 43803310,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'CTA',
       'bal' => 3,
       'bp' => 43803308,
       'br' => 'C',
       'brl' => 1
     },
     'alt' => 'CTATA',
     'alt_cn' => 2,
     'altlen' => 5,
     'chr' => 'chr21',
     'guess' => 'rep',
     'pos' => 43803308,
     'ref' => 'CTA',
     'ref_cn' => 1,
     'reflen' => 3,
     'rep' => 'TA',
     'replen' => 2,
     'sel' => {
       'bc' => [
         [
           'NM_024022.2',
           {
             'bd' => 'B',
             'cpos' => 'c.617-1',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           },
           {
             'bd' => 'B',
             'cpos' => 'c.617-2',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           }
         ]
       ],
       'std' => [
         [
           'NM_024022.2',
           {
             'bd' => 'B',
             'cpos' => 'c.617-2',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'c.617-3',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           }
         ]
       ]
     }
   }
 };
my $splice1_2_varanno={
 'info' => {
     'NM_024022.2' => {
       'bc' => '.',
       'c' => 'c.617-4_617-3delCT',
       'cc' => '.',
       'exin' => 'IVS7',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr21:43803308-43803309',
         'r' => 'chr21:43803312-43803313'
       },
       'func' => 'intron',
       'p' => '.',
       'polar' => '.',
       'r' => 'IC6'
     }
   },
   'var' => {
     'alt' => 'T',
     'altlen' => 1,
     'chr' => 'chr21',
     'guess' => 'del',
     'pos' => 43803309,
     'ref' => 'TAG',
     'reflen' => 3,
     'sel' => {
       'std' => [
         [
           'NM_024022.2',
           {
             'bd' => 0,
             'cpos' => 'c.617-3',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'c.617-4',
             'exin' => 'IVS7',
             'gsym' => 'TMPRSS3',
             'reg' => 'IC6',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $cds_1_varanno={
'info' => {
     'NM_000581.2' => {
       'bc' => '.',
       'c' => 'c.28_38delGCGGCGGCGGCinsGGCGGCGGCGGC',
       'cc' => '.',
       'exin' => 'EX1',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr3:49395672-49395673',
         'r' => 'chr3:49395685-49395686'
       },
       'func' => 'frameshift',
       'p' => 'p.Ala10_Ala13delinsGlyGlyGlyGlyProfs*?',
       'polar' => '.',
       'r' => 'C1'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'CC',
       'bal' => 2,
       'bp' => 49395684,
       'br' => 'C',
       'brl' => 1
     },
     '-' => {
       'ba' => 'GC',
       'bal' => 2,
       'bp' => 49395683,
       'br' => 'G',
       'brl' => 1
     },
     'alt' => 'GGCCGCCGCCGCC',
     'altlen' => 13,
     'chr' => 'chr3',
     'guess' => 'delins',
     'pos' => 49395673,
     'ref' => 'GGCCGCCGCCGC',
     'reflen' => 12,
     'sel' => {
       'std' => [
         [
           'NM_000581.2',
           {
             'bd' => 0,
             'cpos' => 'c.38',
             'exin' => 'EX1',
             'gsym' => 'GPX1',
             'reg' => 'C1',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'c.28',
             'exin' => 'EX1',
             'gsym' => 'GPX1',
             'reg' => 'C1',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $cds_2_varanno={
'info' => {
     'NM_000581.2' => {
       'bc' => '.',
       'c' => 'c.28_34delGCGGCGGinsGG',
       'cc' => '.',
       'exin' => 'EX1',
       'strd' => '-',
       'flanks' => {
         'l' => 'chr3:49395676-49395677',
         'r' => 'chr3:49395685-49395686'
       },
       'func' => 'frameshift',
       'p' => 'p.Ala10_Ala12delinsGlyfs*?',
       'polar' => '.',
       'r' => 'C1'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'C',
       'bal' => 1,
       'bp' => 49395679,
       'br' => 'CGCCGC',
       'brl' => 6
     },
     '-' => {
       'ba' => 'C',
       'bal' => 1,
       'bp' => 49395678,
       'br' => 'CCGCCG',
       'brl' => 6
     },
     'alt' => 'GCC',
     'altlen' => 3,
     'chr' => 'chr3',
     'guess' => 'delins',
     'pos' => 49395677,
     'ref' => 'GCCGCCGC',
     'reflen' => 8,
     'sel' => {
       'std' => [
         [
           'NM_000581.2',
           {
             'bd' => 0,
             'cpos' => 'c.34',
             'exin' => 'EX1',
             'gsym' => 'GPX1',
             'reg' => 'C1',
             'strd' => '-'
           },
           {
             'bd' => 0,
             'cpos' => 'c.28',
             'exin' => 'EX1',
             'gsym' => 'GPX1',
             'reg' => 'C1',
             'strd' => '-'
           }
         ]
       ]
     }
   }
};
my $utr5_1_varanno={
 'info' => {
     'NM_152486.2' => {
       'bc' => '.',
       'c' => 'c.-1T>G',
       'cc' => '.',
       'exin' => 'EX2',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:861319-861320',
         'r' => 'chr1:861322-861323'
       },
       'func' => 'utr-5',
       'p' => '.',
       'polar' => '.',
       'r' => '5U1'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'G',
       'bal' => 1,
       'bp' => 861321,
       'br' => 'T',
       'brl' => 1
     },
     '-' => {
       'ba' => 'G',
       'bal' => 1,
       'bp' => 861321,
       'br' => 'T',
       'brl' => 1
     },
     'alt' => 'G',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'snv',
     'pos' => 861321,
     'ref' => 'T',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'B',
             'cpos' => 'c.-1',
             'exin' => 'EX2',
             'gsym' => 'SAMD11',
             'reg' => '5U1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $utr5_2_varanno={
'info' => {
     'NM_152486.2' => {
       'bc' => '.',
       'c' => 'c.-1_2delTTA',
       'cc' => '.',
       'exin' => 'EX2',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:861319-861320',
         'r' => 'chr1:861324-861325'
       },
       'func' => 'init-loss',
       'p' => 'p.0?',
       'polar' => '.',
       'r' => '5U1-C1'
     }
   },
   'var' => {
     'alt' => 'T',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'del',
     'pos' => 861320,
     'ref' => 'TTTA',
     'reflen' => 4,
     'sel' => {
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'B',
             'cpos' => 'c.-1',
             'exin' => 'EX2',
             'gsym' => 'SAMD11',
             'reg' => '5U1',
             'strd' => '+'
           },
           {
             'bd' => 'b',
             'cpos' => 'c.2',
             'exin' => 'EX2',
             'gsym' => 'SAMD11',
             'reg' => 'C1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };
my $cdna_scoden_1_varanno={
 'info' => {
     'NM_152486.2' => {
       'bc' => '.',
       'c' => 'c.1A>C',
       'cc' => 'ATG=>CTG',
       'exin' => 'EX2',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:861320-861321',
         'r' => 'chr1:861323-861324'
       },
       'func' => 'misstart',
       'p' => 'p.0?',
       'polar' => '.',
       'r' => 'C1'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'C',
       'bal' => 1,
       'bp' => 861322,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'C',
       'bal' => 1,
       'bp' => 861322,
       'br' => 'A',
       'brl' => 1
     },
     'alt' => 'C',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'snv',
     'pos' => 861322,
     'ref' => 'A',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'b',
             'cpos' => 'c.1',
             'exin' => 'EX2',
             'gsym' => 'SAMD11',
             'reg' => 'C1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $cdna_scoden_2_varanno={
   'info' => {
     'NM_152486.2' => {
       'bc' => '.',
       'c' => 'c.1delAinsGC',
       'cc' => '.',
       'exin' => 'EX2',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:861320-861321',
         'r' => 'chr1:861323-861324'
       },
       'func' => 'init-loss',
       'p' => 'p.0?',
       'polar' => '.',
       'r' => 'C1'
     }
   },
   'var' => {
     'alt' => 'GC',
     'altlen' => 2,
     'chr' => 'chr1',
     'guess' => 'delins',
     'pos' => 861322,
     'ref' => 'A',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'b',
             'cpos' => 'c.1',
             'exin' => 'EX2',
             'gsym' => 'SAMD11',
             'reg' => 'C1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $cdna_ecoden_1_varanno={
  'info' => {
     'NM_152486.2' => {
       'bc' => '.',
       'c' => 'c.2046A>T',
       'cc' => 'TGA=>TGT',
       'exin' => 'EX14E',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:879531-879532',
         'r' => 'chr1:879534-879535'
       },
       'func' => 'stop-loss',
       'p' => 'p.*682Cysext*?',
       'polar' => '.=>P0',
       'r' => 'C13E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'T',
       'bal' => 1,
       'bp' => 879533,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'T',
       'bal' => 1,
       'bp' => 879533,
       'br' => 'A',
       'brl' => 1
     },
     'alt' => 'T',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'snv',
     'pos' => 879533,
     'ref' => 'A',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'B',
             'cpos' => 'c.2046',
             'exin' => 'EX14E',
             'gsym' => 'SAMD11',
             'reg' => 'C13E',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $cdna_ecoden_2_varanno={
  'info' => {
     'NM_152486.2' => {
       'bc' => 'c.2046_*1insAAA',
       'c' => 'c.2046A[1>4]',
       'cc' => '.',
       'exin' => 'EX14E',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:879532-879533',
         'r' => 'chr1:879534-879535'
       },
       'func' => 'utr-3',
       'p' => '.',
       'polar' => '.',
       'r' => 'C13E-3U1E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'AAAA',
       'bal' => 4,
       'bp' => 879533,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'GAAA',
       'bal' => 4,
       'bp' => 879532,
       'br' => 'G',
       'brl' => 1
     },
     'alt' => 'GAAAA',
     'alt_cn' => 4,
     'altlen' => 5,
     'chr' => 'chr1',
     'guess' => 'rep',
     'pos' => 879532,
     'ref' => 'GA',
     'ref_cn' => 1,
     'reflen' => 2,
     'rep' => 'A',
     'replen' => 1,
     'sel' => {
       'bc' => [
         [
           'NM_152486.2',
           {
             'bd' => 'B',
             'cpos' => 'c.2046',
             'exin' => 'EX14E',
             'gsym' => 'SAMD11',
             'reg' => 'C13E',
             'strd' => '+'
           },
           {
             'bd' => 'b',
             'cpos' => 'c.*1',
             'exin' => 'EX14E',
             'gsym' => 'SAMD11',
             'reg' => '3U1E',
             'strd' => '+'
           }
         ]
       ],
       'std' => [
         [
           'NM_152486.2',
           {
             'bd' => 'B',
             'cpos' => 'c.2046',
             'exin' => 'EX14E',
             'gsym' => 'SAMD11',
             'reg' => 'C13E',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };
my $intron_1_varanno={
'info' => {
     'NM_000247.1-2' => {
       'bc' => '.',
       'c' => 'c.948+1_949-1delGCTG',
       'cc' => '.',
       'exin' => 'IVS5',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr6_cox_hap2:2892856-2892857',
         'r' => 'chr6_cox_hap2:2892862-2892863'
       },
       'func' => 'splice',
       'p' => '.',
       'polar' => '.',
       'r' => 'IC5'
     }
   },
   'var' => {
     'alt' => 'T',
     'altlen' => 1,
     'chr' => 'chr6_cox_hap2',
     'guess' => 'del',
     'pos' => 2892857,
     'ref' => 'TGCTG',
     'reflen' => 5,
     'sel' => {
       'std' => [
         [
           'NM_000247.1-2',
           {
             'bd' => 'b',
             'cpos' => 'c.948+1',
             'exin' => 'IVS5',
             'gsym' => 'MICA',
             'reg' => 'IC5',
             'strd' => '+'
           },
           {
             'bd' => 'B',
             'cpos' => 'c.949-1',
             'exin' => 'IVS5',
             'gsym' => 'MICA',
             'reg' => 'IC5',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };

my $intron_2_varanno={
 'info' => {
     'NM_000247.1-2' => {
       'bc' => '.',
       'c' => 'c.948_948+1insA',
       'cc' => '.',
       'exin' => 'EX5-IVS5',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr6_cox_hap2:2892856-2892857',
         'r' => 'chr6_cox_hap2:2892858-2892859'
       },
       'func' => 'frameshift',
       'p' => 'p.Ala317fs*?',
       'polar' => '.',
       'r' => 'C5-IC5'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'TA',
       'bal' => 2,
       'bp' => 2892857,
       'br' => 'T',
       'brl' => 1
     },
     '-' => {
       'ba' => 'TA',
       'bal' => 2,
       'bp' => 2892857,
       'br' => 'T',
       'brl' => 1
     },
     'alt' => 'TA',
     'altlen' => 2,
     'chr' => 'chr6_cox_hap2',
     'guess' => 'ins',
     'pos' => 2892857,
     'ref' => 'T',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_000247.1-2',
           {
             'bd' => 'B',
             'cpos' => 'c.948',
             'exin' => 'EX5',
             'gsym' => 'MICA',
             'reg' => 'C5',
             'strd' => '+'
           },
           {
             'bd' => 'b',
             'cpos' => 'c.948+1',
             'exin' => 'IVS5',
             'gsym' => 'MICA',
             'reg' => 'IC5',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };
my $utr5_del_1_varanno ={
   'info' => {
     'NM_005101.3' => {
       'bc' => '.',
       'c' => 'c.-2_3delCCATG',
       'cc' => '.',
       'exin' => 'EX1',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:948950-948951',
         'r' => 'chr1:948957-948958'
       },
       'func' => 'init-loss',
       'p' => 'p.0?',
       'polar' => '.',
       'r' => '5U1E-C1'
     }
   },
   'var' => {
     'alt' => 'G',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'del',
     'pos' => 948951,
     'ref' => 'GCCATG',
     'reflen' => 6,
     'sel' => {
       'std' => [
         [
           'NM_005101.3',
           {
             'bd' => 'B',
             'cpos' => 'c.-2',
             'exin' => 'EX1',
             'gsym' => 'ISG15',
             'reg' => '5U1E',
             'strd' => '+'
           },
           {
             'bd' => 'B',
             'cpos' => 'c.3',
             'exin' => 'EX1',
             'gsym' => 'ISG15',
             'reg' => 'C1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $utr5_del_2_varanno ={     
   'info' => {
     'NM_005101.3' => {
       'bc' => '.',
       'c' => 'c.1_3delATG',
       'cc' => '.',
       'exin' => 'EX1',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr1:948952-948953',
         'r' => 'chr1:948957-948958'
       },
       'func' => 'init-loss',
       'p' => 'p.0?',
       'polar' => '.',
       'r' => 'C1'
     }
   },
   'var' => {
     'alt' => 'C',
     'altlen' => 1,
     'chr' => 'chr1',
     'guess' => 'del',
     'pos' => 948953,
     'ref' => 'CATG',
     'reflen' => 4,
     'sel' => {
       'std' => [
         [
           'NM_005101.3',
           {
             'bd' => 'b',
             'cpos' => 'c.1',
             'exin' => 'EX1',
             'gsym' => 'ISG15',
             'reg' => 'C1',
             'strd' => '+'
           },
           {
             'bd' => 'B',
             'cpos' => 'c.3',
             'exin' => 'EX1',
             'gsym' => 'ISG15',
             'reg' => 'C1',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $ins_intron_exon3_1_varanno={
   'info' => {
     'NM_017780.3' => {
       'bc' => '.',
       'c' => 'c.2377-1_2377insTAA',
       'cc' => '.',
       'exin' => 'IVS5-EX6',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr8:61714085-61714086',
         'r' => 'chr8:61714087-61714088'
       },
       'func' => 'stop-gain',
       'p' => 'p.Gln792_Gln792ins*',
       'polar' => '.',
       'r' => 'IC4-C5'
     }
   },
   'var' => {
     'alt' => 'GTAA',
     'altlen' => 4,
     'chr' => 'chr8',
     'guess' => 'ins',
     'pos' => 61714086,
     'ref' => 'G',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_017780.3',
           {
             'bd' => 'B',
             'cpos' => 'c.2377-1',
             'exin' => 'IVS5',
             'gsym' => 'CHD7',
             'reg' => 'IC4',
             'strd' => '+'
           },
           {
             'bd' => 'b',
             'cpos' => 'c.2377',
             'exin' => 'EX6',
             'gsym' => 'CHD7',
             'reg' => 'C5',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };

my $ins_intron_exon3_2_varanno={
  'info' => {
     'NM_017780.3' => {
       'bc' => '.',
       'c' => 'c.2377_2378insGGA',
       'cc' => '.',
       'exin' => 'EX6',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr8:61714086-61714087',
         'r' => 'chr8:61714088-61714089'
       },
       'func' => 'cds-indel',
       'p' => 'p.Glu793delinsGlyLys',
       'polar' => '.',
       'r' => 'C5'
     }
   },
   'var' => {
     'alt' => 'TGGA',
     'altlen' => 4,
     'chr' => 'chr8',
     'guess' => 'ins',
     'pos' => 61714087,
     'ref' => 'T',
     'reflen' => 1,
     'sel' => {
       'std' => [
         [
           'NM_017780.3',
           {
             'bd' => 'b',
             'cpos' => 'c.2377',
             'exin' => 'EX6',
             'gsym' => 'CHD7',
             'reg' => 'C5',
             'strd' => '+'
           },
           {
             'bd' => 'b',
             'cpos' => 'c.2378',
             'exin' => 'EX6',
             'gsym' => 'CHD7',
             'reg' => 'C5',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $rep_insertion_dup_1_varanno={
   'info' => {
     'NM_001205273.1' => {
       'bc' => 'c.447_++1insAGAG',
       'c' => 'c.446AG[1>3]',
       'cc' => '.',
       'exin' => '.',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr19:32083944-32083945',
         'r' => 'chr19:32083946-32083947'
       },
       'func' => '.',
       'p' => '.',
       'polar' => '.',
       'r' => '.'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'GAGAG',
       'bal' => 5,
       'bp' => 32083945,
       'br' => 'G',
       'brl' => 1
     },
     '-' => {
       'ba' => 'TAGAG',
       'bal' => 5,
       'bp' => 32083943,
       'br' => 'T',
       'brl' => 1
     },
     'alt' => 'TAGAGAG',
     'alt_cn' => 3,
     'altlen' => 7,
     'chr' => 'chr19',
     'guess' => 'rep',
     'pos' => 32083943,
     'ref' => 'TAG',
     'ref_cn' => 1,
     'reflen' => 3,
     'rep' => 'AG',
     'replen' => 2,
     'sel' => {
       'bc' => [
         [
           'NM_001205273.1',
           {
             'bd' => 'B',
             'cpos' => 'c.447',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           },
           {
             'bd' => 'B',
             'cpos' => '++1',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           }
         ]
       ],
       'std' => [
         [
           'NM_001205273.1',
           {
             'bd' => 'B',
             'cpos' => 'c.446',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           }
         ]
       ]
     }
   }
};
my $rep_insertion_dup_2_varanno={
'info' => {
     'NM_001205273.1' => {
       'bc' => 'c.446_447insAAAA',
       'c' => 'c.446A[1>5]',
       'cc' => '.',
       'exin' => 'EX2E',
       'strd' => '+',
       'flanks' => {
         'l' => 'chr19:32083943-32083944',
         'r' => 'chr19:32083945-32083946'
       },
       'func' => 'coding-synon',
       'p' => 'p.=',
       'polar' => '.',
       'r' => 'C2E'
     }
   },
   'var' => {
     '+' => {
       'ba' => 'AAAAA',
       'bal' => 5,
       'bp' => 32083944,
       'br' => 'A',
       'brl' => 1
     },
     '-' => {
       'ba' => 'TAAAA',
       'bal' => 5,
       'bp' => 32083943,
       'br' => 'T',
       'brl' => 1
     },
     'alt' => 'TAAAAA',
     'alt_cn' => 5,
     'altlen' => 6,
     'chr' => 'chr19',
     'guess' => 'rep',
     'pos' => 32083943,
     'ref' => 'TA',
     'ref_cn' => 1,
     'reflen' => 2,
     'rep' => 'A',
     'replen' => 1,
     'sel' => {
       'bc' => [
         [
           'NM_001205273.1',
           {
             'bd' => 'B',
             'cpos' => 'c.446',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           },
           {
             'bd' => 'B',
             'cpos' => 'c.447',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           }
         ]
       ],
       'std' => [
         [
           'NM_001205273.1',
           {
             'bd' => 'B',
             'cpos' => 'c.446',
             'exin' => 'EX2E',
             'gsym' => 'THEG5',
             'reg' => 'C2E',
             'strd' => '+'
           }
         ]
       ]
     }
   }
 };

my $snv_del_individual={
   'NR_024540.1' => {
     'c' => 'n.[1721G>T];[1714_1720delCGATCTA]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14408-14409 chr1:14409-14410',
       'r' => 'chr1:14411-14412 chr1:14418-14419'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
};
my $insert_rep_individual={
   'NR_024540.1' => {
     'c' => 'n.[1718_1720delCTAinsTAGCTATTC];[1718CTA[1>4]]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14409-14410',
       'r' => 'chr1:14414-14415 chr1:14411-14412'
      },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
 };
my $snv_insert_individual={
   'NR_024540.1' => {
     'c' => 'n.[1721G>T];[1718_1720delCTAinsTAGCTATTC]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14408-14409 chr1:14409-14410',
       'r' => 'chr1:14411-14412 chr1:14414-14415'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
 };
my $insert_dup_individual={
   'NR_024540.1' => {
     'c' => 'n.[1718_1720delCTAinsTAGCTATTC];[1719_1720dupTA]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14409-14410',
       'r' => 'chr1:14414-14415 chr1:14411-14412'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
 };
my $del_rep_individual={
 'NR_024540.1' => {
     'c' => 'n.[1714_1720delCGATCTA];[1718CTA[1>4]]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14409-14410',
       'r' => 'chr1:14418-14419 chr1:14411-14412'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
  }
};
my $del_dup_individual={
   'NR_024540.1' => {
     'c' => 'n.[1714_1720delCGATCTA];[1719_1720dupTA]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14409-14410',
       'r' => 'chr1:14418-14419 chr1:14411-14412'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
 };
my $del_delins_individual={
 'NR_024540.1' => {
     'c' => 'n.[1714_1720delCGATCTA];[1717_1721delTCTAGinsAC]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14408-14409',
       'r' => 'chr1:14418-14419 chr1:14415-14416'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   } 
};
my $del_n_individual={
'NR_024540.1' => {
     'c' => 'n.[1714_1720delCGATCTA];[1720_1721insNTAN]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410 chr1:14409-14410',
       'r' => 'chr1:14418-14419 chr1:14411-14412'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
};
my $rep_dup_individual={
 'NR_024540.1' => {
     'c' => 'n.[1718CTA[1>4]];[1719_1720dupTA]',
     'cc' => '.',
     'exin' => 'EX11E',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr1:14409-14410',
       'r' => 'chr1:14411-14412'
     },
     'func' => 'ncRNA',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'R11E'
   }
};
my $splice_individual={
'NM_006158.4' => {
     'c' => 'c.[1408C[7>6]];[1408C[7>6]]',
     'cc' => '.',
     'exin' => 'IVS3',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr8:24811063-24811064',
       'r' => 'chr8:24811066-24811067'
     },
     'func' => 'abnormal-intron',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'IC3'
   }
 };
my $splice1_individual={
'NM_024022.2' => {
     'c' => 'c.[617-3_617-2dupTA];[617-4_617-3delCT]',
     'cc' => '.',
     'exin' => 'IVS7',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr21:43803307-43803308 chr21:43803308-43803309',
       'r' => 'chr21:43803309-43803310 chr21:43803312-43803313'
     },
     'func' => '[splice-3];[intron]',
     'keep' => 1,
     'p' => '.',
     'polar' => '.',
     'r' => 'IC6'
   }
 };
my $cds_individual={
'NM_000581.2' => {
     'c' => 'c.[28_38delGCGGCGGCGGCinsGGCGGCGGCGGC];[28_34delGCGGCGGinsGG]',
     'cc' => '.',
     'exin' => 'EX1',
     'strd' => '-',
     'flanks' => {
       'l' => 'chr3:49395672-49395673 chr3:49395676-49395677',
       'r' => 'chr3:49395685-49395686 chr3:49395685-49395686'
     },
     'func' => 'frameshift',
     'keep' => 1,
     'p' => 'p.[Ala10_Ala13delinsGlyGlyGlyGlyProfs*?];[Ala10_Ala12delinsGlyfs*?]',
     'polar' => '.',
     'r' => 'C1'
   }
 };
my $utr5_individual={
'NM_152486.2' => {
     'c' => 'c.[-1T>G];[-1_2delTTA]',
     'cc' => '.',
     'exin' => 'EX2',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr1:861319-861320 chr1:861319-861320',
       'r' => 'chr1:861322-861323 chr1:861324-861325'
     },
     'func' => '[utr-5];[init-loss]',
     'keep' => 1,
     'p' => '[.];[p.0?]',
     'polar' => '.',
     'r' => '[5U1];[5U1-C1]'
 }
};
my $cdna_scoden_individual={
'NM_152486.2' => {
     'c' => 'c.[1A>C];[1delAinsGC]',
     'cc' => '[ATG=>CTG];[.]',
     'exin' => 'EX2',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr1:861320-861321',
       'r' => 'chr1:861323-861324'
     },
     'func' => '[misstart];[init-loss]',
     'keep' => 1,
     'p' => 'p.[0?];[0?]',
     'polar' => '.',
     'r' => 'C1'
   }
 };
my $cdna_ecoden_individual={
 'NM_152486.2' => {
     'c' => 'c.[2046A>T];[2046A[1>4]]',
     'cc' => '[TGA=>TGT];[.]',
     'exin' => 'EX14E',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr1:879531-879532 chr1:879532-879533',
       'r' => 'chr1:879534-879535 chr1:879534-879535'
     },
     'func' => '[stop-loss];[utr-3]',
     'keep' => 1,
     'p' => '[p.*682Cysext*?];[.]',
     'polar' => '[.=>P0];[.]',
     'r' => '[C13E];[C13E-3U1E]'
   }
 };
my $intron_individual={
'NM_000247.1-2' => {
     'c' => 'c.[948+1_949-1delGCTG];[948_948+1insA]',
     'cc' => '.',
     'exin' => '[IVS5];[EX5-IVS5]',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr6_cox_hap2:2892856-2892857 chr6_cox_hap2:2892856-2892857',
       'r' => 'chr6_cox_hap2:2892862-2892863 chr6_cox_hap2:2892858-2892859'
     },
     'func' => '[splice];[frameshift]',
     'keep' => 1,
     'p' => '[.];[p.Ala317fs*?]',
     'polar' => '.',
     'r' => '[IC5];[C5-IC5]'
   }
 };
my $utr5_other_individual={
'NM_005101.3' => {
     'c' => 'c.[-2_3delCCATG];[1_3delATG]',
     'cc' => '.',
     'exin' => 'EX1',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr1:948950-948951 chr1:948952-948953',
       'r' => 'chr1:948957-948958 chr1:948957-948958'
     },
     'func' => 'init-loss',
     'keep' => 1,
     'p' => 'p.[0?];[0?]',
     'polar' => '.',
     'r' => '[5U1E-C1];[C1]'
   }
 };
my $ins_intron_individual={
'NM_017780.3' => {
     'c' => 'c.[2377-1_2377insTAA];[2377_2378insGGA]',
     'cc' => '.',
     'exin' => '[IVS5-EX6];[EX6]',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr8:61714085-61714086 chr8:61714086-61714087',
       'r' => 'chr8:61714087-61714088 chr8:61714088-61714089'
     },
     'func' => '[stop-gain];[cds-indel]',
     'keep' => 1,
     'p' => 'p.[Gln792_Gln792ins*];[Glu793delinsGlyLys]',
     'polar' => '.',
     'r' => '[IC4-C5];[C5]'
   }
 };
my $rep_insertion_dup_individual={
'NM_001205273.1' => {
     'c' => 'c.[446AG[1>3]];[446A[1>5]]',
     'cc' => '.',
     'exin' => '[.];[EX2E]',
     'strd' => '+',
     'flanks' => {
       'l' => 'chr19:32083944-32083945 chr19:32083943-32083944',
       'r' => 'chr19:32083946-32083947 chr19:32083945-32083946'
     },
     'func' => '[.];[coding-synon]',
     'keep' => 1,
     'p' => '[.];[p.=]',
     'polar' => '.',
     'r' => '[.];[C2E]'
   }
 };

test_anno( $ref_anno,  "NM_006996.2",    "chr1", 169454923, "C",  "C" );
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

test_parse_var($snv_parse,"chr1", 14410, "C", "A");
test_parse_var($insert_parse,"chr1",14410,"CTAG","CGAATAGCTA");
test_parse_var($del_parse,"chr1", 14410, "CTAGATCG", "C");
test_parse_var($rep_parse,"chr1",14410,"CTAGA","CTAGTAGTAGTAGA");
test_parse_var($dup_parse,"chr1",14410,"CTAGA","CTATAGA");
test_parse_var($delins_parse,"chr1",14410,"CTAGA","GT");
test_parse_var($n_parse,"chr1",14410,"C","CNTAN");
test_parse_var($splice_parse_1, "chr8", 24811064, "AGGGGGGG","AGGGGGG" );
test_parse_var($splice_parse_2, "chr8", 24811068, "GGGG","GAGGG");
test_parse_var($splice1_parse_1, "chr21", 43803308, "CTA", "CTATA");
test_parse_var($splice1_parse_2, "chr21", 43803309, "TAG", "T");
test_parse_var($cds_parse_1, "chr3", 49395673, "GGCCGCCGCCGC","GGCCGCCGCCGCC");
test_parse_var($cds_parse_2, "chr3", 49395677, "GCCGCCGC","GCC" );
test_parse_var($utr5_parse_1, "chr1", 861320, "TT", "TG");
test_parse_var($utr5_parse_2, "chr1", 861320, "TTTA", "T");
test_parse_var($cdna_scoden_parse_1, "chr1", 861321, "TA", "TC");
test_parse_var($cdna_scoden_parse_2, "chr1", 861322, "A", "GC");
test_parse_var($cdna_ecoden_parse_1,  "chr1", 879532, "GA", "GT");
test_parse_var($cdna_ecoden_parse_2,  "chr1", 879532, "GA", "GAAAA");
test_parse_var($intron_parse_1,"chr6_cox_hap2", 2892857,"TGCTG", "T" );
test_parse_var($intron_parse_2,"chr6_cox_hap2", 2892857,"TGCTG", "TAGCTG");
test_parse_var($utr5_del_parse_1, "chr1", 948951, "GCCATG","G" );
test_parse_var($utr5_del_parse_2, "chr1", 948953, "CATG", "C");
test_parse_var($ins_intron_exon3_parse_1, "chr8", 61714086, "G","GTAA");
test_parse_var($ins_intron_exon3_parse_2, "chr8", 61714087, "T","TGGA");
test_parse_var($rep_insertion_dup_parse_1, "chr19", 32083943, "TAG", "TAGAGAG");
test_parse_var($rep_insertion_dup_parse_2, "chr19", 32083943, "TA", "TAAAAA");

test_gHGVS($snv_gHGVS,"SNV_gHGVS",$snv_parse);
test_gHGVS($insert_gHGVS,"Insert_gHGVS",$insert_parse);
test_gHGVS($del_gHGVS,"deletion_gHGVS",$del_parse);
test_gHGVS($rep_gHGVS,"repeat_gHGVS",$rep_parse);
test_gHGVS($dup_HGVS,"duplication_gHGVS",$dup_parse);
test_gHGVS($delins_HGVS,"delins_gHGVS",$delins_parse);
test_gHGVS($n_HGVS,"NNN_gHGVS",$n_parse);
test_gHGVS($splice_1_gHGVS,"splice_1_gHGVS",$splice_parse_1);
test_gHGVS($splice_2_gHGVS,"splice_2_gHGVS",$splice_parse_2);
test_gHGVS($splice1_1_gHGVS,"splice1_1_gHGVS",$splice1_parse_1);
test_gHGVS($splice1_2_gHGVS,"splice1_2_gHGVS",$splice1_parse_2);
test_gHGVS($cds_1_gHGVS,"cds_1_gHGVS",$cds_parse_1);
test_gHGVS($cds_2_gHGVS,"cds_2_gHGVS",$cds_parse_2);
test_gHGVS($utr5_1_gHGVS,"utr5_1_gHGVS",$utr5_parse_1);
test_gHGVS($utr5_2_gHGVS,"utr5_2_gHGVS",$utr5_parse_2);
test_gHGVS($cdna_scoden_1_gHGVS,"cdna_1_gHGVS",$cdna_scoden_parse_1);
test_gHGVS($cdna_scoden_2_gHGVS,"cdna_2_gHGVS",$cdna_scoden_parse_2);
test_gHGVS($cdna_ecoden_1_gHGVS,"cdna_ecoden_1_gHGVS",$cdna_ecoden_parse_1);
test_gHGVS($cdna_ecoden_2_gHGVS,"cdna_ecoden_2_gHGVS",$cdna_ecoden_parse_2);
test_gHGVS($intron_1_gHGVS,"intron_1_gHGVS",$intron_parse_1);
test_gHGVS($intron_2_gHGVS,"intron_2_gHGVS",$intron_parse_2);
test_gHGVS($utr5_del_1_gHGVS,"utr5_del_1_gHGVS",$utr5_del_parse_1);
test_gHGVS($utr5_del_2_gHGVS,"utr5_del_2_gHGVS",$utr5_del_parse_2);
test_gHGVS($ins_intron_exon3_1_gHGVS,"ins_intron_exon3_1_gHGVS",$ins_intron_exon3_parse_1);
test_gHGVS($ins_intron_exon3_2_gHGVS,"ins_intron_exon3_2_gHGVS",$ins_intron_exon3_parse_2);
test_gHGVS($rep_insertion_dup_1_gHGVS,"rep_insertion_dup_1_gHGVS",$rep_insertion_dup_parse_1);
test_gHGVS($rep_insertion_dup_2_gHGVS,"rep_insertion_dup_2_gHGVS",$rep_insertion_dup_parse_2);

test_varanno($snv_varanno,"snv_varanno",$snv_parse);
test_varanno($insert_varanno,"insert_varanno",$insert_parse);
test_varanno($del_varanno,"deletion_varanno",$del_parse);
test_varanno($rep_varanno,"repeat_varanno",$rep_parse);
test_varanno($dup_varanno,"duplication_varanno",$dup_parse);
test_varanno($delins_varanno,"delins_varanno",$delins_parse);
test_varanno($n_varanno,"NN_varanno",$n_parse);
test_varanno($splice_1_varanno,"splice_1_varanno",$splice_parse_1);
test_varanno($splice_2_varanno,"splice_2_varanno",$splice_parse_2);
test_varanno($splice1_1_varanno,"splice1_1_varanno",$splice1_parse_1);
test_varanno($splice1_2_varanno,"splice1_2_varanno",$splice1_parse_2);
test_varanno($cds_1_varanno,"cds_1_varanno",$cds_parse_1);
test_varanno($cds_2_varanno,"cds_2_varanno",$cds_parse_2);
test_varanno($utr5_1_varanno,"utr5_1_varanno",$utr5_parse_1);
test_varanno($utr5_2_varanno,"utr5_2_varanno",$utr5_parse_2);
test_varanno($cdna_scoden_1_varanno,"cdna_scoden_1_varanno",$cdna_scoden_parse_1);
test_varanno($cdna_scoden_2_varanno,"cdna_scoden_2_varanno",$cdna_scoden_parse_2);
test_varanno($cdna_ecoden_1_varanno,"cdna_ecoden_1_varanno",$cdna_ecoden_parse_1);
test_varanno($cdna_ecoden_2_varanno,"cdna_ecoden_2_varanno",$cdna_ecoden_parse_2);
test_varanno($intron_1_varanno,"intron_1_varanno",$intron_parse_1);
test_varanno($intron_2_varanno,"intron_2_varanno",$intron_parse_2);
test_varanno($utr5_del_1_varanno,"utr5_del_1_varanno",$utr5_del_parse_1);
test_varanno($utr5_del_2_varanno,"utr5_del_2_varanno",$utr5_del_parse_2);
test_varanno($ins_intron_exon3_1_varanno,"ins_intron_exon3_1_varanno",$ins_intron_exon3_parse_1);
test_varanno($ins_intron_exon3_2_varanno,"ins_intron_exon3_2_varanno",$ins_intron_exon3_parse_2);
test_varanno($rep_insertion_dup_1_varanno,"rep_insertion_dup_1_varanno",$rep_insertion_dup_parse_1);
test_varanno($rep_insertion_dup_2_varanno,"rep_insertion_dup_2_varanno",$rep_insertion_dup_parse_2);


test_individual($snv_del_individual,"snv_Combin_del_individual",$snv_varanno,$del_varanno);
test_individual($insert_rep_individual,"insert_Combin_repeat_individual",$insert_varanno,$rep_varanno);
test_individual($snv_insert_individual,"snv_Combin_insert_individual",$snv_varanno,$insert_varanno);
test_individual($insert_dup_individual,"insert_Combin_duplication_individual",$insert_varanno,$dup_varanno);
test_individual($del_rep_individual,"deletion_Combin_repeat",$del_varanno,$rep_varanno);
test_individual($del_dup_individual,"deletion_Combin_duplication",$del_varanno,$dup_varanno);
test_individual($del_delins_individual,"deletion_Combin_delins",$del_varanno,$delins_varanno);
test_individual($del_n_individual,"deletion_Combin_N",$del_varanno,$n_varanno);
test_individual($rep_dup_individual,"repeat_Combin_duplication",$rep_varanno,$dup_varanno);
test_individual($splice_individual,"splice_individual",$splice_1_varanno,$splice_1_varanno);
test_individual($splice1_individual,"splice1_individual",$splice1_1_varanno,$splice1_2_varanno);
test_individual($cds_individual,"cds_individual",$cds_1_varanno,$cds_2_varanno);
test_individual($utr5_individual,"utr5_individual",$utr5_1_varanno,$utr5_2_varanno);
test_individual($cdna_scoden_individual,"cdna_scoden_individual",$cdna_scoden_1_varanno,$cdna_scoden_2_varanno);
test_individual($cdna_ecoden_individual,"cdna_ecoden_individual",$cdna_ecoden_1_varanno,$cdna_ecoden_2_varanno);
test_individual($intron_individual,"intron_individual",$intron_1_varanno,$intron_2_varanno);
test_individual($utr5_other_individual,"utr5_other_individual",$utr5_del_1_varanno,$utr5_del_2_varanno);
test_individual($ins_intron_individual,"ins_intron_individual",$ins_intron_exon3_1_varanno,$ins_intron_exon3_2_varanno);
test_individual($rep_insertion_dup_individual,"rep_insertion_dup_individual",$rep_insertion_dup_1_varanno,$rep_insertion_dup_2_varanno);
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

sub test_parse_var {
    my ($expect,@fourargs) = @_;
    my $ranno = parse_var(@fourargs);
    if (!is_deeply($ranno, $expect, "for [ @fourargs ]")) {
        explain "The anno infomations are: ", $ranno;
    }
}
sub test_gHGVS {
        my ($expect,$tid,$vara)=@_;
        my $TgHGVS = get_gHGVS($vara);
        if (!is_deeply($TgHGVS,$expect, "for [ $tid]")) {
                explain "The anno infomations are: ",$TgHGVS;
        }
}

sub test_varanno{
        my ($expect,$tid,$vara)=@_;
        my $vAnno = $beda->varanno($vara);
        if (!is_deeply ($vAnno,$expect,"for [$tid]")){
        explain "The anno infomations are: ", $vAnno;
        }
}
sub test_individual{
        my ($expect,$tid,$var1,$var2)=@_;
        my $IndAnno = individual_anno($var1, $var2);
        if (!is_deeply($IndAnno,$expect,"for[$tid]")){
        explain "The anno infomations are: ", $IndAnno;
        }
}
# test if the tabix is available
sub test_tabix {
    my ($file) = @_;
    my $cmd;
    $cmd = "tabix -f -p bed $file";
    system($cmd);
    if ( !is( $?, 0, "tabix .. $cmd" ) ) {
        exit 1;
    }
}

