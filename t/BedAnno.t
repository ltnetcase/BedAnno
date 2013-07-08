# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl BedAnno.t'

our $data;
BEGIN {
  unless(grep /blib/, @INC) {
    chdir 't' if -d 't';
    unshift @INC, '../lib' if -d '../lib';
    $data = '../data';
  }
}
$data ||= "data";

use Test::Most;
BEGIN { use_ok('BedAnno') };

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
    c        => "c.82G>T",
    p        => "p.Glu28*",
    cc       => "GAA=>TAA",
    r        => "C1",
    exin     => "EX1",
    func     => "nonsense",
    polar    => "P-=>.",
    bc	     => '.',
    flanks   => {
        l    => "chr1:169454921-169454922",
        r    => "chr1:169454924-169454925",
        strd => '-'
    }
};

my $del1_anno = {
    c	     => 'c.757delA',
    p	     => 'p.Ile253Serfs*?',
    cc	     => '.',
    r	     => 'C1E',
    exin     => 'EX1E',
    func     => 'frameshift',
    polar    => '.',
    bc	     => '.',
    flanks   => {
	l    => "chr1:248525637-248525638",
	r    => "chr1:248525640-248525641",
	strd => '+'
    }
};

my $splice_anno = {
    c	     => 'c.1408C[7>6]',
    p	     => '.',
    cc	     => '.',
    r	     => 'IC3',
    exin     => 'IVS3',
    func     => 'abnormal-intron',
    polar    => '.',
    bc	     => 'c.1413+1delC',
    flanks   => {
	l    => "chr8:24811063-24811064",
	r    => "chr8:24811066-24811067",
	strd => '-'
    }
};

my $splice1_anno = {
    c	     => 'c.617-3_617-2dupTA',
    p	     => '.',
    cc	     => '.',
    r	     => 'IC6',
    exin     => 'IVS7',
    func     => 'splice-3',
    polar    => '.',
    bc	     => 'c.617-2_617-1insTA',
    flanks   => {
	l    => "chr21:43803307-43803308",
	r    => "chr21:43803309-43803310",
	strd => '-'
    }
};

my $ins_anno = {
    c	     => 'c.338_339insTGTTACGCAGGAGAC',
    p	     => 'p.Thr113_Ala114insValThrGlnGluThr',
    cc	     => '.',
    r	     => 'C2',
    exin     => 'EX2',
    func     => 'cds-indel',
    polar    => '.',
    bc	     => '.',
    flanks   => {
	l    => "chr3:195518111-195518112",
	r    => "chr3:195518113-195518114",
	strd => '-'
    }
};

my $intr_anno = {
    c	     => 'c.83-12787_83-12786insTGTTACGCAGGAGAC',
    p	     => '.',
    cc	     => '.',
    r	     => 'IC1',
    exin     => 'IVS1',
    func     => 'intron',
    polar    => '.',
    bc	     => '.',
    flanks   => {
	l    => "chr3:195518111-195518112",
	r    => "chr3:195518113-195518114",
	strd => '-'
    }
};

my $delins_anno = {
    c	     => 'c.1192_1195delCAACinsC',
    p	     => 'p.Gln398_Arg399delinsArg',
    cc	     => '.',
    r	     => 'C5',
    exin     => 'EX7',
    func     => 'cds-indel',
    polar    => '.',
    bc	     => '.',
    flanks   => {
	l    => "chr22:38119753-38119754",
	r    => "chr22:38119759-38119760",
	strd => '+'
    }
};

my $rep_anno = {
    c	      => 'c.21GGC[6>4]',
    p	      => 'p.Ala11_Ala13delinsAla',
    cc	      => '.',
    r	      => 'C1',
    exin      => 'EX1',
    func      => 'cds-indel',
    polar     => '.',
    bc	      => 'c.33_38delGGCGGC',
    flanks    => {
	l     => "chr3:49395672-49395673",
	r     => "chr3:49395680-49395681",
	strd  => '-'
    }
};

test_anno($snv_anno, "NM_006996.2", "chr1", 169454923, "C", "A");
test_anno($del1_anno, "NM_001004696.1", "chr1", 248525638, "CA", "C");
test_anno($ins_anno, "NM_018406.6", "chr3", 195518112, "T", "TGTCTCCTGCGTAACA");
test_anno($intr_anno, "NM_004532.5", "chr3", 195518112, "T", "TGTCTCCTGCGTAACA");
test_anno($splice_anno, "NM_006158.4", "chr8", 24811064, "AGGGGGGG", "AGGGGGG");
test_anno($splice1_anno, "NM_024022.2", "chr21", 43803308,  "CTA", "CTATA");
test_anno($delins_anno, "NM_001039141.2", "chr22", 38119754, "TCAAC", "TC");
test_anno($rep_anno, "NM_000581.2", "chr3", 49395673, "GGCCGCCGCCGCCGCCGCC", "GGCCGCCGCCGCC");

done_testing();
exit 0;


sub test_anno {
    my ($expect, $tid, @fourargs) = @_;
    my $ranno = $beda->anno(@fourargs);
    if (!is_deeply($$ranno{info}{$tid}, $expect, "for [ $tid,".join(",", @fourargs)." ]")) {
	explain "The anno infomations are: ", $ranno;
    }
}

# test if the tabix is available
sub test_tabix {
    my ($file) = @_;
    my $cmd;
    $cmd = "tabix $file";
    system($cmd);
    if (!is($?, 0, "tabix .. $cmd")) {
	exit 1;
    }
}


