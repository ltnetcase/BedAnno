# NAME

BedAnno - Perl module for annotating variation depend on the BED format database.

## VERSION v1.30

From version 0.32 BedAnno will change to support CG's variant shell list
and use ncbi annotation release 104 as the annotation database
with reformatted database format, and won't give any individual
annotation, so the individual\_anno() is no longer available.
VCF4.1 format variant description (chr, pos, ref, alt) will also
be supported.

# SYNOPSIS

    use BedAnno;
    my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas' );
    my $anno = $beda->anno( 'chr20', 1234567, 1234569, 'AG', 'TGGG' );

# DESCRIPTION

By using this module, we can get variants from whole-genome or exome-capture 
NGS genotyping result annotated. The information contains various possible 
HGVS string together with a most recent strandard HGVS mutation name. 
It can not annotate ambiguous variants (transition, transvertion, or unknown 
break point large deletion and duplication).

_BedAnno_ annotate genomics variations of hg19 by using a BED format database, 
which construct from ncbi anno release 104, combined with tabix index.
This module can directly parse the vcf4.1 format ref and single alt string(no commas in it),
without normalized by vcftools, and can recognize the tandom repeat 
variation and duplication, generate the standard HGVS strings for 
most of complex cases. Also it will ajust the strand of transcript,
and follow the 3' nearest rules to annotate. 

# Methods

## new

- About : Creat a new annotation entry
- Usage :

        my $beda = BedAnno->new( db => "in.bed.gz", tr => 'in.trans.fas.gz', batch => 1 );

- Args    - (all database files should be tabix indexed)
    - Essential Args:
        - _db_ &lt;in.bed.gz>
            - annotation database file. 
        - _tr_ &lt;in.trans.fas.gz>
            - transcript sequence fasta file
        - See ["DATABASE FORMAT"](#database-format) for more infomation. 
    - Optional Args :
        - Common options :
            - _quiet_
                - Suppress warning messege to output.
            - _batch_ \[boolean\]
                - use batch mode annotation, default in daemon mode as an annotation engine.
            - _genome_ \[ "refgenome.fa.gz" \]
                - reference genome fasta, bgzipped and samtools faidxed for use.
            - _genes_ \[ "genes.list" | geneslist\_filehandle \]
                - annotate transcripts for _genes_. e.g. {"ABC" => 1, "DEF" => 1} or "genes.list" 
            - _trans_ \[ "trans.list" | translist\_filehandle \]
                - annotate transcripts in _trans_. e.g. {"NM\_0012.1" => 1, "NM\_0034.2" => 1} or "trans.list" 
            - _mmap_ \[boolean\]
                - allow annotating all other non "BEST" multiple-mapping records, boolen option, default not allowed.
                      e.g. NM\_0123.1 have 3 mapping location or alternate splicing, default only the "BEST" one will be 
                      annotated. See ["DATABASE FORMAT"](#database-format) for the "BEST" definition.
        - Batch mode options :
            - _region_ \[region\_string\]
                - limit to only annotate transcript in _region_. e.g. "chr20:1234567-1234568", prior to _regbed_.
            - _regbed_ \[BED format file\]
                - similar to _region_, with allowing multiple regions. e.g. "in.region.bed". 
        - Notes
            - Batch mode is designed for annotation of a complete list of variants 
                  on same chromosome, read all information of the chr into memory, 
                  and annotate all variants together in the order of chr coordinates.
                  This mode can avoid frequent IO, brought by tabix searching, but need
                  huge memory cost.
- Returns
    - Annotation Engine object entry, please see ["load\_anno"](#load_anno) for more information.

## set methods for properties

    List of Properties:
                        set
    db                  o
    tr                  o
    refbuild            o
    genome              o

    e.g.    : $beda->set_refbuild($refbuild);

## write\_using

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

## get\_cover\_batch

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

## readtr

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

## load\_anno

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

## region\_merge

    About   : merge consecutive same-entries regions
    Usage   : my $rannodb = region_merge($loaded_db);
    Args    : A hash ref of loaded_db.
    Returns : A hash ref of merged db.

## anno

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

## varanno

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
                        p3                => $threeletter_pHGVS,
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
                        primaryTag    => $refstandard_primary_or_not,   # Y/N
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
                    },
                    ...
                }
            }

## P1toP3

    About : Change 1 letter format of pHGVS string to 3 letter format
    Usage : my $p3 = P1toP3($p1);

## finaliseAnno

    About   : finalise the BedAnno::Anno entry by check all tag values,
              and uniform them for AE output usage, query transcript
              oringinated additional resources and add them into the data
              frame.
    Usage   : $beda->finaliseAnno($annEnt);
    Args    : BedAnno entry and a BedAnno::Anno entry
    Returns : A finalised BedAnno::Anno entry

## decide\_major

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

## getTrChange

    About   : Calculate the transcript changes, based on TrPostition
    Usage   : $beda->getTrChange($annoEnt);
    Returns : assign the following tags in annoEnt
                trRef, prot, c, p, cc, polar, func
                prRef, prAlt

## prWalker

    About   : act as trWalker, but on protein sequence
    Usage   : $prVar = $beda->prWalker( $prSeq,  $prBegin, $prEnd, $prRef, $prAlt );
    Args    : prSeq   - whole protein sequence
              prBegin - protein variant Begin position (1 based)
              prEnd   - protein variant End position (1 based)
              prRef   - protein variant reference sequence
              prAlt   - protein variant alt sequence

## trWalker

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
          

## cmpPos

    About   : judge the order for p1 and p2, because the insertion
              will have a reverted order of left & right position
    Usage   : my $cmpRst = BedAnno->cmpPos($p1, $p2);
    Args    : hgvs positio p1 and p2, with out 'c.' or 'n.' flag
    Return  : 0 for same, 1 for normal order, -1 for reversed order.

## \_getCPosFrame\_by\_cdsPos

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

## \_cdsubstr

    About   : substr from transcript seqeunce involving frameshift changing.
    Usage   : my $codonStr = _cdsubstr( $trdbEnt, $start, $length, $replace );
    Args    : trdbEnt is a sub hash in trInfodb which contains csta, csto.
              start is 0 based start position of sub string on transcript seq.
              length is the length of target region on transcript seq.
              replace is the alternative.
    Returns : A substring cut from the transcript sequence, with frameshift involved.
              Behave like function "substr", but for replace mode, 
              it returns the whole transcript seq after replacement.

## \_getPairedCodon

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

## translate

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

## getTrRef

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

## batch\_anno

    About   : The fastest way to annotate multiple snv and 1bp deletion variations,
              indel and other types also can be annotated, but no faster than annotated
              one by one.
    Usage   : $beda = BedAnno->new( db => 'in.bed.gz', tr => 'in.trans.fas', batch => 1);
              $rAnnoRst = $beda->batch_anno($rVars);
    Args    : an array ref of BedAnno::Var entries.
    Returns : an array ref of BedAnno::Anno entries, see varanno().

# BedAnno::Var

    BedAnno::Var sub module

# METHOD

## new

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

## getUnifiedVar

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

## parse\_complex

    About   : parse complex delins variants to recognize 
              repeat and differ strand-pos var.
    Usage   : my $var = $var->parse_complex();
    Args    : variantion entry, which have been uniform to 
              CG's shell list format, with its 'guess':delins.
    Returns : see BedAnno::Var->new()

## guess\_type

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

## get\_internal

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

## get\_gHGVS

    About   : get genomic (chromosomal) HGVS string of variation
    Usage   : my $gHGVS = $var->get_gHGVS();
    Args    : variation entry, after BedAnno::Var->new().
    Returns : chromosomal HGVS string.

# BedAnno::Anno

    BedAnno::Anno sub package

# METHOD

## new

    About   : Create BedAnno::Anno object
    Usage   : my $annoEnt = BedAnno::Anno->new($var);
    Args    : BedAnno::Var entry
    Returns : BedAnno::Anno entry

## getTrPosition

    About   : Assign the BedAnno::Anno obj's trInfo with affected regions
    Usage   : my $AEIndex = $annoEnt->getTrPosition($rannodb, $AEIndex);
    Args    : $rannodb is a chromosome branch of BedAnno object's annodb feature,
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

## cal\_hgvs\_pos

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
                 of it, the bad alignment string start with 'E'.
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

                 c. delins (S/DI) on refSeq (equal/non-equal length)

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

# BedAnno::CNV

    BedAnno::CNV sub package

# METHOD

## new

    About   : Create BedAnno::CNV object
    Usage   : my $cnva = BedAnno::CNV->new( db => $annodb, tr => $trdb, dgv => $dgvdb, cnvPub => $cnvPubdb );
    Args    : required args: db, tr
              Specific args:
                - ovlp_rate   overlapping region rate while searching for hits
                - max_uncov   maximum uncovered region length in querys while searching.
                - dgv         DGV tabix index file as Controls (got from anno_dgv)
                - sfari       SFARI tabix index file, as a case/control mix resource for autism.
                - cnvPub      Well formatted tabix indexed file as a collection of Cases
                - cnvd        CNVD database tabix indexed file.
              Same args in BedAnno method 'new':
                - db, tr, genes, trans, region, regbed, mmap, batch
    Returns : BedAnno::CNV object

## annoCNV

    About   : Annotate single CNV variants
    Usage   : my $rcnvAnno = $cnva->annoCNV( $chr, $start, $end, $copy_number );
    Args    : $chr   - chromosome number
              $start - 0 based start position
              $end   - 1 based end position
              $copy_number - estimated copy number
    Returns : A hash ref of cnv annotation in the following structure
                {
                    cnva_type => $cnv_anno_type,

                    # if hit on transcript
                    anno => {
                        $tid => {
                            gsym   => $gene_symbol,
                            gid    => $gene_id,
                            strd   => $strand,
                            cpos   => $hgvs_position,
                            exin   => $exon_intron_number,
                            regcod => $code_of_region,
                            regtyp => $type_of_region,
                        },
                        ...
                    },

                    # available when resource ok
                    cytoBand  => $cytoBand_info,
                    dgv    => $dgv_sql_rst,
                    sfari  => $sfari_sql_rst,
                    cnvd   => $cnvd_sql_rst,
                    cnvPub => $cnvPub_sql_rst,
                }

## batch\_annoCNV

    About   : Annotate CNV variants in batch mode
    Usage   : my $rcnvAnnos = $cnva->batch_annoCNV( $ref_hash_AllCNV );
    Args    : A hash ref of all cnv annotation in the following structure
                {
                    $chr => {
                        "$start-$stop" => $copy_number, 
                        # start is 1 based, to coordinate with region string format
                        ...
                      },
                      ...
                }
    Returns : A hash ref of annotated cnv variants in the following structure
                {
                    $chr => {
                        "$start-$stop" => {
                           cnva_type => $cnv_anno_type,
                           
                           # if hit on transcript
                           anno => {
                             $tid => {
                               gsym => $gene_symbol,
                               gid  => $gene_id,
                               strd => $strand,
                               cpos => $hgvs_position,
                               exin => $exon_intron_number,
                               regcod => $code_of_region,
                               regtyp => $type_of_region,
                             },
                             ...
                           },

                           # available when resource ok
                           cytoBand  => $cytoBand_info,
                           dgv => $dgv_sql_rst,
                           cnvPub => $cnvPub_sql_rst,
                        },
                        ...
                      },
                      ...
                }

# DATABASE FORMAT

The Format of annotation database is listed as following:

## BED FORMAT DATABASE

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

## TRANSCRIPT FASTA DATABASE

    One-line sequence fasta file
    ============================
    Header format is: ( separate by " ", with "." for unavailable value )

        >rnaAcc.ver rnaLen gene protAcc.ver protLen cdsSta,cdsEnd tags [altStartCodons] [frameshiftSite]
    
     "tags" are a string of multiple properties separated by "|", which are:

         <alignStat>[|altstart][|selenocysteine][|inseqStop][|polyATail]

     when "tags" contain "altstart", alternative start codons will be list in 
     "altStartCodons", separated by ";"

## PRIMARY TAG ASSIGNMENT

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

# SEE ALSO

    HGVS     :  http://www.hgvs.org/mutnomen/recs.html
    Mutalyzer:  https://mutalyzer.nl

# AUTHOR

liutao <liut@geneplus.org.cn>

# COPYRIGHT AND LICENSE

Please check LICENSE file for detail
