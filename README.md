BedAnno
=======

Annotate genomics variations of hg19 by using a BED format database, 
which construct from ncbi anno release 104. This module can directly 
parse the vcf4.1 format ref and single alt string(no commas in it),
without normalized by vcftools, and can recognize the tandom repeat 
variation and duplication, generate the standard HGVS strings for 
most of complex cases. Also it will ajust the strand of transcript,
and follow the 3' nearest rules to annotate. CG's variant shell list
will also be supported.

BED Format for ncbi annotation rel104
----------------------------------------
The start position is in 0 based, stop in 1 based
Tag parsing rules: Entries are separated by "; ", and for tags in entry are separated by "|"

**Tags are:**

1.  Acc.Ver
2.  GeneID
3.  Gene Symbol
4.  Strand
5.  BlockAttr
6.  GenePartsSO
7.  ExIn Num
8.  nHGVS start for block before departing
9.  nHGVS end for block before departing
10. cHGVS start for block before departing
11. cHGVS end for block before departing
12. Length for block before departing
13. MismatchBlock :  $type,$gstart,$gstop,$gseq 
                     ($gseq is in the strand of refseq, '.' for deletion)
14. Primary Tag   :  see [PRIMARY TAG ASSIGNMENT][1] at the bottom.
15. Offset to leftmost of non departing block.

*The BlockAttr, GenePartsSO, and ExIn Num are defined as following:*

     5'==|>>>>|[============]|>>>|>>>|[=========]|>>>|>>>|[============]|>>>>|==3'
     PROM 5U2E D5U1 I5U1 A5U1 5U1 C1  DC1 IC1 AC1 C2E 3U1 D3U1 I3U1 A3U1 3U2E
     167  204  163  447  164  204 316 163 191 164 316 205 163  448  164  205
     .   |EX1 |     IVS1     |  EX2  |   IVS2    |  EX3  |    IVS3      |EX4E|               


**Example:**

    1       155252631       155252633       NM_020897.2|HCN3|57657|+|DC2|163|IVS2|872+1|872+2|708+1|708+2|2||Y|0; NR_073074.1|HCN3|57657|+|DR2|163|IVS2|872+1|872+2|||2||N|0
    1       155252633       155253762       NM_020897.2|HCN3|57657|+|IC2|191|IVS2|872+3|873-3|708+3|709-3|1129||Y|0; NR_073074.1|HCN3|57657|+|IR2|191|IVS2|872+3|873-3|||1129||N|0
    1       155253762       155253764       NM_020897.2|HCN3|57657|+|AC2|164|IVS2|873-2|873-1|709-2|709-1|2||Y|0; NR_073074.1|HCN3|57657|+|AR2|164|IVS2|873-2|873-1|||2||N|0
    1       155253764       155253926       NM_020897.2|HCN3|57657|+|C3|316|EX3|873|1034|709|870|162||Y|0; NR_073074.1|HCN3|57657|+|R3|655|EX3|873|1034|||162||N|0
    1       155253926       155253928       NM_020897.2|HCN3|57657|+|DC3|163|IVS3|1034+1|1034+2|870+1|870+2|2||Y|0; NR_073074.1|HCN3|57657|+|DR3|163|IVS3|1034+1|1034+2|||2||N|0
    1       155253928       155254327       NM_020897.2|HCN3|57657|+|IC3|191|IVS3|1034+3|1035-3|870+3|871-3|399||Y|0; NR_073074.1|HCN3|57657|+|IR3|191|IVS3|1034+3|1035-3|||532||N|0
    1       155254327       155254329       NM_020897.2|HCN3|57657|+|AC3|164|IVS3|1035-2|1035-1|871-2|871-1|2||Y|0; NR_073074.1|HCN3|57657|+|IR3|191|IVS3|1034+3|1035-3|||532||N|399
    1       155254329       155254460       NM_020897.2|HCN3|57657|+|C4|316|EX4|1035|1253|871|1089|219||Y|0; NR_073074.1|HCN3|57657|+|IR3|191|IVS3|1034+3|1035-3|||532||N|401
    1       155254460       155254462       NM_020897.2|HCN3|57657|+|C4|316|EX4|1035|1253|871|1089|219||Y|131; NR_073074.1|HCN3|57657|+|AR3|164|IVS3|1035-2|1035-1|||2||N|0
    1       155254462       155254548       NM_020897.2|HCN3|57657|+|C4|316|EX4|1035|1253|871|1089|219||Y|133; NR_073074.1|HCN3|57657|+|R4|655|EX4|1035|1120|||86||N|0


TRANSCRIPT FASTA DATABASE FORMAT
--------------------------------

   One-line sequence fasta file
   ----------------------------
   Header format is: ( separate by " ", with "." for unavailable value )

       >rnaAcc.ver rnaLen gene protAcc.ver protLen cdsSta,cdsEnd tags [altStartCodons] [frameshiftSite]

   "tags" are a string of multiple properties separated by "|", which are:

       <alignStat>[|altstart][|selenocysteine][|inseqStop][|polyATail]

    when "tags" contain "altstart", alternative start codons will be list in "altStartCodons",
    separated by ";"

PRIMARY TAG ASSIGNMENT
----------------------

Our database need to sort the refSeq record with same Acc,
or with same gene by the following rules:

1. whether the transcript is on the primary assembly.
2. whether the CDS is 3-codons (badcds will be put to the tail)
3. concatenated CDS is longer
4. concatenated Exon is longer
5. union Exon with flank region is longer
6. Chromosome ID number is smaller
7. Position number on forward strand-chromosome is smaller

if LRG\_RefSeqGene file is used, then the primary tag "Y" will
only assign to the reference standard transcripts's first
mapping entry. But for non-RefSeqGene gene/transcript,
assign primary tag "Y" to the first record for same genes' all
transcripts, which is the same with no LRG\_RefSeqGene file case.

For multiple-mapping of a same transcript, add postfix "-N"
(1..n-1) to the other records in the order of sort.

[1]: #primary-tag-assignment
