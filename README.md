BedAnno
=======

Annotate genomics variations of hg19 by using a BED+1 format database, which construct from ucsc hg19 databases.

BED +1 Format for hg19 refseq annotation
----------------------------------------
The start position is in 0 based, stop in 1 based
Tag parsing rules: Entries are separated by "; ", and for tags in entry are separated by "|"

**Tags are:**

- Acc.Ver[-submap]
- Gene Symbol
- Strand
- BlockAttr
- ExIn Num
- n./r. HGVS start for block before departing
- n./r. HGVS end for block before departing
- c. HGVS start for block before departing.
- c. HGVS end for block before departing.
- Original block length
- Primary Tag (If "Y", the primary transcript record in all transcripts of the same gene, otherwise "N")
- Offset to leftmost of non departing block. 

The BlockAttr and the ExIn Num are defined as the following:

    For entry with cmpl 5' cds and cmpl 3' cds:

	       5'UTR   =======================================>   3'UTR

	       5U3E,I5U2,5U2,I5U1,5U1,C1,IC1, C2, IC2,C3E,3U1,I3U1,3U2E
	       EX1  IVS1 EX2 IVS2   EX3  IVS3 EX4 IVS4  EX5   IVS5 EX6E

    For entry with incmpl 5' cds and cmpl 3' cds:

               5'UTR ======> 5'CDS   ===================>   3'UTR

                             C-3P,IC-2,C-2,IC-1,C-1,3U1,I3U1,3U2E
                             EX-4P     EX-3      EX-2        EX-1
                                  IVS-3    IVS-2        IVS-1

               5U2E,I5U1,5U1,C-3E,...............................
               EX-5E       EX-4        EX-3      EX-2        EX-1
                    IVS-4         IVS-3    IVS-2        IVS-1

    For entry with incmpl 3' cds and cmpl 5' cds:

               5'UTR   ===================>   3'CDS ======> 3'UTR

               5U2E,I5U1,5U1,C+1,IC+1,C+2,IC+2,C+3P
               EX+1       EX+2        EX+3     EX+4P
                    IVS+1        IVS+2    IVS+3

               ................................C+3E,3U1,I3U1,3U2E
               EX+1       EX+2        EX+3       EX+4        EX+5E
                    IVS+1        IVS+2    IVS+3         IVS+4

    For entry with incmpl 3' cds and incmpl 5' cds:

              5'UTR ======> 5'CDS   =======>  3'CDS ======> 3'UTR
     
                            C*1P,IC*1,C*2,IC*2,C*3P
                            EX*1P     EX*2     EX*3P
                                 IVS*1    IVS*2
     
              5U2E,I5U1,5U1,C*1, IC*1,C*2,IC*2,C*3E,3U1,I3U1,3U2E
              EX*1       EX*2         EX*3       EX*4        EX*5E
                   IVS*1         IVS*2    IVS*3         IVS*4

    For non-coding RNA

              5' =========> 3'

              R1 IR1 R2 IR2 R3
              EX1    EX2    EX3
                 IVS1   IVS2

**Example entries:**

    chr1    879077  879188  NM_152486.2|SAMD11|+|C12|EX13|1770|1880|1690|1800|111|Y|0
    chr1    879188  879287  NM_152486.2|SAMD11|+|IC12|IVS13|1880|1881|1800|1801|99|Y|0
    chr1    879287  879533  NM_152486.2|SAMD11|+|C13E|EX14E|1881|2126|1801|2046|246|Y|0
    chr1    879533  879582  NM_152486.2|SAMD11|+|3U1E|EX14E|2127|2554|*1|*428|428|Y|0
    chr1    879582  879961  NM_015658.3|NOC2L|-|3U1E|EX19E|2800|2310|*491|*1|491|Y|0; NM_152486.2|SAMD11|+|3U1E|EX14E|2127|2554|*1|*428|428|Y|49
    chr1    879961  880073  NM_015658.3|NOC2L|-|3U1E|EX19E|2800|2310|*491|*1|491|Y|379
    chr1    880073  880180  NM_015658.3|NOC2L|-|C19E|EX19E|2309|2203|2250|2144|107|Y|0
    chr1    880180  880436  NM_015658.3|NOC2L|-|IC18|IVS18|2203|2202|2144|2143|256|Y|0
    chr1    880436  880526  NM_015658.3|NOC2L|-|C18|EX18|2202|2113|2143|2054|90|Y|0

For the primary tag definition, there's a cooresponding sort strategy:

 *Rules are list from prior to minor.*

- whether the transcript is on the primary assembly. "On" is prior.
- concatenated CDS is longer
- concatenated Exon is longer
- union Exon with flank region is longer
- Chromosome ID number is smaller
- Position number on forward strand-chromosome is smaller.

Then assign the 'best' NMid to the first record, and add postfix "-N" (1..n-1) to the following records

