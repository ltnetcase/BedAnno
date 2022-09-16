=head1 BedAnno::Anno

    BedAnno::Anno sub package

=cut

package BedAnno::Anno;
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
    my $var   = shift;
    my $self;
    $self->{var} = $var;
    bless $self, ref($class) || $class;
    return $self;
}

sub TO_JSON {
    return { %{ shift() } };
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
    my ( $annoEnt, $rannodb, $aeIndex ) = @_;

    my $var = $annoEnt->{var};

    # gather the covered entries
    my $new_aeIndex;

    my %tidExblk = ();
    my %rpreLeft = ();

    # debug
    #    print STDERR "DEBUG: getTrPosition for ".Dumper($annoEnt);

    my %hit_badaln_ins = ();
    for ( my $k = $aeIndex ; $k < @$rannodb ; $k++ ) {
        if ( $$rannodb[$k]{sto} < $var->{pos} ) {    # not reach var
            $aeIndex++;
            next;
        }
        elsif ( $$rannodb[$k]{sta} > $var->{end} ) {    # past var
            last;
        }
        else {                                          # covered by var

      #	    # debug info
      #	    print STDERR "blk: $k [ $$rannodb[$k]{sta}, $$rannodb[$k]{sto} ]\n";

            if (    $$rannodb[$k]{sta} <= $var->{pos}
                and $var->{pos} <= $$rannodb[$k]{sto} )
            {
                $new_aeIndex = $k;

                # for edge hit case fetch backward one block.
                $new_aeIndex -= 1
                  if ( $new_aeIndex > 0 and $$rannodb[$k]{sta} == $var->{pos} );
            }

            if ( !exists $$rannodb[$k]{detail} ) {    # parse anno db
                $$rannodb[$k] = BedAnno->assign_detail( $$rannodb[$k] );
            }

            foreach my $tid ( keys %{ $$rannodb[$k]{detail} } ) {
                my $rtidDetail = $$rannodb[$k]{detail}{$tid};

                my $strd = ( $rtidDetail->{strd} eq '+' ) ? 1 : 0;

                my ( $unify_p, $unify_r, $unify_a, $unify_rl, $unify_al ) =
                  $var->getUnifiedVar( $rtidDetail->{strd} );

                # skip non hitted block
                next
                  if ( $unify_p > $$rannodb[$k]{sto}
                    or ( $unify_p + $unify_rl ) < $$rannodb[$k]{sta} );

                if (    $$rannodb[$k]{sta} == $$rannodb[$k]{sto}
                    and $unify_p + $unify_rl == $$rannodb[$k]{sto}
                    and $rtidDetail->{mismatch} =~ /^D/ )
                {
                    $hit_badaln_ins{$tid} = 1;
                }

                my $total_left_ofst =
                  $$rtidDetail{offset} + ( $unify_p - $$rannodb[$k]{sta} );
                my $total_right_ofst = $total_left_ofst + $unify_rl;

                # debug
                #		print STDERR join(" ", "tid: $tid", $rtidDetail->{strd},
                #		    $rtidDetail->{nsta}, $rtidDetail->{nsto},
                #		    $rtidDetail->{mismatch})."\n";
                #		print STDERR "rel: $total_left_ofst $total_right_ofst\n";

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
                            $tmp_trInfo->{trRefComp}->{P0}->[0] =
                              0;
                        }
                        elsif ( $rtidDetail->{exin} =~ /^IVS/ ) {
                            $tmp_trInfo->{trRefComp}->{ $rtidDetail->{exin} }
                              ->[0] = 0;
                        }
                        elsif ( !$strd and $tmp_trInfo->{rnaEnd} =~ /^\+(\d+)/ )
                        {
                            $tmp_trInfo->{trRefComp}->{Z999}->[0] = 0;
                            $tmp_trInfo->{trRefComp}->{Z999}->[1] = $1;
                        }

                        # let right rna positions to record Exon block

                        if ( $total_left_ofst > 0 ) {
                            $rpreLeft{$tid} = $annoEnt->cal_hgvs_pos(
                                tid       => $tid,
                                tidDetail => $rtidDetail,
                                offset    => $total_left_ofst,
                                LR        => 0,                  # preLeft mode
                                noassign  => 1,
                            );
                        }
                        elsif ( !exists $rpreLeft{$tid} ) {

                            # no left block of annotation for this tid
                            # then 5'promoter left or 3'downstream
                            $rpreLeft{$tid} = 0;
                        }

                        unless ( ref( $rpreLeft{$tid} ) ) {

                            # no left block of annotation for this tid
                            # then 5'promoter left or 3'downstream
                            if ( $total_left_ofst == 0 ) {
                                my $preLeft_cDot = "";
                                my $preLeft_nDot = "";
                                my $preLeft_exin = "";
                                my $preLeft_r    = "";
                                if ($strd) {
                                    $preLeft_exin = ".";
                                    $preLeft_r    = "PROM";
                                    if ( $tmp_trInfo->{cdsBegin} ne "" ) {
                                        if ( $tmp_trInfo->{cdsBegin} =~
                                            /^(.*)-u(\d+)$/ )
                                        {
                                            $preLeft_cDot =
                                              $1 . '-u' . ( $2 + 1 );
                                        }
                                        else { # for some sub region anno db entries
                                            $preLeft_cDot = "[uncertain preLeft]:" . $tmp_trInfo->{cdsBegin};
                                        }
                                    }
                                    if ( $tmp_trInfo->{rnaBegin} eq "1" ) {
                                        $preLeft_nDot = -1;
                                    }
                                    elsif ($tmp_trInfo->{rnaBegin} =~ /^\-?\d+$/) {
                                        $preLeft_nDot =
                                          $tmp_trInfo->{rnaBegin} - 1;
                                    }
                                    else { # for some sub region anno db entries
                                        $preLeft_nDot = "[uncertain preLeft]:" . $tmp_trInfo->{rnaBegin};    
                                    }
                                }
                                else {
                                    $preLeft_exin = ".";
                                    $preLeft_r    = "3D";
                                    $preLeft_nDot = '+1';
                                    $preLeft_cDot = '+d1';
                                }
                                $rpreLeft{$tid} = {
                                    nDot => $preLeft_nDot,
                                    cDot => $preLeft_cDot,
                                    exin => $preLeft_exin,
                                    r    => $preLeft_r,
                                };
                            }
                        }

                        if ( ref( $rpreLeft{$tid} ) ) {
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

                if ( $var->{pos} == $var->{end} )
                {    # assume all insertion are on transcript
                    $trinfoEnt->{staInTr} = 1;
                    $trinfoEnt->{stoInTr} = 1;
                }

                # check if hit on transcript region
                if (    $$rannodb[$k]{sta} < $$rannodb[$k]{sto}
                    and $var->{pos} >= $$rannodb[$k]{sta}
                    and $var->{pos} < $$rannodb[$k]{sto} )
                {
                    $trinfoEnt->{staInTr} = 1;
                }
                if (    $$rannodb[$k]{sta} < $$rannodb[$k]{sto}
                    and $var->{end} > $$rannodb[$k]{sta}
                    and $var->{end} <= $$rannodb[$k]{sto} )
                {
                    $trinfoEnt->{stoInTr} = 1;
                }

                # assign right rna positions
                if (    $total_left_ofst == 0
                    and $total_left_ofst == $total_right_ofst
                    and $rtidDetail->{wlen} > 0
                    and $rtidDetail->{mismatch} eq ""
                    and !exists $hit_badaln_ins{$tid} )
                {    # here is the part of insertion in edge,
                        # we don't need to recall cal_hgvs_pos,
                        # just use the already generated infomation
                        # About genepart judgement for this edge insertion case
                        # 1. promoter-any     => promoter
                        # 2. any-3'downstream => 3'downstream
                        # 3. utr-cds          => utr
                        # 4. exon-intron      => exon
                    if ($strd) {
                        $trinfoEnt->{rnaEnd} = $rpreLeft{$tid}->{nDot};
                        $trinfoEnt->{cdsEnd} = $rpreLeft{$tid}->{cDot};
                        $trinfoEnt->{ei_End} = $rpreLeft{$tid}->{exin};
                        $trinfoEnt->{r_End}  = $rpreLeft{$tid}->{r};
                    }
                    else {
                        $trinfoEnt->{rnaBegin} = $rpreLeft{$tid}->{nDot};
                        $trinfoEnt->{cdsBegin} = $rpreLeft{$tid}->{cDot};
                        $trinfoEnt->{ei_Begin} = $rpreLeft{$tid}->{exin};
                        $trinfoEnt->{r_Begin}  = $rpreLeft{$tid}->{r};
                    }

                    # correct genepartSO for different type of edge insertion
                    if (   $trinfoEnt->{r_Begin} eq 'PROM'
                        or $trinfoEnt->{r_End} eq 'PROM' )
                    {
                        $trinfoEnt->{genepartSO} = '167';
                        $trinfoEnt->{trRefComp}->{P0} = [ 0, 0 ];
                        $trinfoEnt->{exin}            = '.';
                        $trinfoEnt->{r}               = 'PROM';
                    }
                    elsif ($trinfoEnt->{r_Begin} eq '3D'
                        or $trinfoEnt->{r_End} eq '3D' )
                    { # 3'downstream insertion will be ignored as intergenic variantion
                        delete $annoEnt->{trInfo}->{$tid};
                        next;
                    }
                    elsif (
                        (
                                $trinfoEnt->{ei_End} =~ /EX/
                            and $trinfoEnt->{ei_Begin} =~ /IVS/
                        )
                        or (    $trinfoEnt->{r_End} =~ /^[53]U/
                            and $trinfoEnt->{r_Begin} =~ /^C/ )
                      )
                    {
                        $trinfoEnt->{genepartSO} =
                          getSOfromR( $trinfoEnt->{r_End} );
                        $trinfoEnt->{trRefComp}->{ $trinfoEnt->{ei_End} } = 0
                          if ( !exists $trinfoEnt->{trRefComp}
                            or !
                            exists $trinfoEnt->{trRefComp}
                            ->{ $trinfoEnt->{ei_End} } );
                        $trinfoEnt->{exin} = $trinfoEnt->{ei_End};
                        $trinfoEnt->{r}    = $trinfoEnt->{r_End};
                    }
                    elsif (
                        (
                                $trinfoEnt->{ei_Begin} =~ /EX/
                            and $trinfoEnt->{ei_End} =~ /IVS/
                        )
                        or (    $trinfoEnt->{r_Begin} =~ /^[53]U/
                            and $trinfoEnt->{r_End} =~ /^C/ )
                      )
                    {
                        $trinfoEnt->{genepartSO} =
                          getSOfromR( $trinfoEnt->{r_Begin} );
                        $trinfoEnt->{trRefComp}->{ $trinfoEnt->{ei_Begin} } = 0
                          if ( !exists $trinfoEnt->{trRefComp}
                            or !
                            exists $trinfoEnt->{trRefComp}
                            ->{ $trinfoEnt->{ei_Begin} } );
                        $trinfoEnt->{exin} = $trinfoEnt->{ei_Begin};
                        $trinfoEnt->{r}    = $trinfoEnt->{r_Begin};
                    }
                    else {    # no genepart change needed
                        if ( $trinfoEnt->{ei_Begin} eq $trinfoEnt->{ei_End} ) {
                            if ( $trinfoEnt->{ei_Begin} =~ /^EX/ ) {
                                $trinfoEnt->{trRefComp}
                                  ->{ $trinfoEnt->{ei_Begin} } = 0
                                  if ( !exists $trinfoEnt->{trRefComp}
                                    or !
                                    exists $trinfoEnt->{trRefComp}
                                    ->{ $trinfoEnt->{ei_Begin} } );

                            }
                            elsif ( $trinfoEnt->{ei_Begin} eq '?' ) {
                                $trinfoEnt->{trRefComp}->{'Q1'} = 0;
                            }
                            else {
                                $trinfoEnt->{trRefComp}
                                  ->{ $trinfoEnt->{ei_Begin} } = [ 0, 0 ];
                            }
                        }
                        $trinfoEnt->{exin} = $trinfoEnt->{ei_Begin};
                        $trinfoEnt->{r}    = $trinfoEnt->{r_Begin};
                    }
                }
                elsif (
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
                            or !
                            exists $tidExblk{$tid}{$tmp_exin}
                            { $rtidDetail->{nsta} . ',' . $rtidDetail->{nsto} }
                          )
                        {    # non Insertion on refseq, and first occured exon
                            $tidExblk{$tid}{$tmp_exin}
                              {     $rtidDetail->{nsta} . ','
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

                            delete $trinfoEnt->{trRefComp}->{Z999}
                              if (  $strd
                                and exists $trinfoEnt->{trRefComp}
                                and exists $trinfoEnt->{trRefComp}->{Z999} );

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
                    elsif ( $rtidDetail->{blka} eq '?' )
                    {    # problem for trRef assign for new version
                        my ( $af_sta, $af_sto );
                        $af_sta =
                          ( $$rannodb[$k]{sta} < $unify_p )
                          ? $unify_p
                          : $$rannodb[$k]{sta};
                        $af_sto =
                          ( $$rannodb[$k]{sto} < $unify_p + $unify_rl )
                          ? $$rannodb[$k]{sto}
                          : ( $unify_p + $unify_rl );
                        $trinfoEnt->{trRefComp}->{'Q1'} += $af_sto - $af_sta;
                    }
                    else {    # for intron region and PROM
                        if ( $rtidDetail->{blka} =~ /^PROM/ ) {
                            $tmp_exin = 'P0';    # for sort
                        }

                        delete $trinfoEnt->{trRefComp}->{Z999}
                          if (  $strd
                            and exists $trinfoEnt->{trRefComp}
                            and exists $trinfoEnt->{trRefComp}->{Z999} );

                        if (   !exists $trinfoEnt->{trRefComp}
                            or !exists $trinfoEnt->{trRefComp}->{$tmp_exin} )
                        {
                            $trinfoEnt->{trRefComp}->{$tmp_exin}->[0] =
                              $$rannodb[$k]{sta} - $unify_p;
                        }
                        my $less =
                            ( $$rannodb[$k]{sto} > ( $unify_p + $unify_rl ) )
                          ? ( $unify_p + $unify_rl )
                          : $$rannodb[$k]{sto};
                        $trinfoEnt->{trRefComp}->{$tmp_exin}->[1] =
                          $less - $unify_p;
                        if ( $tmp_exin eq 'P0' and !$strd ) {
                            $trinfoEnt->{trRefComp}->{$tmp_exin}->[1] =
                              $unify_rl;
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
                if ( ref($rpostRight) ) {
                    if ($strd) {
                        $trinfoEnt->{postEnd} = $rpostRight;
                    }
                    else {
                        $trinfoEnt->{preStart} = $rpostRight;
                    }
                }

                # debug
                #                print STDERR "AnnoEnt Check: "
                #                  . Dumper( $annoEnt->{trInfo}->{$tid} )
                #                  if (  exists $annoEnt->{trInfo}
                #                    and exists $annoEnt->{trInfo}->{$tid} );
            }
        }
    }

    # debug
    #    print STDERR "Annotations before check and uniform:", Dumper($annoEnt);

    # final check and uniform trinfo
    if ( exists $annoEnt->{trInfo} ) {
        foreach my $t ( keys %{ $annoEnt->{trInfo} } ) {
            my $tinfo = $annoEnt->{trInfo}->{$t};
            if ( exists $tinfo->{staInTr} and exists $tinfo->{stoInTr} ) {
                delete $tinfo->{staInTr};
                delete $tinfo->{stoInTr};
            }
            elsif ( !exists $tinfo->{staInTr}
                and !exists $tinfo->{stoInTr}
                and exists $tinfo->{preStart}
                and $tinfo->{preStart}->{r} ne 'PROM' )
            {    # hit 3'downstream only
                delete $annoEnt->{trInfo}->{$t};
                next;
            }
            else {
                delete $tinfo->{staInTr} if ( exists $tinfo->{staInTr} );
                delete $tinfo->{stoInTr} if ( exists $tinfo->{stoInTr} );
            }

            if (
                (
                    $tinfo->{strd} =~ /\+/
                    and
                    ( !exists $tinfo->{rnaEnd} or !defined $tinfo->{rnaEnd} )
                )
                or (
                    $tinfo->{strd} =~ /\-/
                    and (  !exists $tinfo->{rnaBegin}
                        or !defined $tinfo->{rnaBegin} )
                )
              )
            {
                delete $annoEnt->{trInfo}->{$t};
            }
            elsif ( exists $tinfo->{cdsEnd} and $tinfo->{cdsEnd} eq '' )
            {    # Downstream annotation fail
                $tinfo->{cdsBegin} = '';
            }
        }
        if ( 0 == scalar keys %{ $annoEnt->{trInfo} } ) {
            delete $annoEnt->{trInfo};
        }
    }

    #    debug
    #    print STDERR "DEBUG: Done getTrPosition\n";

    $new_aeIndex ||= $aeIndex;
    return $new_aeIndex;
}

sub getSOfromR {
    my $r = shift;
    if ( $r eq '?' ) {
        return 'annotation-fail';
    }
    elsif ( $r eq 'PROM' ) {
        return '167';
    }
    elsif ( $r =~ /^5U/ ) {
        return '204';
    }
    elsif ( $r =~ /^3U/ ) {
        return '205';
    }
    elsif ( $r =~ /^C/ ) {
        return '316';
    }
    elsif ( $r =~ /^R/ ) {
        return '655';
    }
    elsif ( $r =~ /^D/ ) {
        return '163';
    }
    elsif ( $r =~ /^A/ ) {
        return '164';
    }
    elsif ( $r =~ /^I5U/ ) {
        return '447';
    }
    elsif ( $r =~ /^I3U/ ) {
        return '448';
    }
    elsif ( $r =~ /^IC/ or $r =~ /^IR/ ) {
        return '191';
    }
    else {    # default to intergenic_region
        return '605';
    }
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

=cut

sub cal_hgvs_pos {
    my $annoEnt  = shift;
    my %cal_args = @_;
    if ( exists $cal_args{noassign} and $cal_args{noassign} == 0 ) {
        delete $cal_args{noassign};
    }
    my ( $ofst, $tid, $rtidDetail, $lr ) =
      @cal_args{qw(offset tid tidDetail LR)};
    my ( $nDot, $cDot ) = ( '', '' );
    my $strd = ( $$rtidDetail{strd} eq '+' ) ? 1 : 0;
    my $trAlt = $annoEnt->{trInfo}->{$tid}->{trAlt}
      if ( !exists $cal_args{noassign} );

    my $exin = $rtidDetail->{exin};
    my $blka = $rtidDetail->{blka};
    my $gpSO = $rtidDetail->{gpSO};

    my $lofst = $ofst;
    if ($lr) {    # only one time assignment for one var, offset for left
        my $rofst = $$rtidDetail{wlen} - $ofst - 1;
        if ( $lofst < 0 ) {
            if ($strd) {

                # 5'promoter is always available
                $nDot = $$rtidDetail{nsta} + $lofst;
                if ( $$rtidDetail{csta} =~ /^(\S+)\-u(\d+)/ ) {
                    $cDot = $1 . '-u' . ( $2 - $lofst );
                }
            }
            else {    # outside transcript region into 3'downstream
                $nDot = '+' . ( 0 - $lofst );
                if ( $$rtidDetail{csta} =~ /^\*?\d+$/ ) {
                    $cDot = $$rtidDetail{csta} . '+d' . ( 0 - $lofst );
                }

                $exin = '.';
                $blka = '3D';
                $gpSO = '605';
            }
        }
        elsif ( $lofst > $$rtidDetail{wlen} ) {
            if ( !exists $cal_args{noassign} ) {
                confess
                  "Error: exceed the whole length [$lofst>$$rtidDetail{wlen}]";
            }
            else {
                return 0;
            }
        }
        else {    # 0 <= $lofst <= $$rtidDetail{wlen}

            # 1. check if a mismatch block (only exon have this)
            if (    $rtidDetail->{mismatch} ne ""
                and $rtidDetail->{mismatch} ne "?" )
            {
                # for each kind of mismatch block,
                # hit the left edge of var, the refSeq ref will
                # contain start from the start of the mismatch
                # block, case I in database already has
                # a reversed coordinate between left and right.
                # transcript-alt string shoule be assigned to
                # this tid
                $nDot = $rtidDetail->{nsta};
                $cDot = $rtidDetail->{csta};
                if ( !exists $cal_args{noassign} and $lofst > 0 ) {

                    # mismatch info record the transcript-stranded
                    # reference sequence
                    my ( $mType, $mStart, $mEnd, $strand_ref ) =
                      split( /,/, $rtidDetail->{mismatch} );

                    if ( $mType eq 'E' ) {    # edge cross alignment
                        $trAlt = '?';
                    }

                    if ($strd) {              # cut left
                        $trAlt = substr( $strand_ref, 0, $lofst ) . $trAlt
                          if ( $trAlt ne '?' );
                    }
                    else {                    # cut right
                        $trAlt .= substr( $strand_ref, -$lofst, $lofst )
                          if ( $trAlt ne '?' );
                    }
                }
            }

            # 2. check if hit a normal block's right edge
            elsif ( $lofst == $$rtidDetail{wlen} ) {

                # give back the chance for left edge parsing
                if ( !exists $cal_args{noassign} ) {
                    delete $annoEnt->{trInfo}->{$tid};
                }
                return 0;
            }

            # 3. check if annotation-fail
            elsif ( $gpSO eq 'annotation-fail' ) {    # may not exist from 0.73

                #		( $nDot, $cDot ) = ( '?', '?' );
            }

            # 4. check if a promoter
            elsif ( $blka eq 'PROM' ) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $lofst )
                  : ( $rtidDetail->{nsta} - $lofst );
                if ( $rtidDetail->{csta} =~ /^(\S+)\-u(\d+)$/ ) {
                    $cDot =
                        ($strd)
                      ? ( $1 . '-u' . ( $2 - $lofst ) )
                      : ( $1 . '-u' . ( $2 + $lofst ) );
                }
            }

            # 5. check if an exon
            elsif ( $exin =~ /EX/ ) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $lofst )
                  : ( $rtidDetail->{nsta} - $lofst );
                if ( $rtidDetail->{csta} =~ /^(\*?)(\-?\d+)$/ ) {
                    $cDot =
                      ($strd) ? $1 . ( $2 + $lofst ) : $1 . ( $2 - $lofst );
                }
            }

            # 6. intron
            elsif ( $exin =~ /IVS/ ) {
                my $half_length = $rtidDetail->{wlen} / 2;
                if ( $gpSO eq 'abnormal-intron' ) {

                    # for abnormal intron
                    # the length may be less than 2 bp
                    # and then the nsta nsto and csta csto
                    # will point to the neighbor exons' edge
                    if ( $lofst >= $half_length ) {
                        if ( $rtidDetail->{nsto} =~ /\d+/ ) {
                            $nDot =
                              ($strd)
                              ? $rtidDetail->{nsto} . '-' . ( $rofst + 1 )
                              : $rtidDetail->{nsto} . '+' . ( $rofst + 1 );
                        }
                        if ( $rtidDetail->{csto} =~ /\d+/ ) {
                            $cDot =
                              ($strd)
                              ? $rtidDetail->{csto} . '-' . ( $rofst + 1 )
                              : $rtidDetail->{csto} . '+' . ( $rofst + 1 );
                        }
                    }
                    else {
                        if ( $rtidDetail->{nsta} =~ /\d+/ ) {
                            $nDot =
                              ($strd)
                              ? $rtidDetail->{nsta} . '+' . ( $lofst + 1 )
                              : $rtidDetail->{nsta} . '-' . ( $lofst + 1 );
                        }
                        if ( $rtidDetail->{csta} =~ /\d+/ ) {
                            $cDot =
                              ($strd)
                              ? $rtidDetail->{csta} . '+' . ( $lofst + 1 )
                              : $rtidDetail->{csta} . '-' . ( $lofst + 1 );
                        }
                    }
                }
                else {
                    if ( $lofst >= $half_length ) {    # drop into latter part
                        if ( $rtidDetail->{nsto} =~ /^(\d+\+?)(\-?\d+)$/ ) {
                            $nDot =
                                ($strd)
                              ? ( $1 . ( $2 - $rofst ) )
                              : ( $1 . ( $2 + $rofst ) );
                        }
                        if (
                            $rtidDetail->{csto} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/ )
                        {
                            $cDot =
                                ($strd)
                              ? ( $1 . ( $2 - $rofst ) )
                              : ( $1 . ( $2 + $rofst ) );
                        }
                    }
                    else {
                        if ( $rtidDetail->{nsta} =~ /^(\d+\+?)(\-?\d+)$/ ) {
                            $nDot =
                                ($strd)
                              ? ( $1 . ( $2 + $lofst ) )
                              : ( $1 . ( $2 - $lofst ) );
                        }
                        if (
                            $rtidDetail->{csta} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/ )
                        {
                            $cDot =
                                ($strd)
                              ? ( $1 . ( $2 + $lofst ) )
                              : ( $1 . ( $2 - $lofst ) );
                        }
                    }
                }
            }
            else {
                confess "Error: what's this? $rtidDetail->{exin}\n";
            }
        }

        if ( !exists $cal_args{noassign} ) {

           # the Begin End assignment is always from 5'->3' except for insertion
            if ($strd) {
                $annoEnt->{trInfo}->{$tid}->{rnaBegin} = $nDot;
                $annoEnt->{trInfo}->{$tid}->{cdsBegin} = $cDot;
                $annoEnt->{trInfo}->{$tid}->{ei_Begin} = $exin;
                $annoEnt->{trInfo}->{$tid}->{r_Begin}  = $blka;
            }
            else {
                $annoEnt->{trInfo}->{$tid}->{rnaEnd} = $nDot;
                $annoEnt->{trInfo}->{$tid}->{cdsEnd} = $cDot;
                $annoEnt->{trInfo}->{$tid}->{ei_End} = $exin;
                $annoEnt->{trInfo}->{$tid}->{r_End}  = $blka;
            }
            $annoEnt->{trInfo}->{$tid}->{genepartSO} = $gpSO;
            $annoEnt->{trInfo}->{$tid}->{trAlt}      = $trAlt;
        }
    }
    else {    # can be multiple times assignment if large var, offset for right
        my $real_ofst = $lofst - 1;                   # end_ofst change
        my $rofst     = $$rtidDetail{wlen} - $ofst;
        if ( $lofst > $$rtidDetail{wlen} ) {
            my $ex_ofst = ( $lofst - $$rtidDetail{wlen} );
            if ($strd) {
                $nDot = '+' . $ex_ofst;
                if ( $$rtidDetail{csto} =~ /^\*?\d+$/ ) {
                    $cDot = $$rtidDetail{csto} . '+d' . $ex_ofst;
                }
                $exin = '.';
                $blka = '3D';
                $gpSO = '605';
            }
            else {    # outside 5'promoter region ?
                $nDot =
                    ( $$rtidDetail{nsto} =~ /^\-\d+$/ )
                  ? ( $$rtidDetail{nsto} - $ex_ofst )
                  : ( -$ex_ofst );
                if ( $$rtidDetail{csto} =~ /^(\S+)\-u(\d+)/ ) {
                    $cDot = $1 . '-u' . ( $2 + $ex_ofst );
                }
            }
        }
        elsif ( $lofst < 0 ) {    # not proper called
            if ( !exists $cal_args{noassign} ) {
                confess "Error: right preceed the block [$lofst < 0]";
            }
            else {
                return 0;
            }
        }
        else {                    # 0 <= $lofst <= $$rtidDetail{wlen}

            # 1. check if a mismatch block (exon)
            if (    $rtidDetail->{mismatch} ne ""
                and $rtidDetail->{mismatch} ne "?" )
            {
                $nDot = $rtidDetail->{nsto};
                $cDot = $rtidDetail->{csto};
                if ( !exists $cal_args{noassign}
                    and $lofst < $rtidDetail->{wlen} )
                {
                    my ( $mType, $mStart, $mEnd, $strand_ref ) =
                      split( /,/, $rtidDetail->{mismatch} );

                    if ( $mType eq 'E' ) {    # edge cross alignment
                        $trAlt = '?';
                    }

                    my $cut_len = $rtidDetail->{wlen} - $lofst;
                    if ($strd) {              # cut right
                        $trAlt .= substr( $strand_ref, -$cut_len, $cut_len )
                          if ( $trAlt ne '?' );
                    }
                    else {                    # cut left
                        $trAlt = substr( $strand_ref, 0, $cut_len ) . $trAlt
                          if ( $trAlt ne '?' );
                    }
                }
            }

            # 2. check if hit a normal block's left edge
            elsif ( $lofst == 0 ) {
                return 0;
            }

            #	    # 3. check if annotation-fail
            elsif ( $gpSO eq 'annotation-fail' ) {

                #		( $nDot, $cDot ) = ( '?', '?' );
            }

            # 4. check if a promoter
            elsif ( $blka =~ /^PROM/ ) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $real_ofst )
                  : ( $rtidDetail->{nsta} - $real_ofst );
                if ( $rtidDetail->{csta} =~ /^(\S+)\-u(\d+)$/ ) {
                    $cDot =
                        ($strd)
                      ? ( $1 . '-u' . ( $2 - $real_ofst ) )
                      : ( $1 . '-u' . ( $2 + $real_ofst ) );
                }
            }

            # 5. check if an exon
            elsif ( $exin =~ /EX/ ) {
                $nDot =
                    ($strd)
                  ? ( $rtidDetail->{nsta} + $real_ofst )
                  : ( $rtidDetail->{nsta} - $real_ofst );
                if ( $rtidDetail->{csta} =~ /^(\*?)(\-?\d+)$/ ) {
                    $cDot =
                      ($strd)
                      ? $1 . ( $2 + $real_ofst )
                      : $1 . ( $2 - $real_ofst );
                }
            }

            # 6. intron
            elsif ( $exin =~ /IVS/ ) {
                my $half_length = $rtidDetail->{wlen} / 2;

                if ( $gpSO eq 'abnormal-intron' ) {
                    if ( $real_ofst >= $half_length ) {
                        if ( $rtidDetail->{nsto} =~ /\d+/ ) {
                            $nDot =
                              ($strd)
                              ? $rtidDetail->{nsto} . '-' . ( $rofst + 1 )
                              : $rtidDetail->{nsto} . '+' . ( $rofst + 1 );
                        }
                        if ( $rtidDetail->{csto} =~ /\d+/ ) {
                            $cDot =
                              ($strd)
                              ? $rtidDetail->{csto} . '-' . ( $rofst + 1 )
                              : $rtidDetail->{csto} . '+' . ( $rofst + 1 );
                        }
                    }
                    else {
                        if ( $rtidDetail->{nsta} =~ /\d+/ ) {
                            $nDot =
                              ($strd)
                              ? $rtidDetail->{nsta} . '+' . ( $real_ofst + 1 )
                              : $rtidDetail->{nsta} . '-' . ( $real_ofst + 1 );
                        }
                        if ( $rtidDetail->{csta} =~ /\d+/ ) {
                            $cDot =
                              ($strd)
                              ? $rtidDetail->{csta} . '+' . ( $real_ofst + 1 )
                              : $rtidDetail->{csta} . '-' . ( $real_ofst + 1 );
                        }
                    }
                }
                else {
                    if ( $real_ofst >= $half_length ) {  # drop into latter part
                        if ( $rtidDetail->{nsto} =~ /^(\d+\+?)(\-?\d+)$/ ) {
                            $nDot =
                                ($strd)
                              ? ( $1 . ( $2 - $rofst ) )
                              : ( $1 . ( $2 + $rofst ) );
                        }
                        if (
                            $rtidDetail->{csto} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/ )
                        {
                            $cDot =
                                ($strd)
                              ? ( $1 . ( $2 - $rofst ) )
                              : ( $1 . ( $2 + $rofst ) );
                        }
                    }
                    else {
                        if ( $rtidDetail->{nsta} =~ /^(\d+\+?)(\-?\d+)$/ ) {
                            $nDot =
                                ($strd)
                              ? ( $1 . ( $2 + $real_ofst ) )
                              : ( $1 . ( $2 - $real_ofst ) );
                        }
                        if (
                            $rtidDetail->{csta} =~ /^([\*\-]?\d+\+?)(\-?\d+)$/ )
                        {
                            $cDot =
                                ($strd)
                              ? ( $1 . ( $2 + $real_ofst ) )
                              : ( $1 . ( $2 - $real_ofst ) );
                        }
                    }
                }
            }
            else {
                confess "Error: what's this? $rtidDetail->{exin}\n";
            }
        }

        if ( !exists $cal_args{noassign} ) {

           # the Begin End assignment is always from 5'->3' except for insertion
            if ($strd) {
                $annoEnt->{trInfo}->{$tid}->{rnaEnd} = $nDot;
                $annoEnt->{trInfo}->{$tid}->{cdsEnd} = $cDot;
                $annoEnt->{trInfo}->{$tid}->{ei_End} = $exin;
                $annoEnt->{trInfo}->{$tid}->{r_End}  = $blka;
            }
            else {
                $annoEnt->{trInfo}->{$tid}->{rnaBegin} = $nDot;
                $annoEnt->{trInfo}->{$tid}->{cdsBegin} = $cDot;
                $annoEnt->{trInfo}->{$tid}->{ei_Begin} = $exin;
                $annoEnt->{trInfo}->{$tid}->{r_Begin}  = $blka;
            }
            $annoEnt->{trInfo}->{$tid}->{trAlt} = $trAlt;

            if (
                $rtidDetail->{gpSO} eq 'annotation-fail'
                or ( exists $annoEnt->{trInfo}->{$tid}->{genepartSO}
                    and $annoEnt->{trInfo}->{$tid}->{genepartSO} eq
                    'annotation-fail' )
              )
            {
                $annoEnt->{trInfo}->{$tid}->{genepartSO} = 'annotation-fail';
            }
            elsif ( exists $annoEnt->{trInfo}->{$tid}->{genepartSO}
                and $annoEnt->{trInfo}->{$tid}->{genepartSO} ne "span"
                and $rtidDetail->{gpSO} ne
                $annoEnt->{trInfo}->{$tid}->{genepartSO} )
            {
                $annoEnt->{trInfo}->{$tid}->{genepartSO} = "span";
            }
            elsif ( !exists $annoEnt->{trInfo}->{$tid}->{genepartSO} ) {
                carp "Warning: no genepartSO specified in left?";
                $annoEnt->{trInfo}->{$tid}->{genepartSO} = $gpSO;
            }
        }
    }

    if ( exists $cal_args{noassign} ) {
        return {
            nDot => $nDot,
            cDot => $cDot,
            exin => $exin,
            r    => $blka,
        };
    }
    return 1;
}

1;
__END__
