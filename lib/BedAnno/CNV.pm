=head1 BedAnno::CNV

    BedAnno::CNV sub package

=cut

package BedAnno::CNV;
use base qw(BedAnno);
use strict;
use warnings;

use Data::Dumper;
use Carp;
use Tabix;

=head1 METHOD

=head2 new

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
                - db, tr, genes, trans, region, regbed, mmap, batch, cytoBand
    Returns : BedAnno::CNV object

=cut

sub new {
    my ( $class, %args ) = @_;
    my $self;
    $self = $class->SUPER::new(%args);
    bless $self, ref($class) || $class;

    if ( exists $args{dgv} ) {
        $self->set_dgv( $args{dgv} );
    }

    if ( exists $args{sfari} ) {
        $self->set_sfari( $args{sfari} );
    }

    if ( exists $args{cnvPub} ) {
        $self->set_cnvPub( $args{cnvPub} );
    }

    if ( exists $args{cnvd} ) {
        $self->set_cnvd( $args{cnvd} );
    }

    return $self;
}

sub set_dgv {
    my $self  = shift;
    my $dgvdb = shift;
    $self->{dgv} = $dgvdb;
    require GetDGV if ( !exists $self->{dgv_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    $common_opts{ovlp_rate} = $self->{ovlp_rate}
      if ( exists $self->{ovlp_rate} );
    $common_opts{max_uncov} = $self->{max_uncov}
      if ( exists $self->{max_uncov} );
    my $dgv_h = GetDGV->new( db => $dgvdb, %common_opts );
    $self->{dgv_h} = $dgv_h;
    return $self;
}

sub set_sfari {
    my $self    = shift;
    my $sfaridb = shift;
    $self->{sfari} = $sfaridb;
    require GetSFARI if ( !exists $self->{sfari_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    $common_opts{ovlp_rate} = $self->{ovlp_rate}
      if ( exists $self->{ovlp_rate} );
    $common_opts{max_uncov} = $self->{max_uncov}
      if ( exists $self->{max_uncov} );
    my $sfari_h = GetSFARI->new( db => $sfaridb, %common_opts );
    $self->{sfari_h} = $sfari_h;
    return $self;
}

sub set_cnvPub {
    my $self     = shift;
    my $cnvPubdb = shift;
    $self->{cnvPub} = $cnvPubdb;
    require GetCNVPub if ( !exists $self->{cnvPub_h} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    $common_opts{ovlp_rate} = $self->{ovlp_rate}
      if ( exists $self->{ovlp_rate} );
    $common_opts{max_uncov} = $self->{max_uncov}
      if ( exists $self->{max_uncov} );
    my $cnvPub_h = GetCNVPub->new( db => $cnvPubdb, %common_opts );
    $self->{cnvPub_h} = $cnvPub_h;
    return $self;
}

sub set_cnvd {
    my $self  = shift;
    my $cnvdb = shift;
    $self->{cnvd} = $cnvdb;
    require GetCNVD if ( !exists $self->{cnvd} );
    my %common_opts = ();
    $common_opts{quiet} = 1 if ( exists $self->{quiet} );
    $common_opts{ovlp_rate} = $self->{ovlp_rate}
      if ( exists $self->{ovlp_rate} );
    $common_opts{max_uncov} = $self->{max_uncov}
      if ( exists $self->{max_uncov} );
    my $cnvd_h = GetCNVD->new( db => $cnvdb, %common_opts );
    $self->{cnvd_h} = $cnvd_h;
    return $self;
}

=head2 annoCNV

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

=cut

sub annoCNV {
    my $self = shift;
    my ( $chr, $start, $end, $copy_number ) = @_;

    $chr =~ s/^chr//i;
    if ( $chr =~ /^M/i ) {
        $chr = 'MT';
    }

    my $cb = $self->{cytoBand_h}->getCB( $chr, $start, $end )
      if ( defined $self->{cytoBand_h} );
    my $dgv = $self->{dgv_h}->getDGV( $chr, $start, $end, $copy_number )
      if ( defined $self->{dgv_h} );
    my $sfari = $self->{sfari_h}->getSFARI( $chr, $start, $end, $copy_number )
      if ( defined $self->{sfari_h} );
    my $cnvPub =
      $self->{cnvPub_h}->getCNVpub( $chr, $start, $end, $copy_number )
      if ( defined $self->{cnvPub_h} );
    my $cnvd = $self->{cnvd_h}->getCNVD( $chr, $start, $end, $copy_number )
      if ( defined $self->{cnvd_h} );

    my $rAnnos;
    my %open_args;
    my ( $qstart, $qend ) = ( $start, $end );
    if ( $start == $end ) {    #query with 1 bp flank for ins case
        $qstart-- if ( $qstart > 0 );
        $qend++;
    }
    $open_args{region} = $chr . ':' . ($qstart + 1) . '-' . $qend;
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
    my $rcurrent = $self->load_anno(%open_args);
    $rAnnos = $rcurrent->{$chr} if ( exists $rcurrent->{$chr} );

    my %cnv_anno = ();
    $cnv_anno{cytoBand} = $cb     if ( defined $cb );
    $cnv_anno{dgv}      = $dgv    if ( defined $dgv );
    $cnv_anno{cnvPub}   = $cnvPub if ( defined $cnvPub );
    $cnv_anno{sfari}    = $sfari  if ( defined $sfari );
    $cnv_anno{cnvd}     = $cnvd   if ( defined $cnvd );

    if ( !defined $rAnnos ) {
        $self->warn("Warning: no available annotation items in curdb for $chr")
          if ( !exists $self->{quiet} );
        $cnv_anno{cnva_type} = "unknown";
        return \%cnv_anno;
    }

    my $pseudo_var = BedAnno::Var->new( $chr, $start, $end, "=", "?" );
    my $pseudo_anno = BedAnno::Anno->new($pseudo_var);
    $pseudo_anno->getTrPosition( $rAnnos, 0 );
    my $rhitted_blks = $self->get_hitted_blk($pseudo_anno);

    my $ranno_cnv = $self->parse_cnva($rhitted_blks);
    $cnv_anno{anno} = $ranno_cnv->{anno} if ( exists $ranno_cnv->{anno} );
    $cnv_anno{cnva_type} = $ranno_cnv->{cnva_type};

    return \%cnv_anno;
}

=head2 batch_annoCNV

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

=cut

sub batch_annoCNV {
    my $self     = shift;
    my $rCNVvars = shift;
    my %cnvAnnos = ();
    foreach my $chr ( sort keys %$rCNVvars ) {
        my @PosPairs = map { [ split(/\-/) ] } keys %{ $rCNVvars->{$chr} };
        my $rh_cnva = $self->get_cover_batch( $chr, \@PosPairs );
        foreach my $pos_pair ( sort keys %$rh_cnva ) {
            my $cur_anno = $self->parse_cnva( $rh_cnva->{$pos_pair} );

            my $cur_CN = $rCNVvars->{$chr}->{$pos_pair};
            $cur_anno->{cur_CN} = $cur_CN;

            my ( $start, $stop ) = split( /-/, $pos_pair );
            $start -= 1;    # change 1 based pos pair to 0 based start
            $cur_anno->{dgv} =
              $self->{dgv_h}->getDGV( $chr, $start, $stop, $cur_CN )
              if ( defined $self->{dgv_h} );
            $cur_anno->{sfari} =
              $self->{sfari_h}->getSFARI( $chr, $start, $stop, $cur_CN )
              if ( defined $self->{sfari_h} );
            $cur_anno->{cnvPub} =
              $self->{cnvPub_h}->getCNVpub( $chr, $start, $stop, $cur_CN )
              if ( defined $self->{cnvPub_h} );
            $cur_anno->{cytoband} =
              $self->{cytoBand_h}->getCB( $chr, $start, $stop )
              if ( defined $self->{cytoBand_h} );
            $cur_anno->{cnvd} =
              $self->{cnvd_h}->getCNVD( $chr, $start, $stop, $cur_CN )
              if ( defined $self->{cnvd_h} );
            $cnvAnnos{$chr}{$pos_pair} = $cur_anno;
        }
    }
    return \%cnvAnnos;
}

sub parse_cnva {
    my $self         = shift;
    my $rhitted_blks = shift;
    my %cur_anno     = ();
    my $cnv_anno_type;
    foreach my $rTid (@$rhitted_blks) {
        my ( $tid, $rleft, $rright ) = @$rTid;
        my $ref_tida = {};
        $ref_tida = get_region( $rleft, $rright );

        if ( exists $ref_tida->{regtyp} ) {
            if ( !defined $cnv_anno_type ) {
                $cnv_anno_type = $ref_tida->{regtyp};
            }
            elsif ( $cnv_anno_type ne $ref_tida->{regtyp} ) {
                $cnv_anno_type = "Multiple";
            }
        }
        else {
            $self->throw("Error: fail to get region type for $tid\n");
        }

        @$ref_tida{qw(gsym gid strd)} = @$rleft{qw(gsym gid strd)};
        $cur_anno{anno}{$tid} = $ref_tida;
    }

    $cur_anno{cnva_type} =
      ( defined $cnv_anno_type ) ? $cnv_anno_type : "unknown";

    return \%cur_anno;
}

sub get_region {
    my ( $rl, $rr ) = @_;

    confess "Unknown invoking error, no cpos.\n"
      if ( !exists $rl->{cpos} or !exists $rr->{cpos} );

    my %region_anno = ();
    if ( $rl->{cpos} eq '?' or $rr->{cpos} eq '?' ) {
        $region_anno{cpos}   = '?';
        $region_anno{exin}   = 'unknown';
        $region_anno{regcod} = 'unknown';
        $region_anno{regtyp} = 'unknown';
        return \%region_anno;
    }

    my $cmp_stat =
      BedAnno->cmpPos( substr( $rl->{cpos}, 2 ), substr( $rr->{cpos}, 2 ) );

    if ( $cmp_stat == 0 ) {
        $region_anno{cpos} = $rl->{cpos};
    }
    elsif ( $cmp_stat > 0 ) {
        my $latter_part = substr( $rr->{cpos}, 2 );
        $region_anno{cpos} = $rl->{cpos} . '_' . $latter_part;
    }
    else {
        my $latter_part = substr( $rl->{cpos}, 2 );
        $region_anno{cpos} = $rr->{cpos} . '_' . $latter_part;
    }

    if ( $rl->{reg} eq $rr->{reg} ) {
        $region_anno{regcod} = $rl->{reg};
        $region_anno{exin} =
          ( $rl->{reg} eq 'PROM' ) ? "Promoter" : $rl->{exin};
        if ( $region_anno{regcod} =~ /^I|^5U|^3U|^PROM/i ) {
            $region_anno{regtyp} = "LowRisk";
        }
        else {
            $region_anno{regtyp} = "HighRisk";
        }
    }
    else {
        my ( $r5, $r3 ) = ( $rl, $rr );
        if ( $cmp_stat < 0 ) {
            $r5 = $rr;
            $r3 = $rl;
        }

        $region_anno{regcod} = $r5->{reg} . '-' . $r3->{reg};

        my $formor_exin = ( $r5->{reg} eq "PROM" ) ? "Promoter" : $r5->{exin};
        if ( $r3->{reg} eq '3D' ) {
            $region_anno{exin} = $formor_exin . "-DownStream";
        }
        elsif ( $r3->{exin} eq $r5->{exin} ) {    # UTR:CDS or SPLICE:INTRON
            $region_anno{exin} = $r5->{exin};
        }
        else {
            $region_anno{exin} = $formor_exin . '-' . $r3->{exin};
        }

        if ( $r5->{reg} eq 'PROM' and $r3->{reg} eq '3D' ) {
            $region_anno{regtyp} = 'Whole-gene';
        }
        elsif ( ( $r5->{reg} =~ /5U/ or $r5->{reg} eq 'PROM' )
            and ( $r3->{reg} =~ /3U/ or $r3->{reg} eq '3D' ) )
        {
            $region_anno{regtyp} = 'Whole-cds';
        }
        else {
            $region_anno{regtyp} = "HighRisk";
        }
    }

    return \%region_anno;
}

sub DESTROY {
    my $self = shift;
    foreach my $key ( keys %$self ) {
        next if ( $key !~ /_h$/ );
        $self->{$key}->DESTROY()
          if ( defined $self->{$key} and $self->{$key}->can('DESTROY') );
        delete $self->{$key};
    }
    return;
}

1;
__END__

