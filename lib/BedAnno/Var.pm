=head1 BedAnno::Var

    BedAnno::Var sub module

=cut

package BedAnno::Var;
use strict;
use warnings;
use Data::Dumper;
use Carp;

use BedAnno qw($AAcount %AAnumber);
our $MAX_COMPLEX_PARSING = 200;

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
    my ( $chr, $start, $end, $ref, $alt );
    my %var;

    if ( ref( $_[0] ) ) {
        if (   !exists $_[0]->{chr}
            or !exists $_[0]->{begin}
            or ( !exists $_[0]->{end} and !exists $_[0]->{referenceSequence} ) )
        {
            confess "Error: unavailable object. need keys: ",
              "chr, start, alt, ref specified.";
        }

        %var = %{ $_[0] };

        $chr   = $_[0]->{chr};
        $start = $_[0]->{begin};
        $end   = $_[0]->{end} if ( exists $_[0]->{end} );
        if ( exists $_[0]->{variantSequence}
            and defined $_[0]->{variantSequence} )
        {
            $alt = $_[0]->{variantSequence};
        }
        else {
            $alt = "";
        }
        $alt = "" if ( $alt =~ /^null$/i );

        if ( exists $_[0]->{referenceSequence}
            and !defined $_[0]->{referenceSequence} )
        {
            $_[0]->{referenceSequence} = "";
        }
        $ref =
          ( exists $_[0]->{end} and $start eq $end ) ? ""
          : (
            ( exists $_[0]->{referenceSequence} ) ? $_[0]->{referenceSequence}
            : "="
          );
        $ref = "" if ( $ref =~ /^null$/i );

        # clean %var
        delete $var{"chr"};
        delete $var{"begin"};
        delete $var{"end"}               if ( exists $var{end} );
        delete $var{"referenceSequence"} if ( exists $var{referenceSequence} );
        delete $var{"variantSequence"}   if ( exists $var{variantSequence} );
    }
    else {
        confess "Error: not enough args, need at least 4 args." if ( 4 > @_ );
        ( $chr, $start, $end, $ref, $alt ) = @_;
    }

    if ( !defined $end or $end !~ /^\d+$/ ) {    # from VCF v4.1 1-based start
        if ( defined $end ) {
            $alt = $ref;
            $ref = $end;
        }

        $ref = normalise_seq($ref);
        $alt = normalise_seq($alt);
        if ( ( $ref ne $alt or length($ref) > 1 )
            and substr( $ref, 0, 1 ) eq substr( $alt, 0, 1 ) )
        {
            $ref = substr( $ref, 1 );
            $alt = substr( $alt, 1 );
        }
        else {
            $start -= 1;    # change to 0-based start
        }
        my $rl = length($ref);
        $end = $start + $rl;
    }

    my $len_ref = $end - $start;    # chance to annotate long range
    my ( $varType, $implicit_varType, $sm ) =
      guess_type( $len_ref, $ref, $alt );

    $chr = normalise_chr($chr);
    %var = (
        %var,    # remain some extra keys in the var hash, e.g. var_id etc.
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
    $var{altlen} = length($alt) if ( $alt ne '?' );

    my $var;
    $var = {%var};
    bless $var, ref($class) || $class;

    if ( $varType eq 'no-call' or $implicit_varType ne 'delins' ) {
        return $var;
    }

    return $var
      if ( $len_ref > $MAX_COMPLEX_PARSING
        or length($alt) > $MAX_COMPLEX_PARSING );
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
    my $var  = shift;
    my $strd = shift;
    my $norep = shift || 0;

    my $consPos = $$var{pos};
    my $consRef = $$var{ref};
    my $consAlt = $$var{alt};
    my $consRL  = $$var{reflen};

    if ( !exists $$var{altlen} ) {
        return ( $consPos, $consRef, $consAlt, $consRL );
    }

    my $consAL = $$var{altlen};

    if (!$norep and exists $var->{p} ) {    # rep
        $consPos = $$var{p};
        $consRef = $$var{r};
        $consAlt = $$var{a};
        $consRL  = $$var{rl};
        $consAL  = $$var{al};
    }
    elsif ( exists $var->{bp} ) {    # complex bc strand same
        $consPos = $$var{bp};
        $consRef = $$var{br};
        $consAlt = $$var{ba};
        $consRL  = $$var{brl};
        $consAL  = $$var{bal};
    }
    elsif ( exists $var->{$strd} ) {    # complex bc strand diff
        $consPos = $$var{$strd}{bp};
        $consRef = $$var{$strd}{br};
        $consAlt = $$var{$strd}{ba};
        $consRL  = $$var{$strd}{brl};
        $consAL  = $$var{$strd}{bal};
    }

    return ( $consPos, $consRef, $consAlt, $consRL, $consAL );
}

sub normalise_chr {
    my $chr = shift;
    $chr =~ s/^chr//i;
    if ( $chr eq '23' ) {
        $chr = 'X';
    }
    elsif ( $chr eq '24' ) {
        $chr = 'Y';
    }
    elsif ( $chr eq '25' ) {
        $chr = 'MT';
    }
    elsif ( $chr =~ /^[xym]$/i ) {
        $chr = uc($chr);
        $chr .= 'T' if ( $chr eq 'M' );
    }
    elsif ( $chr =~ /^M_/i ) {
        $chr = 'MT';
    }
    return $chr;
}

sub normalise_seq {
    my $seq = shift;
    $seq =~ s/\s+//g;
    $seq = "" if ( $seq eq '.' );
    $seq = uc($seq);
    if ( $seq =~ /[^ACGTN]/ and $seq ne '?' ) {
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
    my ( $ref, $alt, $len_ref, $len_alt ) = @{$var}{qw(ref alt reflen altlen)};

    my $get_rst = get_internal( $ref, $len_ref, $alt, $len_alt );
    if ( $get_rst->{r} != $len_ref ) {
        my ( $imp_guess, $sm ) =
          guess_type_by_length( $get_rst->{r}, $get_rst->{a} );
        $var->{imp} = $imp_guess;
        $var->{sm}  = $sm;
        if ( $get_rst->{'+'} != $get_rst->{'-'} ) {
            for my $strd ( '+', '-' ) {
                $var->{$strd}->{bp}  = $var->{pos} + $get_rst->{$strd};
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
            $var->{br}  = substr( $var->{ref}, $get_rst->{'+'}, $get_rst->{r} );
            $var->{ba}  = substr( $var->{alt}, $get_rst->{'+'}, $get_rst->{a} );
        }
    }

    if ( $var->{sm} == 3 ) {    # equal length long subs
        my @separate_snvs = ();
        for ( my $i = 0 ; $i < $var->{reflen} ; $i++ ) {
            next
              if (
                substr( $var->{ref}, $i, 1 ) eq substr( $var->{alt}, $i, 1 ) );
            push( @separate_snvs, $var->{pos} + $i + 1 );    # 1-based
        }
        $var->{sep_snvs} = [@separate_snvs];
    }

    my $rc_ref = count_content($ref);
    my $rc_alt = count_content($alt);
    my @diff   = map { $$rc_ref[$_] - $$rc_alt[$_] } ( 0 .. $AAcount );

    # check if the sign of all base diff are consistent.
    if ( check_sign( \@diff ) ) {    # possible short tandom repeat variation
        my @absdiff = map { abs } @diff;
        my ( $larger, $smaller, $llen, $slen );
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
        for ( my $rep = $llen ; $rep > 0 ; $rep-- ) {
            while ( $larger =~ /([A-Z]+)(?:\1){$rep}/g ) {
                next if ( exists $has{$1} );

                my $rep_el = $1;
                my $lofs   = length($`);    # $` is the prematched string

                my $cn = check_div( $rep_el, \@absdiff );
                my $lenrep = length($rep_el);
                if ( $cn and check_insrep( $larger, $smaller, $rep_el, $cn ) ) {

                    @$var{qw(p rep replen)} =
                      ( ( $var->{pos} + $lofs ), $rep_el, $lenrep );

                    $rep += 1;              # add the first copy of element

                    my $l_cn       = $rep;
                    my $s_cn       = $rep - $cn;
                    my $l_wholerep = $rep_el x $l_cn;
                    my $s_wholerep = $rep_el x $s_cn;
                    my $l_replen   = $lenrep * $l_cn;
                    my $s_replen   = $lenrep * $s_cn;

                    if ( $llen == $len_ref ) {    # ref is longer
                        @$var{qw(ref_cn alt_cn r a rl al)} = (
                            $l_cn, $s_cn, $l_wholerep, $s_wholerep, $l_replen,
                            $s_replen
                        );
                    }
                    else {                        # alt is longer
                        @$var{qw(alt_cn ref_cn a r al rl)} = (
                            $l_cn, $s_cn, $l_wholerep, $s_wholerep, $l_replen,
                            $s_replen
                        );
                    }

                    $var->{imp} = 'rep';
                    $var->{sm} = ( $var->{rl} == 1 ) ? 1 : 2;

                    return $var;
                }
                else {
                    pos($larger) -= $lenrep * ($rep + 1) - 1;
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
    my ( $rl, $al ) = @_;
    my ( $imp_guess, $sm );
    if ( $rl == 0 ) {
        $imp_guess = 'ins';
        $sm        = 0;
    }
    elsif ( $rl == 1 ) {
        if ( $al == 1 ) {
            $imp_guess = 'snv';
        }
        elsif ( $al == 0 ) {
            $imp_guess = 'del';
        }
        else {
            $imp_guess = 'delins';
        }
        $sm = 1;
    }
    else {
        $imp_guess = ( $al == 0 )   ? 'del' : 'delins';
        $sm        = ( $rl == $al ) ? 3     : 2;
    }

    return ( $imp_guess, $sm );
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
    my ( $len_ref, $ref, $alt ) = @_;

    # imp_varType: implicit variant type is for HGVS naming
    # varType    : to be output as varType name. can be as the key
    #              to get the SO term by Name2SO hash.
    # sm         : single or multiple bases tag
    #              0 - for insertion, 0 base
    #              1 - for single base variants
    #              2 - for non-equal-length multiple bases variants
    #              3 - for equal-length multiple bases delins
    my ( $imp_varType, $varType, $sm );
    if ( $len_ref == 0 ) {
        if ( $ref eq $alt ) {
            $imp_varType = 'ref';
        }
        else {
            $imp_varType = 'ins';
        }
        $sm = 0;
    }
    elsif ( $len_ref == 1 ) {
        if ( $ref eq $alt ) {
            $imp_varType = 'ref';
        }
        elsif ( $alt eq '' ) {
            $imp_varType = 'del';
        }
        elsif ( 1 == length($alt) ) {
            $imp_varType = 'snv';
        }
        else {
            $imp_varType = 'delins';
        }
        $sm = 1;
    }
    elsif ( $len_ref > 1 ) {
        if ( $ref eq $alt ) {
            $imp_varType = 'ref';
        }
        elsif ( $alt eq '' ) {
            $imp_varType = 'del';
        }
        else {
            $imp_varType = 'delins';
        }

        if ( length($ref) != length($alt) ) {
            $sm = 2;    # non-equal-length delins
        }
        else {
            $sm = 3;    # equal-length subs
        }
    }
    $varType = ( $alt eq '?' or $alt eq 'N' ) ? 'no-call' : $imp_varType;
    return ( $varType, $imp_varType, $sm );
}

# check whether the smaller with inserted $cn copies repeat elements
# is the same with larger one
sub check_insrep {
    my ( $larger, $smaller, $repeat, $cn ) = @_;
    my $ind = index( $smaller, $repeat );
    if ( $ind > -1 ) {
        my $temp = $smaller;
        substr( $temp, $ind, 0, ( $repeat x $cn ) );
        if ( $larger eq $temp ) {
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
    my ( $ref, $reflen, $alt, $altlen ) = @_;
    my $shorter = ( $reflen < $altlen ) ? $reflen : $altlen;
    my ( $lgo, $loff, $rgo, $roff ) = ( 1, 0, 1, 0 );
    for ( my $i = 0 ; $i < $shorter ; $i++ ) {
        if ( $lgo and substr( $ref, $i, 1 ) eq substr( $alt, $i, 1 ) ) {
            $loff++;
        }
        else {
            $lgo = 0;
        }

        if ( $rgo
            and substr( $ref, -( $i + 1 ), 1 ) eq
            substr( $alt, -( $i + 1 ), 1 ) )
        {
            $roff++;
        }
        else {
            $rgo = 0;
        }

        last if ( $lgo == 0 and $rgo == 0 );
    }
    my ( $new_ref_len, $new_alt_len );
    if ( $shorter >= $loff + $roff ) {
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
            '-' => ( $shorter - $roff ),
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
    my ( $s, $rc ) = @_;
    my $rcs = count_content($s);
    return 0 unless ( $$rc[0] % $$rcs[0] == 0 );
    my $div = $$rc[0] / $$rcs[0];
    for ( 1 .. $AAcount ) {
        return 0 if ( $$rc[$_] != $$rcs[$_] * $div );
    }
    return $div;
}

sub count_content {
    my $s     = uc(shift);
    my $l     = length($s);
    my @count = ( $l, (0) x $AAcount );
    while ( $s =~ /(.)/g ) {
        confess "unknown base [$1]" if ( !exists $AAnumber{$1} );
        $count[ $AAnumber{$1} ]++;
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
    my $var   = shift;
    my $gHGVS = 'g.';
    if ( $var->{chr} =~ /^M/ ) {    # hit mito chromosome
        $gHGVS = 'm.';
    }

    my $imp = $var->{imp};
    my $sm  = $var->{sm};
    my ( $pos, $ref, $alt, $reflen, $altlen ) = $var->getUnifiedVar('+');

    if ( $imp eq 'snv' ) {
        $gHGVS .= $pos . $ref . '>' . $alt;
    }
    elsif ( $imp eq 'ins' ) {
        $gHGVS .= $pos . '_' . ( $pos + 1 ) . 'ins' . $alt;
    }
    elsif ( $imp eq 'del' or $imp eq 'delins' ) {

        # 1bp del/delins
        $gHGVS .= ( $pos + 1 );
        if ( $sm > 1 ) {
            $gHGVS .= '_' . ( $pos + $reflen );
        }
        $gHGVS .= 'del' . ( ( $ref =~ /^[ACGTN]+$/ ) ? $ref : "" );
        $gHGVS .= 'ins' . $alt if ( $imp =~ /delins/ );
    }
    elsif ( $imp eq 'rep' ) {
        $gHGVS .= ( $pos + 1 );
        if ( $var->{ref_cn} == 1 and $var->{alt_cn} == 2 ) {    # dup
            if ( $sm > 1 ) {
                $gHGVS .= '_' . ( $pos + $var->{replen} );
            }
            $gHGVS .= 'dup' . $var->{rep};
        }
        else {
            $gHGVS .=
              $var->{rep} . '[' . $var->{ref_cn} . '>' . $var->{alt_cn} . ']';
        }
    }
    elsif ( $imp eq 'ref' ) {
        if ($reflen == 1) {
            $gHGVS .= $pos . $ref . '=';
        }
        else {
            $gHGVS .= $pos . '_' . ($pos + $reflen - 1) . '=';
        }
    }
    else {
        confess "Can not recognize type $imp.";
    }
    $var->{gHGVS} = $gHGVS;

    return $var;
}

1;
__END__
