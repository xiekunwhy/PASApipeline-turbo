#!/usr/local/bin/perl

package main;
our $SEE;


package CDNA::PASA_alignment_assembler;

=head1 NAME

CDNA::PASA_alignment_assembler

=cut

=head1 DESCRIPTION

This module is used to assemble compatible cDNA alignments.  The algorithm is as follows:
must describe this here.

=cut

use strict;
use CDNA::CDNA_alignment;
use Data::Dumper;
use Carp;
use FindBin;

## File scoped globals:
my $DELIMETER = "$;,";
our $FUZZLENGTH = 20;
our $PASA_BIN; ## resolve the pasa binary only once per process
our $PASA_BIN_CHECKED = 0; ## whether resolution was attempted (cache negative results too)
my $TMP_COUNTER = 0;


=item new()

=over 4

B<Description:> instantiates a new cDNA assembler obj.

B<Parameters:> none

B<Returns:> $obj_href

$obj_href is the object reference newly instantiated by this new method.

=back

=cut

sub new  {
    my $package_name = shift;
    my $self = {};
    bless ($self, $package_name);
    $self->_init(@_);
    return ($self);
}


sub _init {
    my $self = shift;
    $self->{incoming_alignments} = []; #these are the alignments to be assembled.
    $self->{assemblies} = []; #contains list of all singletons and assemblies.
    $self->{fuzzlength} = $FUZZLENGTH;  #default setting.

    my $pasa_bin = $PASA_BIN;
    unless ($PASA_BIN_CHECKED) {
        $PASA_BIN_CHECKED = 1;
        $pasa_bin = `which pasa 2>/dev/null`;
        $pasa_bin =~ s/\s//g;

        unless ($pasa_bin && -x $pasa_bin) {
            ## fall back to the binary bundled in this PASA tree,
            ## independent of PATH/`which` quirks in batch-job environments.
            my $bundled = "$FindBin::Bin/../bin/pasa";
            $pasa_bin = $bundled if (-x $bundled);
        }

        ## not fatal here anymore: 1- and 2-alignment assemblies are computed
        ## in pure Perl (see pasa_cpp_assemblies); the binary is only required
        ## when actually shelling out (3+ alignments).
        if ($pasa_bin && -x $pasa_bin) {
            $PASA_BIN = $pasa_bin;
        } else {
            $pasa_bin = undef;
        }
    }

    $self->{pasa_bin} = $pasa_bin;

}


=item assemble_alignments()
    
=over 4
    
B<DESCRIPTION:> assembles a series of cDNA aligmnments into one or more cDNA assemblies using a directed acyclic graph.

B<Parameters:> @alignments

@alignments is an array of CDNA::CDNA_alignment objects

B<Returns:> none.

=back

=cut


sub assemble_alignments {
    my $self = shift;
    my @alignments = @_;
    @alignments = sort {$a->{lend}<=>$b->{lend}} @alignments; #keep in order of lend across genomic sequence to provide a layout.
    $self->{incoming_alignments} = [@alignments];
    my %accs;
    my %spliced_orientations; # track so can set later in each assembly based on content.
    my %aligned_orientations;
    my %FL_accs;
    
    foreach my $alignment (@alignments) {
        my $acc = $alignment->get_acc();
        $accs{$acc} = 0;
        my $spliced_orient = $alignment->get_spliced_orientation();
        $spliced_orientations{$acc} = $spliced_orient;
        my $aligned_orient = $alignment->get_orientation() or confess "Error, no orientation set for " . $alignment->toToken();
        $aligned_orientations{$acc} = $aligned_orient;
        $FL_accs{$acc} = $alignment->is_fli();
    }
    $self->{accs_in_assemblies} = \%accs;
    
    my $num_alignments = $#alignments + 1;
    
    $self->force_flexorient('+');
    my @assemblies = $self->pasa_cpp_assemblies('+');
    $self->force_flexorient('-');
    push (@assemblies, $self->pasa_cpp_assemblies('-'));
    
    # sort in order of decreasing score
    @assemblies = reverse sort {$a->{num_contained_aligns}<=>$b->{num_contained_aligns}} @assemblies;  
    
    if ($SEE) {
        print "\n\nScore summary for all assemblies (nr set unchosen):\n\n";
        foreach my $assembly (@assemblies) {
            print "score: " . $assembly->{num_contained_aligns} . ", " . $assembly->toToken . "\n";
        }
        print "\n\n";
    }
    my @report_assemblies;
    my $still_missing = 1;
    foreach my $assembly (@assemblies) {
        print "\nAnalyzing assembly.\n" if $SEE;
        my $contained_aligns_aref = $assembly->{contained_aligns};
        my $have_unseen = 0;
        my $spliced_orient = '?';
        my $is_fli = 0;
        my %aligned_orient_counts;

        foreach my $acc (@$contained_aligns_aref) {
            print "got: $acc\n" if $SEE;
            # check to see if we've encountered this one yet.
            unless ($accs{$acc}) {
                $have_unseen = 1;
            }
            $accs{$acc} = 1;
            my $curr_spliced_orient = $spliced_orientations{$acc};
            print "curr_spliced_orient: $curr_spliced_orient\n" if $SEE;
            if ($curr_spliced_orient ne '?') {
                if ($spliced_orient ne '?' && $spliced_orient ne $curr_spliced_orient) {
                    ## cannot have conflicting spliced orientations in the assembly: corruption.
                    die "Fatal: conflicting spliced orientations in current PASA assembly results (spliced_orient: $spliced_orient, $acc has $curr_spliced_orient).\n";
                }
                $spliced_orient = $curr_spliced_orient; ## retain original spliced orientation.
            }
            if ( (!$is_fli) && $FL_accs{$acc}) {
                $is_fli = 1;
            }
            ## track aligned orientation
            $aligned_orient_counts{ $aligned_orientations{$acc} } ++;

        }
        
        ## set assembly orientation
        $assembly->set_spliced_orientation($spliced_orient);
        if ($spliced_orient eq '?') {
            ## set aligned orientation based on a majority vote
            ## (count desc; ties broken alphabetically so results are
            ##  reproducible - hash key order is randomized per process)
            my @orients = sort { $aligned_orient_counts{$b} <=> $aligned_orient_counts{$a}
                                 || $a cmp $b } keys %aligned_orient_counts;
            my $winning_aligned_orient = shift @orients;
            $assembly->set_orientation($winning_aligned_orient);
        }
        else {
            # got spliced orientation, use it for aligned orientation too.
            $assembly->set_orientation($spliced_orient);
        }
        
        $assembly->set_fli_status($is_fli);
        
        if ($have_unseen) {
            print $assembly->toToken . "\n" if $SEE;
            push (@report_assemblies, $assembly);
        }
        $still_missing = 0;
        foreach my $key (keys %accs) {
            if (! $accs{$key}) {
                $still_missing = 1;
                print "still missing: $key\n" if $SEE;
            }
        }
        if (! $still_missing) {
            last; #got them all.
        }
    }
    
    $self->{assemblies} = \@report_assemblies;
    
    if ($still_missing) {
        die "Didn't obtain assemblies describing all maximal assemblies.\n";
    }
    
    
}


sub pasa_cpp_assemblies {
    my $self = shift;
    my $forced_orient = shift;

    my $prev_input_sep = $/;
    $/ = "\n";

    my $sequence_ref;
    my $incoming_alignments_aref = $self->{incoming_alignments};

    ## Fast path: 1 or 2 alignments are assembled in pure Perl, mirroring the
    ## pasa c++ binary's logic exactly (see _assemble_pair_perl), avoiding two
    ## external process spawns per call.  This matters in the annotation
    ## comparer where pairwise compatibility checks are done per
    ## (gene model x transcript) pair.  Set env var PASA_NO_PERL_PAIR_ASSEMBLY=1
    ## to fall back to the external binary.
    if (@$incoming_alignments_aref <= 2 && ! $ENV{PASA_NO_PERL_PAIR_ASSEMBLY}) {
        my @assemblies = $self->_assemble_pair_perl($forced_orient);
        $/ = $prev_input_sep;
        return (@assemblies);
    }

    # create input file for pasa-cpp implementation:
    ## collision-free temp file token: pid + thread id + counter
    my $thread_id = ($INC{'threads.pm'}) ? threads->tid() : 0;
    $TMP_COUNTER++;
    my $uniq_token = "$$-$thread_id-$TMP_COUNTER";
    
    my $tmpdir = $ENV{TMPDIR};
    unless ($tmpdir) {
        if (-d '/tmp') {
            $tmpdir = "/tmp";
        }
        else {
            $tmpdir = ".";
        }
    }

    my $pasa_input = "$tmpdir/pasa.$uniq_token.$forced_orient.in";
    my $pasa_output = "$tmpdir/pasa.$uniq_token.$forced_orient.out";
    my @assemblies;
  
    open (TMPIN, ">$pasa_input") or die "Can't open file $pasa_input";
    foreach my $alignment (@$incoming_alignments_aref) {
        my $acc = $alignment->get_acc();
        ## commas not allowed in acc name:
        if ($acc =~ /,/) { 
            die "ERROR, $acc accession contains comma(s).  This is not allowed.\n";
        }
        my $orient = $alignment->{fixed_orient};
        my $alignText = "$acc,$orient";
        unless (ref $sequence_ref) {
            $sequence_ref = $alignment->get_genomic_seq_ref();
        }
        foreach my $seg ($alignment->get_alignment_segments()) {
            my ($lend, $rend) = $seg->get_coords();
            $alignText .= ",$lend-$rend";
        }
        print TMPIN $alignText . "\n";
    }
    close TMPIN;
    
    if ($SEE) {
        print "PASA_INPUT ($forced_orient):\n====\n";
        system "cat $pasa_input";
        print "====\n";
    }
    my $pasa_bin = $self->{pasa_bin};
    unless ($pasa_bin && -x $pasa_bin) {
        confess "Error, pasa binary isn't executable or couldn't be found (required for assembling more than two alignments).";
    }
    my $cmd = $pasa_bin . " $pasa_input > $pasa_output";
    my $ret = system $cmd;
    if ($ret) {
        system "mv $pasa_input pasa_killer.input";
        print STDERR "PASA died on input file.  See pasa_killer.input";
        die;
    } else {
        
        # process the output.
        open (TMPOUT, $pasa_output) or die "Can't open $pasa_output";
        while (<TMPOUT>) {
            if (/assembly:\s\(\d+\)\scontains\salignments:\s\[([^\]]+)\]\swith\sstructure\s\[([^\]]+)\]/) {
                print "Extracting assembly output: $_" if $SEE;
                my $acclist = $1;
                my $aligndescript = $2;
                my @x = split (/,/, $aligndescript);
                
                shift @x;
                my $orient = shift @x;
                my @alignSegs;
                my $length = 0;
                foreach my $coordset (@x) {
                    my ($lend, $rend) = sort {$a<=>$b} split (/-/, $coordset);
                    my $seg = new CDNA::Alignment_segment($lend, $rend);
                    $length += ($rend - $lend) + 1;
                    push (@alignSegs, $seg);
                }
                my $assembly = new CDNA::CDNA_alignment($length, \@alignSegs, $sequence_ref);
                
                my @accs = split (/,/, $acclist);
                my $num_accs = $#accs + 1;
                $assembly->{contained_aligns} = [@accs];
                $assembly->{num_contained_aligns} = $num_accs;
                
                $acclist =~ s/,/\//g; #convert list of accessions into a new accession representing a single entry (unity)
                # if we keep the commas, use of this assembly in future PASA runs will break the assembler
                # because of the input file requirements.
                $assembly->set_acc($acclist);
                
                push (@assemblies, $assembly);
            }
            
        }
        close TMPOUT;
        
        if ($SEE) {
            print "PASA_OUTPUT ($forced_orient):\n####\n";
            system "cat $pasa_output";
            print "####\n";
        }
        unlink ($pasa_input, $pasa_output) unless $SEE;
        
    }
    
    $/ = $prev_input_sep; ## restore
    
    return (@assemblies);
}


####
## Pure-Perl reimplementation of the pasa c++ binary for the 1- and
## 2-alignment cases.  Mirrors cdna_alignment_assembler.cpp:
##  - canMerge(): span overlap, equal (fixed) orientation, then a lockstep walk
##    over lend-sorted segments checking splice-junction boundaries with
##    fuzzlength tolerance.
##  - mergeAlignments(): splice-aware coordinate merging.
##  - For a single alignment, the binary simply reports the alignment itself
##    as one assembly.
## Splice-junction flags are derived *structurally* from segment position
## (first segment: right junction; last: left; internal: both; single: none),
## exactly as the c++ side does in CDNA_alignment::refineAlignment().
sub _assemble_pair_perl {
    my $self = shift;
    my $forced_orient = shift;

    my @alignments = @{$self->{incoming_alignments}};
    my $sequence_ref = $alignments[0]->get_genomic_seq_ref();

    if (scalar(@alignments) == 1) {
        my @segs = &_sorted_seg_coords($alignments[0]);
        my $acc = $alignments[0]->get_acc();
        return (&_mk_assembly_from_coords(\@segs, [$acc], $sequence_ref));
    }

    my ($A, $B) = @alignments; # incoming alignments are lend-sorted by assemble_alignments()

    my $compatible = 0;
    if ($A->{fixed_orient} eq $B->{fixed_orient}) {
        my @a_segs = &_sorted_seg_coords($A);
        my @b_segs = &_sorted_seg_coords($B);
        if (&_pair_can_merge(\@a_segs, \@b_segs, $self->{fuzzlength})) {
            $compatible = 1;
            my @merged = &_pair_merge_coords(\@a_segs, \@b_segs);
            my @accs = ($A->get_acc(), $B->get_acc());
            return (&_mk_assembly_from_coords(\@merged, \@accs, $sequence_ref));
        }
    }

    unless ($compatible) {
        ## not mergeable: c++ emits one singleton assembly per alignment,
        ## highest lend index first (its bin walk processes them in reverse)
        my @a_segs = &_sorted_seg_coords($A);
        my @b_segs = &_sorted_seg_coords($B);
        return (&_mk_assembly_from_coords(\@b_segs, [$B->get_acc()], $sequence_ref),
                &_mk_assembly_from_coords(\@a_segs, [$A->get_acc()], $sequence_ref));
    }
}


## mirror of common_subs.cpp overlap(): inclusive coordinate overlap
sub _coords_overlap {
    my ($a_l, $a_r, $b_l, $b_r) = @_;
    return ($a_l <= $b_r && $a_r >= $b_l) ? 1 : 0;
}


## lend-sorted [lend, rend] coordsets for an alignment
## (c++ CDNA_alignment::refineAlignment sorts segments by lend)
sub _sorted_seg_coords {
    my ($alignment) = @_;
    my @segs = map { my ($l, $r) = $_->get_coords(); [$l, $r] } $alignment->get_alignment_segments();
    @segs = sort { $a->[0] <=> $b->[0] } @segs;
    return (@segs);
}


## mirror of cdna_alignment_assembler.cpp canMerge() for two alignments,
## given as lend-sorted segment coordset listrefs.
sub _pair_can_merge {
    my ($a_segs, $b_segs, $fuzzlength) = @_;

    ## span overlap check
    unless (&_coords_overlap($a_segs->[0][0], $a_segs->[-1][1], $b_segs->[0][0], $b_segs->[-1][1])) {
        return (0);
    }

    ## find first overlapping segment pair
    my ($i, $j) = (-1, -1);
  SEG_SEARCH:
    for (my $ii = 0; $ii <= $#$a_segs; $ii++) {
        for (my $jj = 0; $jj <= $#$b_segs; $jj++) {
            if (&_coords_overlap($a_segs->[$ii][0], $a_segs->[$ii][1], $b_segs->[$jj][0], $b_segs->[$jj][1])) {
                ($i, $j) = ($ii, $jj);
                last SEG_SEARCH;
            }
        }
    }

    return (0) if ($i == -1 || $j == -1); # couldn't align two segments of the alignments
    return (0) unless ($i == 0 || $j == 0); # couldn't map first segment of either to the other

    my $n_a = scalar(@$a_segs);
    my $n_b = scalar(@$b_segs);

    ## lockstep comparison of the overlapping segment series
    while ($i < $n_a && $j < $n_b) {
        my ($a_l, $a_r) = @{$a_segs->[$i]};
        my ($b_l, $b_r) = @{$b_segs->[$j]};

        if (&_coords_overlap($a_l, $a_r, $b_l, $b_r)) {
            ## structural splice-junction flags (by segment position, as the c++ side sets them)
            my $a_left_j  = ($i > 0) ? 1 : 0;
            my $a_right_j = ($i < $n_a - 1) ? 1 : 0;
            my $b_left_j  = ($j > 0) ? 1 : 0;
            my $b_right_j = ($j < $n_b - 1) ? 1 : 0;

            if ($a_left_j || $b_left_j) {
                return (0) if ($a_left_j && $b_left_j && $a_l != $b_l); # diff left splice sites
                return (0) if ($a_left_j && ($b_l + $fuzzlength < $a_l)); # not within fuzzdist
                return (0) if ($b_left_j && ($a_l + $fuzzlength < $b_l)); # not within fuzzdist
            }
            if ($a_right_j || $b_right_j) {
                return (0) if ($a_right_j && $b_right_j && $a_r != $b_r); # diff right splice sites
                return (0) if ($a_right_j && ($b_r - $fuzzlength > $a_r)); # not within fuzzdist
                return (0) if ($b_right_j && ($a_r - $fuzzlength > $b_r)); # not within fuzzdist
            }
        } else {
            ## two ordered segments do not overlap each other
            return (0);
        }
        $i++; $j++;
    }

    return (1);
}


## mirror of cdna_alignment_assembler.cpp mergeAlignments() for two alignments,
## given as lend-sorted segment coordset listrefs.  Returns merged, lend-sorted
## coordsets (the c++ constructor re-sorts its merged segments by lend).
sub _pair_merge_coords {
    my ($a_segs, $b_segs) = @_;

    my (%leftsplice, %rightsplice);
    for (my $i = 0; $i <= $#$a_segs; $i++) {
        $leftsplice{$a_segs->[$i][0]} = 1 if ($i > 0);
        $rightsplice{$a_segs->[$i][1]} = 1 if ($i < $#$a_segs);
    }
    for (my $j = 0; $j <= $#$b_segs; $j++) {
        $leftsplice{$b_segs->[$j][0]} = 1 if ($j > 0);
        $rightsplice{$b_segs->[$j][1]} = 1 if ($j < $#$b_segs);
    }

    my @merged;
    foreach my $a_seg (@$a_segs) {
        my ($a_l, $a_r) = @$a_seg;
        my ($merged_lend, $merged_rend) = (-1, -1);
        foreach my $b_seg (@$b_segs) {
            my ($b_l, $b_r) = @$b_seg;
            if (&_coords_overlap($a_l, $a_r, $b_l, $b_r)) {
                $merged_lend = ($leftsplice{$a_l}) ? $a_l
                             : (($leftsplice{$b_l}) ? $b_l
                             : (($a_l < $b_l) ? $a_l : $b_l));
                $merged_rend = ($rightsplice{$a_r}) ? $a_r
                             : (($rightsplice{$b_r}) ? $b_r
                             : (($a_r > $b_r) ? $a_r : $b_r));
                last;
            }
        }
        if ($merged_lend != -1 && $merged_rend != -1) {
            push (@merged, [$merged_lend, $merged_rend]);
        } else {
            ## no overlap; keep the a1 coords.
            push (@merged, [$a_l, $a_r]);
        }
    }

    ## add the unconsumed b coordsets
    foreach my $b_seg (@$b_segs) {
        my ($b_l, $b_r) = @$b_seg;
        my $overlap_flag = 0;
        foreach my $m (@merged) {
            if (&_coords_overlap($b_l, $b_r, $m->[0], $m->[1])) {
                $overlap_flag = 1;
                last;
            }
        }
        unless ($overlap_flag) {
            push (@merged, [$b_l, $b_r]);
        }
    }

    @merged = sort { $a->[0] <=> $b->[0] } @merged;

    return (@merged);
}


## builds a CDNA::CDNA_alignment the same way the pasa output parser does:
## fresh segments without cdna coords/per_id, contained_aligns + acc set.
sub _mk_assembly_from_coords {
    my ($seg_coords_aref, $accs_aref, $sequence_ref) = @_;

    my @alignSegs;
    my $length = 0;
    foreach my $coordset (@$seg_coords_aref) {
        my ($lend, $rend) = @$coordset;
        my $seg = new CDNA::Alignment_segment($lend, $rend);
        $length += ($rend - $lend) + 1;
        push (@alignSegs, $seg);
    }
    my $assembly = new CDNA::CDNA_alignment($length, \@alignSegs, $sequence_ref);

    my @accs = @$accs_aref;
    $assembly->{contained_aligns} = [@accs];
    $assembly->{num_contained_aligns} = scalar(@accs);
    $assembly->set_acc(join("/", @accs));

    return ($assembly);
}


sub force_flexorient {
    my $self = shift;
    my $orient = shift;
    my $alignments_aref = $self->{incoming_alignments};
    my $num_alignments = $#{$alignments_aref} + 1;
    ## Fix orientations of fli and multi-segment alignments:
    for (my $i = 0; $i < $num_alignments; $i++) {
        my $alignment = $alignments_aref->[$i];
        my $num_segments = $alignment->get_num_segments();
        my $spliced_orientation = $alignment->get_spliced_orientation();
        
        if ($spliced_orientation =~ /^[+-]$/) { # set specifically
            $alignment->{fixed_orient} = $spliced_orientation;  ## adding tag, using only in this module.
        } else {
            print "Setting $i to flex orient: $orient.\n" if $SEE;
            $alignment->{fixed_orient} = $orient;
        }
    }
    
}


sub unique_entries {
    my @x = @_;
    my %z;
    foreach my $y (@x) {
        $z{$y}=1;
    }
    return (keys %z);
}


=item get_assemblies()
    
=over 4
    
B<Description:> returns all the alignment assemblies resulting from the assembly procedure.

B<Parameters:> none.

B<Returns:> @assemblies

@assemblies is an array of CDNA::CDNA_alignment objects.

use the get_acc() method of the alignment object to retrieve all the accessions of the cDNAs that were merged into the assembly.

=back

=cut


sub get_assemblies {
    my $self = shift;
    return (@{$self->{assemblies}});
}

=item toAlignIllustration()

=over 4

B<Description:> illustrates the individual cDNAs to be assembled along with the final products.

B<Parameters:> $max_line_chars(optional)

$max_line_chars is an integer representing the maximum number of characters in a single line of output to the terminal.  The default is  100.

B<Returns:> $alignment_illustration_text

$alignment_illustration_text is a string containing a paragraph of text which illustrates the alignments and assemblies. An example is below:

 --->    <-->  <----->     <--->    <----------------	(+)gi|1199466

 --->    <-->  <----->     <--->    <------------   (+)gi|1209702
                        
---->    <-->  <----	(+)AV827070

---->    <-->  <---	(+)AV828861

---->    <-->  <---      (+)AV830936

 --->    <-->  <-	(+)H36350

ASSEMBLIES: (1)

---->    <-->  <----->     <--->    <---------------- (+) gi|1199466, gi|1209702, AV827070, AV828861, AV830936, H36350




=back

=cut

    ;    

sub toAlignIllustration () {
    my $self = shift;
    my $max_line_chars = shift;
    $max_line_chars = ($max_line_chars) ? $max_line_chars : 100; #if not specified, 100  chars / line is default.
    
    ## Get minimum coord for relative positioning.
    my @coords;
    my @alignments = @{$self->{incoming_alignments}};
    foreach my $alignment (@alignments) {
        my @c = $alignment->get_coords();
        push (@coords, @c);
    }
    @coords = sort {$a<=>$b} @coords;
    print "coords: @coords\n" if $::SEE;
    my $min_coord = shift @coords;
    my $max_coord = pop @coords;
    my $rel_max = $max_coord - $min_coord;
    my $alignment_text = "";
    ## print each alignment followed by assemblies:
    my $num_alignments = $#alignments + 1;
    $alignment_text .= "Individual Alignments: ($num_alignments)\n";
    my $i = 0;
    foreach my $alignment (@alignments) {
        $alignment_text .= (sprintf ("%3d ", $i)) . $alignment->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
        $i++;
    }
    
    my @assemblies = @{$self->{assemblies}};
    my $num_assemblies = $#assemblies + 1;
    $alignment_text .= "\n\nASSEMBLIES: ($num_assemblies)\n";
    foreach my $assembly (@assemblies) {
        $alignment_text .= "    " . $assembly->toAlignIllustration($min_coord, $rel_max, $max_line_chars) . "\n";
    }
    
    return ($alignment_text);
}


1;
