#!/usr/local/bin/perl

use strict;
use warnings;
use Digest::MD5 qw(md5_hex);

=head1 NAME

grasta.ph

=cut

=head1 DESCRIPTION

provides a subroutine for performing a grasta pairwise alignment between two sequences, protein or nucleotide.  It returns information for the top match as indicated in the subroutine specification below.

=cut




=item perform_fasta_align()

=over 4

B<Description:> Performs a pairwise alignment between two sequences (protein or nucleotide).

B<Parameters:> $seq1, $seq2

$seq1 and $seq2 are scalar values holding a raw line of text.  Do NOT provide a FASTA-formatted sequence.

B<Returns:> ($per_id, $overlap)

$per_id is the percent identity across the top match.
$overlap is the number of residues in this top alignment.

if the fasta program is not located in your PATH setting, set 
our \$FASTAPATH = /full/path/to/fasta

=back

=cut

    ;


our $FASTAPATH = "";
our $SEE;
our %ALIGN_CACHE;
our $FASTA_PROG;
our $TMP_COUNTER = 0;

####
sub perform_fasta_align {
    my $seq1 = shift;
    my $seq2 = shift;

    ## shortcut: identical sequences need no external alignment
    if ($seq1 eq $seq2) {
        return (100, length($seq1));
    }

    ## cache results; protein comparisons are highly redundant across a locus
    my $cache_key = md5_hex($seq1) . "-" . md5_hex($seq2);
    if (my $cached = $ALIGN_CACHE{$cache_key}) {
        return (@$cached);
    }

    ## resolve the fasta binary only once per process
    unless ($FASTA_PROG) {
        my $prog;
        if ($FASTAPATH) {
            $prog = $FASTAPATH;
        } else {
            $prog = `which fasta`; chomp($prog);
            ## fall back to the fasta bundled with the PASA installation
            if ((!$prog || !-x $prog) && $ENV{PASAHOME} && -x "$ENV{PASAHOME}/bin/fasta") {
                $prog = "$ENV{PASAHOME}/bin/fasta";
            }
        }
        die "Cannot find program fasta (set \$FASTAPATH in fasta.ph or put fasta in PATH)\n" unless $prog && -s $prog && -x $prog;
        $FASTA_PROG = $prog;
    }

    $TMP_COUNTER++;
    my $tmp_token = "$$-" . &hashCode($seq1) . "-" . &hashCode($seq2) . "-" . $TMP_COUNTER;


    my $file1 = "/tmp/$tmp_token.seq1";
    open (SEQ1, ">$file1") or die;
    print SEQ1 ">seq1\n" . &FASTA_format($seq1);
    close SEQ1;
    my $file2 = "/tmp/$tmp_token.seq2";
    open (SEQ2, ">$file2") or die;
    print SEQ2 ">seq2\n" . &FASTA_format ($seq2);
    close SEQ2;

    my $result_file = "/tmp/$tmp_token.grasta_result";

    my $cmd = "$FASTA_PROG -p $file1 $file2 > $result_file";
    my $ret = system ($cmd);
    if ($ret) {
		die "Error: couldn't perform alignment using prog $FASTA_PROG.\n";
    }

    ## Parse the output file
    open (OUTPUT, $result_file) or die "ERROR: Sorry, can't open $result_file";
    my $per_id = 0;
    my $overlap = 0;
    while (<OUTPUT>) {
		#print;

		if (/\s(\S+)% identity .* in (\d+) (nt|aa) overlap/) {
			$per_id = $1;
			$overlap = $2;
			last;
		}

	}
    close OUTPUT;
    print "FASTA: per_id: $per_id, overlap: $overlap\n" if $SEE;
    unlink ($file1, $file2, $result_file);

    $ALIGN_CACHE{$cache_key} = [$per_id, $overlap];

    return ($per_id, $overlap);
}


sub FASTA_format {
    my $seq = shift;
    $seq =~ s/(\w{60})/$1\n/g;
    return ($seq);
}


sub hashCode {
    my ($string) = @_;

    my $hashcode = 7;
    
    foreach my $char (split (//, $string)) {
        $hashcode = 31 * $hashcode + ord($char);
        #print STDERR "String: $string => $hashcode\n";
        $hashcode &= 2**32-1;
    }
    
    $hashcode = sprintf("%x", $hashcode);
    
    # print STDERR "String: $string => $hashcode\n";
    
    return($hashcode);
    
}



1; #end of ph.
