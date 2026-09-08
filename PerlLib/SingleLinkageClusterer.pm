#!/usr/bin/env perl

package main;
our $CLUSTERPATH;


package SingleLinkageClusterer;

## package not to be instantiated.  Just provides a namespace.

## Input: Array containing array-refs of pairs:
##               @_ = ( [1,2], [2,3], [6,7], [7,8], ...)
## Output: Array of all clusters as array-refs.
##              return ([1,2,3] , [6,7,8], ...)

use strict;
use warnings;

__run_test() unless caller;

sub build_clusters {
    my @pairs = @_;

    ## pure-perl union-find; replaces the external slclust invocation,
    ## which cost a fork+exec (and /tmp I/O) per cluster.
    my %parent;
    my $find = sub {
        my ($x) = @_;
        my $root = $x;
        while ($parent{$root} ne $root) {
            $root = $parent{$root};
        }
        ## path compression
        while ($parent{$x} ne $root) {
            my $next = $parent{$x};
            $parent{$x} = $root;
            $x = $next;
        }
        return ($root);
    };

    foreach my $pair (@pairs) {
        my ($a, $b) = @$pair;
        $parent{$a} = $a unless (exists $parent{$a});
        $parent{$b} = $b unless (exists $parent{$b});
        my ($ra, $rb) = ($find->($a), $find->($b));
        $parent{$rb} = $ra if ($ra ne $rb);
    }

    my %clusters;
    foreach my $elem (keys %parent) {
        push (@{$clusters{$find->($elem)}}, $elem);
    }

    return (values %clusters);
}


############
## Testing
###########

sub __run_test {
    
    my @pairs = ( [1,2], [2,3], [4,5] );

    my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);

    use Data::Dumper;
    
    print "Input: " . Dumper(\@pairs);
    print "Output: " . Dumper(\@clusters);

    exit(0);
}


1;
