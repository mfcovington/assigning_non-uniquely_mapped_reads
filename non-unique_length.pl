#!/usr/bin/env perl
# non-unique_length.pl
# Mike Covington
# created: 2013-02-01
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Number::Range;
use List::Util qw(min sum);
use Data::Printer;

# Number::Range prints unnecessary warnings; therefore, turn them off
no warnings 'Number::Range';
local $SIG{__WARN__} = sub { warn $_[0] unless $_[0] =~ m|Use of uninitialized value \$previous in string at .*Number/Range.pm line \d+.|};

# Number::Range stores large ranges in a way that causes problems
# when calculating total lengths, unless max_hash_size is high enough.
my $max_hash_size = 10_000;

# For intial design and testing, will use:
# CLUSTER 281:
# Solyc01g096580.2.1      2215
# Solyc01g096590.2.1      2215
# CLUSTER 1001:
# Solyc05g056050.2.1      3493
# Solyc05g056060.2.1      9140
# Solyc05g056070.2.1      9211
# In the future; however, read in clusters from ID.clusters file and extract all appropriate subclusters from ID.groups file.
# This is the subcluster for CLUSTER 281:
# 2215    Solyc01g096580.2.1|Solyc01g096590.2.1
# These are the subclusters for CLUSTER 1001:
# 5718    Solyc05g056060.2.1|Solyc05g056070.2.1
# 3422    Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1
# 71      Solyc05g056050.2.1|Solyc05g056070.2.1
# Subclusters will be stored in HoA with cluster ID being the hash key and the value being an array of subclusters.

my %clusters;
my $cluster_id = 1001;
my @subclusters = (
    "Solyc05g056060.2.1|Solyc05g056070.2.1",
    "Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1",
    "Solyc05g056050.2.1|Solyc05g056070.2.1"
);
$clusters{$cluster_id} = [@subclusters];
$clusters{281} = ["Solyc01g096580.2.1|Solyc01g096590.2.1"];
# p %clusters;

# build set of genes in ALL clusters
my %gene_set;
for my $cluster_id (keys %clusters) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        $gene_set{$_}++ for split /\|/, $subcluster;
    }
}
# p %gene_set;

# build unique counts hash
my %uniq_counts;
for my $cluster_id (keys %clusters) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        $uniq_counts{$_} = 0 for split /\|/, $subcluster;
    }
}
# p %uniq_counts;

# build subcluster counts hash
my %subcluster_counts;
for my $cluster_id (keys %clusters) {
    $subcluster_counts{$_} = 0 for @{ $clusters{$cluster_id} };
}
# p %subcluster_counts;

# build ranges data structure (HoHoO)
my %ranges;
for my $cluster_id ( keys %clusters ) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        $ranges{$subcluster} = ();
        for (split /\|/, $subcluster) {
            $ranges{$subcluster}{$_} = Number::Range->new();
            $ranges{$subcluster}{$_}->set_max_hash_size($max_hash_size);
        }
    }
}
# p %ranges;

# build unique_ranges hash
my %unique_ranges;
for (keys %gene_set) {
    $unique_ranges{$_} = Number::Range->new();
    $unique_ranges{$_}->set_max_hash_size($max_hash_size);
}
# p %unique_ranges;

# for a cluster, read in seqreads. if a gene matches, consider read. increment uniq_count cluster, if applicable, else deal with multi_read
# deal w/ multi_read = extract subcluster ID to use as hash key for hash of arrays (each element in array is read or at least relevant info of read)
my $sam_filename = "1.1.2_rep1_bwa0.6.2.100.sam";
open my $sam_fh, "<", $sam_filename;

my $gene_regex = join "|", map { quotemeta } keys %gene_set;
# say "regex: $gene_regex";
my %counts;
my $relevant_count;
my $total_count;
my %gene_lengths;
while (<$sam_fh>){
    $total_count++;
    next unless /$gene_regex/;

    # collect gene lengths
    if (/^\@/) {
        next unless /^\@SQ/;
        m/SN:(.+)\tLN:(\d+)$/;
        $gene_lengths{$1} = $2;
        next;
    }
    die "Something may be wrong with the sam file header..."
      unless
      join( "", sort keys %gene_set ) eq join( "", sort keys %gene_lengths );


    my ( $subcluster, $best_count, @best_hits ) = best_hits($_);
    # my @best_hits = split /\|/, $subcluster;

    # for my $gene ( keys %gene_set ){
    #     $counts{$gene}++ if /\Q$gene\E/;
    # }

    for my $gene ( @best_hits ){
        $counts{$gene}++ if /\Q$gene\E/;
    }

    push my @positions, first_pos($_);
    push @positions, other_pos( $_, $best_count ) if $best_count > 1;
    die "Numbers of hits and positions don't match"
      unless scalar @best_hits == scalar @positions;

    # sort best_hits and positions both by sorted best_hits order
    # adapted from: http://www.perlmonks.org/?node_id=720562
    my @sorted_index =
      sort { $best_hits[$a] cmp $best_hits[$b] } 0 .. $#best_hits;
    @best_hits = @best_hits[@sorted_index];
    @positions = @positions[@sorted_index];

    # record counts for unique hits, subclusters hits, and relevant reads
    # also populate %unique_ranges
    if ( $best_count == 1 ) {
        $uniq_counts{ $best_hits[0] }++;
        $unique_ranges{ $best_hits[0] }->addrange( $positions[0] );
        next;
    }
    $subcluster_counts{$subcluster}++;
    $relevant_count++;

    # populate range data structure
    my $index = 0;
    for ( sort keys $ranges{$subcluster} ) {
        $ranges{$subcluster}{$_}->addrange( $positions[$index] );
        $index++;
    }

    # # just for monitoring/testing
    # if ($relevant_count >= 50) {
    #     # p %counts;
    #     say "best for $relevant_count: $subcluster";
    #     say "total: $total_count";
    #     # p %gene_lengths;
    #     # p %subcluster_counts;
    #     # p %ranges;
    #     # p %uniq_counts;
    #     last;
    # }

}
# say "total reads: $total_count";
# say "relevant reads: $relevant_count";
# p %counts;
# p %gene_lengths;
# p %subcluster_counts;
# p %ranges;
# p %uniq_counts;

# build gene_multi_lengths hash
my %gene_multi_lengths;
for my $cluster_id (keys %clusters) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        for ( split /\|/, $subcluster ) {
            $gene_multi_lengths{$_} = Number::Range->new();
            $gene_multi_lengths{$_}->set_max_hash_size($max_hash_size);
        }
    }
}
# p %gene_multi_lengths;

# populate gene_multi_lengths hash
for my $subcluster ( sort keys %subcluster_counts) {
    say "$subcluster: $subcluster_counts{$subcluster} reads";

    for my $gene ( sort keys $ranges{$subcluster} ) {
        $gene_multi_lengths{$gene}->addrange( $ranges{$subcluster}{$gene}->range);
        my $multi_length;
        $multi_length++ for $ranges{$subcluster}{$gene}->range;
        # my $unique_length = $gene_lengths{$gene} - $multi_length;
        # die "Unique length is less than zero... $gene $subcluster" if $unique_length < 0;
        # say "$gene: $multi_length / $unique_length (multi vs unique length)";
        say "$gene: $multi_length (non-unique length)";
    }
}
# p %gene_multi_lengths;

# THIS IS THE CORRECT WAY!!
# build/populate unique_lengths hash
my %unique_lengths;
for my $gene ( keys %gene_set ) {
    $unique_lengths{$gene}++ for $unique_ranges{$gene}->range;
    my $format = $unique_ranges{$gene}->range;
    say "$gene U1-ranges: $format";
}
# p %unique_ranges;
p %unique_lengths;

for my $gene ( keys %gene_multi_lengths ) {
    my $format = $gene_multi_lengths{$gene}->range;
    say "$gene MM-ranges: $format";
}

# NOT THE RIGHT WAY
# my %unique_ranges2;
# $unique_ranges2{$_} = Number::Range->new() for keys %gene_set;
# $unique_ranges2{$_}->set_max_hash_size(3000) for keys %gene_set;
# for my $gene ( keys %gene_lengths ) {
#     my $length = $gene_lengths{$gene};
#     $unique_ranges2{$gene}->addrange("1..$length");
# }
# for my $gene ( sort keys %gene_multi_lengths ) {
#     my $format1 = $gene_multi_lengths{$gene}->range;
#     say "$gene NN-ranges: $format1";
#     $unique_ranges2{$gene}->delrange($gene_multi_lengths{$gene}->range);
#     my $format = $unique_ranges2{$gene}->range;
#     say "$gene U2-ranges: $format";
# }
# # p %unique_ranges2;

# NOT THE RIGHT WAY
# # build unique_lengths hash
# my %uniq_lengths = %gene_lengths;
# for my $gene ( sort keys %gene_multi_lengths ) {
#     $uniq_lengths{$gene}-- for $gene_multi_lengths{$gene}->range;
#     # my $format = $gene_multi_lengths{$gene}->range;
#     # say "$gene U2-ranges: $format";
# }
# p %uniq_lengths;
p %gene_lengths;

# my $format = $ranges{'Solyc05g056060.2.1|Solyc05g056070.2.1'}{'Solyc05g056070.2.1'}->range;
# say $format;

sub best_hits {
    # adapted from harvest_gene_ids.pl
    my $read = shift;
    my ($best_hits) = $read =~ m|X0:i:(\d+)|;
    my @genes = $read =~ m|(Solyc\d{2}g\d{6}\.\d\.\d)|g;

    my $delimiter = "|";
    my $max_best  = 7;
    my @bests = @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
    my $best_hits_string = join $delimiter, sort @bests;
    my $best_hits_count = scalar @bests;
    return ( $best_hits_string, $best_hits_count, @bests );
}

sub first_pos {
    my $read = shift;
    my ( $start, $cigar ) = ( split /\t/, $read )[ 3, 5 ];
    my $end = calc_end( $start, $cigar );

    return ( "$start..$end" );
}

sub other_pos {
    my ( $read, $best_count ) = @_;
    my ($others) = $read =~ m|XA:Z:([^\s]+)|;

    my @other_positions;
    my $best_seen_count = 1;
    for ( split /;/, $others ) {
        $best_seen_count++;
        my ( $start, $cigar ) = m|^[^,]+,[+-](\d+),([^,]+),\d+$|;
        my $end = calc_end( $start, $cigar );
        push @other_positions, "$start..$end";
        last if $best_seen_count == $best_count;
    }

    return ( @other_positions );
}

sub calc_end {
    my ( $start, $cigar ) = @_;

    my $matches    = sum( $cigar =~ m|(\d+)M|g );
    my $insertions = sum( $cigar =~ m|(\d+)I|g ) || 0;
    my $deletions  = sum( $cigar =~ m|(\d+)D|g ) || 0;

    my $length = $matches + $deletions - $insertions;
    my $end    = $start + $length - 1;

    return $end;
}



# my $range = Number::Range->new(""); # convert this into a hash w/ keys == geneIDs in cluster

# $range->addrange("600..450");
# my $format =  $range->range;
# say $format;

# my $i;
# $i++ for ($range->range);
# say $i;

__END__
{
    Solyc01g096580.2.1   496,
    Solyc01g096590.2.1   421,
    Solyc05g056050.2.1   924,
    Solyc05g056060.2.1   1111,
    Solyc05g056070.2.1   431
}
{
    Solyc01g096580.2.1   768,
    Solyc01g096590.2.1   639,
    Solyc05g056050.2.1   1030,
    Solyc05g056060.2.1   2473,
    Solyc05g056070.2.1   1112
}
Solyc01g096580.2.1|Solyc01g096590.2.1: 2215 reads
Solyc01g096580.2.1: 409 (non-unique length)
Solyc01g096590.2.1: 409 (non-unique length)
Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1: 3422 reads
Solyc05g056050.2.1: 529 (non-unique length)
Solyc05g056060.2.1: 529 (non-unique length)
Solyc05g056070.2.1: 529 (non-unique length)
Solyc05g056050.2.1|Solyc05g056070.2.1: 71 reads
Solyc05g056050.2.1: 45 (non-unique length)
Solyc05g056070.2.1: 45 (non-unique length)
Solyc05g056060.2.1|Solyc05g056070.2.1: 5718 reads
Solyc05g056060.2.1: 827 (non-unique length)
Solyc05g056070.2.1: 825 (non-unique length)
Solyc05g056050.2.1 U1-ranges: 14..501,504..679,719..801,807..983
Solyc05g056060.2.1 U1-ranges: 17..615,646..689,1266..1309,1354..1428,1611..1794,1820..1984
Solyc01g096580.2.1 U1-ranges: 28..101,127..248,429..728
Solyc01g096590.2.1 U1-ranges: 17..96,121..242,423..641
Solyc05g056070.2.1 U1-ranges: 218..302,323..406,847..1108
Solyc05g056050.2.1 MM-ranges: 342..404,409..458,460..761,763..849,852..879
Solyc05g056060.2.1 MM-ranges: 1391..1646,1757..1861,1948..2473
Solyc01g096580.2.1 MM-ranges: 60..169,206..471,475..507
Solyc01g096590.2.1 MM-ranges: 54..163,200..465,469..501
Solyc05g056070.2.1 MM-ranges: 6..889
[Finished in 63.2s]
