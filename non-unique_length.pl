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
        $ranges{$subcluster}{$_} = Number::Range->new()
          for split /\|/, $subcluster;
    }
}
# p %ranges;

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
    if ( $best_count == 1 ) {
        $uniq_counts{ $best_hits[0] }++;
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
say "total reads: $total_count";
say "relevant reads: $relevant_count";
p %counts;
p %gene_lengths;
p %subcluster_counts;
# p %ranges;
p %uniq_counts;

my %uniq_lengths = %gene_lengths;
# build gene_multi_lengths hash
my %gene_multi_lengths;
for my $cluster_id (keys %clusters) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        $gene_multi_lengths{$_} = Number::Range->new() for split /\|/, $subcluster;
    }
}
# p %gene_multi_lengths;

for my $subcluster ( sort keys %subcluster_counts) {
    say "$subcluster: $subcluster_counts{$subcluster} reads";

    for my $gene ( sort keys $ranges{$subcluster} ) {
        $gene_multi_lengths{$gene}->addrange( $ranges{$subcluster}{$gene}->range);
        my $multi_length;
        $multi_length++ for $ranges{$subcluster}{$gene}->range;
        # $uniq_lengths{$gene}-- for $ranges{$subcluster}{$gene}->range;
        my $unique_length = $gene_lengths{$gene} - $multi_length;
        die "Unique length is less than zero... $gene $subcluster" if $unique_length < 0;
        say "$gene: $multi_length / $unique_length (multi vs unique length)";
    }
}

for my $gene ( sort keys %gene_multi_lengths ) {
    $uniq_lengths{$gene}-- for $gene_multi_lengths{$gene}->range;
}

# for my $subcluster ( sort keys %subcluster_counts) {
#     say "$subcluster: $subcluster_counts{$subcluster} reads";

#     for my $gene ( sort keys $ranges{$subcluster} ) {
#         $gene_multi_lengths{$gene}->addrange( $ranges{$subcluster}{$gene}->range);
#         # $uniq_lengths{$gene}-- for $ranges{$subcluster}{$gene}->range;
#         # my $unique_length = $gene_lengths{$gene} - $multi_length;
#         # die "Unique length is less than zero... $gene $subcluster" if $unique_length < 0;
#         # say "$gene: $multi_length / $unique_length (multi vs unique length)";
#     }
# }

# p %gene_multi_lengths;
p %uniq_lengths;

my $format = $ranges{'Solyc05g056060.2.1|Solyc05g056070.2.1'}{'Solyc05g056070.2.1'}->range;
say $format;

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
    my $end    = $start + $length;

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
    Solyc01g096580.2.1   2924,
    Solyc01g096590.2.1   2845,
    Solyc05g056050.2.1   5302,
    Solyc05g056060.2.1   9284,
    Solyc05g056070.2.1   10540
}
{
    Solyc01g096580.2.1   768,
    Solyc01g096590.2.1   639,
    Solyc05g056050.2.1   1030,
    Solyc05g056060.2.1   2473,
    Solyc05g056070.2.1   1112
}
{
    Solyc01g096580.2.1|Solyc01g096590.2.1                      2215,
    Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1   3422,
    Solyc05g056050.2.1|Solyc05g056070.2.1                      71,
    Solyc05g056060.2.1|Solyc05g056070.2.1                      5718
}
{
    Solyc01g096580.2.1   709,
    Solyc01g096590.2.1   630,
    Solyc05g056050.2.1   1809,
    Solyc05g056060.2.1   144,
    Solyc05g056070.2.1   1329
}
{
    Solyc01g096580.2.1   356,
    Solyc01g096590.2.1   227,
    Solyc05g056050.2.1   495,
    Solyc05g056060.2.1   1583,
    Solyc05g056070.2.1   227
}
total reads: 5493452
relevant reads: 11426
Solyc01g096580.2.1|Solyc01g096590.2.1: 2215 reads
Solyc01g096580.2.1: 412 / 356 (multi vs unique length)
Solyc01g096590.2.1: 412 / 227 (multi vs unique length)
Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1: 3422 reads
Solyc05g056050.2.1: 534 / 496 (multi vs unique length)
Solyc05g056060.2.1: 534 / 1939 (multi vs unique length)
Solyc05g056070.2.1: 534 / 578 (multi vs unique length)
Solyc05g056050.2.1|Solyc05g056070.2.1: 71 reads
Solyc05g056050.2.1: 46 / 984 (multi vs unique length)
Solyc05g056070.2.1: 46 / 1066 (multi vs unique length)
Solyc05g056060.2.1|Solyc05g056070.2.1: 5718 reads
Solyc05g056060.2.1: 833 / 1640 (multi vs unique length)
Solyc05g056070.2.1: 830 / 282 (multi vs unique length)
6..346,364..641,680..890
[Finished in 60.9s]
