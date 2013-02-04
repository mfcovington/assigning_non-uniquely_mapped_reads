#!/usr/bin/env perl
# FILE_NAME.pl
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
# exit;

# say $_->range for @{ $ranges{"Solyc01g096580.2.1|Solyc01g096590.2.1"} };
# for my $range_ref ( @{ $ranges{"Solyc01g096580.2.1|Solyc01g096590.2.1"} } ) {
#     my %range = %$range_ref;
#     # addrange("300..350")
#     say $range{$_}->addrange("20..25") for keys $range_ref;
#     say $range{$_}->range for keys $range_ref;
# }

# for my $range_ref ( @{ $ranges{"Solyc01g096580.2.1|Solyc01g096590.2.1"} } ) {
#     my %range = %$range_ref;
#     # addrange("300..350")
#     say $range{$_}->addrange("50..55") for keys $range_ref;
#     say $range{$_}->range for keys $range_ref;
# }


# say keys $_ for @{ $ranges{"Solyc01g096580.2.1|Solyc01g096590.2.1"} };


# my $format =  $range->range;

# p %ranges;
# exit;

# for a cluster, read in seqreads. if a gene matches, consider read. increment uniq_count cluster, if applicable, else deal with multi_read
# deal w/ multi_read = extract subcluster ID to use as hash key for hash of arrays (each element in array is read or at least relevant info of read)
my $sam_filename = "1.1.2_rep1_bwa0.6.2.100.sam";
open my $sam_fh, "<", $sam_filename;

my $gene_regex = join "|", map { quotemeta } keys %gene_set;
# say "regex: $gene_regex";
my %counts;
my $count;
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


    my ( $subcluster, $best_count ) = best_hits($_);
    my @best_hits = split /\|/, $subcluster;

    for my $gene ( keys %gene_set ){
        $counts{$gene}++ if /\Q$gene\E/;
    }
    $count++;

    push my @positions, first_pos($_);
    push @positions, other_pos( $_, $best_count ) if $best_count > 1;
    die "Numbers of hits and positions don't match"
      unless scalar @best_hits == scalar @positions;

    # say "POS after others:  @positions";
    # say $subcluster;
    my %hash;
    @hash{@best_hits} = @positions;

    # for ( @positions) {
      # p  $ranges{$subcluster};
    # }
    # say $_ for  keys $ranges{$subcluster};
    # p $ranges{$subcluster};
    # exit;
    next unless $best_count > 1;
    my $index = 0;
    for ( sort keys $ranges{$subcluster}  ) {

        # say $range{$_}->addrange("50..55") for keys $range_ref;
        # say $range{$_}->range for keys $range_ref;
        # my $index = 0;
        # say "v", keys %range;
        # say "b", $range{"Solyc05g056060.2.1"}->range;
        # for (@best_hits) {
            # say "$_ $positions[$index] ";
            # p $range{$_};
            # say $range{$_};
            # $range{$_}->addrange( $positions[$index] );
        # }
        # $range{@best_hits}->addrange(@positions);
        # say "b", $range{$_}->range for keys %range;
        # $$range_ref{$_}->addrange($positions[0]) for keys $range_ref;
        # $range{$_}->addrange($positions[0]) for keys %range;
        # say "dfdf", $_;
        # say $positions[$index];
        $ranges{$subcluster}{$_}->addrange($positions[$index]);
        # exit;
        # my $index = 0;
        # for ( sort keys $range_ref) {
        #     $$range_ref{$_}->addrange($positions[$index]);
        #     say $positions[$index];
        #     say $index;
        #     $index++;
        #     say "sdsf";
        # }
        $index++;

        # $range{'Solyc05g056060.2.1'}->addrange("60..65");
        # say "b", $range{$_}->range for keys %range;
    }
    # use Data::Dumper;
    # say Dumper $ranges{$subcluster};
    # p $ranges{$subcluster};
# exit;

    # just for monitoring/testing
    if ($count == 50) {
        p %counts;
        say "best for $count: $subcluster";
        say "total: $total_count";
        # p %gene_lengths;
        # p %hash;
        p %ranges;
        exit;
    }

}
say "total reads: $total_count";
say "relevant reads: $count";
p %counts;

my $multi_range = Number::Range->new(); # convert this into a hash w/ keys == geneIDs in cluster
my $multi_count;

# loop through reads for sub-cluster
    # extract start position, length, and strand for each gene in subcluster
    # my $start # make hash
    # my $end # make hash
    # add ranges
    # $multi_range->addrange("$start..$end");
    # $multi_count++;


sub best_hits {
    # adapted from harvest_gene_ids.pl
    my $read = shift;
    my ($best_hits) = $read =~ m|X0:i:(\d+)|;
    my @genes = $read =~ m|(Solyc\d{2}g\d{6}\.\d\.\d)|g;

    my $delimiter = "|";
    my $max_best  = 7;
    my @bests = sort @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
    my $best_hits_string = join $delimiter, @bests;
    my $best_hits_count = scalar @bests;
    return ( $best_hits_string, $best_hits_count );
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



my $range = Number::Range->new(""); # convert this into a hash w/ keys == geneIDs in cluster

my $format =  $range->range;
say $format;

$range->addrange("300..350");
$format =  $range->range;
say $format;

$range->addrange("349..355");
$format =  $range->range;
say $format;

$range->addrange("600..450");
$format =  $range->range;
say $format;

# say length($_) for split /,/, $range->range;

my $i;
$i++ for ($range->range);
say $i;

__END__
{
    281    [
        [0] "Solyc01g096580.2.1|Solyc01g096590.2.1"
    ],
    1001   [
        [0] "Solyc05g056060.2.1|Solyc05g056070.2.1",
        [1] "Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1",
        [2] "Solyc05g056050.2.1|Solyc05g056070.2.1"
    ]
}
{
    Solyc01g096580.2.1   1,
    Solyc01g096590.2.1   1,
    Solyc05g056050.2.1   2,
    Solyc05g056060.2.1   2,
    Solyc05g056070.2.1   3
}

012345678910202122232425

012345678910202122232425

012345678910202122232425505152535455

012345678910202122232425505152535455
[Finished in 0.2s]
