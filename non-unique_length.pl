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
use List::Util 'min';
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
p %clusters;

# build set of genes in ALL clusters
my %gene_set;
for my $cluster_id (keys %clusters) {
    for my $subcluster ( @{ $clusters{$cluster_id} } ) {
        $gene_set{$_}++ for split /\|/, $subcluster;
    }
}
p %gene_set;

# for a cluster, read in seqreads. if a gene matches, consider read. increment uniq_count cluster, if applicable, else deal with multi_read
# deal w/ multi_read = extract subcluster ID to use as hash key for hash of arrays (each element in array is read or at least relevant info of read)
my $sam_filename = "1.1.2_rep1_bwa0.6.2.100.sam";
open my $sam_fh, "<", $sam_filename;

my $gene_regex = join "|", map { quotemeta } keys %gene_set;
# say "regex: $gene_regex";
my %counts;
my $count;
my $total_count;
while (<$sam_fh>){
    next if /^@/;
    $total_count++;
    next unless /$gene_regex/;
    my @best_hits = best_hits($_);

    for my $gene ( keys %gene_set ){
        $counts{$gene}++ if /\Q$gene\E/;
    }
    $count++;

    # just for monitoring/testing
    # if ($count == 5) {
    #     p %counts;
    #     say "best for $count: @best_hits";
    #     say "total: $total_count";
    #     exit;
    # }

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
    return join $delimiter,
      sort @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
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
{
    Solyc01g096580.2.1   3323,
    Solyc01g096590.2.1   3245,
    Solyc05g056050.2.1   7589,
    Solyc05g056060.2.1   9854,
    Solyc05g056070.2.1   11071
}
total reads: 5458725
relevant reads: 16047

300..350
300..355
300..355,450..600
207
[Finished in 5.1s]
