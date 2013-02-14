#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-02-12
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Parallel::ForkManager;
use FindBin;
use lib "$FindBin::Bin/../lib";
use Map::Redundant;

my %clusters;
my $cluster_id  = 1001;
my @subclusters = (
    "Solyc05g056060.2.1|Solyc05g056070.2.1",
    "Solyc05g056050.2.1|Solyc05g056060.2.1|Solyc05g056070.2.1",
    "Solyc05g056050.2.1|Solyc05g056070.2.1"
);
$clusters{$cluster_id} = [@subclusters];
$clusters{281} = ["Solyc01g096580.2.1|Solyc01g096590.2.1"];
# use Data::Printer;
# p %clusters;
# exit;
my $sam_dir = "../";
my $out_dir = "../test_out/";

# my @sam_files = "1.1.2_rep1_bwa0.6.2.100.sam";
my @sam_files = "CLUSTER1001.sam";

# my @sam_files = qw(1.1.2_rep1_bwa0.6.2.100.sam 1.1.2_rep1_bwa0.6.2.100.sam 1.1.2_rep1_bwa0.6.2.100.sam);
# my $threads = 3;
# my $pm = new Parallel::ForkManager($threads);

# for my $sam_file (@sam_files){
#     $pm->start and next;
#     # my %clusters = clusters("$sam_dir/$sam_file");
#     coverage("$sam_dir/$sam_file", %clusters);
#     $pm->finish;
# }
# $pm->wait_all_children;

my $sam_file = $sam_dir . "/" . $sam_files[0];
# my $sam_file = $sam_files[0];

my $assign = new Map::Redundant(
    'sam_file'        => $sam_file,
    'out_dir'         => $out_dir,
    'gene_summary'    => 1,
    'count_summary'   => 1,
    'cluster_summary' => 1,
    'size_summary'    => 1,
);

$assign->identify_subclusters;
$assign->summarize_subclusters;
$assign->build_clusters;

# use Data::Printer;
# p $assign->subclusters_hash;
# p $assign->clusters_hash;

# $assign->clusters_hash(\%clusters);
$assign->calculate_coverage();


# say $assign->sam_file;
# # my %hash = %{ $assign->clusters_hash }
# # say "$_: @{ $assign->clusters_hash->{$_} }" for keys  $assign->clusters_hash ;
# use Data::Printer;
# p $assign->clusters_hash;

__END__
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
Solyc05g056050.2.1: 1809 (unique reads)
Solyc05g056050.2.1: 924 (unique length)
Solyc05g056060.2.1: 144 (unique reads)
Solyc05g056060.2.1: 1111 (unique length)
Solyc05g056070.2.1: 1329 (unique reads)
Solyc05g056070.2.1: 431 (unique length)
[Finished in 55.5s]

