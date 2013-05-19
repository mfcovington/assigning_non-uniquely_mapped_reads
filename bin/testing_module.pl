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

my $sam_dir = "../sample_files/";
my $out_dir = "../test_out/";
# my $sam_dir = "../eXpress/bwa/";
# my $out_dir = "../eXpress/bwa/";

# my @sam_files = "1.1.2_rep1_bwa0.6.2.100.sam";
# my @sam_files = "CLUSTER281+1001.sam";
# my @sam_files = "CLUSTER1001.sam";
# my @sam_files = "H1_000_000.sam";
my @sam_files = "H100_000.sam";
# my @sam_files = "IL1.1.2_1.filtered.sam";

my $sam_file = $sam_dir . "/" . $sam_files[0];

my $assign = new Map::Redundant(
    'sam_file'         => $sam_file,
    'out_dir'          => $out_dir,
    'gene_summary'     => 1,
    'count_summary'    => 1,
    'cluster_summary'  => 1,
    'size_summary'     => 1,
    'coverage_summary' => 1,
    'verbose'          => 1,
);
$assign->identify_subclusters;
$assign->summarize_subclusters;
$assign->build_clusters;
$assign->calculate_coverage();

__END__
Identifying Subclusters
Summarizing Subclusters
Building Clusters
Calculating Coverage
[Finished in 57.2s]
