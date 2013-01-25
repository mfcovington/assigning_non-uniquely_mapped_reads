#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-01-24
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use File::Basename;

my $genes_filename = $ARGV[0];
open my $genes_fh, "<", $genes_filename;
my %gene_groups;
while (<$genes_fh>) {
    my ( $read_id, $genes ) = split;
    $gene_groups{$genes}++;
}
close $genes_fh;

my %cluster_sizes;
my $cluster_number;
my $id = fileparse( $genes_filename, ".genes" );
open my $clusters_fh, ">", "$id.clusters";
while (%gene_groups) {
    my ($group) = sort keys %gene_groups;
    my ($seed_gene) = split /\|/, $group;
    my %checked;
    my %cluster;
    while ( defined $seed_gene ) {
        for my $genes ( keys %gene_groups ) {
            next unless $genes =~ m|$seed_gene|;
            my $counts = $gene_groups{$genes};
            delete $gene_groups{$genes};
            my %gene_list;
            $gene_list{$_} = 1 for split /\|/, $genes;
            next if scalar keys %gene_list == 1;
            $cluster{$_} += $counts for keys %gene_list;
        }
        $checked{$seed_gene}++;
        my @remainder = grep { not $checked{$_} } keys %cluster;
        $seed_gene = shift @remainder;
    }
    if ( scalar keys %cluster > 1 ){
        $cluster_number++;
        say $clusters_fh "CLUSTER $cluster_number:";
        say $clusters_fh "$_\t$cluster{$_}" for sort keys %cluster;
        $cluster_sizes{ scalar keys %cluster }++;
    }
}
close $clusters_fh;

open my $sizes_fh, ">", "$id.sizes";
say $sizes_fh "Size\tNumber of clusters";
say $sizes_fh "$_\t$cluster_sizes{$_}"
  for sort { $a <=> $b } keys %cluster_sizes;
close $sizes_fh;

exit;
