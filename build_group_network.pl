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
while (%gene_groups) {
    my ($group) = keys %gene_groups;
    my ($seed_gene) = split /\|/, $group;
    my %checked;
    my %cluster;
    while ( defined $seed_gene ) {
        for my $genes ( keys %gene_groups ) {
            next unless $genes =~ m|$seed_gene|;
            my $counts = $gene_groups{$genes};
            $cluster{$_} += $counts for split /\|/, $genes;
            delete $gene_groups{$genes};
        }
        $checked{$seed_gene}++;
        my @remainder = grep { not $checked{$_} } keys %cluster;
        $seed_gene = shift @remainder;
    }
    $cluster_sizes{ scalar keys %cluster }++;
}

exit;
