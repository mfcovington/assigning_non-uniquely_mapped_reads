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
use Data::Printer;

my $genes_filename = $ARGV[0];
open my $genes_fh, "<", $genes_filename;

my %gene_groups;
while (<$genes_fh>) {
    chomp;
    my ( $read_id, $genes ) = split;
    $gene_groups{$genes}++;
}
close $genes_fh;

my %cluster_sizes;
# put the following into another loop to loop through each cluster
while ( %gene_groups ) {
    # my ( $seed_gene ) = split /\|/, ${ keys %gene_groups }[0];
    my ( $group ) = keys %gene_groups;
    my ( $seed_gene ) = split /\|/, $group;
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
    # p %cluster if scalar keys %cluster == 1; # some clusters have only one gene, because single reads map to multiple positions within the same gene
}

# p %gene_groups;
# p %checked;
# p %cluster;
p %cluster_sizes;

__END__
some clusters have only one gene, because single reads map to multiple positions within the same gene
    1    71,
    2    1117,
    3    234,
    4    88,
    5    40,
    6    17,
    7    22,
    8    9,
    9    7,
    10   3,
    11   2,
    12   2,
    13   1,
    19   1,
    20   3,
    21   1,
    26   1,
    27   1,
    29   1,
    40   1,
    41   1,
    93   1



