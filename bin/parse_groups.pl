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
    chomp;
    my ( $read_id, $genes ) = split;
    $gene_groups{$genes}++;
}
close $genes_fh;

my $id = fileparse( $genes_filename, ".genes" );
open my $groups_fh, ">", "$id.groups";
say $groups_fh "$gene_groups{$_}\t$_"
  for sort { $gene_groups{$b} <=> $gene_groups{$a} } keys %gene_groups;
close $groups_fh;

exit;

__END__
./parse_groups.pl sample_id.genes
