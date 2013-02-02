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
local $SIG{__WARN__} = sub { warn $_[0] unless $_[0] =~ m{Use of uninitialized value \$previous in string at /Users/mfc/perl5/perlbrew/perls/perl-5.16.1/lib/site_perl/5.16.1/Number/Range.pm line 215.|\d+ already in range|\d+ not in range or already removed|previous|\d+ is > \d+|\d+:\d+ is pointless}};


# for a cluster, read in seqreads. if a gene matches, consider read. increment uniq_count cluster, if applicable, else deal with multi_read
# deal w/ multi_read = extract subcluster ID to use as hash key for hash of arrays (each element in array is read or at least relevant info of read)


my $multi_range = Number::Range->new(); # convert this into a hash w/ keys == geneIDs in cluster
my $multi_count;

# loop through reads for sub-cluster
    # extract start position, length, and strand for each gene in subcluster
    # my $start # make hash
    # my $end # make hash
    # add ranges
    # $multi_range->addrange("$start..$end");
    # $multi_count++;






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
annoying warnings:

Use of uninitialized value $previous in string at /Users/mfc/perl5/perlbrew/perls/perl-5.16.1/lib/site_perl/5.16.1/Number/Range.pm line 215.
349 already in range at /Users/mfc/git.repos/assigning_non-uniquely_mapped_reads/non-unique_length.pl line 44.
350 already in range at /Users/mfc/git.repos/assigning_non-uniquely_mapped_reads/non-unique_length.pl line 44.
600 is > 450 at /Users/mfc/git.repos/assigning_non-uniquely_mapped_reads/non-unique_length.pl line 48.


actual output:

300..350
300..355
300..355,450..600
207


