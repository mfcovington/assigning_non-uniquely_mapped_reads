#!/usr/bin/env perl
# Redundant.pl
# Mike Covington
# created: 2013-02-01
#
# Description:
#
package Map::Redundant;
use Moose;
use namespace::autoclean;
use autodie;
use feature 'say';
use Number::Range;
use List::Util qw(min max sum);
use File::Basename;
use File::Path 'make_path';
use IO::Handle;

STDERR->autoflush(1);
STDOUT->autoflush(1);

has 'max_best' => (
    is      => 'rw',
    isa     => 'Num',
    default => 100,
    lazy    => 1,
);

has 'out_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => '../',
    lazy    => 1,
);

has 'sam_dir' => (
    is      => 'rw',
    isa     => 'Str',
    default => '../',
    lazy    => 1,
);

has 'delimiter' => (
    is      => 'rw',
    isa     => 'Str',
    default => '|',
    lazy    => 1,
);

has 'gene_summary' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'count_summary' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'cluster_summary' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'size_summary' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'coverage_summary' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'verbose' => (
    is      => 'rw',
    isa     => 'Bool',
    default => 0,
    lazy    => 1,
);

has 'sam_file' => (
    is  => 'rw',
    isa => 'Str',
);

has 'clusters_hash' => (
    is  => 'rw',
    isa => 'HashRef',
);

has 'subclusters_hash' => (
    is  => 'rw',
    isa => 'HashRef',
);

has 'clustered_genes_hash' => (
    is  => 'rw',
    isa => 'HashRef',
);

has 'id' => (
    is      => 'rw',
    isa     => 'Str',
    default => sub {
        my $self = shift;

        my $id = fileparse( $self->sam_file, ".sam" );
        return $id;
    },
    lazy => 1,
);

sub _make_dir {
    my $self     = shift;
    my $dir_name = shift;

    ( my $filename, $dir_name ) = fileparse( $self->out_file )
      unless defined $dir_name;
    make_path($dir_name) unless -e $dir_name;
}

sub identify_subclusters {
    my $self = shift;

    my $sam_filename  = $self->sam_file;
    my $count_summary = $self->count_summary;
    my $gene_summary  = $self->gene_summary;
    my $verbose       = $self->verbose;
    my $max_best      = $self->max_best;
    my $out_dir       = $self->out_dir;
    my $delimiter     = $self->delimiter;

    say "Identifying Subclusters" if $verbose;

    $self->_make_dir($out_dir);

    my %best_count;
    my %subopt_count;
    my $unmapped_count;
    my $read_count;

    my %subclusters;

    open my $sam_fh, "<", $sam_filename;
    my $id = $self->id;
    open my $genes_fh, ">", "$out_dir/$id.genes" if $gene_summary;

    while (<$sam_fh>) {
        next if $_ =~ m|^@|;

        $read_count++ if $count_summary;

        my ( $read_id, $map_flag ) = split /\t/, $_;
        if ( $map_flag == 4 ) {
            $unmapped_count++ if $count_summary;
            next;
        }

        my ($best_hits)   = $_ =~ m|X0:i:(\d+)|;
        my ($subopt_hits) = $_ =~ m|X1:i:(\d+)|;

        if ($count_summary) {
            $best_count{$best_hits}++     if defined $best_hits;
            $subopt_count{$subopt_hits}++ if defined $subopt_hits;
        }

        next unless defined $best_hits && $best_hits > 1;
        my @genes = $_ =~ m|(Solyc\d{2}g\d{6}\.\d\.\d)|g;
        my $subcluster = join $delimiter,
          sort @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
        $subclusters{$subcluster}++;
        say $genes_fh "$read_id\t$subcluster" if $gene_summary;
    }
    close $sam_fh;
    close $genes_fh if $gene_summary;

    $self->subclusters_hash( \%subclusters );

    if ($count_summary) {
        open my $counts_fh, ">", "$out_dir/$id.counts";

        say $counts_fh "# of reads: $read_count";
        say $counts_fh "# of unmapped: " . ( $unmapped_count // 0 );

        say $counts_fh "\n# of best hits";
        say $counts_fh "$_:\t$best_count{$_}"
          for sort { $a <=> $b } keys %best_count;

        say $counts_fh "\n# of suboptimal hits";
        say $counts_fh "$_:\t$subopt_count{$_}"
          for sort { $a <=> $b } keys %subopt_count;

        close $counts_fh;
    }
}

sub summarize_subclusters {
    my $self = shift;

    my $verbose     = $self->verbose;
    my %subclusters = %{ $self->subclusters_hash };

    say "Summarizing Subclusters" if $verbose;

    my $out_dir = $self->out_dir;
    $self->_make_dir($out_dir);

    my $id = $self->id;

    open my $subclusters_fh, ">", "$out_dir/$id.subclusters";
    say $subclusters_fh "$subclusters{$_}\t$_"
      for sort { $subclusters{$b} <=> $subclusters{$a} } keys %subclusters;
    close $subclusters_fh;
}

sub build_clusters {
    my $self = shift;

    my %subclusters     = %{ $self->subclusters_hash };
    my $id              = $self->id;
    my $cluster_summary = $self->cluster_summary;
    my $size_summary    = $self->size_summary;
    my $out_dir         = $self->out_dir;
    my $verbose         = $self->verbose;

    say "Building Clusters" if $verbose;

    $self->_make_dir($out_dir);

    my %clustered_genes;
    my %cluster_sizes;
    my $cluster_number;
    open my $clusters_fh, ">", "$out_dir/$id.clusters" if $cluster_summary;
    while (%subclusters) {
        my ($group) = sort keys %subclusters;
        my ($seed_gene) = split /\|/, $group;
        my %checked;
        my %current_cluster;
        while ( defined $seed_gene ) {
            for my $genes ( keys %subclusters ) {
                next unless $genes =~ m|$seed_gene|;
                my $counts = $subclusters{$genes};
                delete $subclusters{$genes};
                my %gene_list;
                $gene_list{$_} = 1 for split /\|/, $genes;
                next if scalar keys %gene_list == 1;
                $current_cluster{$_} += $counts for keys %gene_list;
            }
            $checked{$seed_gene}++;
            my @remainder = grep { not $checked{$_} } keys %current_cluster;
            $seed_gene = shift @remainder;
        }
        if ( scalar keys %current_cluster > 1 ) {
            $cluster_number++;
            if ($cluster_summary) {
                say $clusters_fh "CLUSTER $cluster_number:";
                say $clusters_fh "$_\t$current_cluster{$_}"
                  for sort keys %current_cluster;
            }
            $cluster_sizes{ scalar keys %current_cluster }++;
            $clustered_genes{$cluster_number} = [ sort keys %current_cluster ];
        }
    }
    close $clusters_fh if $cluster_summary;

    if ($size_summary) {
        open my $sizes_fh, ">", "$out_dir/$id.sizes";
        say $sizes_fh "Size\tNumber of clusters";
        say $sizes_fh "$_\t$cluster_sizes{$_}"
          for sort { $a <=> $b } keys %cluster_sizes;
        close $sizes_fh;
    }

    %subclusters = %{ $self->subclusters_hash };

    my %clustered_subclusters;
    for my $cluster_id ( keys %clustered_genes ) {
        my $regex = join "|",
          map { quotemeta } @{ $clustered_genes{$cluster_id} };
        $clustered_subclusters{$cluster_id} = [];
        for ( sort keys %subclusters ) {
            push $clustered_subclusters{$cluster_id}, $_ if /$regex/;
        }
    }

    $self->clustered_genes_hash( \%clustered_genes );
    $self->clusters_hash( \%clustered_subclusters );
}

sub calculate_coverage {    # adapted from non-unique_length.pl
    my $self = shift;

    my $sam_filename     = $self->sam_file;
    my %clustered_genes  = %{ $self->clustered_genes_hash };
    my %clusters         = %{ $self->clusters_hash };
    my $max_best         = $self->max_best;
    my $delimiter        = $self->delimiter;
    my $verbose          = $self->verbose;
    my $id               = $self->id;
    my $coverage_summary = $self->coverage_summary;
    my $out_dir          = $self->out_dir;
    $self->_make_dir($out_dir);

    say "Calculating Coverage" if $verbose;

    # Number::Range prints unnecessary warnings; therefore, turn them off
    no warnings 'Number::Range';
    local $SIG{__WARN__} = sub {
        warn $_[0]
          unless $_[0] =~
m|Use of uninitialized value \$previous in string at .*Number/Range.pm line \d+.|;
    };

    # Number::Range stores large ranges in a way that causes problems
    # when calculating total lengths, unless max_hash_size is high enough.
    # Since there is a (mis-annotated!) gene of 219kb in my dataset:
    my $max_hash_size = 220_000;
    # my $max_hash_size = 10_000;
    # my $max_hash_size = 52;   # 52 and lower cause problems

    my %gene_set;
    my %unique_counts;
    my %ranges;
    my %subcluster_counts;
    for my $cluster_id ( keys %clusters ) {
        for my $subcluster ( @{ $clusters{$cluster_id} } ) {

            # build set of genes in ALL clusters
            $gene_set{$_}++ for split /\|/, $subcluster;

            # build unique counts hash
            $unique_counts{$_} = 0 for split /\|/, $subcluster;

            # build ranges data structure (HoHoO)
            $ranges{$subcluster}{$_} = {} for split /\|/, $subcluster;

            # build subcluster counts hash
            $subcluster_counts{$_} = 0 for @{ $clusters{$cluster_id} };
        }
    }

    # build unique_ranges hash
    my %unique_ranges;
    $unique_ranges{$_} = {} for keys %gene_set;

    # for a cluster, read in seqreads. if a gene matches, consider read. increment unique_count. cluster, if applicable, else deal with multi_read
    # deal w/ multi_read = extract subcluster ID to use as hash key for hash of arrays (each element in array is read or at least relevant info of read)
    open my $sam_fh, "<", $sam_filename;

    my $gene_regex = join "|", map { quotemeta } keys %gene_set;
    my %counts;
    my $relevant_count;
    my $total_count;
    my %gene_lengths;

    while (<$sam_fh>) {
        $total_count++;
        say "Processing read # $total_count"
          if $verbose && $total_count % 250_000 == 0;
        next unless /$gene_regex/;

        # collect gene lengths
        if (/^\@/) {
            next unless /^\@SQ/;
            m/SN:(.+)\tLN:(\d+)$/;
            $gene_lengths{$1} = $2;
            next;
        }

        my ( $subcluster, $best_count, @best_hits ) =
          best_hits( $_, $max_best, $delimiter );

        # bypass instances where genes of interest are only suboptimal hits
        next unless join( "-", @best_hits ) =~ /$gene_regex/;

        for my $gene (@best_hits) {
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
        # also populate %unique_ranges
        if ( $best_count == 1 ) {
            $unique_counts{ $best_hits[0] }++;
            my ( $start, $end ) = split /\.{2}/, $positions[0];
            add_range( $start, $end, $unique_ranges{ $best_hits[0] });
            next;
        }
        $subcluster_counts{$subcluster}++;
        $relevant_count++;

        # populate range data structure
        my $index = 0;
        for ( keys %{ $ranges{$subcluster} } ) {
            my ( $start, $end ) = split /\.{2}/, $positions[$index];
            add_range( $start, $end, $ranges{$subcluster}{$_} );
            $index++;
        }
    }

    for ( keys %unique_ranges ) {
        collapse_ranges( $unique_ranges{ $_ } );
    }

    for my $subcluster ( keys %ranges ) {
        for my $gene ( keys $ranges{$subcluster} ) {
            collapse_ranges( $ranges{$subcluster}{$gene} );
        }
    }

    die "Something may be wrong with the sam file header..."
      unless
      join( "", sort keys %gene_set ) eq join( "", sort keys %gene_lengths );

    my %gene_multi_lengths;
    my %unique_lengths;
    open my $coverage_fh, ">", "$out_dir/$id.coverage" if $coverage_summary;
    for my $cluster_id ( sort { $a <=> $b } keys %clusters ) {
        say $coverage_fh "CLUSTER: $cluster_id";
        for my $subcluster ( @{ $clusters{$cluster_id} } ) {

            # build gene_multi_lengths hash
            for ( split /\|/, $subcluster ) {
                $gene_multi_lengths{$_} = Number::Range->new();
                $gene_multi_lengths{$_}->set_max_hash_size($max_hash_size);
            }

            # populate gene_multi_lengths hash
            say $coverage_fh
              "$subcluster: $subcluster_counts{$subcluster} reads";

            for my $gene ( sort keys $ranges{$subcluster} ) {
                $gene_multi_lengths{$gene}
                  ->addrange( $ranges{$subcluster}{$gene}->range );
                my $multi_length;
                $multi_length = range_length( $ranges{$subcluster}{$gene} );
                say $coverage_fh "$gene: $multi_length (non-unique length)";
            }
        }

        # build/populate unique_lengths hash
        say $coverage_fh "-----";
        for my $gene ( @{ $clustered_genes{$cluster_id} } ) {
            $unique_lengths{$gene} = 0;
            $unique_lengths{$gene} = range_length( $unique_ranges{$gene} );
            say $coverage_fh "$gene: $unique_counts{$gene} (unique reads)";
            say $coverage_fh "$gene: $unique_lengths{$gene} (unique length)";
        }
        say $coverage_fh "";
    }
    close $coverage_fh;

    # compare max gene length to max_hash_size parameter from Number::Range
    my @lengths;
    push @lengths, $gene_lengths{$_} for keys %gene_lengths;
    my $max_length = max @lengths;
    say <<WARNING if $max_length > $max_hash_size;
The maximum gene length ($max_length) is greater than
the maximum hash size ($max_hash_size).
You should probably increase \$max_hash_size, just in case.
WARNING
}

sub best_hits {    # adapted from harvest_gene_ids.pl
    my ( $read, $max_best, $delimiter ) = @_;
    my ($best_hits) = $read =~ m|X0:i:(\d+)|;
    my @genes = $read =~ m|(Solyc\d{2}g\d{6}\.\d\.\d)|g;

    my @bests = @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
    my $best_hits_string = join $delimiter, sort @bests;
    my $best_hits_count = scalar @bests;
    return ( $best_hits_string, $best_hits_count, @bests );
}

sub first_pos {
    my $read = shift;
    my ( $start, $cigar ) = ( split /\t/, $read )[ 3, 5 ];
    my $end = calc_end( $start, $cigar );

    return ("$start..$end");
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

    return (@other_positions);
}

sub calc_end {
    my ( $start, $cigar ) = @_;

    my $matches    = sum( $cigar =~ m|(\d+)M|g );
    my $insertions = sum( $cigar =~ m|(\d+)I|g ) || 0;
    my $deletions  = sum( $cigar =~ m|(\d+)D|g ) || 0;

    my $length = $matches + $deletions - $insertions;
    my $end    = $start + $length - 1;

    return $end;
}

sub add_range {
    my ( $start, $end, $range_ref ) = @_;
    if ( exists $range_ref->{$start} ) {
        $range_ref->{$start} = max( $end, $range_ref->{$start} );
    }
    else {
        $range_ref->{$start} = $end;
    }
}

sub collapse_ranges {
    my $range_ref = shift;
    return if scalar keys %$range_ref == 0;

    my @cur_interval;
    my @result;
    my %temp_ranges;

    for my $start ( sort { $a <=> $b } keys %$range_ref ) {
        unless (@cur_interval) {
            @cur_interval = ( $start, $range_ref->{$start} );
            next;
        }
        my ( $cstart, $cend ) = @cur_interval;
        if ( $start <= $cend + 1 ) {
            @cur_interval = ( $cstart, max( $range_ref->{$start}, $cend ) );
        }
        else {
            push @result, @cur_interval;
            $temp_ranges{ $cur_interval[0] } = $cur_interval[1];
            @cur_interval = ( $start, $range_ref->{$start} );
        }
    }
    push @result, @cur_interval;
    $temp_ranges{ $cur_interval[0] } = $cur_interval[1];
    %$range_ref = %temp_ranges;
}

sub range_length {
    my $range_ref = shift;
    my $length = 0;
    for ( keys %$range_ref ) {
        $length += $range_ref->{$_} - $_ + 1 if $range_ref->{$_};
    }
    return $length;
}

__PACKAGE__->meta->make_immutable;
