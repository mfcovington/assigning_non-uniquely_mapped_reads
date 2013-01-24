#!/usr/bin/env perl
# FILE_NAME.pl
# Mike Covington
# created: 2013-01-23
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use List::Util 'min';

my ( $sam_filename, $count_summary, $gene_summary );
my $max_best   = 7;
# my $max_subopt = 0;
my $options    = GetOptions(
    "sam_filename=s" => \$sam_filename,
    "count_summary"  => \$count_summary,
    "gene_summary"   => \$gene_summary,
    "max_best=i"     => \$max_best,
    # "max_subopt=i"   => \$max_subopt,
);

open my $sam_fh, "<", $sam_filename;

my %best_count;
my %subopt_count;
my $unmapped_count;
my $read_count;

for (<$sam_fh>) {
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
    print "$read_id";
    print "\t$_"
      for @genes[ 0 .. min( $#genes, $best_hits - 1, $max_best - 1 ) ];
    say "\n";

}

if ($count_summary) {
    say "# of reads: $read_count";
    say "# of unmapped: " . ( $unmapped_count // 0 );

    say "\n# of best hits";
    say "$_:\t$best_count{$_}" for sort { $a <=> $b } keys %best_count;

    say "\n# of suboptimal hits";
    say "$_:\t$subopt_count{$_}" for sort { $a <=> $b } keys %subopt_count;
}

__END__
Sample of what data look like:

HWI-ST611_0181:6:1101:6600:2218#0       0       Solyc12g056220.1.1      265     23      44M     *       0       0       CAAGGAGTTGCTTGGGCTTTTGGTGGTATGATCTTTGCTTTGGT    eeeefeegffghhfffihiihag^ffefcdacehfhiiiihifa    XT:A:U  NM:i:0  X0:i:1  X1:i:1  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc08g081190.2.1,+415,44M,1;
HWI-ST611_0181:6:1101:7891:2170#0       16      Solyc05g056070.2.1      803     0       44M     *       0       0       ATTAAATTGTACTATTTATCATATAGTATATTAGCCAAAAACAC    iiiiiiiihhhhiihihiiiihiiihgfhhiiiiiiigggggce    XT:A:R  NM:i:0  X0:i:2  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc05g056060.2.1,-2387,44M,0;
HWI-ST611_0181:6:1101:9025:2242#0       0       Solyc02g085940.2.1      77      0       44M     *       0       0       TTATGGTTCTCACGGTACACAAATCCGTGCTCAGTCTCGAATTC    degfggbefdfhfhh`bdghhhhhhhfbegffhhfghdghhbgg    XT:A:R  NM:i:0  X0:i:2  X1:i:0  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc02g085950.2.1,-335,44M,0;
HWI-ST611_0181:6:1101:10019:2156#0      16      Solyc02g071030.1.1      460     0       44M     *       0       0       GCACAAAGCATCTTGGCCATCTGGGCTTGCCAAGTTGTGTTGAT    a^X^gchfgfaaahhhf_[fb_egfdab_dfffgfddggecac^    XT:A:R  NM:i:0  X0:i:2  X1:i:3  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc02g071020.2.1,+373,44M,0;Solyc02g071010.1.1,-460,44M,1;Solyc02g070980.1.1,-460,44M,1;Solyc02g070990.1.1,-460,44M,1;
HWI-ST611_0181:6:1101:10539:2186#0      16      Solyc02g071020.2.1      575     0       44M     *       0       0       GGATCACCTCAAGTTCACGGTTCTTGGCAAAAGTTTCAGGGTCT    ihihgfhhhiiiiiiihiiiiiihiiihiihiiiiiigggggee    XT:A:R  NM:i:0  X0:i:2  X1:i:1  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc02g071030.1.1,+258,44M,0;Solyc02g070990.1.1,+258,44M,1;
HWI-ST611_0181:6:1101:11639:2219#0      16      Solyc04g058150.2.1      239     23      44M     *       0       0       AGCAGCTGGAGAAGGATGCAAATGTGGATCAAACTGCACTTGTG    hffhiihhiiiiigiihhiiihfiihhihiiiihhhfgggggee    XT:A:U  NM:i:0  X0:i:1  X1:i:1  XM:i:0  XO:i:0  XG:i:0  MD:Z:44 XA:Z:Solyc04g058100.2.1,-240,44M,1;
