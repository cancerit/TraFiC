package Sanger::CGP::TraFiC::Cluster;

##########LICENCE##########
# Copyright (c) 2015 Genome Research Ltd.
#
# Author: Cancer Genome Project cgpit@sanger.ac.uk
#
# This file is part of TraFiC.
#
# TraFiC is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
# details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained within
# this distribution should be interpreted as being equivalent to a list of years
# including the first and last year specified and all consecutive years between
# them. For example, a copyright statement that reads 'Copyright (c) 2005, 2007-
# 2009, 2011-2012' should be interpreted as being identical to a statement that
# reads 'Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012' and a copyright
# statement that reads "Copyright (c) 2005-2012' should be interpreted as being
# identical to a statement that reads 'Copyright (c) 2005, 2006, 2007, 2008,
# 2009, 2010, 2011, 2012'."
##########LICENCE##########


use strict;
use autodie qw(:all);
use English qw( -no_match_vars );
use warnings FATAL => 'all';
use Carp qw(croak);
use Const::Fast qw(const);
use File::Which qw(which);

use Sanger::CGP::TraFiC;

const my $TAB => "\t";

const my @CLASSES => qw(SINGLE_END INTER_CHROM ABERRANT);
const my @HEADER_ROOT => qw(CHR FAMILY);
const my @HEADER_MAIN => qw(L_POS R_POS TOTAL_READS SINGLE_END_COUNT INTER_CHROM_COUNT ABERRANT_COUNT SINGLE_END_READS INTER_CHROM_READS ABERRANT_READS);

sub new {
    my ($class) = @_;
    croak "Unable to find standard unix 'sort' in path." unless(defined which('sort'));
    my $self = { };
    bless $self, $class;
    $self->{'min_reads'} = 5;
    return $self;
}

sub reciprocal_clusters {
  my $self = shift;
  # results are small so load into hash structures to simplify searching by family/chr
  my $pos = $self->clustered_struct('+');
  my $neg = $self->clustered_struct('-');
  $self->reciprocal_hits($pos, $neg);
  return 1;
}

sub reciprocal_hits {
  my ($self, $pos, $neg) = @_;
  # keyed as family->chr
  # if family in pos not found in neg then can skip the whole block
  # if no chr for a family then you can skip that too...
  my @reciprocal_clusters;
  my @families = sort keys %{$pos};
  for my $family(@families) {
    unless(exists $neg->{$family}) {
      delete $pos->{$family};
      next;
    }
    my @chrs = sort keys %{$pos->{$family}};
    for my $chr(@chrs) {
      unless(exists $neg->{$family}->{$chr}) {
        delete $pos->{$family}->{$chr};
        next;
      }
      push @reciprocal_clusters, @{$self->compare_hits($chr, $family, $pos->{$family}->{$chr}, $neg->{$family}->{$chr})};
    }
  }
  $self->output_reciprocal(\@reciprocal_clusters);
  return 1;
}

sub output_reciprocal {
  my ($self, $clusters) = @_;

  # output clusters
  open my $RECIP, '>', "$self->{reciprocal_clusters}.raw";
  for my $cluster(@{$clusters}) {
    my ($chr, $family, $pos, $neg) = @{$cluster};
    print $RECIP (join $TAB,
                      $chr, $family, # shared
                      $pos->{'l_pos'}, $pos->{'r_pos'}, @{$pos->{'detail'}}, # pos cluster info
                      $neg->{'l_pos'}, $neg->{'r_pos'}, @{$neg->{'detail'}}) # neg cluster info
                  ,$RS;
  }
  close $RECIP;

  # sort clusters
  # safest to use unix command as will be fast and lean
  my $sort_command = sprintf "sort -k 1,1 -k 3,4n -k 2,2 -o %s %s",
                              "$self->{reciprocal_clusters}.sort", # output file
                              "$self->{reciprocal_clusters}.raw"; # input file
  croak "Failed to sort file using command: $sort_command" unless(system($sort_command) == 0);

  # write header
  open $RECIP, '>', $self->{'reciprocal_clusters'};
  print $RECIP '#',(join $TAB, @HEADER_ROOT);
  for(@HEADER_MAIN) {
    print $RECIP $TAB,'P_',$_;
  }
  for(@HEADER_MAIN) {
    print $RECIP $TAB,'N_',$_;
  }
  print $RECIP $RS;

  # append the sorted results
  open my $SORT, '<', "$self->{reciprocal_clusters}.sort";
  while(<$SORT>) { print $RECIP $_; }
  close $SORT;
  close $RECIP;

  # cleanup intermediate files
  unlink "$self->{reciprocal_clusters}.raw";
  unlink "$self->{reciprocal_clusters}.sort";
  return 1;
}

sub compare_hits {
  my ($self, $chr, $family, $pos_hits, $neg_hits) = @_;
  my $neg_tot = scalar @{$neg_hits};
  my $neg_i = 0;
  my @reciprocal;
  PRIMARY: for my $pos(@{$pos_hits}) {
    # if pos l_pos is > neg r_pos+200
    while($neg_i < $neg_tot) {
      if($neg_hits->[$neg_i]->{'l_pos'} < $pos->{'r_pos'}) {
        $neg_i++;
        next;
      }
      if($neg_hits->[$neg_i]->{'r_pos'}+200 < $pos->{'l_pos'}) {
        # neg data is prior to pos data so catch up
        $neg_i++;
        next;
      }
      if($neg_hits->[$neg_i]->{'l_pos'}-200 > $pos->{'r_pos'}) {
        # neg is ahead of pos so catchup pos and re-check this entry
        next PRIMARY;
      }
      push @reciprocal, [$chr, $family, $pos, $neg_hits->[$neg_i]];
      $neg_i++;
    }
  }
  return \@reciprocal;
}

sub clustered_struct {
  my ($self, $strand) = @_;
  my %clusters;
  open my $CLUSTERS, '<', $self->{$strand.'_clusters'};
  while (my $cluster = <$CLUSTERS>) {
    next if(substr($cluster, 0,1) eq '#');
    chomp $cluster;
    my ($chr, $family, $l_pos, $r_pos, @detail) = split /\t/, $cluster;
    push @{$clusters{$family}{$chr}}, { 'l_pos' => int $l_pos,
                                        'r_pos' => int $r_pos,
                                        'detail' => \@detail,
                                        };
  }
  close $CLUSTERS;
  return \%clusters;
}

sub cluster_hits {
  my ($self) = @_;
  $self->cluster_file('+');
  $self->cluster_file('-');
  return 1;
}

sub cluster_file {
  my ($self, $strand) = @_;
  open my $HITS, '<', $self->{$strand.'_sorted_hits'};
  my ($p_chr, $p_pos, $p_fam) = (q{}, 0, q{});
  my $cluster = [];
  open my $CLUSTERED, '>', $self->{$strand.'_clusters'};
  print $CLUSTERED '#',(join $TAB, @HEADER_ROOT, @HEADER_MAIN),$RS;
  while (my $line = <$HITS>) {
    chomp $line;
    my ($rn, $chr, $pos, $fam, $class) = split $TAB, $line;
    next if($fam eq 'Alu' && $class ne 'INTER_CHROM');

    if($p_chr eq $chr && $p_fam eq $fam && ($pos - $p_pos) <= 200) {
      push @{$cluster}, [$rn, $chr, $pos, $class];
    }
    else {
      $cluster = $self->finalise_cluster($cluster, $p_fam, $strand, $CLUSTERED);
      push @{$cluster}, [$rn, $chr, $pos, $class];
    }
    ($p_chr, $p_pos, $p_fam) = ($chr, $pos, $fam);
  }
  $self->finalise_cluster($cluster, $p_fam, $strand, $CLUSTERED);
  close $CLUSTERED;
  close $HITS;
  return 1;
}

sub finalise_cluster {
  my ($self, $cluster, $fam, $strand, $CLUSTERED) = @_;
  if(scalar @{$cluster} >= $self->{'min_reads'}) {
    my %classes = map { $_ => [] } @CLASSES;
    for my $item(@{$cluster}) {
      push @{$classes{$item->[3]}}, $item->[0];
    }
    # chr, family, l_pos, r_pos, read_count_*, read_names_*
    print $CLUSTERED (join $TAB, $cluster->[0]->[1],
                                $fam,
                                $cluster->[0]->[2],
                                $cluster->[-1]->[2],
                                scalar @{$classes{'SINGLE_END'}} + scalar @{$classes{'INTER_CHROM'}} + scalar @{$classes{'ABERRANT'}},
                                scalar @{$classes{'SINGLE_END'}},
                                scalar @{$classes{'INTER_CHROM'}},
                                scalar @{$classes{'ABERRANT'}},
                                ${_format_readnames($classes{'SINGLE_END'})},
                                ${_format_readnames($classes{'INTER_CHROM'})},
                                ${_format_readnames($classes{'ABERRANT'})},
                      ),$RS;
  }
  $cluster = [];
  return $cluster;
}

sub _format_readnames {
  my $reads = shift;
  my $ret_val = q{.};
  if(defined $reads && @{$reads} > 0) {
    $ret_val = join ',', @{$reads};
  }
  return \$ret_val
}

sub set_min_reads {
  my ($self, $min_reads) = @_;
  $self->{'min_reads'} = $min_reads || 5;
  return 1;
}

sub set_input_output {
  my ($self, $root_path) = @_;
  croak "$root_path does not exist." unless(-e $root_path);
  croak "$root_path is not a directory." unless(-d $root_path);
  croak "$root_path cannot be written to." unless(-w $root_path);

  $self->{'+_sorted_hits'} = "$root_path/pos_hits.txt";
  $self->{'-_sorted_hits'} = "$root_path/neg_hits.txt";

  $self->{'+_clusters'} = "$root_path/pos_clusters.txt";
  $self->{'-_clusters'} = "$root_path/neg_clusters.txt";
  $self->{'reciprocal_clusters'} = "$root_path/reciprocal_clusters.txt";

  return 1;
}


1;

__END__

=head1 NAME

Sanger::CGP::TraFiC::Cluster

=head1 SYNOPSIS

 my $cluster = Sanger::CGP::TraFiC::Cluster->new;

 $cluster->set_input_output($output_folder);
    # which should already contain data from Formatter.pm

 $cluster->set_min_reads($minimum_reads_required);
 $cluster->cluster_hits;
 $cluster->reciprocal_clusters;

=head1 GENERAL

Take the merged repeat-masker *.out and fasta data generated by Formatter.pm and generate clusters
of nearby/overlapping events.

Generates three files:

=over

=item pos_clusters.txt

Data for reads anchored on positive strand that are able to be clustered with greater or equal to
the value of I<min_reads>.

=item neg_clusters.txt

Data for reads anchored on negative strand that are able to be clustered with greater or equal to
the value of I<min_reads>.

=item reciprocal_clusters.txt

Where entries in pos/neg_clusters.txt support each other they are additionally output here.  The
pos_cluster is always the primary entry.

=back

=head2 Format of output

The output files all follow the same format

=head3 Core fields

These are present at the start of each line in every file.

=over

=item CHR

Chromosome or sequence identifier.

=item FAMILY

Masking identified reads as hitting repeat of this type.

=back

=head3 Data fields

These are present in the pos/neg_clusters.txt files as described here.

In the reciprocal_clusters.txt file these are presented twice, the positive cluster first (prefixed B<P_>)
and then the negative cluster data (prefixed B<N_>).

=over

=item L_POS

Left most position of cluster.

=item R_POS

Left most position of cluster.

=item TOTAL_READS

Total number of reads that support this cluster

=item SINGLE_END_COUNT

Number of single end mapped reads that contribute to this cluster.

=item INTER_CHROM_COUNT

Number of inter-chromosomal mapped reads that contribute to this cluster.

=item ABERRANT_COUNT

Number of aberrantly paired reads that contribute to this cluster.

=item SINGLE_END_READS

Names of reads contributing to I<SINGLE_END_COUNT>.

=item INTER_CHROM_READS

Names of reads contributing to I<INTER_CHROM_COUNT>.

=item ABERRANT_READS

Names of reads contributing to I<ABERRANT_COUNT>.

=back

=head1 METHODS

=head2 Constructor/configuration

=head3 new

No options, sets up object with default values.

=head3 set_input_output

=over

=item Description

Define the path to the input files (generated by L<Sanger::CGP::TraFiC::Formatter>).  The output will
be written to this area as it is expected for these file to be generated in a single step.

=item Args

Path to output of L<Sanger::CGP::TraFiC::Formatter>

=back

=head3 set_min_reads

=over

=item Description

Define the minimum number of reads required to generate a cluster.  Defaults to 5.

=item Args

Integer, number of reads that must exist within a cluster.

=back

=head2 Processing

=head3 cluster_hits

=over

=item Description

Load the hit files, cluster and output in L<cluster format|Format of output>.

=item Note

No arguments as paths are determined at higher level.

=back

=head3 reciprocal_clusters

=over

=item Description

Compare the clustered data for pos and neg data identifying reciprocal entities.  A reciprocal output
file is generated.  See L<formatting section|Format of output>.

=item Note

No arguments as paths are determined at higher level.

=back
