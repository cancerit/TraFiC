package Sanger::CGP::TraFiC::Filter;

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
use File::Path qw(make_path);
use File::Temp;
use List::Util qw(min max);
use File::Which qw(which);

use Const::Fast qw(const);

use Sanger::CGP::TraFiC;

const my $TAB => "\t";

const my $POS_L_POS => 0;
const my $POS_R_POS => 1;
const my $NEG_L_POS => 9;
const my $NEG_R_POS => 10;

sub new {
    my ($class) = @_;
    my $self = { };
    bless $self, $class;
    croak "Unable to find standard unix 'sort' in path." unless(defined which('sort'));
    $self->{'min_filter_reads'} = 2;
    $self->{'min_germline_reads'} = 5;
    $self->{'divergence'} = 200;
    $self->{'readlength'} = 100;
    $self->{'filters'} = [];
    $self->{'assume_sorted'} = 0;
    return $self;
}

sub assume_sorted {
  shift->{'assume_sorted'} = 1;
}

sub set_output {
  my ($self, $root_dir) = @_;
  make_path($root_dir) unless(-e $root_dir);
  $self->{'out_loc'} = $root_dir;
  return 1;
}

sub set_clusters {
  my ($self, $clusters) = @_;

  croak "$clusters is not a valid file\n" unless(-e $clusters && -f $clusters);

  $self->{'clusters'} = $clusters;

  # grab the header line for output purposes
  open my $CL, '<', $clusters;
  $self->{'header'} = <$CL>;
  close $CL;

  return 1;
}

sub filter {
  my ($self) = @_;
  my $cluster_set = $self->_load_clusters($self->{'clusters'});

  if($cluster_set->{'cluster_count'} > 0) {
    for my $filter(@{$self->{'filters'}}){
      $self->_filter( $cluster_set,
                      $self->_load_clusters($filter, $self->{'min_filter_reads'}, $cluster_set->{'clusters'}));
    }
  }
  else {
    warn "No clusters in tumour input file, expect only a header line in output file.\n";
  }

  if($cluster_set->{'cluster_count'} > 0 && defined $self->{'germline'}) {
     $self->_filter($cluster_set,
                    $self->_load_clusters($self->{'germline'}, $self->{'min_germline_reads'}, $cluster_set->{'clusters'}));
  }

  if($cluster_set->{'cluster_count'} > 0 && defined $self->{'trans_elements'}) {
     $self->_filter($cluster_set,
                    $self->_load_clusters($self->{'trans_elements'}, $self->{'divergence'}, $cluster_set->{'clusters'}, 'divergence'),
                    1);
  }

  $self->output($cluster_set->{'clusters'});
  return 1;
}

sub output {
  my ($self, $cluster_set) = @_;
  open my $OUT, '>', "$self->{out_loc}/filtered.txt.tmp";


  for my $family(keys %{$cluster_set}) {
    for my $chr(keys %{$cluster_set->{$family}}) {

      for my $cluster(@{$cluster_set->{$family}->{$chr}}) {
        print $OUT join $TAB, $chr,
                              $family,
                              @{$cluster};
        print $OUT $RS;
      }
    }
  }
  close $OUT;

  my $sort_command = sprintf "sort -k 1,1 -k 3,4n -k 2,2 -T %s -S 200M %s",
                              $self->{'out_loc'},
                              "$self->{out_loc}/filtered.txt.tmp";

  open my $FINAL, '>', "$self->{out_loc}/filtered.txt";
  print $FINAL $self->{'header'};

  open my $proc, '-|', $sort_command;
  while(my $line = <$proc>) {
    print $FINAL $line;
  }
  close $proc;

  close $FINAL;

  unlink "$self->{out_loc}/filtered.txt.tmp";
  return 1;
}

sub _filter {
  my ($self, $clusters, $filters, $is_te) = @_;
  my $query = $clusters->{'clusters'}; # we want these
  my $filter = $filters->{'clusters'}; # we want to remove these
  my @families = sort keys %{$query};

  my $started_with = $clusters->{'cluster_count'};

  PRIMARY: for my $family(@families) {
    # if this family doesn't exist in the filter skip it
    next unless(exists $filter->{$family});
    my @chrs = sort keys %{$query->{$family}};
    for my $chr(@chrs) {
      # skip this chromosome if not in filters
      next unless(exists $filter->{$family}->{$chr});

      my $initial_clusters = scalar @{$query->{$family}->{$chr}};


      if($is_te) {
        $query->{$family}->{$chr} = $self->_filter_tes($query->{$family}->{$chr}, $filter->{$family}->{$chr});
      }
      else {
        $query->{$family}->{$chr} = $self->_filter_set($query->{$family}->{$chr}, $filter->{$family}->{$chr});
      }
      my $remaining_set_clusters = scalar @{$query->{$family}->{$chr}};
      $clusters->{'cluster_count'} -= $initial_clusters - $remaining_set_clusters;

      # as filter data is only loaded if there are corresponding
      # family/chrs to filter against ensure that empty keys
      # are removed from structure
      if($remaining_set_clusters == 0) {
        delete $query->{$family}->{$chr};
        my @remaining_keys = keys %{$query->{$family}};
        delete $query->{$family} if(scalar @remaining_keys == 0);
      }

      # no point filtering when nothing left
      last PRIMARY if($clusters->{'cluster_count'} == 0);

      croak "Cluster count has become negative, this suggests a programming error" if($clusters->{'cluster_count'} < 0);
    }
  }

  my $removed = $started_with - $clusters->{'cluster_count'};

  print "$removed clusters removed by $filters->{file}\n";
  return 1;
}

sub _filter_tes {
  my ($self, $queries, $tes) = @_;
  my $filter_tot = scalar @{$tes};
  my $filter_i = 0;

  my $read_len = $self->{'readlength'};

  my @l_queries = @{$queries};

  my @new_clusters;
  PRIMARY: while(scalar @l_queries > 0) {
    last if($filter_i == $filter_tot);
    my $query = shift @l_queries;

    my $q_l_pos = (min $query->[$POS_R_POS], $query->[$NEG_R_POS]) + $read_len;
    my $q_r_pos = max $query->[$POS_L_POS], $query->[$NEG_L_POS];

    while($filter_i < $filter_tot) {
      if($tes->[$filter_i]->[$POS_R_POS] < $q_l_pos) {
        # filter data is prior to query data so catch up
        $filter_i++;
        next;
      }
      if($tes->[$filter_i]->[$POS_L_POS] > $q_r_pos) {
        # filter is ahead of query so query must be valid
        push @new_clusters, $query;
        # move on to next query
        next PRIMARY;
      }
      # if we get here the filter must overlap the query
      next PRIMARY;
    }
    push @new_clusters, $query;
  }
  push @new_clusters, @l_queries if(scalar @l_queries > 0);

  return \@new_clusters;
}

sub _filter_set {
  my ($self, $queries, $filters) = @_;

  my $filter_tot = scalar @{$filters};
  my $filter_i = 0;

  my $double_read_len = $self->{'readlength'} * 2;

  my @l_queries = @{$queries};

  my @new_clusters;
  PRIMARY: while(scalar @l_queries > 0) {
    last if($filter_i == $filter_tot);
    my $query = shift @l_queries;

    my $q_l_pos = $query->[$POS_L_POS];
    my $q_r_pos = $query->[$POS_R_POS];

    while($filter_i < $filter_tot) {
      if($filters->[$filter_i]->[$POS_R_POS]+$double_read_len < $q_l_pos) {
        # filter data is prior to query data so catch up
        $filter_i++;
        next;
      }
      if($filters->[$filter_i]->[$POS_L_POS]-$double_read_len > $q_r_pos) {
        # filter is ahead of query so query must be valid
        push @new_clusters, $query;
        # move on to next query
        next PRIMARY;
      }
      # if we get here the filter must overlap the query
      next PRIMARY;
    }
    push @new_clusters, $query;
  }
  push @new_clusters, @l_queries if(scalar @l_queries > 0);

  return \@new_clusters;
}

sub add_filter_file {
  my ($self, $filter) = @_;
  if(defined $filter) {
    if(ref $filter eq 'ARRAY') {
      for(@{$filter}) {
        croak "$_ is not a valid file\n" unless(-e $_ && -f $_);
        push @{$self->{'filters'}}, $_;
      }
    }
    else {
      croak "$filter is not a valid file\n" unless(-e $filter && -f $filter);
      push @{$self->{'filters'}}, $filter;
    }
  }
  return scalar @{$self->{'filters'}};
}

sub _load_clusters {
  my ($self, $input_set, $limit, $query_clusters, $divergence) = @_;

  my $cluster_file = $input_set;

  my $dir;

  unless($self->{'assume_sorted'}) {
    $dir = File::Temp->newdir('traficFilterXXXX', (DIR => $self->{'out_loc'}));
    # safest to use unix command as will be fast and lean
    my $sort_command = sprintf "sort -k 1,1 -k 3,4n -k 2,2 -T %s -S 200M -o %s %s",
                                $dir, # temp area for sorting
                                "$dir/sorted", # output file
                                $input_set; # input file
    croak "Failed to sort file using command: $sort_command" unless(system($sort_command) == 0);
    $cluster_file = "$dir/sorted";
  }

  my %data_set;
  my $cluster_count = 0;
  open my $CLUSTERS, '<', $cluster_file;
  while (my $line = <$CLUSTERS>) {
    next if($line =~ m/^#/);
    chomp $line;

    my ($chr, $family, @detail) = split $TAB, $line;
    if(defined $limit) {
      my $score = $detail[2]; # can be reads or divergence score
      if(defined $divergence) {
        next unless($limit > $score); # keep the record in the set
      }
      else {
        # this is only applied to reciprocal clustered data
        $score += $detail[11] if(scalar @detail > 10);
        next if($limit > $score); # discard the record
      }
    }

    if(defined $query_clusters) {
      next unless(exists $query_clusters->{$family} && exists $query_clusters->{$family}->{$chr});
    }

    push @{$data_set{$family}{$chr}}, \@detail;

    $cluster_count++;
  }
  close $CLUSTERS;

  return {'cluster_count' => $cluster_count,
          'clusters'      => \%data_set,
          'file'          => $input_set} ;
}

sub set_min_filter_reads {
  my ($self, $min_reads) = @_;
  $self->{'min_filter_reads'} = $min_reads if(defined $min_reads);
  return 1;
}

sub set_min_germline_reads {
  my ($self, $min_reads) = @_;
  $self->{'min_germline_reads'} = $min_reads if(defined $min_reads);
  return 1;
}

sub set_divergence {
  my ($self, $divergence) = @_;
  $self->{'divergence'} = $divergence if(defined $divergence);
  return 1;
}

sub set_readlength {
  my ($self, $readlength) = @_;
  $self->{'readlength'} = $readlength if(defined $readlength);
  return 1;
}

sub set_germline_file {
  my ($self, $file) = @_;

  croak "$file is not a valid file\n" unless(-e $file && -f $file);

  $self->{'germline'} = $file;
  return 1;
}

sub set_te_file {
  my ($self, $file) = @_;
  croak "$file is not a valid file\n" unless(-e $file && -f $file);
  $self->{'trans_elements'} = $file;
  return 1;
}

1;

__END__

=head1 NAME

Sanger::CGP::TraFiC::Filter

=head1 SYNOPSIS

 use Sanger::CGP::TraFiC::Filter;
 my $filter = Sanger::CGP::TraFiC::Filter->new;
 $filter->set_output($output_dir);
 $filter->set_clusters($sample_clusters);

 $filter->add_filter_file($filter_on_this);
 $filter->add_filter_file($filter_on_this_too);
  ...
  # OR
 $filter->add_filter_file(\@list_of_filter_files);

 $filter->filter;

=head1 GENERAL

Provides bulk filtering of a set of clusters.  Will allow more arbitrary processing of data after
main processing is completed.

=head1 METHODS

=head2 User functions

=head3 set_output

Set the output location for output files.

=head3 set_clusters

Specify file which contains the clusters that are of interest.

=head3 add_filter_file

Add files that are used to filter the content of C<set_clusters>.  Accepts both a simple scalar for
one file or an array reference for several files.

It is recommended that the first file added is the matched normal sample as this will remove the
most noise from the data.

=head3 filter

Filter the content of C<set_clusters> against each of the items in C<add_filter_file>.  Entries in
C<add_filter_file> are only sorted/loaded if data remains in the data set loaded from C<set_cluster>.

=head3 output

Output filtered data.  Includes sorting.

=head3 set_min_reads

Set the minimum number of reads that must support a filtering record  for it to be used.

Defaults to 5.

=head2 Internal functions

=head3 _load_clusters

Sort and then load cluster information into a data structure.

Sort is always invoked to ensure that data is consistent.  Most of the input should be relatively
small.
