package Sanger::CGP::TraFiC::Formatter;

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
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempdir);
use File::Which qw(which);
use Const::Fast qw(const);

use Sanger::CGP::TraFiC;

const my $TAB => "\t";

sub new {
    my ($class) = @_;
    my $self = { };
    bless $self, $class;
    croak "Unable to find standard unix 'sort' in path." unless(defined which('sort'));
    $self->{'silent'} = 0;
    return $self;
}

sub set_output {
  my ($self, $root_path) = @_;
  croak "$root_path does not exist." unless(-e $root_path);
  croak "$root_path is not a directory." unless(-d $root_path);
  croak "$root_path cannot be written to." unless(-w $root_path);
  $self->{'temp'} = tempdir ( 'TraFiC_XXXX', DIR => $root_path );
  make_path($self->{'temp'});

  $self->{'+_tmp'} = "$self->{temp}/pos_hits.txt";
  $self->{'-_tmp'} = "$self->{temp}/neg_hits.txt";

  open my $pos, '>', $self->{'+_tmp'};
  open my $neg, '>', $self->{'-_tmp'};

  $self->{'+'} = $pos;
  $self->{'-'} = $neg;

  $self->{'+_final'} = "$root_path/pos_hits.txt";
  $self->{'-_final'} = "$root_path/neg_hits.txt";
  return 1;
}

sub format_rm_results {
  my ($self, $out_list, $fa_list) = @_;
  my $rm_hits = $self->load_rm_output($out_list);
  $self->output_matched_hits($rm_hits, $fa_list);
  $self->sort_hits;
  return 1;
}

sub sort_hits {
  my ($self) = @_;
  close $self->{'+'};
  close $self->{'-'};
  # remove the fh entries so destroy doesn't try again.
  delete $self->{'+'};
  delete $self->{'-'};
  $self->sort_file('+');
  $self->sort_file('-');
  return 1;
}

sub sort_file {
  my ($self, $strand) = @_;
  # safest to use unix command as will be fast and lean
  # sort -k 2,2 -k 3,3n -o file -T $self->{temp} -S BUFFER_SIZE
  my $sort_command = sprintf "sort -k 4,4 -k 2,2 -k 3,3n -T %s -S 200M -o %s %s",
                              $self->{temp},
                              $self->{$strand.'_final'}, # output file
                              $self->{$strand.'_tmp'}; # input file
  croak "Failed to sort file using command: $sort_command" unless(system($sort_command) == 0);
  return 1;
}

sub output_matched_hits {
  my ($self, $rm_hits, $fa_list) = @_;

  my ($readname, $chr, $pos, $class, $strand);
  my $remaining_hits = scalar(keys %{$rm_hits});
  print "Correlating $remaining_hits records\n" unless($self->{'silent'});
  SHORT: for my $in_fa(@{$fa_list}) {
    open my $FA_FILE, '<', $in_fa;
    while (my $line = <$FA_FILE>) {
      <$FA_FILE>; # discard the sequence line
      chomp $line;
      ($readname, $chr, $pos, $class, $strand) = split /\t/, substr($line, 1); # ditch leading '>'
      if(exists $rm_hits->{$readname}) {
        print {$self->{$strand}} join $TAB, ($readname,
                                          $chr,
                                          $pos,
                                          $rm_hits->{$readname},
                                          $class,
                                          $strand);
        print {$self->{$strand}} $RS;
        delete $rm_hits->{$readname}; # removing things from the hash will speed up subsequent lookups
        if(--$remaining_hits == 0) {
          close $FA_FILE;
          print "All records correlated\n" unless($self->{'silent'});
          last SHORT;
        }
        print "\t$remaining_hits remaining\n" if($self->{'silent'} != 1 && $remaining_hits % 5000 == 0);
      }
    }
    close $FA_FILE;
  }
  return 1;
}

sub load_rm_output {
  my ($self, $out_list) = @_;
  my %rm_hits;
  my $last_score = 0;
  print "Loading repeat masker output\n" unless($self->{'silent'});
  for my $rm_out(@{$out_list}) {
    print "\t$rm_out\n" unless($self->{'silent'});
    open my $RM_FILE, '<', $rm_out;
    while (my $line = <$RM_FILE>) {
      chomp $line;
      $line =~ s/^\s+//;
      $line =~ s/\s+$//;

      # best hits always end with an id
      # alternate hits (shorter matches) are terminated with
      # an additional (unheaded) column containing '*'
      # we only want the best hits
      next unless($line =~ m/\d$/);

      my @bits = split /[\s\/]+/, $line;

      # all reads are output adjacent, however not ordered by score.
      if(exists $rm_hits{$bits[4]}) {
        # if this read is already stored but this record has a better score
        # then replace it, otherwise skip
        $rm_hits{$bits[4]} = $bits[11] if($bits[0] > $last_score);
      }
      else {
        $rm_hits{$bits[4]} = $bits[11]; # readname and family are only bits needed
      }
      $last_score = $bits[0];
    }
    close $RM_FILE;
  }
  print "Loading complete\n" unless($self->{'silent'});
  return \%rm_hits;
}

sub silence {
  shift->{'silent'} = 1;
}

sub unsilence {
  shift->{'silent'} = 0;
}

sub DESTROY {
  my $self = shift;
  close $self->{'+'} if(defined $self->{'+'});
  close $self->{'-'} if(defined $self->{'-'});

  remove_tree($self->{'temp'}) if(defined $self->{'temp'});;
  return 1;
}


1;

__END__

=head1 NAME

Sanger::CGP::TraFiC::Formatter

=head1 SYNOPSIS

 use Sanger::CGP::TraFiC::Formatter;
 my $formatter = Sanger::CGP::TraFiC::Formatter->new;
 $formatter->set_output($output_dir);
 $formatter->format_rm_results($repeatmasker_files, $fasta_files);

=head1 GENERAL

Designed to take one or more repeatmasker *.out files along with the input files to generate two
output files (positive anchors, negative anchors) containing the full information about the mapped
end of those reads that were successfully processed by repeat masker.

All the data required is encoded in the fasta header for each read:

=over

=item 1. readname

=item 2. chr

=item 3. pos

=item 4. repeat family

=item 5. strand (+/-)

=back

2,3 and 5 are based on the data for the mapped end of the pair.

This means we only need to store the readname and family of repeat in memory from the rm output as
the rest can be read back from the rm input files with no memory overhead.

=head1 METHODS

=head2 Constructor/configuration

=head3 new

=over

=item Description

Initialises the object and ensures that unix sort is available in path.

 my $formatter = Sanger::CGP::TraFiC::Formatter->new;

=item Errors

Will C<croak> if unable to find C<sort> in path with this message:

 Unable to find standard unix 'sort' in path

=back

=head3 set_output

=over

=item Description

Specify the output folder for collated results.

 $formatter->set_output($output_path);

=back

=head3 silence

=over

=item Description

Allows user to silence internal messages.  Mainly added for testing system

 $formatter->silence;

=back

=head3 unsilence

=over

=item Decription

Allows user resume internal messages.

 $formatter->unsilence;

=back

=head2 Processing

=head3 format_rm_results

=over

=item Description

Main processing function.  C<set_output> must be called prior to this.

 $formatter->format_rm_results(\@rm_out_files, \@rm_in_fasta);

=item Args

 rm_out_files - \@ of repeat masker output files (*.out).
 rm_in_fasta  - \@ reference of the fasta files presented to repeat masker.

=item Returns

Nothing is returned, results are written to the specified output location as pos/neg_hits.txt

=back

=head3 load_rm_output

=over

=item Description

Load the repeat-masker output file data.  Only best/longest hits are retained.  Not really intended
for external use, but no reason you can't use it if you want to parse a RM file for the basic data
retrieved.

=item Args

 rm_out_files - \@ of repeat masker output files (*.out).

=item Returns

 \% where key is readname and value is the family that masked this read.

=back


