#!/usr/bin/perl

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
use Try::Tiny qw(try catch);
use Getopt::Long qw(:config no_auto_abbrev no_ignore_case);
use Pod::Usage;
use File::Glob ':glob';
use File::Path qw(make_path);
use File::Copy qw(move);
use Carp;
$Carp::Verbose = 1;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use Sanger::CGP::TraFiC::Formatter;
use Sanger::CGP::TraFiC::Cluster;

try {
  my $options = option_builder();
  my $formatter = Sanger::CGP::TraFiC::Formatter->new;
  $formatter->set_output($options->{'o'});
  $formatter->format_rm_results($options->{'r'}, $options->{'f'});
  my $cluster = Sanger::CGP::TraFiC::Cluster->new;
  $cluster->set_input_output($options->{'o'});
  $cluster->set_min_reads($options->{'m'});
  warn "Clustering\n";
  $cluster->cluster_hits;
  warn "Finding reciprocal clusters\n";
  $cluster->reciprocal_clusters;
}
catch {
  croak $_;
};

sub get_input {
  #  Expand file lists with % wildcards
  my ($list) = @_;
  my @files;
  foreach my $i_val(@{$list}) {
    if($i_val =~ m/\%/) {
      $i_val =~ s/\%/*/g;
      push @files, bsd_glob($i_val);
    }
    else {
      push @files, $i_val;
    }
  }
  return \@files;
}


sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'     => \$opts{'h'},
		'o|output=s' => \$opts{'o'},
		'r|rm=s@'    => \$opts{'r'},
		'f|fa=s@'    => \$opts{'f'},
		'm|min=i'    => \$opts{'m'},
	);

	pod2usage(0) if($opts{'h'});

  pod2usage(q{Input ('-f') must be defined at least once}) unless($opts{'f'});
  pod2usage(q{Input ('-r') must be defined at least once}) unless($opts{'r'});

  if($opts{'o'}) {
    $opts{'o'} =~ s{[/\\]$}{};
    unless(-e $opts{'o'}) {
      make_path($opts{'o'}) or croak "Unable to create $opts{o}";
    }
  }
  else {
    warn qq{Output will be to current directory\n};
    $opts{'o'} = '.';
  }

  $opts{'f'} = get_input($opts{'f'});
  $opts{'r'} = get_input($opts{'r'});

  pod2usage("Unexpected data remains on command after parsing:\n\t".(join ' ', @ARGV)) if(@ARGV > 0);

	return \%opts;
}

__END__

=head1 NAME

TraFic_cluster.pl - Generate files suitable for input to cluster discovery

=head1 SYNOPSIS

TraFic_cluster.pl [-h] -o /some/path/ -r %.out -f %.fa

  Required options:

    --rm      (-r)  Repeat masker *.out files
    --fa      (-f)  Modified fasta format originally supplied to repeat masker

      NOTE: Both can be defined once with/without wildcard % or multiple times

  Other options:

    --output  (-o)  Output written to specified folder (current working directory otherwise):
                      - pos_hits.txt
                      - neg_hits.txt
                      - pos_clusters.txt
                      - neg_clusters.txt

    --min     (-m)  Minimum number of reads required to emit a cluster. [5]

    --help    (-h)  This message

  Examples:
    TraFic_cluster.pl -o $HOME/trafic -r COLO-829.%.fa.out -f COLO-829.fa

=cut
