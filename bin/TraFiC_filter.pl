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

use Sanger::CGP::TraFiC::Filter;

try {
  my $options = option_builder();
  my $filter = Sanger::CGP::TraFiC::Filter->new;
  $filter->set_output($options->{'o'});

  $filter->set_clusters($options->{'i'});
  $filter->add_filter_file($options->{'f'});
  $filter->set_germline_file($options->{'g'});
  $filter->set_te_file($options->{'t'});

  $filter->assume_sorted if($options->{'s'});
  $filter->set_divergence($options->{'d'});
  $filter->set_readlength($options->{'r'});
  $filter->set_min_filter_reads($options->{'fm'});
  $filter->set_min_germline_reads($options->{'gm'});

  $filter->filter;
}
catch {
  croak $_;
};


sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'      => \$opts{'h'},
		'o|output=s'  => \$opts{'o'},
		'i|input=s'   => \$opts{'i'},
		'f|filter=s@' => \$opts{'f'},
		'fm|filtmin=i'=> \$opts{'fm'},
		'g|germ=s'   => \$opts{'g'},
		'gm|germmin=i'=> \$opts{'gm'},
		't|tes=s'     => \$opts{'t'},
		'd|div=i'     => \$opts{'d'},
		'r|readlen=i' => \$opts{'r'},
		's|sorted'    => \$opts{'s'},
	);

	pod2usage(0) if($opts{'h'});

  pod2usage(q{Option '-i' must be defined}) unless(defined $opts{'i'});
  pod2usage(q{Option '-f'  or '-t' must be defined at least once}) unless(defined $opts{'f'} || defined $opts{'t'} );

  pod2usage(q{Option '-d' has no effect unless '-t' is defined.}) if($opts{'d'} && !defined $opts{'t'});
  pod2usage(q{Option '-r' has no effect unless '-t' is defined.}) if($opts{'r'} && !defined $opts{'t'});

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

  pod2usage("Unexpected data remains on command after parsing:\n\t".(join ' ', @ARGV)) if(@ARGV > 0);

	return \%opts;
}

__END__

=head1 NAME

TraFic_filter.pl - Filter clusters

=head1 SYNOPSIS

TraFic_filter.pl [-h] -o /some/path/ -r %.out -f %.fa

  Required options:

    --input   (-i)  Single data file to be filtered

      NOTE: All files must be formatted to include the first 5 fields as described in Cluster.pm
            see end of this message.

  Other options:

    --filter  (-f)  One or more files to be used as a filter.

    --filtmin (-fm) Minimum number of reads required in wt pos/neg cluster data to exclude a cluster. [2]

    --germ    (-g)  Pre-calculated Germline events file.

    --germmin (-gm) Minimum number of reads required in germline cluster data to exclude a cluster. [5]

    --output  (-o)  Output written to specified folder as filtered.txt. [.]

    --tes     (-t)  Transposable element regions from genome, format mimics clustered file core
                    fields, see below.  Exception is TOTAL_READS is replaced with the divergence value.

    --div     (-d)  Only user records from -t where divergence is < this value. [200]

    --readlen (-r)  Approximation of read length in input data, as clusters are merged sensible
                    to use max value. [100]

    --sorted  (-s)  Assume all inputs are sorted

  Examples:
    TraFic_filter -o tumour/filtered/ -i tumour/reciprocal_clusters.txt -f normal/reciprocal_clusters.txt \
     -f ref/known_false_pos.txt -f ref/normal_panel.txt \
     -t hg19_TEs.txt -d 200

=head2 Input format

=over

=item CHR

Chromosome or sequence identifier.

=item FAMILY

Masking identified reads as hitting repeat of this type.

=item L_POS

Left most position of cluster.

=item R_POS

Left most position of cluster.

=item TOTAL_READS

Total number of reads that support this cluster

=back
