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
use File::Path qw(make_path);
use File::Copy qw(move);
use Carp;
$Carp::Verbose = 1;

use FindBin qw($Bin);
use lib "$Bin/../modules";

use Sanger::CGP::TraFiC::Parser;

try {
  my $options = option_builder();
  my $parser = Sanger::CGP::TraFiC::Parser->new;
  $parser->set_output($options->{'o'});
  $parser->process_files($options->{'i'}, $options->{'r'});
}
catch {
  croak $_;
};

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'        => \$opts{'h'},
		'o|output=s'    => \$opts{'o'},
		'i|input=s@'    => \$opts{'i'},
		'r|readgroup=s' => \$opts{'r'},
	);

	pod2usage(0) if($opts{'h'});

  pod2usage(q{Input '-i' must be defined at least once}) unless($opts{'i'});
  pod2usage(q{Only 1 input file should be defined with '-r'}) if($opts{'r'} && @{$opts{'i'}} > 1);
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

candidateReads.pl - Find candidate reads for transposon or viral detection

=head1 SYNOPSIS

candidateReads.pl [-h] -o /some/path/ -i sample.bam -i lane.ump.bam ...

  Required options:

    --input   (-i)    Define at least once

  Other options:

    --output  (-o)    Output written to specified folder as candidates.fa [.]

    --readgroup (-r)  Readgroup id for parallel processing.
                        - only one input file should be specified

    --help      (-h)  This message

  Examples:
    perl ~/candidateReads.pl -o $HOME/trafic/COLO-829 -i COLO-829.bam

=cut
