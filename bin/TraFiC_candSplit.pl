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
use Carp;
use File::Path qw(make_path);
$Carp::Verbose = 1;

use FindBin qw($Bin);
use lib "$Bin/../modules";

try {
  my $options = option_builder();
  my $record_count = 0;
  my $file_count = 0;
  my $OUT_FH;
  while (my $in = <>) {
    if($record_count == 0) {
      close $OUT_FH if(defined $OUT_FH);
      my $file_name = $options->{'o'}.++$file_count;
      warn "creating $file_name\n";
      open $OUT_FH, '>', $file_name;
    }
    print $OUT_FH $in;
    $in = <>;
    print $OUT_FH $in;
    $record_count = 0 if(++$record_count == $options->{'r'});
  }
  close $OUT_FH if(defined $OUT_FH);
  print "$file_count split files generated\n";
}
catch {
  croak "An error occurred while splitting files\n$_" if($_);
};

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'        => \$opts{'h'},
		'o|output=s'    => \$opts{'o'},
		'r|records=n'   => \$opts{'r'},
	);

	pod2usage(0) if($opts{'h'});

	pod2usage(q{Output prefix '-o' must be defined.}) unless($opts{'o'});

	$opts{'r'} = 350_000 unless($opts{'r'});

  pod2usage("Unexpected data remains on command after parsing:\n\t".(join ' ', @ARGV)) if(@ARGV > 0);

	return \%opts;
}

__END__

=head1 NAME

splitCandidates.pl - Split candidate fa files into chunks with non-padded suffix (easier to farm distribute)

=head1 SYNOPSIS

splitCandidates.pl [-h] -o /some/path/prefix. < candidates.fa

  Required options:

    --output   (-o)   Output prefix, directories must be present.

  Other options:

    --records  (-r)   Number of candidate reads to be output per file. [350000]

    --help      (-h)  This message.

  Examples:
    perl splitCandidates.pl -o $HOME/trafic/split/COLO-829. < candidates.fa

=cut
