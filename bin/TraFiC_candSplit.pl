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
use lib "$Bin/../lib";

try {
  my $options = option_builder();
  input_setup($options);

use Data::Dumper;

  my $record_count = 0;
  my $file_count = 0;
  my $OUT_FH;
  for my $in_file(@{$options->{'i'}}) {
    open my $INF, '<', $in_file;
    my $file_name;
    while (my $in = <$INF>) {
      if($record_count == 0) {
        close $OUT_FH if(defined $OUT_FH);
        $file_name = $options->{'o'}.++$file_count;
        warn "creating $file_name\n";
        open $OUT_FH, '>', $file_name;
      }
      print $OUT_FH $in or die "Failed to write to $file_name";
      $in = <$INF>;
      print $OUT_FH $in or die "Failed to write to $file_name";
      $record_count = 0 if(++$record_count == $options->{'r'});
    }
    close $INF;
  }
  close $OUT_FH if(defined $OUT_FH);
  print "$file_count split files generated\n";
}
catch {
  croak "An error occurred while splitting files\n$_" if($_);
};

sub input_setup {
  my $options = shift;
  my @in_files;
  die "Options '-i' ($options->{i}) does not exist as file or directory." unless(-e $options->{'i'});
  if(-d _) {
    my $indir = $options->{'i'};
    opendir(my $dh, $indir) || die "Unable to read from $indir: $OS_ERROR";
    while(my $thing = readdir $dh) {
      next if($thing =~ m/^[.]/);
      push @in_files, "$indir/$thing" if($thing =~ m/[.]fa$/);
    }
    closedir($dh);
  }
  else {
    push @in_files, $options->{'i'};
  }
  $options->{'i'} = \@in_files;
  return 1;
}

sub option_builder {
	my ($factory) = @_;

	my %opts;

	&GetOptions (
		'h|help'        => \$opts{'h'},
		'i|input=s'       => \$opts{'i'},
		'o|output=s'    => \$opts{'o'},
		'r|records=n'   => \$opts{'r'},
	);

	pod2usage(0) if($opts{'h'});

  pod2usage(q{Input file/folder '-i' must be defined.}) unless($opts{'i'});
	pod2usage(q{Output prefix '-o' must be defined.}) unless($opts{'o'});

	$opts{'r'} = 350_000 unless($opts{'r'});

  pod2usage("Unexpected data remains on command after parsing:\n\t".(join ' ', @ARGV)) if(@ARGV > 0);

	return \%opts;
}

__END__

=head1 NAME

TraFic_candSplit.pl - Split candidate fa files into chunks with non-padded suffix (easier to farm distribute)

=head1 SYNOPSIS

TraFic_candSplit.pl [-h] -o /some/path/prefix. < candidates.fa

  Required options:

    --input    (-i)   Input file or directory containing *.fa
    --output   (-o)   Output prefix, directories must be present.

  Other options:

    --records  (-r)   Number of candidate reads to be output per file. [350000]

    --help      (-h)  This message.

  Examples:
    TraFic_candSplit.pl -o $HOME/trafic/split/COLO-829. < candidates.fa

=cut
