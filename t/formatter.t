#!usr/bin/perl

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
use autodie;
use Test::More;
use Try::Tiny qw(try catch finally);
use File::Temp;
use Const::Fast qw(const);

use FindBin qw($Bin);
my $test_data = "$Bin/data";

const my $MODULE => 'Sanger::CGP::TraFiC::Formatter';

# these values are for checking loading of records is correct
const my $RM_LINE_COUNT => 4182;
const my $RM_HEADER_LINES => 3;
const my $RM_SUB_OPT => 8;
const my $RM_ALT_FAM => 3;

# these values are to check correct repeat is selected for reads with multiple answers (not marked as sub optimal '*')
const my %RM_CORRECT_FAM => ('first' => ['HS17_08691:6:2114:19099:58478', 'ERV'],
                                'second' => ['HS17_08691:8:2210:16491:31098', 'L1'], );

const my $POS_EXPECTED_HITS => 2041;
const my $NEG_EXPECTED_HITS => 2127;

my $mod;

subtest 'Instatiation checks' => sub {
  use_ok($MODULE);
  try {
    note('Clear path and confirm that errors occur due to inability to find sort');
    local $ENV{'PATH'} = q{};
    my $bad_mod = $MODULE->new;
    fail('Formatter new should fail as sort is not found in path');
  }
  catch {
    pass('Formatter new fails as sort is not found in path');
  };

  $mod = new_ok($MODULE);
  $mod->silence;
};

my $dir = File::Temp->newdir('TraFiC_tests_XXXX', CLEANUP => 1);
ok($mod->set_output($dir), "Set output path to $dir");

my $rm_data;
subtest 'format_rm_results check' => sub {
  subtest 'load_rm_output check' => sub {
    $rm_data = $mod->load_rm_output(["$test_data/rm.out"]);
    is(ref $rm_data, 'HASH', 'Expected hash ref');

    my $expected_hits = $RM_LINE_COUNT - $RM_HEADER_LINES - $RM_SUB_OPT - $RM_ALT_FAM;

    is(scalar(keys %{$rm_data}), $expected_hits, qq{Expect $expected_hits repeat-masker hits.});

    is($rm_data->{ $RM_CORRECT_FAM{'first'}->[0] }, $RM_CORRECT_FAM{'first'}->[1], qq{First hit selected});

    is($rm_data->{ $RM_CORRECT_FAM{'second'}->[0] }, $RM_CORRECT_FAM{'second'}->[1], qq{Second hit selected});
  };

  ok($mod->output_matched_hits($rm_data, ["$test_data/candidate.fa"]), 'output_matched_hits check');

  subtest 'sort_hits check' => sub {
    ok($mod->sort_hits, "Sorting/moving hit files");
    is(count_records("$dir/pos_hits.txt"), $POS_EXPECTED_HITS, "Expecting $POS_EXPECTED_HITS hits on positive strand");
    is(count_records("$dir/neg_hits.txt"), $NEG_EXPECTED_HITS, "Expecting $NEG_EXPECTED_HITS hits on negative strand");
  };
};

done_testing();


sub count_records {
  my $file = shift;
  my $lines;
  open my $LCOUNT, $file;
  while(<$LCOUNT>){$lines++ unless($_ =~ m/^#/);};
  close $LCOUNT;
  return $lines;
}
