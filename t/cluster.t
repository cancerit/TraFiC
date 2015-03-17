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
use File::Temp;
use File::Copy qw(copy);
use Const::Fast qw(const);
use FindBin qw($Bin);
my $data_path = "$Bin/data";

const my $MODULE => 'Sanger::CGP::TraFiC::Cluster';
const my $MIN_READ_HIGH => 6;
const my $MIN_READ_DEFAULT => 5;
const my $POS_CLUSTER_DEFAULT => 3;
const my $POS_CLUSTER_HIGH => 2;
const my $NEG_CLUSTER_DEFAULT => 5;
const my $NEG_CLUSTER_HIGH => 3;

const my $EXPECTED_RECIP => 1;

const my $READ_LIST => [qw(HS17_08691:7:1313:4088:84581 HS20_08443:7:1103:9634:78733
                              HS25_09078:8:2203:8897:14960 HS20_08443:7:1102:6388:97835)];

const my $READ_STR => 'HS17_08691:7:1313:4088:84581,HS20_08443:7:1103:9634:78733,HS25_09078:8:2203:8897:14960,HS20_08443:7:1102:6388:97835';

const my $READ_NULL_LIST => q{.};

const my $POS_LARGE_FAMILY => {'family' => 'L1', 'chr' => 12, 'members' => 1};

my $tmp_dir = &setup;

use_ok($MODULE);

my $mod = new_ok($MODULE);

ok($mod->set_input_output($tmp_dir), "Set input/output path to $tmp_dir");

subtest 'min_reads value checks' => sub {
  ok($mod->set_min_reads($MIN_READ_HIGH), "Setting min_reads is successful");
  is($mod->{'min_reads'}, $MIN_READ_HIGH, "Correct value in min_reads");
  ok($mod->set_min_reads(), "set_min_reads with undef successful");
  is($mod->{'min_reads'}, $MIN_READ_DEFAULT, "Correct value in min_reads when set with undef");
};

subtest 'readname formatting' => sub {
  is(${Sanger::CGP::TraFiC::Cluster::_format_readnames([])}, $READ_NULL_LIST, 'Empty list generates correct readname string');
  is(${Sanger::CGP::TraFiC::Cluster::_format_readnames($READ_LIST)}, $READ_STR, 'Populated list generates correct readname string');
};

subtest 'clustering' => sub {
  $mod->set_min_reads($MIN_READ_HIGH);
  ok($mod->cluster_file('+'), 'Clustered pos data, min_reads='.$MIN_READ_HIGH);
  is(count_records("$tmp_dir/pos_clusters.txt"), $POS_CLUSTER_HIGH, 'Expected pos clusters for min_reads='.$MIN_READ_HIGH);
  ok($mod->cluster_file('-'), 'Clustered neg data, min_reads='.$MIN_READ_HIGH);
  is(count_records("$tmp_dir/neg_clusters.txt"), $NEG_CLUSTER_HIGH, 'Expected neg clusters for min_reads='.$MIN_READ_HIGH);

  $mod->set_min_reads;
  ok($mod->cluster_hits, 'Full clustering with defaults');
  is(count_records("$tmp_dir/pos_clusters.txt"), $POS_CLUSTER_DEFAULT, 'Expected pos clusters for min_reads='.$MIN_READ_DEFAULT);
  is(count_records("$tmp_dir/neg_clusters.txt"), $NEG_CLUSTER_DEFAULT, 'Expected neg clusters for min_reads='.$MIN_READ_DEFAULT);
};

subtest 'reciprocals' => sub {
  my ($pos_clusters, $neg_clusters);
  ok($pos_clusters = $mod->clustered_struct('+'), 'Generating pos cluster structure');

  ok(exists $pos_clusters->{ $POS_LARGE_FAMILY->{'family'} }, 'Family exists (pos)');

  ok(exists $pos_clusters->{ $POS_LARGE_FAMILY->{'family'} }->{ $POS_LARGE_FAMILY->{'chr'} }, 'Chr exists (pos)');
  is(scalar @{$pos_clusters->{ $POS_LARGE_FAMILY->{'family'} }->{ $POS_LARGE_FAMILY->{'chr'} } },
      $POS_LARGE_FAMILY->{'members'}, 'Expected number of entries in family/chr');

  ok($neg_clusters = $mod->clustered_struct('-'), 'Generating neg cluster structure');
  ok($mod->reciprocal_hits($pos_clusters, $neg_clusters), 'Reciprocal hits generated');

  is(count_records("$tmp_dir/reciprocal_clusters.txt"), $EXPECTED_RECIP, 'Expected reciprocal clusters');
};


done_testing();


sub setup {
  my $dir = File::Temp->newdir('TraFiC_tests_XXXX', CLEANUP => 1);
  copy "$data_path/pos_hits.txt", "$dir/pos_hits.txt";
  copy "$data_path/neg_hits.txt", "$dir/neg_hits.txt";
  return $dir;
}

sub count_records {
  my $file = shift;
  my $lines = 0;
  open my $LCOUNT, $file;
  while(<$LCOUNT>){$lines++ unless($_ =~ m/^#/);};
  close $LCOUNT;
  return $lines;
}
