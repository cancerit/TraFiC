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
use Test::More;
use Try::Tiny qw(try catch finally);
use File::Temp;

use constant MODULE => 'Sanger::CGP::TraFiC::Filter';

use_ok(MODULE);

my $mod;

try {
  note('Clear path and confirm that errors occur due to inability to find sort');
  local $ENV{'PATH'} = q{};
  my $bad_mod = MODULE->new;
  fail('Formatter new should fail as sort is not found in path');
}
catch {
  pass('Formatter new fails as sort is not found in path');
};

$mod = new_ok(MODULE);

my $dir = File::Temp->newdir('TraFiC_tests_XXXX', CLEANUP => 1);
ok($mod->set_output($dir), "Set output path to $dir");

#fail('need output from cluster test in testData before I can complete this');

done_testing();
