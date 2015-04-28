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

use constant MODULE => 'Sanger::CGP::TraFiC::Parser';

# dropped optional tags and minimised seq/qual where not important to test
use constant PROPER_PAIR =>
        ( ['HS17_08691:6:1210:6217:36448',99,'1','10002','0','100M','=','10010','108','ACGT','EEEE'],
          ['HS17_08691:6:1210:6217:36448',147,'1','10010','25','100M','=','10002','-108','ACGT','EEEE'] );

use constant UNMAPPED_PAIR =>
        ( ['HS20_08443:7:1101:10000:18717',77,'*','0','0','*','*','0','0','ACGT','EEEE'],
          ['HS20_08443:7:1101:10000:18717',141,'*','0','0','*','*','0','0','ACGT','EEEE'] );

use constant SINGLE_END_MAPPED =>
        ( ['HS25_09078:8:1112:19570:49276',89,'1','9992',37,'2S3M1D95M','=','9992','0','ACGT','EEEE'],
          ['HS25_09078:8:1112:19570:49276',181,'1','9992','0','*','=','9992','0','ACGT','EEEE'] );


use constant INTER_CHROM_MAPPED =>
        ( ['HS25_09078:8:2110:14096:65629',65,'1','213325515',37,'63M1D37M','2','126897','0','ACGT','EEEE'],
          ['HS25_09078:8:2110:14096:65629',129,'2','126897',23,'6M1D94M','1','213325515','0','ACGT','EEEE'] );

use constant ABERRANT_MAPPED =>
        ( ['HS17_08691:6:1209:11505:95018',145,'1','16387',25,'100M','=','16393','-94','ACGT','EEEE'],
          ['HS17_08691:6:1209:11505:95018',97,'1','16393',37,'94M2D3M3S','=','16387','94','ACGT','EEEE'] );

use constant MAPQ_ORDER_GOOD =>
        ( ['HS17_08691:6:1209:11505:95018',81,'1','16387',37,'100M','=','16393','-94','ACGT','EEEE'],
          ['HS17_08691:6:1209:11505:95018',161,'1','16393',0,'94M2D3M3S','=','16387','94','TGCA','EEEE'] );

use constant MAPQ_ORDER_GOOD_OUT => qq{>HS17_08691:6:1209:11505:95018\t1\t16387\tABERRANT\t-\nTGCA\n};

use constant MAPQ_ORDER_SWAP =>
        ( ['HS17_08691:6:1209:11505:95018',81,'1','16387',0,'100M','=','16393','-94','ACGT','EEEE'],
          ['HS17_08691:6:1209:11505:95018',161,'1','16393',37,'94M2D3M3S','=','16387','94','TGCA','EEEE'] );

use constant MAPQ_ORDER_SWAP_OUT => qq{>HS17_08691:6:1209:11505:95018\t1\t16393\tABERRANT\t+\nACGT\n};


use_ok(MODULE);

my $mod = new_ok(MODULE);

is(scalar @{$mod->pair_to_candidate(PROPER_PAIR)}, 0, 'Proper pairs discarded');
is(scalar @{$mod->pair_to_candidate(UNMAPPED_PAIR)}, 0, 'Unmapped pairs discarded');

is(scalar @{$mod->pair_to_candidate(SINGLE_END_MAPPED)}, 2, 'Retained single end mapped');
is($mod->pair_to_candidate(SINGLE_END_MAPPED)->[0], 'SINGLE_END', 'Single end mapped, correct class');

is(scalar @{$mod->pair_to_candidate(INTER_CHROM_MAPPED)}, 2, 'Retained inter-chromosomal mapped pair');
is($mod->pair_to_candidate(INTER_CHROM_MAPPED)->[0], 'INTER_CHROM', 'Inter-chromosomal mapped, correct class');

is(scalar @{$mod->pair_to_candidate(ABERRANT_MAPPED)}, 2, 'Retained aberrant pair');
is($mod->pair_to_candidate(ABERRANT_MAPPED)->[0], 'ABERRANT', 'Aberrant pair, correct class');

is(${$mod->pair_to_candidate(MAPQ_ORDER_GOOD)->[1]}, MAPQ_ORDER_GOOD_OUT, 'Check order is determined by MAPQ when both ends mapped (1)');

is(${$mod->pair_to_candidate(MAPQ_ORDER_SWAP)->[1]}, MAPQ_ORDER_SWAP_OUT, 'Check order is determined by MAPQ when both ends mapped (2)');

done_testing();
