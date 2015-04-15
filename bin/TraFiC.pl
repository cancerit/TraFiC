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
use Pod::Usage;
use Cwd qw(abs_path);
use Getopt::Long;
use File::Path qw(make_path remove_tree);
use File::Copy qw(move);
use Const::Fast qw(const);
use Carp;
$Carp::Verbose = 1;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use PCAP::Cli;
use PCAP::Threaded;
use Sanger::CGP::TraFiC;
use Sanger::CGP::TraFiC::Implement;

const my @VALID_PROCESS => qw(select orphans split mask cluster filter);
my %index_max = ( 'select'  => -1,
                  'orphans'  => 2,
                  'split'   => 1,
                  'mask'    => -1,
                  'cluster' => 1,
                  'filter'  => 1,);


{
  my $options = setup();
  Sanger::CGP::TraFiC::Implement::prepare($options);
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});

  # register any process that can run in parallel here
  $threads->add_function('select', \&Sanger::CGP::TraFiC::Implement::select_candidate_reads);
  $threads->add_function('orphans', \&Sanger::CGP::TraFiC::Implement::merge_select_orphans);
  $threads->add_function('split', \&Sanger::CGP::TraFiC::Implement::split_candidate_reads);
  $threads->add_function('mask', \&Sanger::CGP::TraFiC::Implement::repeat_mask_te);
  $threads->add_function('cluster', \&Sanger::CGP::TraFiC::Implement::cluster_te);

  # start processes here (in correct order obviously), add conditions for skipping based on 'process' option
  $threads->run($options->{'ref_count'}*2, 'select', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'select');
  $threads->run(2, 'orphans', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'orphans');
  $threads->run(2, 'split', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'split');

  if(!exists $options->{'process'} || $options->{'process'} eq 'mask') {
    unless($options->{'split_count'}) { # duplicated here incase this wasn't available at startup.
      $options->{'split_files'} = Sanger::CGP::TraFiC::Implement::split_files($options->{'tmp'});
      $options->{'split_count'} = scalar @{$options->{'split_files'}};
    }
    $threads->run($options->{'split_count'}, 'mask', $options);
  }

  $threads->run(2, 'cluster', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'cluster');

  if(!exists $options->{'process'} || $options->{'process'} eq 'filter') {
    Sanger::CGP::TraFiC::Implement::filter_results($options);

    cleanup($options);
  }
}

sub cleanup {
  my $options = shift;
  my $exit = system("tar -C $options->{tmp} -cZf $options->{outdir}/clusters.tar.gz clusters");
  die "Failed to archive clusters" if($exit);
  if(-d "$options->{outdir}/logs") {
    remove_tree("$options->{outdir}/logs_old") if(-e "$options->{outdir}/logs_old");
    move "$options->{outdir}/logs", "$options->{outdir}/logs_old";
  }
  move "$options->{tmp}/logs", "$options->{outdir}/logs";
  remove_tree $options->{'tmp'};
}

sub setup {
  my %opts;
  pod2usage(-msg  => "\nERROR: Option must be defined.\n", -verbose => 1,  -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
  GetOptions( 'h|help' => \$opts{'h'},
              'm|man' => \$opts{'m'},
              'c|cpus=i' => \$opts{'threads'},
              'o|outdir=s' => \$opts{'outdir'},
              't|tumour=s' => \$opts{'tumour'},
              'n|normal=s' => \$opts{'normal'},
              'td|tdiv=i' => \$opts{'tdiv'},
              'nd|ndiv=i' => \$opts{'ndiv'},
              'd|maskdb=s' => \$opts{'maskdb'},
              'e|engine=s' => \$opts{'engine'},
              'r|rm=s' => \$opts{'rm'},
              's|support=i' => \$opts{'support'},
              'k|known=s' => \$opts{'known'},
              'di|div=i' => \$opts{'div'},
              'g|germline=s' => \$opts{'germline'},
              'fm|filtmin=i' => \$opts{'filtmin'},
              'gm|germmin=i' => \$opts{'germmin'},
              'rl|length=i' => \$opts{'length'},
              'p|process=s' => \$opts{'process'},
              'i|index=i' => \$opts{'index'},
              'l|limit=i' => \$opts{'limit'},
              'v|version' => \$opts{'version'},
              'j|jobs' => \$opts{'jobs'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});

  if($opts{'version'}) {
    print 'Version: ',Sanger::CGP::TraFiC::Implement->VERSION,"\n";
    exit 0;
  }

  PCAP::Cli::file_for_reading('tumour', $opts{'tumour'});
  PCAP::Cli::file_for_reading('normal', $opts{'normal'});
  PCAP::Cli::file_for_reading('maskdb', $opts{'maskdb'});
  PCAP::Cli::file_for_reading('rm', $opts{'rm'});
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

  # defaults
  $opts{'tdiv'} = 5 unless(defined $opts{'tdiv'});
  $opts{'ndiv'} = 10 unless(defined $opts{'ndiv'});
  $opts{'support'} = 5 unless(defined $opts{'support'});
  $opts{'div'} = 200 unless(defined $opts{'div'});
  $opts{'filtmin'} = 2 unless(defined $opts{'filtmin'});
  $opts{'germmin'} = 5 unless(defined $opts{'germmin'});
  $opts{'length'} = 100 unless(defined $opts{'length'});

  for my $optional(qw(known germline)){
    if(defined $opts{$optional}) {
      PCAP::Cli::file_for_reading($optional, $opts{$optional});
    }
    else {
      delete $opts{$optional};
    }
  }


  $opts{'outdir'} = abs_path($opts{'outdir'});

  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
  delete $opts{'limit'} unless(defined $opts{'limit'});

  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpTraFiC');
  make_path($tmpdir) unless(-d $tmpdir);
  my $progress = File::Spec->catdir($tmpdir, 'progress');
  make_path($progress) unless(-d $progress);
  my $logs = File::Spec->catdir($tmpdir, 'logs');
  make_path($logs) unless(-d $logs);

  $opts{'tmp'} = $tmpdir;

  $opts{'ordered_seqs'} = Sanger::CGP::TraFiC::Implement::get_sequences($opts{'tumour'});
  $opts{'ref_count'} = scalar @{$opts{'ordered_seqs'}};
  $opts{'split_files'} = Sanger::CGP::TraFiC::Implement::split_files($opts{'tmp'});
  $opts{'split_count'} = scalar @{$opts{'split_files'}};

  $opts{'engine'} = 'crossmatch' unless(defined $opts{'engine'});

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    my $jobs = $index_max{$opts{'process'}};
    if($jobs == -1) {
      if($opts{'process'} eq 'select') {
        $jobs = $opts{'ref_count'}*2;
      }
      elsif($opts{'process'} eq 'mask') {
        die "ERROR: Previous steps don't appear to have been completed\n" unless($opts{'split_count'});
        $jobs = $opts{'split_count'};
      }
    }

    croak "I don't know how to determine the jobs for this step..." if($jobs == -1);

    if(exists $opts{'jobs'} && defined $opts{'jobs'}) {
      print "Jobs to complete process '$opts{process}' = $jobs\n";
      exit;
    }

    if(exists $opts{'index'}) {
      my $max = $index_max{$opts{'process'}};
      $max = $jobs if($max==-1 && $jobs);
      croak "'limit' processing is not configured" if(exists $opts{'limit'});
      if($max==-1){
        croak "I don't know how to determine the max indicies for this step";
      }

#      die "ERROR: based on reference and exclude option index must be between 1 and $refs\n" if($opts{'index'} < 1 || $opts{'index'} > $max);
      die "ERROR: Need to figure out how to determine number of jobs for RM step" if($opts{'index'} < 1 || $opts{'index'} > $max);
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);

      die "No max has been defined for this process type\n" if($max == 0);

      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max, 1);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  # now safe to apply defaults
  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  return \%opts;
}

__END__

=head1 TraFiC.pl

Reference implementation of Cancer Genome Project transposable element detection.

=head1 SYNOPSIS

TraFiC.pl [options]

  Required parameters:
    -outdir    -o   Folder to output result to.
    -tumour    -t   Tumour BAM file
    -normal    -n   Normal BAM file
    -maskdb    -d   RepeatMasker database

  RepeatMasker
    -engine    -e   RepeatMasker engine [crossmatch]
    -rm        -r   Path to RepeatMasker binary [search $PATH]
    -tdiv      -td  Tumour max divergence from RepeatMasker database [5]
    -ndiv      -nd  Normal max divergence from RepeatMasker database [10]

  Clustering
    -support   -s   Minimum supporting reads required to emit a cluster. [5]

  Filtering
    -known     -k   Known genomic events (e.g. TE regions)
    -div       -di  Only user records from -k where divergence is < this value. [200]
    -germline  -g   Pre-calculated germline events (aka Normal Panel)
    -filtmin   -fm  Minimum number of reads required in wt pos/neg cluster data to exclude a cluster. [2]
    -germmin   -gm  Minimum number of reads required in germline cluster data to exclude a cluster. [5]
    -length    -rl  Max read length in input data [100]

  General
    -cpus      -c   Number of cores to use. [1]
                     - recommend max 4 during 'input' process.
    -limit     -l   When defined with '-cpus' internally thread concurrent processes.
                     - requires '-p', specifically for pindel/pin2vcf steps

  Targeted processing (further detail under OPTIONS):
    -process   -p   Only process this step then exit, optionally set -index
    -index     -i   Optionally restrict '-p' to single job
    -jobs      -j  Declare with -p to determine how many jobs are needed for this step

  Other:
    -help      -h   Brief help message.
    -man       -m   Full documentation.
    -version   -v   Version

  File list can be full file names or wildcard, e.g.
    pindel.pl -c 4 -r some/genome.fa[.gz] -o myout -t tumour.bam -n normal.bam

  Run with '-m' for possible input file types.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  select
  split
  mask
  cluster
  filter

=item B<-index>

Possible index ranges for processes above are:

  select   = 1..2
  split  = 1
  mask = 1..<n>
  cluster   = 1
  filter    = 1

If you want STDOUT/ERR to screen ensure index is set even for single job steps.

=back
