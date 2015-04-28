package Sanger::CGP::TraFiC::Implement;

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
use warnings FATAL => 'all';
use English qw( -no_match_vars );
use autodie qw(:all);
use Const::Fast qw(const);
use File::Which qw(which);
use File::Path qw(make_path remove_tree);
use File::Temp qw(tempfile);
use FindBin qw($Bin);
use File::Basename qw(fileparse);

use Sanger::CGP::TraFiC;

use PCAP::Threaded;
use PCAP::Bam;

const my $CAND_READS => q{ -c %s -o %s -i %s};
const my $CAND_ORPH => q{ -i %s};
const my $CAND_SPLIT => q{ -o %s -i %s};
const my $CLUSTER => q{ -r %s/rm/%%/%s.%%.out -f %s/split/%s.%% -m %s -o %s};

### the core work

sub select_candidate_reads {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    for my $chr(@{$options->{'ordered_seqs'}}) {
      next if($iter++ != $index); # skip to the relevant input in the list

      ## build command for this index
      #

      my $sample = sanitised_sample_from_bam($input);
      my $output = File::Spec->catdir($tmp, 'candidates', $sample);
      make_path($output) unless(-e $output);

      my $command = "$^X ";
      $command .= _which('TraFiC_candReads.pl');
      $command .= sprintf $CAND_READS, $chr, $output, $input;
      if(exists $options->{'badloci'}) {
        $command .= ' -e '.$options->{'badloci'};
      }

      PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

      #
      ## The rest is auto-magical
      PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
    }
  }
  return 1;
}

sub merge_select_orphans {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    ## build command for this index
    #

    my $sample = sanitised_sample_from_bam($input);
    my $inputdir = File::Spec->catdir($tmp, 'candidates', $sample);

    my $command = "$^X ";
    $command .= _which('TraFiC_candOrph.pl');
    $command .= sprintf $CAND_ORPH, $inputdir;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub split_candidate_reads {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my @inputs = ($options->{'tumour'}, $options->{'normal'});
  my $iter = 1;
  for my $input(@inputs) {
    next if($iter++ != $index); # skip to the relevant input in the list

    ## build command for this index
    #

    my $sample = sanitised_sample_from_bam($input);
    my $inputdir = File::Spec->catdir($tmp, 'candidates', $sample);
    my $outdir = File::Spec->catdir($tmp, 'split');
    make_path($outdir) unless(-e $outdir);

    my $command = "$^X ";
    $command .= _which('TraFiC_candSplit.pl');
    $command .= sprintf $CAND_SPLIT, File::Spec->catfile($outdir, $sample.'.'), $inputdir;

    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);

    #
    ## The rest is auto-magical
    PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  }
  return 1;
}

sub repeat_mask_te {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my $infile = $options->{'split_files'}->[$index-1];

  # build command

 	my $rm_outdir = File::Spec->catdir($tmp, 'rm', $index);
	if(-e $rm_outdir){
		remove_tree($rm_outdir); # always cleanup any failed attempt first
	}
	make_path($rm_outdir);

  my $n_name = $options->{'normal_name'};
  my $t_name = $options->{'tumour_name'};
  my $file_sample = fileparse($infile);

  my $divergence;
  # files are name sorted so index won't line up
	if($file_sample =~ m/^$n_name[.][[:digit:]]+$/) {
		$divergence = $options->{'ndiv'};
	}
	elsif($file_sample =~ m/^$t_name[.][[:digit:]]+$/) {
	  $divergence = $options->{'tdiv'};
	}
	else {
	  die "Unexpected file name: $file_sample\n";
	}

	my $repeatmasker;
	if(defined $options->{'rm'}) {
	  $repeatmasker = $options->{'rm'};
	}
	else {
	  $repeatmasker = _which('RepeatMasker');
	}
	die "Unable to find RepeatMasker binary, try setting '-rm'" unless(defined $repeatmasker);

  my $command = "cd $rm_outdir;"; # you can set -dir, but tmp files written to cwd
	$command .= $repeatmasker;
	$command .= " -no_is -nolow";
	$command .= " -lib $options->{maskdb}";
	$command .= " -e $options->{engine}";
	$command .= " -norna $options->{accuracy}";
	$command .= " -frag 2500000";
	$command .= " -div $divergence"; # should this be configurable?
	$command .= " -dir $rm_outdir";
	$command .= " $infile";

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);

  my $basefile = "$rm_outdir/$file_sample.*";
  system('rm', '-f', "$basefile.cat.gz", "$basefile.log" , "$basefile.masked", "$basefile.tbl");

  return 1;
}

sub cluster_te {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

  my $sample;
  if($index == 1) {
    $sample = $options->{'tumour_name'};
  }
  elsif($index == 2) {
    $sample = $options->{'normal_name'};
  }

  my $cluster_outdir = File::Spec->catdir($tmp, 'clusters', $sample);
	if(-e $cluster_outdir){
		remove_tree($cluster_outdir); # always cleanup any failed attempt first
	}
	make_path($cluster_outdir);

  my $command = "$^X ";
  $command .= _which('TraFiC_cluster.pl');
  $command .= sprintf $CLUSTER,
                      $tmp, $sample, # repeatmaster.out
                      $tmp, $sample, # repeatmasker in modified fasta
                      $options->{'support'},
                      $cluster_outdir;
  $command .= ' -p' if($sample eq $options->{'tumour_name'});


  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
  return 1;
}

sub filter_results {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $clusters_tumour = File::Spec->catdir($tmp, 'clusters', $options->{'tumour_name'});
  my $clusters_normal = File::Spec->catdir($tmp, 'clusters', $options->{'normal_name'});

  my $np;
  if(exists $options->{'germline'}) {
    $np = $options->{'germline'}
  }
  else {
    $np = File::Spec->catfile($tmp, 'faked_np.txt');
    open my $NPFH, '>', $np;
    print $NPFH "#\n";
    close $NPFH;
  }

 	my $command = "$^X ";
  $command .= _which('TraFiC_filter.pl');
	$command .= " -o $options->{outdir}";
	$command .= ' -i '.File::Spec->catfile($clusters_tumour, 'reciprocal_clusters.txt');
	$command .= ' -f '.File::Spec->catfile($clusters_normal, 'pos_clusters.txt');
	$command .= ' -f '.File::Spec->catfile($clusters_normal, 'neg_clusters.txt');
	$command .= " -g $np";
	if(exists $options->{'known'}) {
	  $command .= " -t $options->{known}";
	  $command .= " -d $options->{div}";
  }
	$command .= " -fm $options->{filtmin}";
	$command .= " -gm $options->{germmin}";

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	unless(exists $options->{'germline'}) {
	  unlink $np;
	}

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return 1;
}

### Utility methods

sub prepare {
  my $options = shift;
  $options->{'tumour_name'} = (PCAP::Bam::sample_name($options->{'tumour'}))[0];
  $options->{'normal_name'} = (PCAP::Bam::sample_name($options->{'normal'}))[0];
  return 1;
}

sub sanitised_sample_from_bam {
  my $sample = (PCAP::Bam::sample_name(shift))[0];
  $sample =~ s/[^a-z0-9_-]/_/ig; # sanitise sample name
  return $sample;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  return $path;
}

sub get_sequences {
  my $bam_file = shift;
  my $bam = PCAP::Bam->new($bam_file);
  my @ordered_sq = sort @{$bam->header_sq};
  return \@ordered_sq;
}

sub split_files {
  my $tmp_folder = shift;
  my @files;
  my $indir = File::Spec->catdir($tmp_folder, 'split');
  if(-e $indir) {
    opendir(my $dh, $indir) || die "Unable to read from $indir: $OS_ERROR";
    while(my $thing = readdir $dh) {
      next if($thing =~ m/^[.]/);
      push @files, "$indir/$thing" if($thing =~ m/.*[.][[:digit:]]+$/);
    }
    closedir($dh);
  }
  @files = sort @files; # make sure these are always in same order
  return \@files;
}

1;
