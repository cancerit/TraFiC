package Sanger::CGP::TraFiC::Parser;

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
use Carp qw(croak);
use List::Util qw(first);
use File::Which qw(which);
use Capture::Tiny qw(capture);

use Const::Fast qw(const);

use Sanger::CGP::TraFiC;

const my $SINGLE_END  => 'SINGLE_END';
const my $INTER_CHROM => 'INTER_CHROM'; # between 2 chromosomes
const my $ABERRANT    => 'ABERRANT'; # within a chromosome INTRA (not using as name as easy to miss similarity)

# sam field indexes
const my $QNAME => 0;
const my $FLAG  => 1;
const my $RNAME => 2;
const my $POS   => 3;
const my $MAPQ  => 4;
const my $RNEXT => 6;
const my $SEQ   => 9;

const my $N_LIMIT => q{NNNNN};

# sam flag values
const my $PROPER_PAIR   => 0x0002;
const my $READ_UNMAPPED => 0x0004;
const my $MATE_UNMAPPED => 0x0008;
const my $READ_REVCOMP  => 0x0010; # will only be testing the mapped read
const my $READ_ONE      => 0x0040; # assuming paired reads only

sub new {
    my ($class) = @_;
    my $self = { };
    bless $self, $class;

    return $self;
}

sub grab_bam_header {
  my $bam = shift;
  my $command = which('samtools');
  die "Can't find 'samtools' in path\n" unless($command);
  $command .= ' view -H '.$bam;
  my ($sti,$sto,$err) = capture { system($command); };
  my $cleaned = q{};
  for my $line(split /\n/, $sti) {
    $cleaned .= $line."\n" if($line =~ m/^[@](HD|SQ)/);
  }
  return $cleaned;
}

sub set_output {
  my ($self, $root_path, $bam, $chr) = @_;
  croak "$root_path does not exist." unless(-e $root_path);
  croak "$root_path is not a directory." unless(-d $root_path);
  croak "$root_path cannot be written to." unless(-w $root_path);
  my $outfile;
  my $tmpfile;
  if(defined $chr) {
    if($chr eq '*') {
      $outfile = "$root_path/candidates.orphans.fa";
      $tmpfile = "$root_path/biobambamtmp.orphans.tmp";
      $self->{'merged_orph'} = "$root_path/merged_orphans.bam";
    }
    else {
      $outfile = "$root_path/$chr.candidates.fa";
      $tmpfile = "$root_path/$chr.biobambamtmp.tmp";
    }
    if(defined $bam) {
    my $orphans = "$root_path/$chr.o.sam";
      open my $orph, '>', $orphans;
      print $orph grab_bam_header($bam);
      $self->{'ORPHAN_FH'} = $orph;
    }
  }
  else {
    $outfile = "$root_path/candidates.fa";
    $tmpfile = "$root_path/biobambamtmp.tmp";
  }
  open my $se, '>', $outfile;
  $self->{'CANDIDATE_FH'} = $se;

  $self->{'biobambam_tmpfile'} = $tmpfile;
  return 1;
}

sub process_orphans {
  my ($self, $indir) = @_;
  my $file_list = $self->find_orphan_files($indir);
  my $command = sprintf 'bamsort level=1 disablevalidation=0 inputformat=sam tmpfile=%s O=%s I=',
                        $self->{'biobambam_tmpfile'},
                        $self->{'merged_orph'};
  my $full_command = $command.join(q{ I=}, @{$file_list} );
  system($full_command) == 0 or die "Processing failed: $full_command\n";

  $self->process_files([$self->{'merged_orph'}]);
  unlink $self->{'merged_orph'}; # contained within this process so clean it up
  return 1;
}

sub find_orphan_files {
  my ($self, $indir) = @_;
  my @orphans;
  opendir(my $dh, $indir) || die "Unable to read from $indir: $OS_ERROR";
  while(my $thing = readdir $dh) {
    next if($thing =~ m/^[.]/);
    push @orphans, "$indir/$thing" if($thing =~ m/.*[.]o[.]sam$/);
  }
  closedir($dh);
  return \@orphans;
}

sub process_files {
  my ($self, $file_list, $chr) = @_;
  my $out_fh = $self->{'CANDIDATE_FH'};
  for my $in(@{$file_list}) {
    warn "Processing: $in\n";
    my ($file_type) = $in =~ m/[.](sam|bam|cram)$/;
    my $command = sprintf 'bamcollate2 exclude=SECONDARY,QCFAIL,DUP,SUPPLEMENTARY collate=1 outputformat=sam T=%s inputformat=%s filename=%s classes=F,F2%s%s',
                          $self->{'biobambam_tmpfile'},
                          $file_type,
                          $in,
                          ((defined $chr) ? ',O,O2' : q{}),
                          ((defined $chr) ? " ranges=$chr" : q{});

    open my $process, '-|', $command;
    my ($read1, $read2);
    MAIN: while($read1 = <$process>) {
      next if((substr $read1,0,1) eq '@');
      $read2 = <$process>;
      chomp ($read1, $read2);
      my @r1 = split /\t/, $read1;
      my @r2 = split /\t/, $read2;
      while($r1[0] ne $r2[0]) {
        print {$self->{'ORPHAN_FH'}} $read1,"\n";
        $read1 = $read2;
        @r1 = @r2;
        $read2 = <$process>;
        unless($read2) {
          # the old read2 has not been output yet just moved to $read1
          print {$self->{'ORPHAN_FH'}} $read1,"\n";
          last MAIN;
        }
        chomp $read2;
        @r2 = split /\t/, $read2;
      }
      my ($class, $fasta) = @{$self->pair_to_candidate(\@r1, \@r2)};
      next unless(defined $class);
      print $out_fh ${$fasta};
#      last if($. > 1_000_000);
    }
    close $process;
  }
  return 1;
}

sub pair_to_candidate{
  my ($self, $read1, $read2) = @_;
  # make assumptions that flags are correct
  my @r1 = @{$read1};
  my @r2 = @{$read2};
  # both ends unmapped discard
  return [] if(($r1[$FLAG] & $READ_UNMAPPED) && ($r2[$FLAG] & $READ_UNMAPPED));
  # if proper pair return, paranoid version so check each end is mapped
  return [] if(($r1[$FLAG] & $PROPER_PAIR) && !($r1[$FLAG] & $READ_UNMAPPED) && !($r2[$FLAG] & $READ_UNMAPPED));

  my $swap = 0; # want r1 to always be the mapped read for ease of use
  my $class;
  if((!($r1[$FLAG] & $READ_UNMAPPED) && ($r2[$FLAG] & $READ_UNMAPPED)) || (!($r2[$FLAG] & $READ_UNMAPPED) && ($r1[$FLAG] & $READ_UNMAPPED)) ) {
    $class = $SINGLE_END;
    $swap = 1 if($r1[$FLAG] & $READ_UNMAPPED);
  }
  else { # both ends mapped
    return [] if(($r1[$MAPQ] == 0 && $r2[$MAPQ] == 0) || $r1[$MAPQ] == $r2[$MAPQ]); # can't determine which end should be anchor
    $swap = 1 if($r1[$MAPQ] < $r2[$MAPQ]);
    if($r1[$RNAME] eq $r2[$RNAME]) {
      $class = $ABERRANT;
    }
    else {
      $class = $INTER_CHROM;
    }
  }

  my $fasta = $self->record_out(\@r1, \@r2, $swap, $class);
  return [] unless(defined $fasta);
  return [$class, $fasta];
}

sub record_out {
  my ($self, $r1, $r2, $swap, $class) = @_;
  # returns undef on fail and scalar ref on pass
  ($r2, $r1) = ($r1,$r2) if($swap == 1);
  my $fasta;
  return $fasta if($r1->[$MAPQ] == 0); # reads with non-specific anchors are of little use
  return $fasta unless(index($r2->[$SEQ], $N_LIMIT) == -1);

  $fasta = join "\t", ( ">$r1->[$QNAME]",
                        $r1->[$RNAME],
                        $r1->[$POS],
                        $class,
                        ($r1->[$FLAG] & $READ_REVCOMP) ? '-' : '+');
  $fasta .= "$RS$r2->[$SEQ]$RS";
  return \$fasta;
}

sub DESTROY {
  my $self = shift;
  close $self->{'CANDIDATE_FH'} if(defined $self->{'CANDIDATE_FH'});
  close $self->{'ORPHAN_FH'} if(defined $self->{'ORPHAN_FH'});
  return 1;
}

1;

__END__

=head1 NAME

Sanger::CGP::TraFiC::Parser

=head1 SYNOPSIS

 use Sanger::CGP::TraFiC::Parser;
 my $parser = Sanger::CGP::TraFiC::Parser->new;
 $parser->set_output('output/folder);
 $parser->process_files($list_of_bams);

=head2 Constants

=over

=item SINGLE_END

Only one end of the pair was mapped

=item INTER_CHROM

Both ends of pair where mapped to different chromosomes

=item ABERRANT

Both ends mapped to same chromosome but not as a proper pair (i.e. intRA chromosomal)

=back

=head1 METHODS

=head3 set_output

=over

=item Description

Set the root path that output should be sent to.  Generates file-handles
to the 6 files required for processing.

 $parser->set_output('/some/existing/path');

=back

=head3 process_files

=over

=item Description

Find and output candidate reads in the supplied files.

 $parser->process_files([filelist, ...]);
  or
 $parser->process_files([singlefile], readgroup_id);

=item Notes

When I<readgroup_id> is supplied only reads within the group will be
handled.  Using this inappropriately will result in no output for files
that don't have the specified group.  There will be no warnings when this occurs.

=back

=head3 pair_to_candidate

=over

=item Description

Takes a pair of reads and determines if they are suitable
for processing as aberrant.

 my ($class, $detail, $fasta) = $parser->pair_to_candidate($sam_r1, $sam_r2);
 next unless(defined $class);

=item Returns

\@ of the following values

 $class : Type of aberrant read FastA record (see L</Constants>)
 $fasta : FastA record ready for printing

=item Notes

When the pair of reads input into this function are not found to be candidates
an empty array ref is returned.

=back
