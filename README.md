LICENCE
=======
Copyright (c) 2015 Genome Research Ltd.

Author: Cancer Genome Project cgpit@sanger.ac.uk

This file is part of TraFiC.

TraFiC is free software: you can redistribute it and/or modify it under
the terms of the GNU Affero General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

1. The usage of a range of years within a copyright statement contained within
this distribution should be interpreted as being equivalent to a list of years
including the first and last year specified and all consecutive years between
them. For example, a copyright statement that reads 'Copyright (c) 2005, 2007-
2009, 2011-2012' should be interpreted as being identical to a statement that
reads 'Copyright (c) 2005, 2007, 2008, 2009, 2011, 2012' and a copyright
statement that reads "Copyright (c) 2005-2012' should be interpreted as being
identical to a statement that reads 'Copyright (c) 2005, 2006, 2007, 2008,
2009, 2010, 2011, 2012'."

TraFiC
======

[![DOI](https://zenodo.org/badge/10817/cancerit/TraFiC.svg)](http://dx.doi.org/10.5281/zenodo.17793)

TraFiC (Transposon Finder in Cancer) uses paired-end sequencing data for the identification of
somatic insertion of transposable elements.  The identification of somatic TEs is performed in
four steps:

1. Selection of candidate reads.
2. Transposable element masking.
3. Clustering and prediction of insertions sites.
4. Filtering of germline events.

## Interpretation of results

The final results are found in the `filtered.txt` file.  Columns are defined as follows:

* **CHR** - Chromosome of somatic insertion
* **FAMILY** - Transposon type
* **P\_L\_POS** - Left position of positive cluster supporting somatic insertion
* **P\_R\_POS** - Right position of positive cluster supporting somatic insertion
* **P\_TOTAL\_READS** - Total number of reads supporting somatic insertion in positive cluster
* **P\_SINGLE\_END\_COUNT** - Number of Single-end reads supporting somatic insertion in positive cluster
* **P\_INTER\_CHROM\_COUNT** - Number of Inter-chromosomal reads supporting somatic insertion in positive cluster
* **P\_ABERRANT\_COUNT** - Number of Intra-chromosomal reads supporting somatic insertion in positive cluster
* **P\_SINGLE\_END\_READS** - Names of Single-end reads supporting somatic insertion in positive cluster
* **P\_INTER\_CHROM\_READS** - Names of Inter-chromosomal reads supporting somatic insertion in positive cluster
* **P\_ABERRANT\_READS** - Names of Intra-chromosomal reads supporting somatic insertion in positive cluster
* **N\_L\_POS** - Left position of negative cluster supporting somatic insertion
* **N\_R\_POS** - Right position of negative cluster supporting somatic insertion
* **N\_TOTAL\_READS** - Total number of reads supporting somatic insertion in negative cluster
* **N\_SINGLE\_END\_COUNT** - Number of Single-end reads supporting somatic insertion in negative cluster
* **N\_INTER\_CHROM\_COUNT** - Number of Inter-chromosomal reads supporting somatic insertion in negative cluster
* **N\_ABERRANT\_COUNT** - Number of Intra-chromosomal reads supporting somatic insertion in negative cluster
* **N\_SINGLE\_END\_READS** - Names of Single-end reads supporting somatic insertion in negative cluster
* **N\_INTER\_CHROM\_READS** - Names of Inter-chromosomal reads supporting somatic insertion in negative cluster
* **N\_ABERRANT\_READS** - Names of Intra-chromosomal reads supporting somatic insertion in negative cluster

Events supported by at minimum of 5 reads in the both the positive and negative clusters are considered highly accurate.

---

###Dependencies/Install

Please install the following first:

* [PCAP-core](http://github.com/ICGC-TCGA-PanCancer/PCAP-core/releases)
* [RepeatMasker](http://www.repeatmasker.org/RMDownload.html)
    * We recommend the [RMBlast](http://www.repeatmasker.org/RMBlast.html) engine as it is easy to obtain (crossmatch is slightly faster)

Please see these for any child dependencies.

Once complete please run:

```
./setup.sh /some/install/location
```

This installs [Bedtools](http://bedtools.readthedocs.org/en/latest/) and the TraFiC codebase.

---

##Creating a release
####Preparation
* Commit/push all relevant changes.
* Pull a clean version of the repo and use this for the following steps.

####Cutting the release
**This should be carried out using hubflow tools**
1. Update `perl/lib/Sanger/CGP/TraFiC.pm` to the correct version.
2. Update `Changes` to show major items.
3. Run `./prerelease.sh`
4. Check all tests and coverage reports are acceptable.
5. Commit the updated docs tree and updated module/version.
6. Push commits.
7. Use the GitHub tools to add details to the generated tag.
8. Add relevant *.tar.gz for reference data to release.
