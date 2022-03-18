[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"


[![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/)
&nbsp;
[![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses)


<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [CONTENT](#content)
   - [HOW TO RUN](#how-to-run)
   - [OUTPUT](#output)
   - [VERSIONS](#versions)
   - [LICENCE](#licence)
   - [CITATION](#citation)
   - [CREDITS](#credits)
   - [ACKNOWLEDGEMENTS](#Acknowledgements)
   - [WHAT'S NEW IN](#what's-new-in)


<br /><br />
## AIM

Return homopolymers info per DNA sequence in a batch of DNA sequences, as well as statistics about homopolymers for this batch.


<br /><br />
## CONTENT

| File or folder | Description |
| --- | --- |
| **main.nf** | file that can be executed using a CLI (command line interface)
| **nextflow.config** | parameter settings for the main.nf file |
| **dataset** | Folder containing some datasets than can be used as examples |
| **example_of_result** | Folder containing examples of result obtained with the dataset |


<br /><br />
## HOW TO RUN

See Protocol 136 (ask me).


### If error message

If an error message appears, like:
```
Unknown error accessing project `gmillot/homopolymer` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/homopolymer
```
Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```


### Using the committed version on gitlab:

1) Create the scm file:

```bash
providers {
    pasteur {
        server = 'https://gitlab.pasteur.fr'
        platform = 'gitlab'
    }
}
```

And save it as 'scm' in the .nextflow folder. For instance in:
\\wsl$\Ubuntu-20.04\home\gael\.nextflow

Warning: ssh key must be set for gitlab, to be able to use this procedure (see protocol 44).


2) Mount a server if required:

```bash
DRIVE="C"
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
```

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like
```
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
```


3) Then run the following command from here \\wsl$\Ubuntu-20.04\home\gael:

```bash
nextflow run -hub pasteur gmillot/homopolymer -r v1.0.0
```

If an error message appears, like:
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Fhomopolymer
```
Make the distant repo public

If an error message appears, like:

```
permission denied
```

See chmod in protocol 44.


### Using a cluster

Start with:

```bash
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/homopolymer" # where the bin folder of the main.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export SINGU_CONF=singularity/3.8.3
export SINGU_CONF_AFTER=bin/singularity # on maestro
export GIT_CONF=git/2.25.0
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${SINGU_CONF}/${SINGU_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER},${CONF_BEFORE}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER}"
# cd ${EXEC_PATH} # not required when using the gitlab repo to run the script
# chmod 755 ${EXEC_PATH}/bin/*.* # not required when using the gitlab repo to run the script
module load ${JAVA_CONF} ${SINGU_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF}

```

Then run:

```bash
# distant main.nf file
HOME="$ZEUSHOME/homopolymer/" ; nextflow run --modules ${MODULES} -hub pasteur gmillot/homopolymer -r v7.10.0 -c $HOME/nextflow.config ; HOME="/pasteur/appa/homes/gmillot/"

# local main.nf file ($HOME changed to allow the creation of .nextflow into /$ZEUSHOME/homopolymer/. See NFX_HOME in the nextflow soft script)
HOME="$ZEUSHOME/homopolymer/" ; nextflow run --modules ${MODULES} main.nf ; HOME="/pasteur/appa/homes/gmillot/"
```


## OUTPUT


**report.html** report of the analysis

**reports** folder containing all the reports of the different processes as well as the **nextflow.config** file used

**figures** folder containing all the figures in the **report.html** in the .png format

**files** folder containing the following tables

<br /><br />
*<FILE_NAME>*_homopol_summary.tsv

| Column | Description |
| --- | --- |
| **name** | name of the sequence |
| **seq_length** | nb of bases in the sequence |
| **nucleotide** | nucleotide of the homopolymer |
| **starting_position** | position of the first nucleotide of the homopolymer of max size |
| **relative_position** | relative position of the starting_position value when the first base of the sequence is 0 and the last one is 1. The formula used is y = (starting_position - 1) / (seq_length - max_size) to get 0 <= y <= 1) |
| **max_size** | number of times the nucleotide is repeated in the homopolymer (homopolymer length) |
| **nb** | number homopolymers in the sequence (including homopolyers of size 1) |
| **mean_size** | average homopolymer size in the sequence (including homopolyers of size 1) |
| **homopol_obs_distrib** | number of homopol of size 1, 2, ..., n (semi-colon separator) |
| **homopol_theo_distrib** | number of homopol of size 1, 2, ..., n (semi-colon separator) |

If several longest homopolymers in a sequence, results are semi-colon separated in each cell.

<br /><br />
barplot_stat.tsv

| Column | Description |
| --- | --- |
| **length** | homopolymer length |
| **freq** | frequency |
| **kind** | observed or random (theoretical) homopolymers |

<br /><br />
scatterplot_stat.tsv

| Column | Description |
| --- | --- |
| **length** | homopolymer length |
| **kind** | observed or random (theoretical) homopolymers |
| **mean** | frequency mean along all the sequences of the batch |
| **sd** | frequency standard deviation along all the sequences of the batch |
| **CI95.inf** | 95% lower Confidence Interval of the mean |
| **CI95.sup** | 5% upper Confidence Interval of the mean |


<br /><br />
## VERSIONS


The different releases are tagged [here](https://gitlab.pasteur.fr/gmillot/homopolymer/-/tags)

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses.

<br /><br />
## CITATION


Not yet published

<br /><br />
## CREDITS


[Gael A. Millot](https://gitlab.pasteur.fr/gmillot), Hub-CBD, Institut Pasteur, Paris, France

<br /><br />
## ACKNOWLEDGEMENTS


Yoann Dufresne, Hub-CBD, Institut Pasteur, Paris

The mentioned softwares and packages developers & maintainers

Gitlab developers

<br /><br />
## WHAT'S NEW IN

### v3.0

1) Completely modified. Now the file is a nextflow and outputs include tables ,graphs and stats.


### v2.0

1) New features included in the result table


### v1.0

1) Everything



