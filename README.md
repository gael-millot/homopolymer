[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"


[![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/)
&nbsp;
[![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses)


<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNING](#warning)
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
## WARNING

The algorithm works as if it splits the input sequence according to homopolymers and then returns the info.

Example with the input sequence ATTTAAGCGGG:
<br />
A
<br />
TTT
<br />
AA
<br />
G
<br />
C
<br />
GGG


<br /><br />
## CONTENT

| File or folder | Description |
| --- | --- |
| **homopolymer.nf** | File that can be executed using a CLI (command line interface)
| **nextflow.config** | Parameter settings for the homopolymer.nf file |
| **dataset** | Folder containing some datasets than can be used as examples |
| **example_of_result** | Folder containing examples of result obtained with the dataset |


<br /><br />
## INPUT

A fasta file


<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL2](https://github.com/nextflow-io/nextflow)<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu<br />
[Singularity/apptainer](https://github.com/apptainer/apptainer)<br />

<br /><br />
### 2. Local running (personal computer)


#### 2.1. homopolymer.nf file in the personal computer

- Mount a server if required:

<pre>
DRIVE="Z"
sudo mkdir /mnt/z
sudo mount -t drvfs $DRIVE: /mnt/z
</pre>

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `homopolymer.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the homopolymer.nf and nextflow.config files are (example: \\wsl$\Ubuntu-20.04\home\gael):

<pre>
nextflow run homopolymer.nf -c nextflow.config
</pre>

with -c to specify the name of the config file used.

<br /><br />
#### 2.3. homopolymer.nf file in the public gitlab repository

Run the following command from where you want the results:

<pre>
nextflow run -hub pasteur gmillot/homopolymer -r v1.0.0
</pre>

<br /><br />
### 3. Distant running (example with the Pasteur cluster)

#### 3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

<pre>
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/homopolymer" # where the bin folder of the homopolymer.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export SINGU_CONF=apptainer/1.1.5
export SINGU_CONF_AFTER=bin/singularity # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${SINGU_CONF}/${SINGU_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER}"
cd ${EXEC_PATH}
# chmod 755 ${EXEC_PATH}/bin/*.* # not required if no bin folder
module load ${JAVA_CONF} ${SINGU_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF}
</pre>

<br /><br />
#### 3.2. homopolymer.nf file in a cluster folder

Modify the second line of the code below, and run from where the homopolymer.nf and nextflow.config files are (which has been set thanks to the EXEC_PATH variable above):

<pre>
HOME_INI=$HOME
HOME="${ZEUSHOME}/homopolymer/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/homopolymer/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} homopolymer.nf -c nextflow.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
#### 3.3. homopolymer.nf file in the public gitlab repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

<pre>
VERSION="v1.0"
HOME_INI=$HOME
HOME="${ZEUSHOME}/homopolymer/" # $HOME changed to allow the creation of .nextflow into /$ZEUSHOME/homopolymer/, for instance. See NFX_HOME in the nextflow software script
trap '' SIGINT
nextflow run --modules ${MODULES} -hub pasteur gmillot/homopolymer -r $VERSION -c $HOME/nextflow.config
HOME=$HOME_INI
trap SIGINT
</pre>

<br /><br />
### 4. Error messages and solutions

#### Message 1
```
Unknown error accessing project `gmillot/homopolymer` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/homopolymer
```

Purge using:
<pre>
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
</pre>

#### Message 2
```
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Fhomopolymer
```

Contact Gael Millot (distant repository is not public).

#### Message 3

```
permission denied
```

Use chmod to change the user rights.


<br /><br />
## OUTPUT


**report.html** report of the analysis

**reports** folder containing all the reports of the different processes as well as the **nextflow.config** file used

**figures** folder containing all the figures in the **report.html** in the .png format

**files** folder containing the following files:

- homopol_summary.tsv

Warning: columns takes into account the **min_length** parameter in the **nextflow.config** file, meaning that the results are only for homopolymer lengths equal or above **min_length**. But, the two **homopol_obs_distrib** and **homopol_theo_distrib** columns provides the distribution of all the polymers, whatever **min_length**.
<br /><br />
| Column | Description |
| --- | --- |
| **name** | name of the input sequence of the batch |
| **seq_length** | nb of bases in the input sequence |
| **max_homopol_size** | number of times the nucleotide is repeated in the longest homopolymer (i.e., length). If several longest homopolymers in the input sequence, results are semi-colon separated in each cell |
| **nucleotide** | nucleotide of the homopolymer of **max_homopol_size**. If several longest homopolymers in the input sequence, results are semi-colon separated in each cell |
| **starting_position** | position of the first nucleotide of the homopolymer of **max_homopol_size**. If several longest homopolymers in the input sequence, results are semi-colon separated in each cell |
| **relative_position** | relative position of the **starting_position** value in the input sequence when the first base of the sequence is 0 and the last one is 1. The formula used is y = 0 if **seq_length** is 1 and y = (starting_position - 1) / seq_length otherwise, to get y between 0 and (starting_position - 1) / seq_length (i.e., not 1). If several longest homopolymers in the input sequence, results are semi-colon separated in each cell |
| **nb** | number of consecutive homopolymers in the input sequence |
| **mean_size** | average homopolymer size among the number of consecutive homopolymers in the input sequence (not considering the homopolymers below the **min_length** parameter in the **nextflow.config** file, meaning that the mean is computed only on the length of the considered homopolymers, not using the whole input sequence length) |
| **homopol_obs_distrib** | number of homopol of size 1, 2, ..., n (semi-colon separator). Warning: the **min_length** parameter in the **nextflow.config** is ignored |
| **homopol_theo_distrib** | number of homopol of size 1, 2, ..., n (semi-colon separator). Warning: the **min_length** parameter in the **nextflow.config** is ignored |

<br /><br />

- boxplot_stat.tsv

From **homopol_obs_distrib** and **homopol_theo_distrib** columns of the homopol_summary.tsv file

| Column | Description |
| --- | --- |
| **length** | homopolymer length |
| **freq** | frequency (sum of all the homopolymer numbers of size **length** in all the input sequences |
| **kind** | observed or random (theoretical) homopolymers |

<br /><br />

- scatterplot_stat.tsv

From **homopol_obs_distrib** and **homopol_theo_distrib** columns of the homopol_summary.tsv file

| Column | Description |
| --- | --- |
| **length** | homopolymer length |
| **kind** | observed or random (theoretical) homopolymers |
| **mean** | frequency mean along all the sequences of the batch (sum of the number of homopolymers of size **categ** in each input sequence) / (number of sequences) |
| **sd** | frequency standard deviation along all the sequences of the batch (sd of the corresponding  **mean**) |
| **CI95.inf** | 95% lower Confidence Interval of the **mean**, according to the normal law(**mean**, **sd**) |
| **CI95.sup** | 5% upper Confidence Interval of the **mean**, according to the normal law(**mean**, **sd**) |

<br /><br />

- t_test.tsv

the t test table displayed in the report.html file

| Column | Description |
| --- | --- |
| **length** | homopolymer length |
| **obs.mean** | mean of the observed homopolymers in the batch of sequences |
| **theo.mean** | mean of the random homopolymers |
| **obs.sd** | standard deviation of the observed homopolymers in the batch of sequences |
| **theo.sd** | standard deviation of the random homopolymers |
| **df** | degree of freedom of the t test |
| **t** | t test statistics |
| **p.value** | p value |
| **BH.adj.p.value** | Benjamini Hochberg adjusted p values along all the t tests performed |

<br /><br />

- graph_stat.RData

.RData file containing all the objects used to make the report.html, that can be reused if necessary.

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

### v4.1

Boxplots modified


### v4.0

Completemy rewritten


### v3.2

1) Minimum length of homopolymer added as parameter, among other things


### v3.1

1) Many things improved


### v3.0

1) Completely modified. Now the file is a nextflow and outputs include tables ,graphs and stats.


### v2.0

1) New features included in the result table


### v1.0

1) Everything



