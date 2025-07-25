| Usage | Requirement |
| :--- | :--- |
| [![Nextflow](https://img.shields.io/badge/code-Nextflow-blue?style=plastic)](https://www.nextflow.io/) | [![Dependencies: Nextflow Version](https://img.shields.io/badge/Nextflow-v23.04.4.5881-blue?style=plastic)](https://github.com/nextflow-io/nextflow) |
| [![License: GPL-3.0](https://img.shields.io/badge/licence-GPL%20(%3E%3D3)-green?style=plastic)](https://www.gnu.org/licenses) | [![Dependencies: Apptainer Version](https://img.shields.io/badge/Apptainer-v1.2.3-blue?style=plastic)](https://github.com/apptainer/apptainer) |
| | [![Dependencies: Graphviz Version](https://img.shields.io/badge/Graphviz-v2.42.2-blue?style=plastic)](https://www.graphviz.org/download/) |

<br /><br />
## TABLE OF CONTENTS


   - [AIM](#aim)
   - [WARNINGS](#warnings)
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
## WARNINGS

The algorithm works as if it splits the input sequence according to homopolymers and then returns the info.

With the input sequence **ATTTAAGCGGG**, the homopolymers are:
<br />
**A**
<br />
**TTT**
<br />
**AA**
<br />
**G**
<br />
**C**
<br />
**GGG**


<br /><br />
## CONTENT
<br />

| Files and folder | Description |
| :--- | :--- |
| **main.nf** | File that can be executed using a linux terminal, a MacOS terminal or Windows 10 WSL2. |
| **nextflow.config** | Parameter settings for the *main.nf* file. Users have to open this file, set the desired settings and save these modifications before execution. |
| **bin folder** | Contains files required by the *main.nf* file. |
| **Licence.txt** | Licence of the release. |

<br /><br />
## INPUT
<br />

| Required files |
| :--- |
| A fasta file. |

<br />

The dataset used in the *nextflow.config* file, as example, is available at https://zenodo.org/records/10681460.

<br />

| File name | Description |
| :--- | :--- |
| **test.fasta** | Fasta file . Available [here](https://zenodo.org/records/10681460/files/test.fasta). |


<br /><br />
## HOW TO RUN

### 1. Prerequisite

Installation of:<br />
[nextflow DSL2](https://gael-millot.github.io/protocols/docs/Protocol%20152-rev0%20DSL2.html#_Toc159933761)<br />
[Graphviz](https://www.graphviz.org/download/), `sudo apt install graphviz` for Linux ubuntu<br />
[Apptainer](https://gael-millot.github.io/protocols/docs/Protocol%20135-rev0%20APPTAINER.html#_Toc160091693)<br />


### 2. Local running (personal computer)


####	2.1. *main.nf* file in the personal computer

- Mount a server if required:

```
DRIVE="Z" # change the letter to fit the correct drive
sudo mkdir /mnt/share
sudo mount -t drvfs $DRIVE: /mnt/share
```

Warning: if no mounting, it is possible that nextflow does nothing, or displays a message like:
<pre>
Launching `main.nf` [loving_morse] - revision: d5aabe528b
/mnt/share/Users
</pre>

- Run the following command from where the *main.nf* and *nextflow.config* files are (example: \\wsl$\Ubuntu-20.04\home\gael):

```
nextflow run main.nf -c nextflow.config
```

with -c to specify the name of the config file used.


#### 2.2.	*main.nf* file in the public git repository

Run the following command from where you want the results:

```
nextflow run gael-millot/homopolymer # github, or nextflow run http://github.com/gael-millot/homopolymer
nextflow run -hub pasteur gmillot/homopolymer -r v1.0.0 # gitlab
```


### 3. Distant running (example with the Pasteur cluster)

####	3.1. Pre-execution

Copy-paste this after having modified the EXEC_PATH variable:

```
EXEC_PATH="/pasteur/zeus/projets/p01/BioIT/gmillot/homopolymer" # where the bin folder of the main.nf script is located
export CONF_BEFORE=/opt/gensoft/exe # on maestro

export JAVA_CONF=java/13.0.2
export JAVA_CONF_AFTER=bin/java # on maestro
export APP_CONF=apptainer/1.3.5
export APP_CONF_AFTER=bin/apptainer # on maestro
export GIT_CONF=git/2.39.1
export GIT_CONF_AFTER=bin/git # on maestro
export GRAPHVIZ_CONF=graphviz/2.42.3
export GRAPHVIZ_CONF_AFTER=bin/graphviz # on maestro

MODULES="${CONF_BEFORE}/${JAVA_CONF}/${JAVA_CONF_AFTER},${CONF_BEFORE}/${APP_CONF}/${APP_CONF_AFTER},${CONF_BEFORE}/${GIT_CONF}/${GIT_CONF_AFTER}/${GRAPHVIZ_CONF}/${GRAPHVIZ_CONF_AFTER}"
cd ${EXEC_PATH}
chmod 755 ${EXEC_PATH}/bin/*.*
module load ${JAVA_CONF} ${APP_CONF} ${GIT_CONF} ${GRAPHVIZ_CONF}
```


####	3.2. *main.nf* file in a cluster folder

Modify the second line of the code below, and run from where the *main.nf* and *nextflow.config* files are (which has been set thanks to the EXEC_PATH variable above):

```
HOME_INI=$HOME
HOME="${HELIXHOME}/homopolymer/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/homopolymer/, for instance. See NFX_HOME in the nextflow software script
nextflow run main.nf -c nextflow.config
HOME=$HOME_INI
```


####	3.3. *main.nf* file in the public git repository

Modify the first and third lines of the code below, and run (results will be where the EXEC_PATH variable has been set above):

```
VERSION="v1.0"
HOME_INI=$HOME
HOME="${HELIXHOME}/homopolymer/" # $HOME changed to allow the creation of .nextflow into /$HELIXHOME/homopolymer/, for instance. See NFX_HOME in the nextflow software script
nextflow run gael-millot/homopolymer -r $VERSION -c $HOME/nextflow.config #github, or nextflow run http://github.com/gael-millot/homopolymer -r $VERSION -c $HOME/nextflow.config
nextflow run -hub pasteur gmillot/homopolymer -r $VERSION -c $HOME/nextflow.config # gitlab
HOME=$HOME_INI
```


### 4. Error messages and solutions

####	Message 1
<pre>
Unknown error accessing project `gmillot/homopolymer` -- Repository may be corrupted: /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot/homopolymer
</pre>

Purge using:
```
rm -rf /pasteur/sonic/homes/gmillot/.nextflow/assets/gmillot*
```

####	Message 2
<pre>
WARN: Cannot read project manifest -- Cause: Remote resource not found: https://gitlab.pasteur.fr/api/v4/projects/gmillot%2Fhomopolymer
</pre>

Contact Gael Millot (distant repository is not public).

####	Message 3
<pre>
permission denied
</pre>

Use chmod to change the user rights. Example linked to files in the bin folder: 
```
chmod 755 bin/*.*
```


<br /><br />
## OUTPUT


An example of results is present at this address: https://zenodo.org/records/10681595/files/homopolymer_2_1708386120.zip.

| homopolymer_<br />\<HOMOPOLYMER_MIN_LENGTH\>_<br />\<UNIQUE_ID\><br />folder | Description |
| :--- | :--- |
| **report.html** | Report of the analysis. |
| **reports** | Folder containing all the reports of the different processes, including the *nextflow.config* file used. |
| **figures** | Folder containing the graphs in png format that are used in the *report.html* file, as well as the corresponding svg vectorial files if needed. |
| **files** | Folder containing the following files:<br /><ul><li>**homopol_summary.tsv**<br />Column description: <br /><ul><li>Warning: columns takes into account the *min_length* parameter in the *nextflow.config* file, meaning that the results are only for homopolymer lengths equal or above *min_length*. But, the two *homopol_obs_distrib* and *homopol_theo_distrib* columns provides the distribution of all the polymers, whatever *min_length*.<br /></li><li>**name**: name of the input sequence of the batch.<br /></li><li>**seq_length**: nb of bases in the input sequence.<br /></li><li>**max_homopol_size**: number of times the nucleotide is repeated in the longest homopolymer (i.e., length). If several longest homopolymers in the input sequence, results are semi-colon separated in each cell.<br /></li><li>**nucleotide**: nucleotide of the homopolymer of *max_homopol_size*. If several longest homopolymers in the input sequence, results are semi-colon separated in each cell.<br /></li><li>**starting_position**: position of the first nucleotide of the homopolymer of *max_homopol_size*. If several longest homopolymers in the input sequence, results are semi-colon separated in each cell.<br /></li><li>**relative_position**: relative position of the *starting_position* value in the input sequence when the first base of the sequence is 0 and the last one is 1. The formula used is y = 0 if *seq_length* is 1 and y = (starting_position - 1) / seq_length otherwise, to get y between 0 and (starting_position - 1) / seq_length (i.e., not 1). If several longest homopolymers in the input sequence, results are semi-colon separated in each cell.<br /></li><li>**nb**: number of consecutive homopolymers in the input sequence.<br /></li><li>**mean_size**: average homopolymer size among the number of consecutive homopolymers in the input sequence (not considering the homopolymers below the *min_length* parameter in the *nextflow.config* file, meaning that the mean is computed only on the length of the considered homopolymers, not using the whole input sequence length).<br /></li><li>**homopol_obs_distrib**: number of homopol of size 1, 2, ..., n (semi-colon separator). Warning: the *min_length* parameter in the *nextflow.config* is ignored.<br /></li><li>**homopol_theo_distrib**: number of homopol of size 1, 2, ..., n (semi-colon separator). Warning: the *min_length* parameter in the *nextflow.config* is ignored.<br /></ul><li>**boxplot_stat.tsv**<br />Column description: <br /><ul><li>From *homopol_obs_distrib* and *homopol_theo_distrib* columns of the homopol_summary.tsv file.<br /></li><li>**length**: homopolymer length.<br /></li><li>**freq**: frequency (sum of all the homopolymer numbers of size *length* in all the input sequences.<br /></li><li>**kind**: observed or random (theoretical) homopolymers.<br /></ul><li>**scatterplot_stat.tsv**<br />Column description: <br /><ul><li>From *homopol_obs_distrib* and *homopol_theo_distrib* columns of the homopol_summary.tsv file.<br /></li><li>**length**: homopolymer length.<br /></li><li>**kind**: observed or random (theoretical) homopolymers.<br /></li><li>**mean**: frequency mean along all the sequences of the batch (sum of the number of homopolymers of size *categ* in each input sequence) / (number of sequences).<br /></li><li>**sd**: frequency standard deviation along all the sequences of the batch (sd of the corresponding  *mean*).<br /></li><li>**CI95.inf**: 95% lower Confidence Interval of the *mean*, according to the normal law(*mean*, *sd*).<br /></li><li>**CI95.sup**: 5% upper Confidence Interval of the *mean*, according to the normal law(*mean*, *sd*).<br /></ul><li>**plot_raw_values.tsv**<br />Column description: <br /><ul><li>From *homopol_obs_distrib* and *homopol_theo_distrib* columns of the homopol_summary.tsv file. <br /></li><li>Each line is a proportion of one polymer length in one sequence. The total nomber of rows should be: number of different homopolymer lengths x number of sequences x number of kind.<br /></li><li>**prop**: proportion of one polymer length in one sequence.<br /></li><li>**length**: homopolymer length.<br /></li><li>**gene_name**: name of the sequence.<br /></li><li>**kind**: observed or random (theoretical) homopolymers.<br /></ul><li>**t_test.tsv**<br />Column description: <br /><ul><li>The t test table displayed in the report.html file.<br /></li><li>**length**: homopolymer length.<br /></li><li>**obs.mean**: mean of the observed homopolymers in the batch of sequences.<br /></li><li>**theo.mean**: mean of the random homopolymers.<br /></li><li>**obs.sd**: standard deviation of the observed homopolymers in the batch of sequences.<br /></li><li>**theo.sd**: standard deviation of the random homopolymers.<br /></li><li>**df**: degree of freedom of the t test.<br /></li><li>**t**: t test statistics.<br /></li><li>**p.value**: p value.<br /></li><li>**BH.adj.p.value**: Benjamini Hochberg adjusted p values along all the t tests performed.<br /></ul><li>**graph_stat.RData**<br />.RData file containing all the objects used to make the report.html, that can be reused if necessary. |

<br /><br />
## VERSIONS


The different releases are tagged [here](https://github.com/gael-millot/homopolymer/tags).

<br /><br />
## LICENCE


This package of scripts can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
Distributed in the hope that it will be useful, but without any warranty; without even the implied warranty of merchandability or fitness for a particular purpose.
See the GNU General Public License for more details at https://www.gnu.org/licenses or in the Licence.txt attached file.

<br /><br />
## CITATION


Not yet published.

<br /><br />
## CREDITS


[Gael A. Millot](https://github.com/gael-millot), Hub, Institut Pasteur, Paris, France.

<br /><br />
## ACKNOWLEDGEMENTS


The developers & maintainers of the mentioned softwares and packages, including:

- [R](https://www.r-project.org/)
- [ggplot2](https://ggplot2.tidyverse.org/)
- [Nextflow](https://www.nextflow.io/)
- [Apptainer](https://apptainer.org/)
- [Docker](https://www.docker.com/)
- [Git](https://git-scm.com/)
- [Github](https://github.com/)
- [Gitlab](https://about.gitlab.com/)
- [Bash](https://www.gnu.org/software/bash/)
- [Ubuntu](https://ubuntu.com/)

Special acknowledgement to [Yoann Dufresne](https://github.com/yoann-dufresne), Hub, Institut Pasteur, Paris, France.

<br /><br />
## WHAT'S NEW IN

### 5.4

- In the nextflow.config file, upgrade singularity -> apptainer (real one).


### 5.3

- In the nextflow.config file, upgrade singularity -> apptainer.


### v5.2

- Bugs fixed in report.
- README improved for github.
- nextflow.config input from zenodo.


### v5.1

Bug fixed in report.


### v5.0

Nextflow DSL1 -> DSL2.


### v4.3

Plot_raw_values.tsv file added and boxplot_stat_log.tsv file modified.


### v4.2

Plot and html report modified.


### v4.1

Boxplots modified.


### v4.0

Completemy rewritten.


### v3.2

Minimum length of homopolymer added as parameter, among other things.


### v3.1

Many things improved.


### v3.0

Completely modified. Now the file is a nextflow and outputs include tables ,graphs and stats.


### v2.0

New features included in the result table.


### v1.0

Everything.



