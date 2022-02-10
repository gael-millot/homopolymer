[//]: # "#to make links in gitlab: example with racon https://github.com/isovic/racon"
[//]: # "tricks in markdown: https://openclassrooms.com/fr/courses/1304236-redigez-en-markdown"


[![python](https://img.shields.io/badge/code-Python-blue?style=plastic)](https://www.python.org/)
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

Return homopolymers info, including the largest homopolymer, per DNA sequence in a batch of DNA sequences


<br /><br />
## CONTENT

| File or folder | Description |
| --- | --- |
| **homopolymer.py** | file that can be executed using a CLI (command line interface) |
| **dataset** | Folder containing some datasets than can be used as examples |
| **example_of_result** | Folder containing examples of result obtained with the dataset |


<br /><br />
## HOW TO RUN

### local terminal


`  python3 homopolymer.py ./dataset/integrases.fasta ./result.tsv  `

arg0: .py script<br />
arg1: fasta file input<br />
arg2: tsv file output<br />
<br /><br />

### Using a cluster

Start with:

`  alias python3='module load Python/3.6.0 ; python3'  `

Then run as for the local terminal

<br /><br />
## OUTPUT

A table with the following columns<br /><br />
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

If several longest homopolymers in a sequence, results are semi-colon separated in each cell.

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

### v2.0

1) New features included in the result table


### v1.0

1) Everything



