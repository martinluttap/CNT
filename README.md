# CNT: A Semi-Automatic Translator from CWL to Nextflow for Genomics Workflows. 

CNT is a semi-automatic translator which performs both tool-level and workflow-level translation of CWL files into Nextflow. 
CNT utilizes [https://github.com/common-workflow-language/cwl-utils](cwl-utils), a library of classes for manipulating CWL files, to transform them into an intermediate representation (e.g. JSON) then proceed to rearrange it according to [https://nextflow-io.github.io/patterns/](patterns) devised by Nextflow's community.

## Installation

## Usage

## Development

## Release

## Citation

CNT was published in IEEE BIBE 2023. If you use CNT in any of your work, please consider citing the publication! This will allow us to track how useful CNT is, and also identify potential improvements.   
```
@inproceedings{cnt.bibe23,
 author = {Martin L. Putra and
            In Kee Kim and
            Robert L. Grossman and
            Haryadi S. Gunawi},
 booktitle = {The 23rd IEEE International Conference on BioInformatics and BioEngineering (BIBE)},
 title = {CNT: Semi-Automatic Translation from CWL to Nextflow for Genomic Workflows},
 year = {2023}
}
```