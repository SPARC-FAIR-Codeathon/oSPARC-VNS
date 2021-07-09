
# o²S²PARC: Reusable Models of Visceral Nerve Stimulation

## Project Aims

SPARC has produced datasets of significant interest (e.g. Quantified Morphology of the [Human](https://sparc.science/datasets/65?type=dataset) and 
[Rat Vagus](https://sparc.science/datasets/60?type=dataset) Nerves) for the development and translation of electroceuticals from animal models to patients. 
These datasets use [standardized file formats](https://www.biorxiv.org/content/10.1101/2020.09.22.306670v1.full) ([MBF-XML](https://www.incf.org/mbf-file-format-v-40), NeuroLucida, MBF Biosciences) which 
can be re-used to build [high-fidelity models for visceral nerve stimulation and recording](https://gitlab.unimelb.edu.au/lab-keast-osborne-release/ViNERS) to predict how device designs and electrical stimulation parameters 
need to be scaled from pre-clinical rodent models to human patients. 

Using [o²S²PARC](https://osparc.io/), this project aims to build tools and example workflows to enable interoperability and re-usability of these and other datasets to
predict nerve stimulation thresholds and responses for measured nerve geometries, enabling in-silico evaluation of new devices and stimulation protocols and assessment of how responses change across species. 
This project entails:
- deploying existing model modules, 
- developing clear interfaces for stimulus and electrode array specification, 
- developing clear visualizations for recruitment of different fibre classes
- providing tools are usable by clinicians without specialized programming training. 
- expand the library of neural interface device architectures offered to include other device architectures developed by SPARC-funded researchers. (e.g. OT2OD024908). This project will advance capabilities with using the cloud deployment and containerization tools, enabling more rapid deployment of improved models in the future, which is the current bottleneck for my activities aimed at deploying this work on o²S²PARC. 


## Relevence
This project entails deploying visceral nerve modelling work so that other researchers can utilize these models (and components from these models) to 
predict performance for new nerves and neural interfaces. The delivered model on o²S²PARC will be modular, not monolithic: 
the components delivered will enable common, discrete modelling tasks (such as finite-element mesh generation from nerve cross-sections) to be 
accomplished by other researchers, accomplishing significant time savings and reduction of reduplication of effort. Clear, well-documented graphical 
user interfaces will enable non-experts to construct models for new devices and nerves based on data uploaded to the SPARC DRC, expanding the
accessibility of modelling work. 

