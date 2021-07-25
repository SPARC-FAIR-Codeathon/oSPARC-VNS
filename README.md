
# oSPARC-VNS
o²S²PARC: Reusable Models of Visceral Nerve Stimulation

## What is Visceral Nerve Stimulation (VNS)?

Stimulating visceral nerves can be utalized as a theraputic neuromodulation thearpy for a wide range of medical conditions. Here we are specifially focusing on the Vagus nerve in the pelvis.
What's the problem or gap with VNS?
- There is poor interface design, stimulation selectivity, and low recording signal-to-noise ratio. We only have rat/animal models to go off of currently. We are trying to expand the animal model to a human model.
How do we fix this?
- New therapties are needed to advance the development of VNS in bioelectronic medicine. We have created a computational pipeline for modeling Visceral Nerve Ensemble Recroding & Stimulation (ViNERS) to accelerate the development and applications of VNS technologies.
ViNERS Work Flow
![image](https://user-images.githubusercontent.com/73490478/126045056-db34f078-50d8-4ab2-a7a0-b8a3fe68808a.png)
How to use ViNERS
1. we will step this through as we go
Installation & Dependencies
- i'll sparce these two out as we go and see where they fit best. But I'm just adding them here for now
Output
- we will email you your results!
Additonal / planned features
- all of the example readmes had these at the bottom. IDK what to do with this yet.
Refrence
- Calvin
Team
- All of us (edited) 


o²S²PARC: Reusable Models of Visceral Nerve Stimulation
Project Aims
SPARC has produced datasets of significant interest (e.g. Quantified Morphology of the Human and Rat Vagus Nerves) for the development and translation of electroceuticals from animal models to patients. These datasets use standardized file formats (MBF-XML, NeuroLucida, MBF Biosciences) which can be re-used to build high-fidelity models for visceral nerve stimulation and recording to predict how device designs and electrical stimulation parameters need to be scaled from pre-clinical rodent models to human patients.

Using o²S²PARC, this project aims to build tools and example workflows to enable interoperability and re-usability of these and other datasets to predict nerve stimulation thresholds and responses for measured nerve geometries, enabling in-silico evaluation of new devices and stimulation protocols and assessment of how responses change across species. This project entails:

deploying existing model modules,
developing clear interfaces for stimulus and electrode array specification,
developing clear visualizations for recruitment of different fibre classes
providing tools are usable by clinicians without specialized programming training.
expand the library of neural interface device architectures offered to include other device architectures developed by SPARC-funded researchers. (e.g. OT2OD024908). This project will advance capabilities with using the cloud deployment and containerization tools, enabling more rapid deployment of improved models in the future, which is the current bottleneck for my activities aimed at deploying this work on o²S²PARC.
Relevence
This project entails deploying visceral nerve modelling work so that other researchers can utilize these models (and components from these models) to predict performance for new nerves and neural interfaces. The delivered model on o²S²PARC will be modular, not monolithic: the components delivered will enable common, discrete modelling tasks (such as finite-element mesh generation from nerve cross-sections) to be accomplished by other researchers, accomplishing significant time savings and reduction of reduplication of effort. Clear, well-documented graphical user interfaces will enable non-experts to construct models for new devices and nerves based on data uploaded to the SPARC DRC, expanding the accessibility of modelling work.

