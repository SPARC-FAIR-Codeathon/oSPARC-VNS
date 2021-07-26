
# o²S²PARC: Reusable Models of Visceral Nerve Stimulation

### What is Visceral Nerve Stimulation (VNS)?
Electrical stimulation of nerves innervating the viscera (organs of the body) can be utilized as a neuromodulation therapy for a wide range of medical conditions affecting various organs. The vagus nerve (10th cranial nerve) in particular has been the target for therapies for epilepsy, inflammation, heart failure, and more. 

### What's the problem or gap with VNS?
Many VNS therapies are not currently sufficiently reliable enough to justify widespread clinical adoption: experimental therapeutic approaches which are extremely effective in some patients have no or even negative effect in other patients. Related, many promising therapeutic approaches in animal models fail to translate to human neuroanatomies. Finally, each VNS therapy must be customised to each patient, which is time-consuming for clinicians and can be disheartening for patients when invasive procedures fail to show immediate dramatic effect. 

### How are we addressing this problem?
Computational modelling of the nerve-electrode interface can shed insight into how different approaches will translate between animal models and human patients, as well as to predict how variations in individual patients’ neuroanatomy will influence outcomes. 
This project consisted of deploying an existing MATLAB-based nerve modelling pipeline, [ViNERS](https://gitlab.unimelb.edu.au/lab-keast-osborne-release/ViNERS/-/wikis/home) on oSPARC and building a web front-end interface to visualise model inputs and outputs. 

### Project aims


SPARC has produced datasets of significant interest (e.g. Quantified Morphology of the [Human](https://sparc.science/datasets/65?type=dataset) and 
[Rat Vagus](https://sparc.science/datasets/60?type=dataset) Nerves) for the development and translation of electroceuticals from animal models to patients. 
These datasets use [standardized file formats](https://www.biorxiv.org/content/10.1101/2020.09.22.306670v1.full) ([MBF-XML](https://www.incf.org/mbf-file-format-v-40), NeuroLucida, MBF Biosciences) which 
can be re-used to build [high-fidelity models for visceral nerve stimulation and recording](https://gitlab.unimelb.edu.au/lab-keast-osborne-release/ViNERS) to predict how device designs and electrical stimulation parameters 
need to be scaled from pre-clinical rodent models to human patients. 

Using [o²S²PARC](https://osparc.io/), this project aims to build tools and example workflows to enable interoperability and re-usability of these and other datasets to predict nerve stimulation thresholds and responses for measured nerve geometries, enabling in-silico evaluation of new devices and stimulation protocols and assessment of how responses change across species. This project entails:

- deploying existing model modules (4 modules deployed)
- developing clear interfaces for nerve and electrode array specification,
- providing tools are usable by clinicians without specialized programming training, 
- and expand the library of neural interface device architectures to include human-relevant devices

### Conceptual workflow
![workflow](https://user-images.githubusercontent.com/63089004/126910461-b3c2c36a-7a74-410c-bb86-5435dcb72baa.png)

### For Scientists
- [Using the web interface](https://github.com/SPARC-FAIR-Codeathon/oSPARC-VNS/wiki/How-to-use-the-web-interface)
- [Simulation inputs](link)
- [Simulation outputs](link)

### For Developers
- [frontend](link)
- [services](link)

### Outcomes

This project entailed deploying visceral nerve modelling work so that other researchers can utilize these models (and components from these models) to predict performance for new nerves and neural interfaces. The delivered model on o²S²PARC is modular, not monolithic: each component completes a common, discrete modelling task such as finite-element mesh generation from nerve cross-sections. By modularising this code, we hope to help other researchers save time and reduction of reduplication of effort. Our web front-end enables non-experts to construct models for new devices and nerves based on data uploaded to the SPARC DRC, expanding the accessibility of modelling work.

For a detailed list of outcomes, see [outcomes](-/wiki/outcomes). 

