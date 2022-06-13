# Mutational signatures

This is a semester project carried out at the Fellay lab at EPFL about mutational signatures.

## Authors

 Ester Simkova, under the supervision of Konstantin Popadin and Jacques Fellay
 
## The project

The aim of the project was to extend the mutagens tracing done by the software [SigProfilerExtractor:](https://github.com/AlexandrovLab/SigProfilerExtractor) to non-human and non-nuclear genomes. A pipeline was developed and assessed using an in silico mutagenesis experiment on the SARS-cov2 genome.
 
## Repository Organisation

- [source:](https://github.com/estersimkova/mutational_signatures/tree/main/source) folder which contains the source code
  - [in_silico_mutagenesis_positive_control.ipynb:](https://github.com/estersimkova/mutational_signatures/blob/main/source/in_silico_mutagenesis_positive_control.ipynb) the jupyter notebook explaining the whole pipeline with an example and instructions on how to run the code
  - [in_silico_mutagenesis.py:](https://github.com/estersimkova/mutational_signatures/blob/main/source/in_silico_mutagenesis.py) the python script to run the in silico mutagenesis experiment
  - [cosine_similarity_through_time.ipynb:](https://github.com/estersimkova/mutational_signatures/blob/main/source/cosine_similarity_through_time.ipynb) the jupyter notebook showing the effect of averaging of the trinucleotide content of the SARS-cov2 genome on the cosine similarity with the original signature used in the in silico mutagenesis experiment
- [results:](https://github.com/estersimkova/mutational_signatures/tree/main/results) folder which contains subfolders with the results of running the pipeline with different signatures and the results of the software SigProExtractor on them
- [data:](https://github.com/estersimkova/mutational_signatures/tree/main/data) folder with the input data to the pipeline; the SARS-cov2 reference sequence, the COSMIC signatures, the trinucleotide counts...
  
