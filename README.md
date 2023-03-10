# PROtein-binDIng-enerGY-predicts-antibody-antigen-interaction-interfaces
## Project: “Characterizing therapeutic antibody interactions”.
### - Mined interaction descriptors:
 + Used Biopython to download, clean up the pdbs and separated the download files into com.pdb (complex), rec.pdb (antibody), and lig.pdb.
 + Counted frequency of amino acids type in each of the CDR sequences.
 + Used Prodigy to estimate binding energy (dG) between antibody and antigen for all the given pdbID complexes
### - Expanded the interaction networks used in a pipeline of a graphML-based prediction models.

## What is PRODIGY / Binding Affinity Prediction Repository?

Understanding the structural characteristics of protein-protein interactions is crucial for developing therapies and for understanding biological processes and disorders. Accurately predicting the binding strength for a certain protein-protein complex is a crucial component of this. In this project, we used PROtein binDIng enerGY prediction (PRODIGY)  a web service that determines the binding affinity of protein-protein complexes based on their three-dimensional (3D) structure. Our straightforward yet incredibly accurate prediction model, based on intermolecular interactions and characteristics extracted from non-interface surfaces, is implemented by the PRODIGY server.
