# LigandNet workflow
This work is tested on python3.6

* User can this workflow generate the machine learning model to predict the state of a ligand either active (that binds to a particular protein) or decoy (inactive).
* Users can use their own actives and decoys (in SMILES or sdf file) or will be asked uniprot id of a protein to extract actives from pharos database and search for decoys in our already generated decoys sets. If decoys are not found, they can be generated from DUD--E website (or using decoys finder 2) and come back for model generation.


- actives folder contains the (user's sample) actives for a protein in .sdf and .txt files
- decoys folder contains the (user's sample) decoys for a protein in .sdf and .txt files
- decoys_database contains the sample decoys if the user does not have the decoys in sdf or txt files

- actives and decoys for a protein in text file should contain the header with SMILES for smiles column.

## Usage:
	./run.sh
