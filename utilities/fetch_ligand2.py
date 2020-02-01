#**********************************************************#
# Author: Govinda KC                                       #
# UTEP, Compuational Science                               #
# Last modified: 1/15/2020                                 #
#***********************************************************
# A script to retrieve ligand structures (as SMILES format) for
# a given target.
# Usage: python fetch_ligand.py TARGET
# TARGET -> gene (ITK) or UniProt Accession id (Q08881)

import csv
import json
from tqdm import tqdm
from urllib.request import urlopen
import sys
import os
#LIGANDS_DIR = './'+sys.argv[1]

#if not os.path.isdir(LIGANDS_DIR): os.mkdir(LIGANDS_DIR)

pharos = "https://pharos.nih.gov/idg/api/v1/targets/"

class Pharos_Data:
    def __init__(self, target, act_path):
        self.target = target
        self.act_path = act_path

    def fetch_ligand(self):
        # List for storing ligand data with header for csv file
        lig_data = [("SMILES", "name", "assay_type", "activity_value")]
        # Fetch target
        req = urlopen(pharos + self.target)
        target_dict = json.loads(req.read())
        # Retrieve ligand links for this target
        print("Fetching ligands for target: ", self.target)
        req = urlopen(pharos + '{0}'.format(target_dict['id']) +
                      "/links(kind=ix.idg.models.Ligand)")
        link = json.loads(req.read())
        if type(link) is dict:
            link = [link]  # Fixes bug when target has only one ligand.
        for l in tqdm(link):
            name = ""
            assay = "N/A"
            value = "N/A"
            # Extract ligand name
            for p in l['properties']:
                if p['label'] == "IDG Ligand" or p['label'] == "IDG Drug":
                    name = p['term']
                    break
            # Extract ligand activity
            for p in l['properties']:
                if p['label'] in {"Ki", "Kd", "IC50", "EC50", "AC50", "Potency"}:
                    assay = p['label']
                    value = p['numval']
                    break
            # Extract SMILES string
            req = urlopen(l['href'] + "/properties(label=CHEMBL Canonical SMILES)")
            try:
                ligand = json.loads(req.read())[0]
            except:
                print('not downloaded')
            if not ligand:
            # Skip molecule if no SMILES was found
                continue
            # Add information to list
            
            lig_data.append((ligand['text'], name, assay, value))
        # Save all data to a CSV file
        #with open(LIGANDS_DIR + "/" + target + ".txt", 'w') as out_file:
        with open(self.act_path+'/'+self.target + ".txt", 'w') as out_file:
            writer = csv.writer(out_file, delimiter=',')
            writer.writerows(lig_data)


if __name__ == "__main__":

    #LIGANDS_DIR = './'+sys.argv[1]

    #if not os.path.isdir(LIGANDS_DIR): os.mkdir(LIGANDS_DIR)

    pharos = "https://pharos.nih.gov/idg/api/v1/targets/"

    target = sys.argv[1]
    act_path = sys.argv[2]
    fl = Pharos_Data(target, act_path)
    fl.fetch_ligand()
