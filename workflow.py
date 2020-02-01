#!/usr/bin/env python
# coding: utf-8

#  Ligandnet workflow

#**************************************
# Govinda KC                          #
# UTEP, Computational Science         #
# Last modified: 1/25/20              #
# *************************************

# Import libraries

import warnings
import os, sys, json, glob
sys.path.append('utilities')
from train2 import Train
from fetch_ligand2 import Pharos_Data
from utility import FeatureGenerator # for features generation of  txt file
from utility2 import FeatureGenerator2 # for features generation of sdf file
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from sklearn import metrics
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.utils.class_weight import compute_class_weight
import joblib
from sklearn.neural_network import MLPClassifier

class Run_Workflow:
    def __init__(self, actives, decoys):
        self.actives = actives
        self.decoys = decoys
        self.results = dict()
        
    def get_fingerprints(self,smiles):
        try:
            fg = FeatureGenerator(smiles)
            features = fg.toTPATF()
            return features
        except Exception as e: print(e)
    
    def get_models(self):
        # Get features at first!
        if not self.fp_generation():
            print('Error: features extraction failed!')
            return
        try:
            t = Train(self.actives_x, self.decoys_x)
            t.train_models()
        except Exception as e: print(e)

    def fp_generation(self):
        # Fingerprint generation
        print('Pleae wait! Fingerprints are getting generated......')
        if self.decoys[-4:] == '.sdf' and self.actives[-4:] == '.sdf':
            # Get fingerprints for actives
            self.actives_x = self.sdf_fp_active()
            # Get fingerprints for decoys
            self.decoys_x = self.sdf_fp_decoy()
            return True
        elif self.decoys[-4:] == '.sdf':
            df = pd.read_csv(self.actives)
#             df = pd.read_csv(open(self.actives,'rU'))#, encoding='utf-8', engine='c')
            # Get fingerprints for actives
            df['tpatf'] = df.SMILES.apply(self.get_fingerprints)
            self.actives_x = np.array([f for f in df.tpatf.values], dtype = np.float32)
            # Get fingerprints for decoys
            self.decoys_x = self.sdf_fp_decoy()
            return True
        else:
            df = pd.read_csv(self.actives)
            df2 = pd.read_csv(self.decoys)
#             df = pd.read_csv(open(self.actives,'rU'))#, encoding='utf-8', engine='c')
#             df2 = pd.read_csv(open(self.decoys, 'rU'))#, encoding='utf-8', engine='c')
            # Get fingerprints for actives
            df['tpatf'] = df.SMILES.apply(self.get_fingerprints)
            # Get fingerprints for decoys
            df2['tpatf'] = df2.SMILES.apply(self.get_fingerprints)
            # numpy arrays
            self.actives_x = np.array([f for f in df.tpatf.values], dtype = np.float32)
            self.decoys_x = np.array([f for f in df2.tpatf.values], dtype = np.float32)
            return True
        return False

    def sdf_fp_decoy(self):
        try:
            fg2 = FeatureGenerator2(self.decoys)
            feat_decoy = fg2.sepTPATF()
            return feat_decoy
        except Exception as e: print(e) 

    def sdf_fp_active(self):
        try:
            fg2 = FeatureGenerator2(self.actives)
            feat_active = fg2.sepTPATF()
            return feat_active
        except Exception as e: print(e)


# If users have their own actives and decoys
def actives_decoys():
    active_file = input("Uniprot id of the file? Example: P07948  \n")   
    active_file = active_file.strip()
    print('Looking for active and decoy files....')
    
    # active in .txt
    actives = main_path+'actives/'+active_file+'.txt'
    
    if not os.path.isfile(actives):
        # active in .sdf
        actives = main_path+'actives/'+active_file+'.sdf'
    
    # decoy in .txt..
    decoys = main_path+'decoys/'+"decoys_" + active_file +".txt"
    
    if not os.path.isfile(decoys):
        # decoy in .sdf..
        decoys = main_path+'decoys/'+ "decoys_" +active_file+".sdf"
    if os.path.isfile(actives) and os.path.isfile(decoys):
        print('Actives and Decoys are found!')
    return actives, decoys

# Searches decoys in our database for give active file (Uniprot id)
def actives_bt_not_decoys():
    
    active_file = input("Uniprot id of the file? Example: P07948 \n")
    active_file = active_file.strip()
    actives = main_path+'actives/'+active_file+'.txt'
    if not os.path.isfile(actives):
        actives = main_path+'actives/'+active_file+'.sdf'
    
    # Path for decoys database
    decoys_database = '../decoys_database'

    print('Searching decoys .....')
    if not os.path.isfile(os.path.join(decoys_database, active_file+".sdf")):
        print("Decoys are not found, exiting! Look for decoys in DUDE website and come back!")
        sys.exit(1)

#     decoys = os.path.join(decoys_database, active_file+".txt")
    decoys = os.path.join(decoys_database, "decoys_" +active_file+".sdf")
    if os.path.isfile(actives) and os.path.isfile(decoys):
        print('Actives and decoys are extracted!')
    return actives, decoys
    

# when user does not have both active and decoys
# actives can be downloaded from pharos.nih.gov and decoys from our database.
def no_actives_and_decoys():
    
    active_file = input("Uniprot id of the file? Example: P07948 \n")
    active_file = active_file.strip()
    active_dir = main_path+'/'+ "actives"
    pdata = Pharos_Data(active_file, active_dir )
    
    print('Actives for a given protein are getting downloaded from Pharos website!')
    pdata.fetch_ligand()
    
    actives = main_path+'actives/'+active_file+'.txt'
    
    print('Searching decoys .....')
    
    decoys_database = '../decoys_database/'
    
    if not os.path.isfile(os.path.join(decoys_database, "decoys_" +active_file+".sdf")):
        print("Decoys are not found, exiting! Look for decoys in DUDE website and come back!")
        sys.exit(1)
    decoys = os.path.join(decoys_database, active_file+".sdf")
    
    if os.path.isfile(actives) and os.path.isfile(decoys):
        print('Actives and decoys are extracted!')
    return actives, decoys


# Start here
def start_workflow():
    
    print('Actives and decoys should either be in sdf file or text file (with header "SMILES" for txt files!)')
    print('ACTIVES AND DECOYS FILE NAMES SHOULD BE LIKE THAT: P07948.txt(or .sdf) and decoys_P07948.txt (or .sdf) ')
    print('PLEASE, MAKE SURE YOU HAVE FOLDERS "actives" and "decoys"')
    print('DO YOU HAVE "actives" and "decoys" FOLDERS? Type y for Yes and n for No!')
    check = input()
    if check != 'y':
        print('Exiting...')
        sys.exit(1)
    print("Do you have actives? Please type y for Yes and n for No !")
    answer1 = input()
          
    print("Do you have decoys? Please type y for Yes and n for No !")
    answer2 = input()

    if answer1 == 'y' and answer2 == 'y':
        actives, decoys = actives_decoys()
        rw = Run_Workflow(actives, decoys)
        rw.get_models()
    elif answer1 == 'y' and answer2 == 'n':
        actives, decoys = actives_bt_not_decoys()
        rw = Run_Workflow(actives, decoys)
        rw.get_models()
    elif answer1 == 'n' and answer2 == 'n':
        actives, decoys = no_actives_and_decoys()
        rw = Run_Workflow(actives, decoys)
        rw.get_models()
    else:
        print('Please provide the right information!. Exiting!')
        sys.exit(1)
        
if __name__ == '__main__':
    # Path for working directory
    print("Please, provide the path for working directory. Example: /Users/gvin/ligandnet_workflow/test_ligandnet/ \n")
    main_path = input()
    main_path = main_path.strip()
    os.chdir(main_path) 
    dirs = ["actives", "decoys"]
    for _dir in dirs:
        if not os.path.isdir(_dir): os.makedirs(_dir)
    if main_path[-1]!='/':
        main_path = main_path+'/'
    # Start Function
    start_workflow()

