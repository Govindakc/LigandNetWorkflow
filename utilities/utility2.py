
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import tempfile
import shutil
import sys
import pandas as pd
import numpy as np
class FeatureGenerator2:
    
    def __init__(self, sdf):
        self.sdf = sdf
        self.temp_dir = tempfile.mkdtemp()
    
    def sepTPATF(self):
        #sdf_file = open(main_path+pname+'.sdf')
        sdf_file = open(self.sdf)
        
        #output2 = open(self.temp_dir+'/'+'temp.csv', 'w')
        w = open(self.temp_dir+'/'+'temp.sdf', 'w') 

        
        feat_list = []     
        for line in sdf_file:
            w.write(line)
            if "$$$$" in line:
                w.close()
                if not os.path.isfile(os.path.join(self.temp_dir, "temp.sdf")): return None
                    
                try:
                    feat = self.toTPATF()
                    feat_list.append(feat)
                except:
                    pass
                w = open(self.temp_dir+'/'+'temp.sdf', 'w')
                
        decoys_fp = np.array(feat_list, dtype = np.float32)
        return decoys_fp
        self._cleanup()
        print('Temp File deleted')

    def toTPATF(self):
        features = []
        script_path = "../mayachemtools/bin/TopologicalPharmacophoreAtomTripletsFingerprints.pl"
        
        if not os.path.isfile(script_path): print(" Mayachemtools does not exist")        

        # Generate the sdf file
        #self.toSDF()
        # Now generate the TPATF features
        # Check if the sdf file exists
        if not os.path.isfile(os.path.join(self.temp_dir, "temp.sdf")):
            
            return None
        command = "perl " + script_path + " -r " + os.path.join(self.temp_dir, "temp") + " --AtomTripletsSetSizeToUse FixedSize -v ValuesString -o " + os.path.join(self.temp_dir, "temp.sdf")
        os.system(command)
        with open(os.path.join(self.temp_dir, "temp.csv"), 'r') as f:
            for line in f.readlines():
                if "Cmpd" in line:
                    line = line.split(';')[5].replace('"','')
                    #features = ','.join(line.split(' '))
                    features = [int(i) for i in line.split(" ")]

        # Clean up the temporary files
        #self._cleanup()
        return features
       
    def _cleanup(self):
        shutil.rmtree(self.temp_dir)


if __name__=='__main__':
    sdf_path = sys.argv[1]
    fgs = FeatureGenerator2(sdf_path)
    fgs.sepTPATF()

