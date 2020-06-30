from os import listdir, environ
import pandas as pd
OUT_DIR = environ['OUT_DIR']
# merge all the dosage files
merged_file = None
for dosage_file in listdir(OUT_DIR+"/condout/ldclump_dosage"):
    if "dosage_" in dosage_file:
        dosage = pd.read_csv(OUT_DIR+"/condout/ldclump_dosage/"+dosage_file)
        if merged_file is None:
            merged_file = dosage
        else: merged_file = pd.concat([merged_file, dosage])

merged_file.to_csv(OUT_DIR+"/condout/ldclump_dosage/dosage.ld", index=None)
