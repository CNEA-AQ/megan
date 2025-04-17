import pandas as pd
import numpy as np

# Set file paths
base_dir = "./EF_input_250205/"
ef_file = base_dir + "EFv250205.csv"

# Read original files
EF = pd.read_csv(ef_file)

SpCrop = pd.read_csv(base_dir + "SpeciationCrop241129.csv")
SpHerb = pd.read_csv(base_dir + "SpeciationHerb241129.csv")
SpShru = pd.read_csv(base_dir + "SpeciationShrub241129.csv")
SpTree = pd.read_csv(base_dir + "SpeciationTree241129.csv").iloc[:, :3]

# Add vegetation type labels
SpCrop["Gtyp"] = "Crop"
SpHerb["Gtyp"] = "Herb"
SpShru["Gtyp"] = "Shrub"
SpTree["Gtyp"] = "Tree"

SpColumnNames = ["EcotypeID", "VegID", "SpFrac", "Gtyp"]
SpCrop.columns = SpColumnNames
SpHerb.columns = SpColumnNames
SpShru.columns = SpColumnNames
SpTree.columns = SpColumnNames

## Merge GrowthForm from EF with Tree data
#df = pd.merge(EF[["VegID", "GrowthForm"]], SpTree, on="VegID", how="right")
#
## Process needleleaf trees (Ntr)
#SpNtr = df[df["GrowthForm"] == "Ntr"].copy()
#NtEtFracs = SpNtr.groupby("EcotypeID")["SpFrac"].sum().reset_index(name="total")
#for _, row in NtEtFracs.iterrows():
#    k = max(row["total"], 1e-23)
#    et = row["EcotypeID"]
#    SpNtr.loc[SpNtr["EcotypeID"] == et, "SpFrac"] /= k
#
## Process broadleaf trees (Btr)
#SpBtr = df[df["GrowthForm"] == "Btr"].copy()
#BtEtFracs = SpBtr.groupby("EcotypeID")["SpFrac"].sum().reset_index(name="total")
#for _, row in BtEtFracs.iterrows():
#    k = max(row["total"], 1e-23)
#    et = row["EcotypeID"]
#    SpBtr.loc[SpBtr["EcotypeID"] == et, "SpFrac"] /= k
#
## Prepare output for Ntr and Btr
#SpNtr = SpNtr[SpColumnNames]
#SpNtr["Gtyp"] = "Ntr"
#SpBtr = SpBtr[SpColumnNames]
#SpBtr["Gtyp"] = "Btr"

# Combine all speciation data
SpDF = pd.concat([SpCrop, SpHerb, SpShru, SpTree], ignore_index=True)
SpDF.to_csv("SpeciationAll.csv", index=False)

# ====================================
# Calculate emission factors (EF) grouped by growth type (Gtyp) and EcotypeID

# Reload input files
EF = pd.read_csv(ef_file)
SP = pd.read_csv("SpeciationAll.csv")

# Merge EF and SP
M = pd.merge(EF, SP, on="VegID")

# Drop unused columns
drop_cols = ["References", "Comment", "X", "Family", "GenusGroup", "CommonName", "Type", "GrowthForm"]
M = M.drop(columns=[col for col in drop_cols if col in M.columns])

# Define EF and LDF columns
efs = [f"EF{i}" for i in range(1, 20)] + [f"LDF{i}" for i in range(3, 7)]

# Multiply by speciation fraction
for col in efs:
    if col in M.columns:
        M[col] = M[col] * M["SpFrac"]

# Group and sum by Gtyp and EcotypeID
DF = M.groupby(["Gtyp", "EcotypeID"])[efs].sum().reset_index()
#DF = M.groupby(["EcotypeID"])[efs].sum().reset_index()
DF = DF.sort_values("Gtyp")
DF = DF.round(3)

# Write output file (no header or index)
DF.to_csv("GtEFbyEcotype_growthform.csv", sep=",", index=False, header=False, float_format="%.3f")

