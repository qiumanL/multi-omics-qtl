import sys
sys.path.append("/scratch/cgg/jiany17/software/MetaXcan/software/")

import metax
from metax import PredictionModel
model = PredictionModel.load_model("/scratch/cgg/jiany17/psychENCODE/predixcan.Dan/output/WGS.db")
# dir(model)    # view items in the model
# model.weights.head()
f = open("/scratch/cgg/jiany17/psychENCODE/predixcan.Dan/output/WGS.db.weights","w")
f.write("rsid\tgene\tweight\teffect_allele\tnon_effect_allele\n")
for i in range(len(model.weights)):
    f.write("\t".join(list(map(str,model.weights.loc[i]))))
    f.write("\n")

f.close()

