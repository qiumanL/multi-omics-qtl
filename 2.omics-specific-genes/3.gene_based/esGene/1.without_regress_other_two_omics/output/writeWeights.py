import sys
sys.path.append("MetaXcan-0.6.11/software/")

import metax
### new add
sys.path.append("MetaXcan-0.6.11/software/metax")

from metax import PredictionModel
model = PredictionModel.load_model("braingvex.db")
# dir(model)    # view items in the model
# model.weights.head()
f = open("braingvex.db.weights","w")
f.write("rsid\tgene\tweight\teffect_allele\tnon_effect_allele\n")
for i in range(len(model.weights)):
    f.write("\t".join(list(map(str,model.weights.loc[i]))))
    f.write("\n")
f.close()

f = open("braingvex.db.extra","w")
f.write("gene\tgene_name\tn_snps_in_model\tpred_perf_r2\tpred_perf_pval\tpred_perf_qval\n")
for i in range(len(model.extra)):
    f.write("\t".join(list(map(str,model.extra.loc[i]))))
    f.write("\n") 
f.close() 


