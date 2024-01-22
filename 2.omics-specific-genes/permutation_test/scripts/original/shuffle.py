import random
import sys
import numpy as np
import pandas as pd
seed_num=int(sys.argv[1])
exp_raw=open('braingvex.predict_predicted_expression.txt2')
genelist=exp_raw.readline()
#exp=np.genfromtxt('braingvex.predict_predicted_expression.txt2',delimiter="\t",skip_header=1,usecols=np.arange(0,200))
exp=pd.read_csv('braingvex.predict_predicted_expression.txt2',sep="\t",header=0,index_col=0)
samplelist=exp.index.values
print(samplelist)
np.random.seed(seed_num)
r=np.random.permutation(185)
exp=exp.astype('str')
####
exp_shuffled=open('braingvex.predict_predicted_expression.txt2.shuffled','w')
sample_order=open('sample.order','w')
shuffled_sample_order=open('sample.order.shuffled','w')
exp_shuffled.write(genelist)
b=0
for i in r:
	print(i-1)
	exp_shuffled.write(samplelist[b]+'\t'+'\t'.join(exp.iloc[i-1])+'\n')
	sample_order.write(samplelist[b]+'\n')
	shuffled_sample_order.write(str(i)+'\t'+samplelist[i-1]+'\n')
	b=b+1
