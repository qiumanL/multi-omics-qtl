import sys
gene=sys.argv[1]
qtl=sys.argv[2]
qtl=open(gene+'.'+qtl)
lines=qtl.readlines()
dic={}
for index, eachline in enumerate(lines):
	eachlist=eachline.strip().split('\t')
	snp=eachlist[1]
	dic[(index+1)]=snp
#print(dic[0])
cs_data=open('cs.summary.txt')
cs_data.readline()
cs_dic={}
for line in cs_data:
	cs_lists=line.split(' ')
	cs_index=cs_lists[0]
	cs_snp_index=cs_lists[4].split(',')
	cs_dic[cs_index]=cs_snp_index
output=open('cs.snp.txt','w')
for each_i in cs_dic:
	for each in cs_dic[each_i]:
		output.write(dic[int(each)]+'\n')

#for line in qtl:
#	list=line.split('\t')
	
