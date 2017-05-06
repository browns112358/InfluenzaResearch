import numpy as np
import Bio
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
import glob
from pprint import pprint
from Bio import SeqUtils
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.model_selection import cross_val_score
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
import time
import mutate
from Bio.SeqRecord import SeqRecord
from scipy.sparse import lil_matrix

startTime = time.time()
File = "bucket/new.fasta"
#importing our data
mySeq = SeqIO.parse(File, "fasta")
with open('sample.txt') as m:
	records = motifs.parse(m, 'dreme')

#mutation
Mset =set()
for m in records:
	Mset.add(m.instances[0])
for ii in range(8):
	Mset = mutate.mutate(Mset, 0.3)
#Mset = mutate.RandMotif(10000, 8)
print "Number of motifs: ", len(Mset)
Mlist =list(Mset)
#initialize matrixs
#452009 in sequences.fasta, 590527 in complete.fasta
#LENGTH=364944
LENGTH =100000
#feature_matrix = np.zeros([LENGTH, len(Mset)], dtype = np.int8)
feature_matrix = lil_matrix((LENGTH, len(Mset)), dtype = np.int8)
label = np.zeros([LENGTH,], dtype = np.int8)
allSpecies ={}
#look for motifs in each sequence and grab the labels
s=-1
big=0
for seq in mySeq:
	k = seq.description.rfind(":")
	species = seq.description[k+1:].lower()
	#testing
	m = species.rfind("|")
	if ((species[m+1:] == "no genes for sequence") or (len(seq.seq)<10)) or (s >=(LENGTH-1)):
		#Dont use this data
		continue
	if species not in allSpecies:
		allSpecies[species]=0
	allSpecies[species] = allSpecies[species]+1
	s += 1
	if (species =='human'):
		label[s]=1
	else:
		label[s]=-1
	m = -1
	for mymotif in Mlist:
		m += 1
		test = SeqUtils.nt_rna_search(repr(seq.seq), mymotif)
		if (len(test)>2):
			big +=1
		if (len(test) > 1):
			feature_matrix[s,m] = test[1]+1
		else:
			feature_matrix[s,m]=0
print 'number with more than 1 match ', big
print 'number of sequences ', s+1
#Save matrixs
#np.savetxt('features.out', feature_matrix, delimiter=',', fmt = '%1.0f')
#np.savetxt('labels.out', label, delimiter=',', fmt = '%1.0f')

#Learn an AdaBoost model and get accuracy results
Ada = AdaBoostClassifier(n_estimators=500)
Ada.fit(feature_matrix, label)
score = Ada.score(feature_matrix, label)
print "AdaBoost accuracy: ", score*100, "%"
#scores = cross_val_score(Ada, feature_matrix, label, cv=5)
#print "AdaBoost accuracy: ", np.mean(scores)*100, "%"

#showing off stats
imp = Ada.feature_importances_
ind = np.argsort(imp)[::1]
print "------------------feature ranking--------------------"
for f in range(1000):
	matrix = feature_matrix.getcol(ind[f]).todense()
	if (matrix.sum()>0):
		ave_pos = np.true_divide(matrix.sum(),(matrix!=0).sum())
		stdev = np.nanstd(np.where(matrix!=0,matrix,np.nan))
		print("%d. feature %d (Average position = %f +/- %f)      %s" % (f+1, ind[f], ave_pos, stdev, Mlist[ind[f]]))
pred = Ada.predict(feature_matrix)
eval =  np.logical_and(np.not_equal(pred,label), label ==1)
print "-----------------Misclassified data----------------"
gen = SeqIO.parse(File, 'fasta')
ii =0
fraqs=[]
for g in gen:
	if (ii >= LENGTH):
		break
	if eval[ii] :
		record = SeqRecord(g.seq,g.description, '', '')
		fraqs.append(record)
		print("Description : ", g.description)
		print("Sequence : ", repr(g._seq))
		print ("========================================")
	ii += 1
SeqIO.write(fraqs, "misclassified.fasta", "fasta")

sum=0
print "-------------------Data used---------------------"
for k, v in allSpecies.items():
	print k, ":", v
	sum += v


print 'function took: ', (time.time()-startTime)/60, 'minutes'
