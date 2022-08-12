from Bio import AlignIO
import pandas as pd
import os
import glob

align = AlignIO.read('aligned.fa', 'fasta')
name = []
description = []
for record in align:
	name.append(record.name)
	description.append(record.description)

align2 = pd.DataFrame(align)

cuts = [0, 255, 804, 2718, 8553, 10053, 10971, 11841, 12090, 
	    12684, 13023, 13440, 16235, 18038, 19619, 20657, 21551, 
	    25383, 26219, 26471, 27190, 27386, 27758, 27886, 28258, 
	    29532, 29673]


#creates variables hhm inputs (1-27) (hmminput1 - hmminput27)
i = 1
dicts = {}
keys = range(len(cuts))
for i in keys:
	if cuts[i] !=29673:
		dicts[i] = pd.DataFrame(align2.iloc[:, cuts[i]:cuts[i+1]])
	if cuts[i] == 29673:
		dicts[i] = pd.DataFrame(align2.iloc[:, cuts[i]:])


for key in dicts:
	hmminput = dicts[key]
	seqs = []
	for row in hmminput.iterrows():
		with open('hmminput' + str(key) +'.fasta', 'a') as f:
			f.write('> ' + name[row[0]] + ' ' + description[row[0]] + '\n')
		row = row[1]
		seq = row.to_string(header=None, index=False).strip('\n')
		seq = seq.replace('\n', '')
		#seq = [','.join(ele.split()) for ele in seq]
		seqs.append(seq)
		with open('hmminput' + str(key) +'.fasta', 'a') as f:
			f.write(seq + '\n')

for filename in glob.glob('*.fasta'):
    cmd = 'hmmbuild out_'+ filename + '.hmm ' + filename
    print('building ' + filename + ' hmm')
    print(cmd)
    os.system(cmd)