"""
name: Sheethal Umesh Nagalakshmi
lang: python3
file: nucleosome_types.py
desc: Divides nucleosomal DNA sequence patterns into 4 types

"""

import pandas as pd
import numpy as np

def type_nuc(infile,outfile1,outfile2,outfile3,outfile4):

	"""
		:param infile1: Input file containing ww/ss DNA sequence patterns 
		:param outfile1: Output file containing type1 nuclesome DNA sequence patterns
		:param outfile2: Output file containing type2 nuclesome DNA sequence patterns
		:param outfile3: Output file containing type3 nuclesome DNA sequence patterns
		:param outfile4: Output file containing type4 nuclesome DNA sequence patterns
		:return: dataframe containing type1, type2, type3 and type4 nuclesome DNA sequence 
		patterns
	
	"""
	# Reading a file containing the nucleosomal ww/ss DNA sequence patterns 
	data = pd.read_csv(infile,sep = '\t')
	data.drop(data.columns[0], axis=1,inplace = True)
	
	# Assigning the co-efficient value that refers to the dinucleotide positions in the 
	# minor and major groove binding sites 
	coef = 1.125
	
	# Creating a numpy array to store the 4 types of nucleosomal DNA sequence patterns
	type1 = ([])
	type2 = ([])
	type3 = ([])
	type4 = ([])
	data2 = pd.DataFrame()
	
	# Assigning the dinucleotide positions in minor-groove binding sites
	minor = data.T.iloc[pd.np.r_[4:8,14:18,25:29,36:40,46:50,56:60]]
	# Assigning the dinucleotide positions in major-groove binding sites
	major = data.T.iloc[pd.np.r_[9:14,20:23,31:34,41:45,51:55,61:65]]
	
	# Finding the sum of ww and ss dinucleotides in nucleosomes
	data2['minor_ww'] = minor.isin(['ww']).sum(axis=0)
	data2['minor_ss'] = minor.isin(['ss']).sum(axis=0)
	data2['major_ww'] = major.isin(['ww']).sum(axis=0)
	data2['major_ss'] = major.isin(['ss']).sum(axis=0)
	
	# Assigning the nucleosomes as type 1, type 2, type 3 and type 4 nucleosomes based on  
	# ww and ss dinucleotides in minor and major-groove binding sites
	 
	type1 = data[(data2['minor_ww'] >= (data2['major_ww']*coef)) & 
	(data2['minor_ss'] <= (data2['major_ss']*coef))]
	
	type2 = data[(data2['minor_ww'] >= (data2['major_ww']*coef)) & 
	(data2['minor_ss'] > (data2['major_ss']*coef))]
	
	type3 = data[(data2['minor_ww'] < (data2['major_ww']*coef)) & 
	(data2['minor_ss'] <= (data2['major_ss']*coef))]
	
	type4 = data[(data2['minor_ww'] < (data2['major_ww']*coef)) & 
	(data2['minor_ss'] > (data2['major_ss']*coef))]
	
	# Writing the 4 types of nucleosome DNA sequence patterns to output files
	pd.DataFrame(type1).to_csv(outfile1,sep = '\t',index = False)
	pd.DataFrame(type2).to_csv(outfile2,sep = '\t',index = False)
	pd.DataFrame(type3).to_csv(outfile3,sep = '\t',index = False)
	pd.DataFrame(type4).to_csv(outfile4,sep = '\t',index = False)

def main():
	infiles = ['1/all_combined_ww_ss.txt','2/all_combined_ww_ss.txt',]
	outfiles1 = ['1/type1.txt','2/type1.txt']
	outfiles2 = ['1/type2.txt','2/type2.txt']
	outfiles3 = ['1/type3.txt','2/type3.txt']
	outfiles4 = ['1/type4.txt','2/type4.txt']
	l = len(infiles)
	for i in range(l):
		infile = infiles[i]
		outfile1 = outfiles1[i]
		outfile2 = outfiles2[i]
		outfile3 = outfiles3[i]
		outfile4 = outfiles4[i]
		type_nuc(infile,outfile1,outfile2,outfile3,outfile4)
	
		
if __name__ == '__main__':
	main()
	



	
