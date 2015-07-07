#! /usr/bin python

def editGTFscaffold(gtf_file):
	with open(gtf_file, 'r') as file:
		for line in file:
			editMe = line.split('\t')[0]
			num = editMe.split('_')[1]
			nZeros = 5 - len(num) 
			pasteMe = 'scaffold' + '0'*nZeros + num
			line.split('\t')[0] = pasteMe