import glob
import os
import sys
import csv
from collections import Counter

import pandas as pd
index_df = pd.read_csv(r'10ntR1_10ntR2.csv')
#print(index_df)

index_list = index_df['index'].tolist()
#print(index_list)

files_list = []
for file in glob.glob("*.fastq"):
	files_list.append(file)

for x in range(0,len(files_list)):
	file1 = files_list[x]
	name = file1[0:5]

	header1 = []
	read1_first = []
	read1_second = []
	qual1_first = []
	qual1_second = []
	header2 = []
	read2_first = []
	read2_second = []
	qual2_first = []
	qual2_second = []

	counter = 0
	with open(file1,'r') as tsvin:
		for row in tsvin:
			counter = counter + 1
			if counter == 1:
				header1.append(row[:-1])
			if counter == 2:
				read1_first.append(row[0:10])
				read1_second.append(row[30:40])
			if counter == 3:
				next
			if counter == 4:
				qual1_first.append(row[0:10])
				qual1_second.append(row[30:40])
				counter = 0

	output_file = open(name+".fastq", "w")

	for x in range(0,len(header1)):
		if 1 == 1: #header1[x][:-15] == header2[x][:-15]:
			#print(read1[x]+"CCTGTTCGTGTGGTATCGGT"+read2[x])
			if read1_first[x] in index_list:
				index1 = read1_first[x]
				qual1 = qual1_first[x]
			else:
				index1 = read1_first[x]
				qual1 = qual1_first[x]
			if read1_second[x] in index_list:
				index2 = read1_second[x]
				qual2 = qual1_second[x]
			else:
				index2 = read1_second[x]
				qual2 = qual1_second[x]

			output_file.write("%s\n%s%s\n+\n%s%s\n" % (header1[x],index1,index2,qual1,qual2))
	output_file.close()