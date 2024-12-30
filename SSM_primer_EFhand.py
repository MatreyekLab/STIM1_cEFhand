## Primer design for site saturation mutagenesis of cEF hand of STIM1
import Bio
from Bio.Seq import Seq
from datetime import datetime

## Sequence before cEF hand.
before = "aagcccctgtgtcacagtgaggatgagaaa"


## Sequence of ORF to mutagenize
stim_cEF = "ctcagcttcgaggcagtccgtaacatccacaaactgatggacgatgatgccaatggtgatgtggatgtggaagaaagtgatgagttcctgagggaagacctcaattac"

## Sequence after cEF hand.
after = "catgacccaacagtgaaacacagcaccttc"

## Degenerate codon of choice
deg = "NNK"

## Concatenate all of the above
fullstring = Seq(before + stim_cEF + after)
#print(fullstring)

## Numbers for For Loop
startnum = 0

## Number of residues * 3 + 3 (since not inclusive)
endnum = 36 * 3 + 3

plist = list()
plistf = list()
plistrc = list()
text_file = open(str(datetime.now()) + "SSM_primer_EFhand_Output.txt", "w")
before_offset = len(before)

for x in range(startnum, endnum, 3):
	pos = fullstring[(before_offset + x):(before_offset + x + 3)]
	plist.append(str(pos))
	# posf adds two strings - one before NNK and one after the 3 letters NNK replaces
	posf = fullstring[(before_offset + x - 18):(before_offset + x)] + deg + fullstring[(before_offset + x + 3):(before_offset + x + 3 + 15)]
	plistf.append(str(posf))
	posr = fullstring[(before_offset + x - 25):(before_offset + x)]
	posrc = posr.reverse_complement()
	plistrc.append(str(posrc))

ntcount = 0

##forward primers
for x in range(0, 36):
	print("EFhand_f_%s\t%s" % (x+2, plistf[x]))
	text_file.write("EFhand_f_%s\t%s\n" % (x+2, plistf[x]))
	ntcount = ntcount + len(plistf[x])

##reverse primers
for x in range(0, 36):
	print("EFhand_r_%s\t%s" % (x+2, plistrc[x]))
	text_file.write("EFhand_r_%s\t%s\n" % (x+2, plistrc[x]))
	ntcount = ntcount + len(plistrc[x])

print("total cost: $ %d" % (ntcount*0.18))
text_file.close()