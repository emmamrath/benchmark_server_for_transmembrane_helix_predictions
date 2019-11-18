#!/usr/bin/env python

# python parse_mempype_ensemble_topograpies.py alpha_mempype_tCagy082.txt alpha_mempype_tCagy082.noblanks.mempype

# python parse_mempype_ensemble_topograpies.py [input-file] [output-file]

import sys
import re


# input file :

# Sequence name:
# 
# 1A91_A
# 
# Prediction status: completed
# Prediction summary:
# 
# Organelle membranes, 2 Transmembrane helices
# Detailed prediction results
# 
# Predicted sequence features:
# Prediction	Presence	Start	End	Detail
# Signal peptide (SPEP):	NO	-	-	not present
# Non cytoplasmic region (ENSEMBLE):	-	1	6	view on sequence
# Transmembrane region (ENSEMBLE):	YES	7	30	view on sequence
# Cytoplasmic region (ENSEMBLE):	-	31	46	view on sequence
# Transmembrane region (ENSEMBLE):	YES	47	77	view on sequence
# Non cytoplasmic region (ENSEMBLE):	-	78	79	view on sequence
# GPI anchor (PredGPI):	NO	-	-	not present
# Sequence name:
# 
# 1A0S_P
# 
# Prediction status: completed
# Prediction summary:
# 
# Internal membranes, Globular
# Detailed prediction results
# 
# Predicted sequence features:
# Prediction	Presence	Start	End	Detail
# Signal peptide (SPEP):	NO	-	-	not present
# Transmembrane region (ENSEMBLE):	NO	-	-	not present
# GPI anchor (PredGPI):	NO	-	-	not present
# Sequence name:
# 
# 4UBP_B
# 
# Prediction status: in processing queue. Position = 8
# Sequence name:
# 
# 5CSM_A
# 
# Prediction status: completed
# Prediction summary:
# 
# Internal membranes, Globular
# Detailed prediction results
# 
# Predicted sequence features:
# Prediction	Presence	Start	End	Detail
# Signal peptide (SPEP):	NO	-	-	not present
# Transmembrane region (ENSEMBLE):	NO	-	-	not present
# GPI anchor (PredGPI):	NO	-	-	not present
# Global Prediction Counts
# Cell Membrane	:	



# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	i = 0
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if ( inline.find('Prediction status:') > -1 ):

				if (curr_id != ''):
					if (curr_seq == ''):
						curr_seq = '-'
					output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
					outfile.write( output_line )

				j = i - 2
				curr_id = inlines1[j]
				#curr_id = curr_id[0:6]
				this_pdbid = curr_id[0:4]
				this_pdbid = this_pdbid.lower()
				this_chain = curr_id[5:6]
				curr_id = this_pdbid + '_' + this_chain
				curr_seq = ''

			elif ( inline.find("Transmembrane region (ENSEMBLE):\tYES") > -1 ):
				bits = inline.split("\t")
				this_from = bits[2]
				this_to = bits[3]
				curr_seq = curr_seq + this_from + '-' + this_to + ','

		i = i + 1

	if (curr_id != ''):
		if (curr_seq == ''):
			curr_seq = '-'
		output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

