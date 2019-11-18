#!/usr/bin/env python

# python parse_tmloop_output.py alpha_polytopic_bitopic_noblanks.fasta_TMLOOP_single_summary.txt alpha_polytopic_biotopic_noblanks.tmloop

# python parse_tmloop_output.py [input-file] [output-file]

import sys
import re


# input file :
# TMLOOP Prediction Analysis
# TMLOOP Prediction Parameters
# Prediction method : Single Motif Method
# Minimum Interloop Length (I) = 30
# Results
# 
# Protein Matches:	50 / 481
# *!* Loop(s) found:
# 	>2b6o_A
# 	Loop found
# 		Loop type:	loop-in-turn-helix-out loop 1 & 2 aquaglycerolporin like
# 		Loop score:	1.000	(1 / 1)
# 		Average Pattern starting point:	 62
# 		Average Pattern length:	  8
# 		Estimated domain boundaries:	63	to	77
# 	Loop found
# 		Loop type:	loop-in-turn-helix-out loop 1 & 2 aquaglycerolporin like
# 		Loop score:	1.000	(1 / 1)
# 		Average Pattern starting point:	178
# 		Average Pattern length:	  8
# 		Estimated domain boundaries:	179	to	193
# 	Loop found
# 		Loop type:	loop-in-turn-helix-out loop 1 aquaglycerolporin like
# 		Loop score:	1.000	(1 / 1)
# 		Average Pattern starting point:	 62
# 		Average Pattern length:	 10
# 		Estimated domain boundaries:	63	to	77
# *!* Loop(s) found:
# 	>1h6i_A
# 	Loop found
# 		Loop type:	loop-in-turn-helix-out loop 1 & 2 aquaglycerolporin like
# 		Loop score:	1.000	(1 / 1)
# 		Average Pattern starting point:	 70
# 		Average Pattern length:	  8
# 		Estimated domain boundaries:	71	to	85
# 	Loop found
# 		Loop type:	loop-in-turn-helix-out loop 1 & 2 aquaglycerolporin like
# 		Loop score:	1.000	(1 / 1)
# 		Average Pattern starting point:	186
# 		Average Pattern length:	  8
# 		Estimated domain boundaries:	187	to	201


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

	curr_id = ''
	curr_seq = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):

			# 	>2b6o_A
			if (inline[0:1] == '>'):

				# finish processing previous id
				if (curr_id != ''):
					if (curr_seq == ''):
						curr_seq = '-'
					output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
					outfile.write( output_line )

				# start processing this new id
				curr_id = inline[1:7]
				curr_seq = ''

			# 		Estimated domain boundaries:	71	to	85
			if (inline.find( 'Estimated domain boundaries:' ) > -1):
				bits = inline.split( 'Estimated domain boundaries:' )
				bit1 = bits[1]
				bits2 = bit1.rsplit()
				this_from = bits2[0]
				this_to = bits2[2]
				curr_seq = curr_seq + this_from + '-' + this_to + ','

	# finish processing last id
	if (curr_id != ''):
		if (curr_seq == ''):
			curr_seq = '-'
		output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

