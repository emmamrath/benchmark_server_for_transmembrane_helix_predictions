#!/usr/bin/env python

# python parse_hmmtop_topology.py alpha_polytopic_bitopic_noblanks.fasta.hmmtop alpha_polytopic_bitopic_noblanks.hmmtop

# python parse_hmmtop_topology.py [input-file] [output-file]

import sys
import re


# input file :
# 2xzb_A:::::
# <HTML><HEAD>
# <TITLE>Server</TITLE>
# </HEAD><BODY>
# <HTML><HEAD><TITLE>Results of the prediction</TITLE></HEAD><BODY><PRE>
# Protein: noname
# Length:  1033
# N-terminus:  IN 
# Number of transmembrane helices: 7
# Transmembrane helices: 107-126 139-158 303-322 800-819 857-876 926-945 964-983 
# 
# Total entropy of the model:  17.0252
# Entropy of the best path:  17.0296
# 
# The best path:
# 
#      seq  GKAENYELYQ VELGPGPSGD MAAKMSKKKA GRGGGKRKEK LENMKKEMEI    50
#      pred IIIIIIIIII IIIIIIIIII IIIIIIIIII IIIIIIIIII IIIIIIIIII 
# 
#      seq  NDHQLSVAEL EQKYQTSATK GLSASLAAEL LLRDGPNALR PPRGTPEYVK   100
#      pred IIIIIIIIII IIIIIIIIII IIIIIIIIII IIIIIIIIII IIIIIIIiii 
# 
#      seq  FARQLAGGLQ CLMWVAAAIC LIAFAIQASE GDLTTDDNLY LALALIAVVV   150
#      pred iiiiiiHHHH HHHHHHHHHH HHHHHHoooo ooooooooHH HHHHHHHHHH 
# 
#      seq  LYARSYQMQD HYRQATGSEN LYFQ  524
#      pred HHHiiiiiii iiiiiiiiII IIII
# 
# </PRE></BODY></HTML>
# <BR><BR><FONT SIZE="-2">If you are going to use these results in your work, please cite:<UL>
# <LI>G.E Tusn&aacute;dy and I. Simon (1998) 
# <em>J. Mol. Biol.</em> <B>283</B>, 489-506.
# <LI>G.E Tusn&aacute;dy and I. Simon (2001) 
# <em>Bioinformatics</em> <b>17</b>, 849-850.</UL></FONT>


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
	curr_topology = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			new_id = ''
			if (len(inline) >= 11):
				if (inline[6:11] == ':::::'):
					new_id = inline[0:6]
			if (new_id != ''):

				# write out previous id
				if (curr_id != ''):
					if (curr_topology == ''):
						curr_topology = '-'
					output_line = curr_id + ':::::' + curr_topology + ':::::' + "\r\n"
					outfile.write( output_line )

				# start the new sequence
				curr_id = new_id
				curr_topology = ''

			if (len(inline) >= 5):
				if (inline[0:5] == 'pred '):
					new_topology = inline[5:]
					new_topology = new_topology.replace( ' ', '' )
					curr_topology = curr_topology + new_topology

	# write out last id
	if (curr_id != ''):
		if (curr_topology == ''):
			curr_topology = '-'
		output_line = curr_id + ':::::' + curr_topology + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

