#!/usr/bin/env python

# python parse_tmpred_web_page.py alpha_polytopic_bitopic_noblanks.fasta.tmpred alpha_polytopic_bitopic_noblanks.tmpred

# python parse_tmpred_web_page.py [input-file] [output-file]

import sys
import re


# input file :
# 2xzb_A:::::
# min=17 | max=33 | html |  | plain_text | # GKAENYELYQVELGPGPSGDMAAKMSKKKAGRGGGKRKEKLENMKKEMEINDHQLSVAELEQKYQTSATKGLSASLAAELLLRDGPNALRPPRGTPEYVKFARQLAGGLQCLMWVAAAICLIAFAIQASEGDLTTDDNLYLALALIAVVVVTGCFGYYQEFKSTNIIASFKNLVPQQATVIRDGDKFQINADQLVV
# QDSYGQEWTFGQRLYQQYTCYTVFFISIEMCQIADVLIRKTRRLSAFQQGFFRNRILVIAIVFQVCIGCFLCYCPGMPNIFNFMPIRFQWWLVPMPFGLLIFVYDEIRKLGVRCCPGSWWDQELYY <p>
# <TITLE>TMpred output for unknown </TITLE><DOCUMENT><body bgcolor='#ddffdd'><H1>TMpred output for unknown </H1><B>[ISREC-Server]</B> Date: Wed Mar 14 11:13:10 2012
#  <P>
# <HR>tmpred -par=matrix.tab -html -min=17 -max=33 -def -in=wwwtmp/.TMPRED.13015.9253.seq -out=wwwtmp/.TMPRED.13015.9253.out -out2=wwwtmp/.TMPRED.13015.9253.out2 -out3=wwwtmp/.TMPRED.13015.9253.txt # >wwwtmp/.TMPRED.13015.9253.err<BR>Sequenc
# e: GKA...LYY,   length:    1033<BR>
# Prediction parameters: TM-helix length between 17 and 33
# <HR>
#  
# <H2> 1.) Possible transmembrane helices </H2>
# The sequence positions in brackets denominate the core region.<BR>
# Only scores above  500 are considered significant.<BR>
# <PRE>
# Inside to outside helices :  10 found
#       from        to    score center
#  108 ( 108) 126 ( 126)   2419    118
#  139 ( 139) 158 ( 156)   2217    148
#  303 ( 306) 322 ( 322)   2650    314
# 
# Outside to inside helices :  12 found
#       from        to    score center
#  105 ( 107) 126 ( 126)   1847    115
#  139 ( 139) 158 ( 158)   2439    149
#  252 ( 252) 271 ( 269)     65    261
# 
# </PRE> <HR>
# <H2> 2.) Table of correspondences </H2>
# Here is shown, which of the inside->outside helices correspond
# to which of the outside->inside helices.<BR>
# <BLOCKQUOTE>
# Helices shown in brackets are considered insignificant.<BR>
# A "+"-symbol indicates a preference of this orientation.<BR>
# A "++"-symbol indicates a strong preference of this orientation.<BR>
# </BLOCKQUOTE>
# <PRE>
#  
#            inside->outside | outside->inside
#    108- 126 (19) 2419 ++   |   105- 126 (22) 1847      
#    139- 158 (20) 2217      |   139- 158 (20) 2439 ++   
#                           |(  252- 271 (20)   65 ++ ) 
#    303- 322 (20) 2650      |   303- 324 (22) 2606      
#                           |   326- 355 (30) 1581 ++   
#                           |(  550- 569 (20)  374 ++ ) 
#    584- 601 (18)  586  +   |(  584- 602 (19)  432    ) 
# (  623- 647 (25)  120 ++ ) |
#    801- 819 (19) 1263      |   796- 817 (22) 1768 ++   
#    860- 882 (23)  963 ++   |(  866- 895 (30)  263    ) 
# (  925- 945 (21)  466    ) |   925- 944 (20)  677 ++   
#    963- 982 (20) 2591 ++   |   963- 981 (19) 2145      
#    993-1010 (18) 1315      |   993-1012 (20) 1360      
# 
# </PRE><HR>
# <H2> 3.) Suggested models for transmembrane topology</H2>
# These suggestions are purely speculative and should be used with
# <B>extreme caution</B> since they are based on the assumption that
# all transmembrane helices have been found.<BR>
# In most cases, the Correspondence Table shown above or the
# prediction plot that is also created should be used for the
# topology assignment of unknown proteins.<BR>
# <PRE>
# 
# 2 possible models considered, only significant TM-segments used
# 
# *** the models differ in the number of TM-helices ! ***
# 
# -----> STRONGLY prefered model: N-terminus inside
# 10 strong transmembrane helices, total score : 17034
#  # from   to length score orientation
#  1  108  126 (19)    2419 i-o
#  2  139  158 (20)    2439 o-i
#  3  303  322 (20)    2650 i-o
# 
# ------> alternative model
#  9 strong transmembrane helices, total score : 14615
#  # from   to length score orientation
#  1  105  126 (22)    1847 o-i
#  2  139  158 (20)    2217 i-o
#  3  303  324 (22)    2606 o-i


# 2 possible models considered, only significant TM-segments used
# 
# -----> slightly prefered model: N-terminus inside
#  1 strong transmembrane helices, total score : 2288
#  # from   to length score orientation
#  1   24   41 (18)    2288 i-o
# 
# ------> alternative model
#  1 strong transmembrane helices, total score : 2250
#  # from   to length score orientation
#  1   24   41 (18)    2250 o-i


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
	in_an_id = 0
	in_prefered_model = 0
	seen_prefered_model = 0
	seq_length = 0
	curr_topology_seq = ''
	prev_topology_char = '?'
	prev_from = 1
	for inline in inlines1:
		#inline = inline.strip()
		#if (inline != ''):

		# 2bg9_B:::::
		is_id_line = 0
		if (len(inline) >= 11):
			if (inline[6:11] == ':::::'):
				is_id_line = 1

		# e: GKA...LYY,   length:    1033<BR>
		is_length_line = 0
		if ( inline.find('length:') != -1 ):
			is_length_line = 1

		if (is_id_line == 1):
			# output previous sequence
			if (in_an_id == 1):
				if (curr_topology_seq == ''):
					curr_topology_seq = '-'
				else:
					if (seq_length == 0):
						print curr_id, ': did not find length of sequence for it'
					if (prev_topology_char == '?'):
						print curr_id, ': did not find topology char for it'
					for i in range( prev_from, (seq_length + 1) ):
						curr_topology_seq = curr_topology_seq + prev_topology_char
				output_line = curr_id + ':::::' + curr_topology_seq + ':::::' + "\r\n"
				outfile.write( output_line )

			curr_id = inline[0:6]
			in_an_id = 1
			in_prefered_model = 0
			seen_prefered_model = 0
			seq_length = 0
			curr_topology_seq = ''
			prev_topology_char = '?'
			prev_from = 1

		elif (is_length_line == 1):
			bits = inline.split( 'length:' )
			bit1 = bits[1]
			bits2 = bit1.split( '<BR>' )
			bit2 = bits2[0]
			seq_length = bit2.strip()
			seq_length = int(seq_length)

		# -----> STRONGLY prefered model: N-terminus inside
		# -----> slightly prefered model: N-terminus inside

		elif (inline.find('prefered model') != -1):
			if (seen_prefered_model == 0):
				in_prefered_model = 1

		elif (in_prefered_model > 0):
			in_prefered_model = in_prefered_model + 1
			if (in_prefered_model >= 4):
				inline_strip = inline.strip()
				if (inline_strip == ''):
					in_prefered_model = 0
				else:
					# 1  108  126 (19)    2419 i-o
					# 2  139  158 (20)    2439 o-i
					# 3  303  322 (20)    2650 i-o
					# 4  326  355 (30)    1581 o-i
					# 5  584  601 (18)     586 i-o
					# 6  796  817 (22)    1768 o-i
					# 7  860  882 (23)     963 i-o
					# 8  925  944 (20)     677 o-i
					# 9  963  982 (20)    2591 i-o
					#10  993 1012 (20)    1360 o-i
					from_resi = inline[2:7]
					from_resi = from_resi.strip()
					from_resi = int(from_resi)
					to_resi = inline[7:12]
					to_resi = to_resi.strip()
					to_resi = int(to_resi)
					from_topology_char = inline[26:27]
					to_topology_char = inline[28:29]
					if (prev_topology_char != '?'):
						if (prev_topology_char != from_topology_char):
							print curr_id, 'has inconsistent topology :', inline
					for i in range( prev_from, from_resi ):
						curr_topology_seq = curr_topology_seq + from_topology_char
					for i in range( from_resi, (to_resi + 1) ):
						curr_topology_seq = curr_topology_seq + 'H'
					prev_from = to_resi + 1
					prev_topology_char = to_topology_char

	# output last sequence
	if (in_an_id == 1):
		if (curr_topology_seq == ''):
			curr_topology_seq = '-'
		else:
			if (seq_length == 0):
				print curr_id, ': did not find length of sequence for it'
			if (prev_topology_char == '?'):
				print curr_id, ': did not find topology char for it'
			for i in range( prev_from, (seq_length + 1) ):
				curr_topology_seq = curr_topology_seq + prev_topology_char
		output_line = curr_id + ':::::' + curr_topology_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

