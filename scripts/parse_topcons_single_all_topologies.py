#!/usr/bin/env python

# python parse_topcons_single_all_topologies.py alpha_polytopic_bitopic_noblanks__single.topcons.net__all_topologies.txt alpha_polytopic_bitopic.scampi_seq alpha_polytopic_bitopic.stmhmm alpha_polytopic_bitopic.hmmtop alpha_polytopic_bitopic.memsat

# python parse_topcons_single_all_topologies.py [input-file] [output-file-scampi-seq] [output-file-stmhmm] [output-file-hmmtop] [output-file-memsat]

import sys
import re


# input file :
# >2zt9_G
# scampi-seq
# ooooMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiii
# 
# stmhmm
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# hmmtop
# ooooMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiii
# 
# memsat
# oooooooMMMMMMMMMMMMMMMMMMMiiiiiiiiiii
# 
# >3din_E
# scampi-seq
# oMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMoo
# 
# stmhmm
# iiMMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMi
# 
# hmmtop
# iiiiiMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMii
# 
# memsat
# iiiiiiMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMiiiiiii


# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name_scampi_seq = sys.argv[2]
	output_file_name_stmhmm = sys.argv[3]
	output_file_name_hmmtop = sys.argv[4]
	output_file_name_memsat = sys.argv[5]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile_scampi_seq = open( output_file_name_scampi_seq, "w" )
	outfile_stmhmm = open( output_file_name_stmhmm, "w" )
	outfile_hmmtop = open( output_file_name_hmmtop, "w" )
	outfile_memsat = open( output_file_name_memsat, "w" )

	curr_id = ''
	next_is_scampi_seq = 0
	next_is_stmhmm = 0
	next_is_hmmtop = 0
	next_is_memsat = 0
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				curr_id = inline[1:7]
			elif (inline == 'scampi-seq'):
				next_is_scampi_seq = 1
			elif (inline == 'stmhmm'):
				next_is_stmhmm = 1
			elif (inline == 'hmmtop'):
				next_is_hmmtop = 1
			elif (inline == 'memsat'):
				next_is_memsat = 1
			elif (next_is_scampi_seq == 1):
				output_line = curr_id + ':::::' + inline + ':::::' + "\r\n"
				outfile_scampi_seq.write( output_line )
				next_is_scampi_seq = 0
			elif (next_is_stmhmm == 1):
				output_line = curr_id + ':::::' + inline + ':::::' + "\r\n"
				outfile_stmhmm.write( output_line )
				next_is_stmhmm = 0
			elif (next_is_hmmtop == 1):
				output_line = curr_id + ':::::' + inline + ':::::' + "\r\n"
				outfile_hmmtop.write( output_line )
				next_is_hmmtop = 0
			elif (next_is_memsat == 1):
				output_line = curr_id + ':::::' + inline + ':::::' + "\r\n"
				outfile_memsat.write( output_line )
				next_is_memsat = 0
	outfile_scampi_seq.close()
	outfile_stmhmm.close()
	outfile_hmmtop.close()
	outfile_memsat.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

