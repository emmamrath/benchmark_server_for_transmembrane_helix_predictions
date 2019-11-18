#!/usr/bin/env python

# python parse_phobius_output.py alpha_polytopic_biotopic_noblanks.fasta.phobius.txt alpha_polytopic_biotopic_noblanks.phobius

# python parse_phobius_output.py [input-file] [output-file]

import sys
import re


# input file :
# <!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
# <html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en-US" lang="en-US"><head>
# <title>Phobius prediction</title>
# <meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">
# </head>
# <body>
# <h2>Phobius prediction</h2><hr><pre>ID   2xzb_A
# FT   TOPO_DOM      1    104       CYTOPLASMIC.
# FT   TRANSMEM    105    127       
# FT   TOPO_DOM    128    138       NON CYTOPLASMIC.
# FT   TRANSMEM    139    158       
# FT   TOPO_DOM    159    299       CYTOPLASMIC.
# FT   TRANSMEM    300    326       
# FT   TOPO_DOM    327    331       NON CYTOPLASMIC.
# FT   TRANSMEM    332    353       
# FT   TOPO_DOM    354    798       CYTOPLASMIC.
# FT   TRANSMEM    799    819       
# FT   TOPO_DOM    820    921       NON CYTOPLASMIC.
# FT   TRANSMEM    922    944       
# FT   TOPO_DOM    945    963       CYTOPLASMIC.
# FT   TRANSMEM    964    983       
# FT   TOPO_DOM    984    994       NON CYTOPLASMIC.
# FT   TRANSMEM    995   1011       
# FT   TOPO_DOM   1012   1033       CYTOPLASMIC.
# //
# ID   2b6o_A
# FT   TOPO_DOM      1     11       CYTOPLASMIC.
# FT   TRANSMEM     12     35       
# FT   TOPO_DOM     36     40       NON CYTOPLASMIC.
# FT   TRANSMEM     41     62       
# FT   TOPO_DOM     63     82       CYTOPLASMIC.

# ID   MTH_DROMEa signal peptide
# FT   SIGNAL        1     24       
# FT   REGION        1      3       N-REGION.
# FT   REGION        4     19       H-REGION.
# FT   REGION       20     24       C-REGION.
# FT   TOPO_DOM     25    218       NON CYTOPLASMIC.
# FT   TRANSMEM    219    238       
# FT   TOPO_DOM    239    249       CYTOPLASMIC.
# FT   TRANSMEM    250    269       


# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	inlines1b = []
	for inline in inlines1:
		if ( inline.find('<pre>ID') != -1 ):
			bits = inline.split( '<pre>ID' )
			bit0 = bits[0] + '<pre>'
			bit1 = 'ID' + bits[1]
			inlines1b.append( bit0 )
			inlines1b.append( bit1 )
		else:
			inlines1b.append( inline )

	outfile = open( output_file_name, "w" )

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines1b:
		#inline = inline.strip()
		#if (inline != ''):
		if (len(inline) >= 2):
			if (inline[0:2] == 'ID'):

				if (in_a_seq == 1):
					if (curr_seq == ''):
						curr_seq = '-'
					output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
					outfile.write( output_line )

				curr_id = inline[5:11]
				curr_seq = ''
				in_a_seq = 1
				prev_to = 0
			elif (inline[0:2] == 'FT'):
				this_from = inline[13:20]
				this_from = this_from.strip()
				this_from = int(this_from)
				this_to = inline[20:27]
				this_to = this_to.strip()
				this_to = int(this_to)
				this_what = inline[5:13]
				this_where = inline[34:]
				this_where = this_where.strip()
				this_char = '_'
				if (this_what == 'SIGNAL  '):
					this_char = 'S'
					for i in range( this_from, (this_to + 1) ):
						curr_seq = curr_seq + this_char
					preV_to = this_to
				elif (this_what == 'TRANSMEM'):
					this_char = 'H'
					for i in range( this_from, (this_to + 1) ):
						curr_seq = curr_seq + this_char
					preV_to = this_to
				elif (this_what == 'TOPO_DOM'):
					if (this_where == 'CYTOPLASMIC.'):
						this_char = 'i'
					elif (this_where == 'NON CYTOPLASMIC.'):
						this_char = 'o'
					else:
						print curr_id, ': has unexpected TOPO_DOM value :', this_where
					for i in range( this_from, (this_to + 1) ):
						curr_seq = curr_seq + this_char
					prev_to = this_to
				elif (this_what == 'REGION  '):
					do_nothing = 1
				else:
					print curr_id, ': has unexpected value :', this_what

	if (in_a_seq == 1):
		if (curr_seq == ''):
			curr_seq = '-'
		output_line = curr_id + ':::::' + curr_seq + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

