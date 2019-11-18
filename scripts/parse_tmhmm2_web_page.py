#!/usr/bin/env python

# python parse_tmhmm2_web_page.py alpha_polytopic_bitopic_noblanks.fasta.tmhmm alpha_polytopic_bitopic_noblanks.tmhmm

# python parse_tmhmm2_web_page.py [input-file] [output-file]

import sys
import re


# input file :
# 2xzb_A	len=1033	ExpAA=169.10	First60=0.00	PredHel=8	Topology=i107-129o139-158i304-326o331-353i799-821o921-943i963-984o994-1011i
# 2b6o_A	len=263	ExpAA=135.47	First60=40.27	PredHel=6	Topology=i13-35o40-62i83-105o125-147i160-182o202-224i
# 2zt9_G	len=37	ExpAA=21.82	First60=21.82	PredHel=1	Topology=o4-26i
# 1l0l_D	len=241	ExpAA=9.39	First60=0.00	PredHel=0	Topology=o


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
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			curr_id = inline[0:6]
			bits1 = inline.split( 'len=' )
			bits1a = bits1[1]
			bits1b = bits1a.split( "\t" )
			this_seq_length = bits1b[0]
			this_seq_length = int(this_seq_length)
			bits2 = inline.split( 'Topology=' )
			topology_string = bits2[1]
			prev_resi = 1
			prev_location = ''
			curr_location = ''
			curr_resi_string = ''
			output_string = ''
			for i in range( 0, len(topology_string) ):
				j = i + 1
				this_char = topology_string[i:j]

				# when see 'i' or 'o', finish off the preceding MMM sequence
				# need to use prev_resi and this just-finished resi
				# set this just-finished resi + 1 as the prev_resi

				if ((this_char == 'i') or (this_char == 'o')):
					if (curr_resi_string != ''):
						new_resi = int(curr_resi_string)
						for k in range( prev_resi, (new_resi + 1) ):
							output_string = output_string + 'M'
						prev_resi = new_resi + 1
					prev_location = curr_location
					curr_location = this_char
					curr_resi_string = ''

				# when see '-', now have the full preceding number so finish off the preceding i/o before it
				# need to use prev_resi and this just-finished resi
				# set this just-finished resi as the prev_resi

				elif (this_char == '-'):
					new_resi = int(curr_resi_string)
					if (curr_resi_string != ''):
						if (prev_location != ''):
							for k in range( prev_resi, new_resi ):
								output_string = output_string + prev_location
					prev_resi = new_resi
					curr_resi_string = ''

				# when see the first number after seeing i/o

				elif ((curr_location == 'i') or (curr_location == 'o')):
					prev_location = curr_location
					curr_location = 'M'
					curr_resi_string = this_char

				# when see a consecutive number

				else:
					curr_resi_string = curr_resi_string + this_char

			# at the end of the topology info for this sequence, finish off the output

			if ((curr_location == 'i') or (curr_location == 'o')):
				for k in range( prev_resi, (this_seq_length + 1) ):
					output_string = output_string + curr_location

			if (output_string == ''):
				output_string = '-'

			output_line = curr_id + ':::::' + output_string + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

