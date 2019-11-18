#!/usr/bin/env python

# python ../list_ids_having_duplicate_field2_entries.py alpha_polytopic_bitopic_peptide.fix_opm_posn

# python list_ids_having_duplicate_field2_entries.py [input-file]

import sys
import re


# global variables


######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	#output_file_name = sys.argv[2]

	#outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_field2 = bits[1]
			save_field2 = []
			this_id_has_duplicate_field2 = 0
			bits2 = this_field2.rsplit( ';;;;;' )
			for bit2 in bits2:
				if (bit2 != ''):
					bits3 = bit2.rsplit( ',,,,,' )
					this_field2_name = bits3[0]
					for this_save_field2 in save_field2:
						if (this_save_field2 == this_field2_name):
							this_id_has_duplicate_field2 = 1
					save_field2.append( this_field2_name )
			if (this_id_has_duplicate_field2 == 1):
				#output_line = inline + "\r\n"
				#outfile.write( output_line )
				print this_id
	#outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

