#!/usr/bin/env python

# python read_opm_txt_all_types.py All_TM_subunits_in_OPM.txt opm_feb2012.all_opm.TMH

# python read_opm_txt_all_types.py [input-file] [output-file]

import sys
import re


# global variables


######################################################
def reformat_input_line( inline ):

	formatted_line = ''

	this_pdbid = inline[0:4] + '_' + inline[5:6]

	bits1 = inline.rsplit( 'Segments: ' )
	segments = bits1[1]
	each_segment = segments.rsplit( ', ' )
	segments_string = ''
	for this_segment in each_segment:
		bits2 = this_segment.rsplit( '(' )
		bit3 = bits2[1]
		bits4 = bit3.rsplit( ')' )
		final_segment = bits4[0]
		final_segment = re.sub(r'\s', '', final_segment) # remove all white space
		segments_string = segments_string + final_segment + ','

	segments_string = 'MEMBRANE-ALPHA-HELIX,,,,,' + segments_string + ';;;;;'
	formatted_line = this_pdbid + ':::::' + segments_string + ':::::'

	return formatted_line

# 1ar1 A - Tilt: 7 - Segments: 1(31-55), 2(90-113), 3(129-151), 4(178-202), 5(219-247), 6(269-294), 7(305-326), 8(338-363), 9(371-394), 10(407-429), 11(443-464), 12(488-512)
# 1ar1 B - Tilt: 14 - Segments: 1(37-60), 2(75-94)
# 1b9u A - Tilt: 41 - Segments: 1( 4- 30)
# 1c17 A - Tilt: 3 - Segments: 1(8-32), 2(52-75)

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			output_line = reformat_input_line( inline )
			output_line = output_line + "\r\n"
			outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

