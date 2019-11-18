#!/usr/bin/env python

# python ../extract_pdbid_from_opm_classes_webpage.py OPM_Alpha_helical_polytopic_superfamilies.html OPM_Alpha_helical_polytopic_superfamilies.pdbid

# python extract_pdbid_from_opm_classes_webpage.py [input-file] [output-file]

import sys
import re

# input file :
# ...
# var FiltersEnabled = 0
# Text[0]=["2l6x &raquo; Green-light absorbing proteorhodopsin","<img src=phpthumb/phpThumb.php?src=../images/png/2l6x.png&amp;w=200&amp;h=200>"]
# Text[1]=["1h68 &raquo; Sensory rhodopsin II, monomer","<img src=phpthumb/phpThumb.php?src=../images/png/1h68.png&amp;w=200&amp;h=200>"]
# ...
# Text[339]=["3v5u &raquo; Sodium/calcium exchanger","<img src=phpthumb/phpThumb.php?src=../images/png/3v5u.png&amp;w=200&amp;h=200>"]
# Text[340]=["3tij &raquo; Concentrative nucleoside transporter","<img src=phpthumb/phpThumb.php?src=../images/png/3tij.png&amp;w=200&amp;h=200>"]
# Style[0]=["white","#000","#0066cc","#ffffff","","","center","","","","","","","",100,"",2,2,10,10,51,1,0,"",""]
# applyCssFilter()
# ...

# global variables

######################################################
def main():

	global outfile
	input_file = sys.argv[1]
	output_file_name = sys.argv[2]

	outfile = open( output_file_name, "w" )

	infile = open( input_file, "r" )
	inlines = infile.readlines()
	infile.close()

	in_list_of_pdbids = 0
	found_end_of_pdbids = 0
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			if ((in_list_of_pdbids == 0) and (found_end_of_pdbids == 0)):
				if (inline == 'var FiltersEnabled = 0'):
					in_list_of_pdbids = 1
			else:
				if ((in_list_of_pdbids == 1) and (found_end_of_pdbids == 0)):
					bit1 = inline[0:9]
					if (bit1 == 'Style[0]='):
						found_end_of_pdbids = 1
					else:
						bits = inline.split( '=["' )
						bit1 = bits[1]
						pdbid = bit1[0:4]
						output_line = pdbid + "\r\n"
						outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

