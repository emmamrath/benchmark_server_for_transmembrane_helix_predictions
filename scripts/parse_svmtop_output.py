#!/usr/bin/env python

# python parse_svmtop_output.py alpha_polytopic_bitopic_noblanks.svmtop_download_results alpha_polytopic_biotopic_noblanks.svmtop

# python parse_svmtop_output.py [input-file] [output-file]

import sys
import re


# input file :
#Start Time: 2012-03-23 12:10:19
#Finish Time: 2012-03-24 00:38:34
#Time Elapsed: 12:28:15
#
#>1s5h_C
#MAPMLSGLLARLVKLLLGRHGSALHWRAAGAATVLLVIVLLAGSYLAVLAERGAPGAQLITYPRALWWSVETATCVGYGDLYPVTLWGRLVAVVVMVAGITSFGLVTAALATWFVGREQERRGH
#iiiiiiiiiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHHHHHHoooooooooooooooooooooooooooooooooooooHHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiii
#9988777788889999999999999999999999999999999999999999999999999999999999999998999999999999999999999999999999999999999999999999
#>3din_E
#MKTFFLIVHTIISVALIYMVQVQMSKFSELGGAFGSGGLHTVFGRRKGLDTGGKITLVLSVLFFVSCVVTAFVLTR
#oHHHHHHHHHHHHHHHHHHHHHHHiiiiiiiiiiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHHHHHHHHHHHHHo
#9945899999999999987641789988887613523546788999999998771478899999999999999799


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
		if ( inline.find( 'is a probable soluble protein.' ) == -1 ):
			inlines1b.append( inline )
		else:
			bits = inline.split( 'is a probable soluble protein.' )
			bit0 = 'is a probable soluble protein.'
			inlines1b.append( bit0 )
			inlines1b.append( bits[1] )

	outfile = open( output_file_name, "w" )

	i = 0
	for inline in inlines1b:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				this_id = inline[1:7]
				this_seq = inlines1b[ i + 2 ]
				this_seq = this_seq.strip()
				if ( this_seq.find( 'is a probable soluble protein.' ) != -1 ):
					this_seq = '-'
				output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
				outfile.write( output_line )
		i = i + 1

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

