#!/usr/bin/env python

# python sort_blastall_results.py pdb_select25_expect0.eukaryota.1line.blastall pdb_select25_expect0.eukaryota.1line.blastall.sorted

# python sort_blastall_results.py [input-file-1] [output-file]

# 1dkc_A:::::1dkc_A	3ug9_A	62.50	8	3	0	3	10	43	50	   12	18.1:::::
# 2pjh_A:::::2pjh_A	3rze_A	36.11	36	17	2	10	39	342	377	0.93	21.6:::::
# 1zza_A:::::1zza_A	3chx_A	50.00	18	8	1	42	59	186	202	3.9	19.6:::::

# The fields for tabular BLAST output are:
# 1	Query	The query sequence id
# 2	Subject	The matching subject sequence id
# 3	% id
# 4	alignment length
# 5	mistmatches
# 6	gap openings
# 7	q.start
# 8	q.end
# 9	s.start
# 10	s.end
# 11	e-value
# 12	bit score


import sys
import re


# global variables
save_pdbid = []
save_id = []
output_id = []
output_pdbid = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	sortlines = []

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.split( ':::::' )
			this_id = bits[0]
			this_blastall = bits[1]
			bits2 = this_blastall.split( "\t" )
			this_percent_identity = float(bits2[2])
			this_formatted = "%3.2f" % this_percent_identity
			this_formatted = this_formatted.zfill( 6 )
			sortline = this_formatted + ':::::' + this_id
			sortlines.append( sortline )

	sortlines.sort()

	outfile = open( output_file_name, "w" )

	for sortline in sortlines:
		bits = sortline.split( ':::::' )
		this_id = bits[1]
		this_formatted = bits[0]
		output_line = this_id + ':::::' + this_formatted + ':::::' + "\r\n"
		outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

