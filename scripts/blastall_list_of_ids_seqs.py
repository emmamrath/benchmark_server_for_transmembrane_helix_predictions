#!/usr/bin/env python

# python blastall_list_of_ids_seqs.py pdb_select25_expect0.eukaryota.1line alpha_polytopic_bitopic.fasta pdb_select25_expect0.eukaryota.1line.blastall

# python blastall_list_of_ids_seqs.py [input-file-1] [input-file-2] [output-file]

# input-file-2 has had the blast 'formatdb' run on it, eg. :
#
# formatdb -i alpha_polytopic_bitopic.fasta
#
# which produces 3 output files, eg. :
#
# alpha_polytopic_bitopic.fasta.psq
# alpha_polytopic_bitopic.fasta.pin
# alpha_polytopic_bitopic.fasta.phr


import sys
import re
import os


# global variables


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_fastaid = bits[1]
			this_seq = bits[2]

			tempfile_name = 'python_blastall_temp.fasta'
			tempfile = open( tempfile_name, "w" )
			output_line = this_fastaid + "\r\n" + this_seq + "\r\n"
			tempfile.write( output_line )
			tempfile.close()
			temp2file_name = 'python_blastall_temp2.out'

			this_command = 'blastall -p blastp -d ' + input_file_2 + ' -i ' + tempfile_name + ' -m 8 -v 1 -b 1 -o ' + temp2file_name + ' -e 50'

			#subprocess.call(["ls", "-l"])
			this_result = os.system( this_command )

			print 'Result for', this_id, 'is', this_result

			infile2 = open( temp2file_name, "r" )
			inlines2 = infile2.readlines()
			infile2.close()

			i = 0
			for inline2 in inlines2:
				inline2 = inline2.strip()
				if (inline2 != ''):
					if (i == 0):
						#print inline2
						output_line = this_id + ':::::' + inline2 + ':::::' + "\r\n"
						outfile.write( output_line )
				i = i + 1

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

