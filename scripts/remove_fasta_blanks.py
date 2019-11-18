#!/usr/bin/env python

# python remove_fasta_blanks.py alpha_polytopic_bitopic.fasta alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic.leading_blanks_removed

# python remove_fasta_blanks.py [input-file] [output-file] [output-file-2]

import sys
import re


# global variables
save_id = []
save_seq = []
save_num_leading_blanks = []
resn_list_3char = ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']
resn_list_1char = ['R',  'H',  'K',  'D',  'E',  'S',  'T',  'N',  'Q',  'C',  'G',  'P',  'A',  'V',  'I',  'L',  'M',  'F',  'Y',  'W']


######################################################
def is_valid_resn( this_char ):
	is_valid = 0
	for this_resn in resn_list_1char:
		if (this_resn == this_char):
			is_valid = 1
	return is_valid

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name = sys.argv[2]
	output_file_name_2 = sys.argv[3]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# parse the sequences because they are over one line

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				if (in_a_seq == 1):
					save_id.append( curr_id )
					save_seq.append( curr_seq )
					save_num_leading_blanks.append( 0 )
				in_a_seq = 1
				curr_id = inline[1:7]
				curr_seq = ''
			else:
				curr_seq = curr_seq + inline
	if (in_a_seq == 1):
		save_id.append( curr_id )
		save_seq.append( curr_seq )
		save_num_leading_blanks.append( 0 )

	# remove leading blanks and replace mid-sequence blanks by 'G'

	replacement_char = 'G'
	i = 0
	for this_id in save_id:
		this_seq = save_seq[i]
		seen_start_of_valid_resns = 0
		new_seq = ''
		skipped_resns = 0
		for j in range( 0, len(this_seq) ):
			k = j + 1
			this_char = this_seq[j:k]
			is_valid = is_valid_resn( this_char )
			if (is_valid == 1):
				new_seq = new_seq + this_char
				seen_start_of_valid_resns = 1
			else:
				if (seen_start_of_valid_resns == 1):
					new_seq = new_seq + replacement_char
				else:
					skipped_resns = skipped_resns + 1
		save_seq[i] = new_seq
		save_num_leading_blanks[i] = skipped_resns
		i = i + 1

	# write out the new sequence and record how many leading blanks were skipped

	outfile = open( output_file_name, "w" )
	outfile2 = open( output_file_name_2, "w" )

	max_line_length = 60
	i = 0
	for this_id in save_id:
		this_seq = save_seq[i]
		this_num_leading_blanks = save_num_leading_blanks[i]

		output_line = ">" + this_id + "\r\n"
		outfile.write( output_line )
		end_seq = 0
		remaining_seq = this_seq
		while (end_seq == 0):
			length_seq = len(remaining_seq)
			if (length_seq <= max_line_length):
				write_seq = remaining_seq[ 0: length_seq ]
				output_line = write_seq + "\r\n"
				outfile.write( output_line )
				end_seq = 1
			else:
				write_seq = remaining_seq[ 0: max_line_length ]
				remaining_seq = remaining_seq[ max_line_length: ]
				output_line = write_seq + "\r\n"
				outfile.write( output_line )

		output_line2 = this_id + ":::::" + str(this_num_leading_blanks) + ":::::" + "\r\n"
		outfile2.write( output_line2 )
		i = i + 1

	outfile.close()
	outfile2.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

