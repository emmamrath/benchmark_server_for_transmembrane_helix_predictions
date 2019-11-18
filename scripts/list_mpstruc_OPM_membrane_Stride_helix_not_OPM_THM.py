#!/usr/bin/env python

# python list_mpstruc_OPM_membrane_Stride_helix_not_OPM_THM.py alpha_polytopic_bitopic.mpstruc_category.sorted alpha_polytopic_bitopic.opm_tm_segments alpha_polytopic_bitopic.stride_struc2D

# python list_mpstruc_OPM_membrane_Stride_helix_not_OPM_THM.py [input-file-mpstruc_category] [input-file-opm_tm_segments] [input-file-stride_struc2D]

import sys
import re


# global variables
save_id = []
save_mpstruc_category = []
save_opm_tm_segments = []
save_stride_struc2D = []

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	input_file_3 = sys.argv[3]

	infile = open( input_file_1, "r" )
	inlines = infile.readlines()
	infile.close()

	# 1ijd_A:::::_iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMmoommommo_______:::::
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_id.append( this_id )
			save_mpstruc_category.append( this_seq )
			save_opm_tm_segments.append( '' )
			save_stride_struc2D.append( '' )

	infile = open( input_file_2, "r" )
	inlines = infile.readlines()
	infile.close()

	# 1iij_A:::::CCHHHHHHHHHHHHHHHHHHHHHHTTTTHHHHHTT:::::
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			i = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					save_opm_tm_segments[i] = this_seq
				i = i + 1

	infile = open( input_file_3, "r" )
	inlines = infile.readlines()
	infile.close()

	# 1iij_A:::::CCHHHHHHHHHHHHHHHHHHHHHHTTTTHHHHHTT:::::
	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			i = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					save_stride_struc2D[i] = this_seq
				i = i + 1

	#print len(save_id)
	#print len(save_mpstruc_category)
	#print len(save_opm_tm_segments)
	#print len(save_stride_struc2D)

	i = 0
	for this_id in save_id:
		this_mpstruc_category = save_mpstruc_category[i]
		this_opm_tm_segments = save_opm_tm_segments[i]
		this_stride_struc2D = save_stride_struc2D[i]
		num_consecutive_mismatch = 0
		this_id_got_string_of_mismatch = 0
		#print this_id
		#print this_mpstruc_category
		#print this_opm_tm_segments
		#print this_stride_struc2D
		#print 
		for j in range( 0, len(this_opm_tm_segments) ):
			k = j + 1
			char1 = this_opm_tm_segments[j:k]
			char2 = this_stride_struc2D[j:k]
			if ((char1 == 'm') and (char2 == 'H')):
				num_consecutive_mismatch = num_consecutive_mismatch + 1
			else:
				num_consecutive_mismatch = 0
			if (num_consecutive_mismatch >= 6):
				this_id_got_string_of_mismatch = 1
		if (this_id_got_string_of_mismatch == 1):
			print this_id, ':', this_mpstruc_category
			print this_opm_tm_segments
			print this_stride_struc2D
			print 
		i = i + 1

	#outfile = open( output_file_name, "w" )
	#output_line = output_line + "\r\n"
	#outfile.write( output_line )
	#outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

