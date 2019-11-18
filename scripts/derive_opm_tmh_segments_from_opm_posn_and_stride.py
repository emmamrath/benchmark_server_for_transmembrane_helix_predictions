#!/usr/bin/env python

# python ../derive_opm_tmh_segments_from_opm_posn_and_stride.py alpha_polytopic_bitopic_peptide.fix_opm_posn alpha_polytopic_bitopic_peptide.fix_stride_struc2D alpha_polytopic_bitopic_peptide.fix_opm_stride_tmh_segments

# python derive_opm_tmh_segments_from_opm_posn_and_stride.py [input-file-2] [input-file-3] [output-file]

# file 1 = opm_posn :
# 1ifp_A:::::MEMBRANE-SEGMENTS,,,,,13-35,;;;;;INSIDE-SEGMENTS,,,,,36-44,;;;;;OUTSIDE-SEGMENTS,,,,,1-12,;;;;;:::::
# 1n7l_A:::::MEMBRANE-SEGMENTS,,,,,2-2,4-5,8-8,31-53,;;;;;INSIDE-SEGMENTS,,,,,1-1,3-3,6-7,9-30,;;;;;:::::
# 1oqw_A:::::MEMBRANE-SEGMENTS,,,,,1-24,;;;;;OUTSIDE-SEGMENTS,,,,,25-144,;;;;;:::::

# file 2 = stride_struc2D :
# 1ifp_A:::::_HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH_:::::
# 1ijd_A:::::_TTGGGGGTS_HHHHHHHHHHHHHHHHHHHHHHHHHS_SHHHHHT________:::::
# 1ijd_B:::::_____HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTS____:::::
# 1grm_A:::::_EEEEEEEEEEEEE__:::::
# 1jb0_K:::::____________________HHHHHHHHHHH_____SS__S_______________TTGGGHHHHHHHHHHHHTT________:::::
# the chars are : b,B,E,G,H,I,S,T,C,_

# output file :
# 1a11_A:::::iiiiHHHHHHHHHHHHHHHHHHooo:::::
# 1a91_A:::::iiiiiiiHHHHHHHHHHHHHHHHHHHHHHHHHHHHHooooooooooooooHHHHHHHHHHHHHHHHHHHHHHHHHiiii:::::


import sys
import re


# global variables
min_num_tmh_residues = 7
look_for_lost_tmh_by_counting_min_num_tmh_residues = 0
tmh_residue_types = [ 'H', 'I']
save_id = []
save_opm_posn_MEMBRANE = []
save_opm_posn_INSIDE = []
save_opm_posn_OUTSIDE = []
save_stride_struc2D = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	# 1ifp_A:::::MEMBRANE-SEGMENTS,,,,,13-35,;;;;;INSIDE-SEGMENTS,,,,,36-44,;;;;;OUTSIDE-SEGMENTS,,,,,1-12,;;;;;:::::
	# 1n7l_A:::::MEMBRANE-SEGMENTS,,,,,2-2,4-5,8-8,31-53,;;;;;INSIDE-SEGMENTS,,,,,1-1,3-3,6-7,9-30,;;;;;:::::

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()
	this_i = 0
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			save_id.append( this_id )
			save_opm_posn_MEMBRANE.append( '' )
			save_opm_posn_INSIDE.append( '' )
			save_opm_posn_OUTSIDE.append( '' )
			save_stride_struc2D.append( '' )
			bit1 = bits[1]
			bits2 = bit1.rsplit( ';;;;;' )
			for segment_posn in bits2:
				if (segment_posn != ''):
					bits3 = segment_posn.rsplit( ',,,,,' )
					this_opm_posn_type = bits3[0]
					this_opm_posns = bits3[1]
					if (this_opm_posn_type == 'MEMBRANE-SEGMENTS'):
						save_opm_posn_MEMBRANE[this_i] = this_opm_posns
					else:
						if (this_opm_posn_type == 'INSIDE-SEGMENTS'):
							save_opm_posn_INSIDE[this_i] = this_opm_posns
						else:
							if (this_opm_posn_type == 'OUTSIDE-SEGMENTS'):
								save_opm_posn_OUTSIDE[this_i] = this_opm_posns
							else:
								print this_id, 'didnt have valid opm_posn', this_opm_posn_type
			this_i = this_i + 1

	# 1grm_A:::::_EEEEEEEEEEEEE__:::::
	# 1jb0_K:::::____________________HHHHHHHHHHH_____SS__S_______________TTGGGHHHHHHHHHHHHTT________:::::
	# the chars are : B,E,G,H,I,S,T,_

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_i = -1
			for i in range( 0, len(save_id) ):
				if (this_id == save_id[i]):
					this_i = i
			if (this_i == -1):
				print this_id, 'in stride_struc2D was not found'
			else:
				this_stride_struc2D = bits[1]
				save_stride_struc2D[this_i] = this_stride_struc2D

	# 1jb0_K:::::iiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHoooooooooooooooooooooooooooHHHHHHHHHHHHHHHiiiiiiiiii:::::

	outfile = open( output_file_name, "w" )
	for ix in range( 0, len(save_id) ):
		this_id = save_id[ix]
		this_opm_posn_MEMBRANE = save_opm_posn_MEMBRANE[ix]
		this_opm_posn_INSIDE = save_opm_posn_INSIDE[ix]
		this_opm_posn_OUTSIDE = save_opm_posn_OUTSIDE[ix]
		this_stride_struc2D = save_stride_struc2D[ix]
		this_seq_length = len(this_stride_struc2D)

		output_opm_tmh_segments_array = []
		for i in range( 0, this_seq_length ):
			output_opm_tmh_segments_array.append( '_' )

		if (this_opm_posn_INSIDE != ''):
			opm_posn_INSIDE_array = this_opm_posn_INSIDE.rsplit( ',' )
			for this_opm_posn_INSIDE in opm_posn_INSIDE_array:
				if (this_opm_posn_INSIDE != ''):
					this_from_is_below_zero = 0
					if (this_opm_posn_INSIDE[0:1] == '-'):
						this_opm_posn_INSIDE = this_opm_posn_INSIDE[1:]
						this_from_is_below_zero = 1
					bits = re.split( '-', this_opm_posn_INSIDE, maxsplit=1 )
					this_from = bits[0]
					this_to = bits[1]
					this_from = this_from.strip()
					this_from_minus_1 = int(this_from) - 1
					if (this_from_is_below_zero == 1):
						this_from_minus_1 = 0
					this_to = this_to.strip()
					this_to = int(this_to)
					if ( this_to > len(output_opm_tmh_segments_array) ):
						print_line = this_id + ' has opm_posn_INSIDE residues from ' + this_from + ' that are outside the length of the sequence'
						print print_line
					else:
						for i in range( this_from_minus_1, this_to ):
							output_opm_tmh_segments_array[i] = 'i'

		if (this_opm_posn_OUTSIDE != ''):
			opm_posn_OUTSIDE_array = this_opm_posn_OUTSIDE.rsplit( ',' )
			for this_opm_posn_OUTSIDE in opm_posn_OUTSIDE_array:
				if (this_opm_posn_OUTSIDE != ''):
					this_from_is_below_zero = 0
					if (this_opm_posn_OUTSIDE[0:1] == '-'):
						this_opm_posn_OUTSIDE = this_opm_posn_OUTSIDE[1:]
						this_from_is_below_zero = 1
					bits = re.split( '-', this_opm_posn_OUTSIDE, maxsplit=1 )
					this_from = bits[0]
					this_to = bits[1]
					this_from = this_from.strip()
					this_from_minus_1 = int(this_from) - 1
					if (this_from_is_below_zero == 1):
						this_from_minus_1 = 0
					this_to = this_to.strip()
					this_to = int(this_to)
					if ( this_to > len(output_opm_tmh_segments_array) ):
						print_line = this_id + ' has opm_posn_OUTSIDE residues from ' + this_from + ' that are outside the length of the sequence'
						print print_line
					else:
						for i in range( this_from_minus_1, this_to ):
							output_opm_tmh_segments_array[i] = 'o'

		if (this_opm_posn_MEMBRANE != ''):
			opm_posn_MEMBRANE_array = this_opm_posn_MEMBRANE.rsplit( ',' )
			for this_opm_posn_MEMBRANE in opm_posn_MEMBRANE_array:
				if (this_opm_posn_MEMBRANE != ''):
					this_from_is_below_zero = 0
					if (this_opm_posn_MEMBRANE[0:1] == '-'):
						this_opm_posn_MEMBRANE = this_opm_posn_MEMBRANE[1:]
						this_from_is_below_zero = 1
					bits = re.split( '-', this_opm_posn_MEMBRANE, maxsplit=1 )
					this_from = bits[0]
					this_to = bits[1]
					this_from = this_from.strip()
					this_from_minus_1 = int(this_from) - 1
					if (this_from_is_below_zero == 1):
						this_from_minus_1 = 0
					this_to = this_to.strip()
					this_to = int(this_to)
					if ( this_to > len(output_opm_tmh_segments_array) ):
						print_line = this_id + ' has opm_posn_MEMBRANE residues from ' + this_from + ' that are outside the length of the sequence'
						print print_line
					else:
						for i in range( this_from_minus_1, this_to ):
							output_opm_tmh_segments_array[i] = 'm'


		current_count = 0
		for i in range( 0, len(this_stride_struc2D) ):
			this_char = this_stride_struc2D[i]
			is_tmh = 0
			for tmh_char in tmh_residue_types:
				if (this_char == tmh_char):
					is_tmh = 1
			if ((is_tmh == 1) and (output_opm_tmh_segments_array[i] == 'm')):
				if (current_count > 0):
					current_count = current_count + 1
				else:
					current_count = 1
			else:
				if (look_for_lost_tmh_by_counting_min_num_tmh_residues == 1):
					if (current_count >= min_num_tmh_residues):
						this_from = i - current_count
						this_to = i
						for j in range( this_from, this_to ):
							output_opm_tmh_segments_array[j] = 'M'
					else:
						if ((current_count == (min_num_tmh_residues - 1)) and (this_char == 'T')):
							this_from = i - current_count
							this_to = i
							for j in range( this_from, this_to ):
								output_opm_tmh_segments_array[j] = 'M'
				current_count = 0
		if (look_for_lost_tmh_by_counting_min_num_tmh_residues == 1):
			if (current_count >= min_num_tmh_residues):
				this_from = len(this_stride_struc2D) - current_count
				this_to = current_count
				for j in range( this_from, this_to ):
					output_opm_tmh_segments_array[i] = 'M'

		output_opm_tmh_segments = ''
		for i in range( 0, len(this_stride_struc2D) ):
			output_opm_tmh_segments = output_opm_tmh_segments + output_opm_tmh_segments_array[i]
		output_line = this_id + ':::::' + output_opm_tmh_segments + ':::::' + "\r\n"
		outfile.write( output_line )

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

