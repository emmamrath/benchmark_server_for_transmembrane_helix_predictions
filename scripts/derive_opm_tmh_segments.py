#!/usr/bin/env python

# python ../derive_opm_tmh_segments.py alpha_polytopic_bitopic_peptide.fix_opm_TMsubunit alpha_polytopic_bitopic_peptide.fix_opm_posn alpha_polytopic_bitopic_peptide.fix_pdb_struc2D alpha_polytopic_bitopic_peptide.fix_pdb_seq alpha_polytopic_bitopic_peptide.fix_opm_tmh_segments

# python derive_opm_tmh_segments.py [input-file-1] [input-file-2] [input-file-3] [input-file-4] [output-file]

# file 1 = opm_TMsubunit :
# 1ar1_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,31-55,90-113,129-151,178-202,219-247,269-294,305-326,338-363,371-394,407-429,443-464,488-512,;;;;;:::::
# 1ar1_B:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,37-60,75-94,;;;;;:::::
# 1b9u_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,4-30,;;;;;:::::

# file 2 = opm_posn :
# 1ifp_A:::::MEMBRANE-SEGMENTS,,,,,13-35,;;;;;INSIDE-SEGMENTS,,,,,36-44,;;;;;OUTSIDE-SEGMENTS,,,,,1-12,;;;;;:::::
# 1n7l_A:::::MEMBRANE-SEGMENTS,,,,,2-2,4-5,8-8,31-53,;;;;;INSIDE-SEGMENTS,,,,,1-1,3-3,6-7,9-30,;;;;;:::::
# 1oqw_A:::::MEMBRANE-SEGMENTS,,,,,1-24,;;;;;OUTSIDE-SEGMENTS,,,,,25-144,;;;;;:::::

# file 3 = pdb_struc2D :
# 1ifp_A:::::_HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH_:::::
# 1ijd_A:::::_TTGGGGGTS_HHHHHHHHHHHHHHHHHHHHHHHHHS_SHHHHHT________:::::
# 1ijd_B:::::_____HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHTS____:::::
# 1grm_A:::::_EEEEEEEEEEEEE__:::::
# 1jb0_K:::::____________________HHHHHHHHHHH_____SS__S_______________TTGGGHHHHHHHHHHHHTT________:::::
# the chars are : B,E,G,H,I,S,T,_

# file 4 = pdb_seq :
# 1a11_A:::::GSEKMSTAISVLLAQAVFLLLTSQR:::::
# 1a91_A:::::MENLNMDLLYMAAAVMMGLAAIGAAIGIGILGGKFLEGAARQPDLIPLLRTQFFIVMGLVDAIPMIAVGLGLYVMFAVA:::::
# 1ehk_C:::::EEKPKGALAVILVLTLTILVFWLGVYAVFFARG:::::

# file 5 = pdb_seq :
# 1a11_A:::::____OOOOOOOOOOIIIIIIII___:::::
# 1a91_A:::::_______OOOOOOOOOOOOIOIIIIIIIIIIIII______________I_IIIIIIIIIIIIOOOOOOOOOOOOO____:::::

# output file :
# 1a11_A:::::iiiiHHHHHHHHHHHHHHHHHHooo:::::
# 1a91_A:::::iiiiiiiHHHHHHHHHHHHHHHHHHHHHHHHHHHHHooooooooooooooHHHHHHHHHHHHHHHHHHHHHHHHHiiii:::::


import sys
import re


# global variables
save_id = []
save_opm_TMsubunit = []
save_opm_posn_MEMBRANE = []
save_opm_posn_INSIDE = []
save_opm_posn_OUTSIDE = []
save_pdb_struc2D = []
save_pdb_seq = []
save_membrane_near = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	input_file_3 = sys.argv[3]
	input_file_4 = sys.argv[4]
	output_file_name = sys.argv[5]

	# 1ar1_A:::::OPM-TRANSMEMBRANE-SUBUNIT,,,,,31-55,90-113,129-151,178-202,219-247,269-294,305-326,338-363,371-394,407-429,443-464,488-512,;;;;;:::::

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			bit1 = bits[1]
			bits2 = bit1.rsplit( ';;;;;' )
			bit2 = bits2[0]
			bits3 = bit2.rsplit( ',,,,,' )
			this_opm_TMsubunit = bits3[1]
			save_id.append( this_id )
			save_opm_TMsubunit.append( this_opm_TMsubunit )
			save_opm_posn_MEMBRANE.append( '' )
			save_opm_posn_INSIDE.append( '' )
			save_opm_posn_OUTSIDE.append( '' )
			save_pdb_struc2D.append( '' )
			save_pdb_seq.append( '' )
			save_membrane_near.append( '' )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# 1ifp_A:::::MEMBRANE-SEGMENTS,,,,,13-35,;;;;;INSIDE-SEGMENTS,,,,,36-44,;;;;;OUTSIDE-SEGMENTS,,,,,1-12,;;;;;:::::
	# 1n7l_A:::::MEMBRANE-SEGMENTS,,,,,2-2,4-5,8-8,31-53,;;;;;INSIDE-SEGMENTS,,,,,1-1,3-3,6-7,9-30,;;;;;:::::

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
				print this_id, 'in opm_posn was not found'
			else:
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

	# 1grm_A:::::_EEEEEEEEEEEEE__:::::
	# 1jb0_K:::::____________________HHHHHHHHHHH_____SS__S_______________TTGGGHHHHHHHHHHHHTT________:::::
	# the chars are : B,E,G,H,I,S,T,_

	infile3 = open( input_file_3, "r" )
	inlines3 = infile3.readlines()
	infile3.close()
	for inline in inlines3:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_i = -1
			for i in range( 0, len(save_id) ):
				if (this_id == save_id[i]):
					this_i = i
			if (this_i == -1):
				print this_id, 'in pdb_struc2D was not found'
			else:
				this_pdb_struc2D = bits[1]
				save_pdb_struc2D[this_i] = this_pdb_struc2D

	# 1a11_A:::::GSEKMSTAISVLLAQAVFLLLTSQR:::::

	infile4 = open( input_file_4, "r" )
	inlines4 = infile4.readlines()
	infile4.close()
	for inline in inlines4:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_i = -1
			for i in range( 0, len(save_id) ):
				if (this_id == save_id[i]):
					this_i = i
			if (this_i == -1):
				print this_id, 'in pdb_seq was not found'
			else:
				this_pdb_seq = bits[1]
				save_pdb_seq[this_i] = this_pdb_seq

	# 1jb0_K:::::iiiiiiiiiiiiiiiiiiiiHHHHHHHHHHHoooooooooooooooooooooooooooHHHHHHHHHHHHHHHiiiiiiiiii:::::

	outfile = open( output_file_name, "w" )
	for ix in range( 0, len(save_id) ):
		this_id = save_id[ix]
		this_opm_TMsubunit = save_opm_TMsubunit[ix]
		this_opm_posn_MEMBRANE = save_opm_posn_MEMBRANE[ix]
		this_opm_posn_INSIDE = save_opm_posn_INSIDE[ix]
		this_opm_posn_OUTSIDE = save_opm_posn_OUTSIDE[ix]
		this_pdb_struc2D = save_pdb_struc2D[ix]
		this_pdb_seq = save_pdb_seq[ix]
		this_membrane_near = save_membrane_near[ix]
		this_seq_length = len(this_pdb_seq)
		#print this_id
		#print this_pdb_seq
		#print this_opm_TMsubunit
		#print this_opm_posn_MEMBRANE
		#print this_pdb_struc2D

		output_opm_tmh_segments_array = []
		for i in range( 0, this_seq_length ):
			output_opm_tmh_segments_array.append( '_' )
		#print 'length pdb_seq = ', len(this_pdb_seq)

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

		if (this_opm_TMsubunit != ''):
			opm_TMsubunit_array = this_opm_TMsubunit.rsplit( ',' )
			for this_opm_TMsubunit in opm_TMsubunit_array:
				if (this_opm_TMsubunit != ''):
					this_from_is_below_zero = 0
					if (this_opm_TMsubunit[0:1] == '-'):
						this_opm_TMsubunit = this_opm_TMsubunit[1:]
						this_from_is_below_zero = 1
					bits = re.split( '-', this_opm_TMsubunit, maxsplit=1 )
					this_from = bits[0]
					this_to = bits[1]
					this_from = this_from.strip()
					this_from_minus_1 = int(this_from) - 1
					if (this_from_is_below_zero == 1):
						this_from_minus_1 = 0
					this_to = this_to.strip()
					this_to = int(this_to)
					num_B = 0
					num_E = 0
					num_G = 0
					num_H = 0
					num_I = 0
					num_S = 0
					num_T = 0
					num__ = 0
					for i in range( this_from_minus_1, this_to ):
						j = i + 1
						this_char = this_pdb_struc2D[i:j]
						if (this_char != ''):
							if (this_char == 'B'):
								num_B = num_B + 1
							else:
								if (this_char == 'E'):
									num_E = num_E + 1
								else:
									if (this_char == 'G'):
										num_G = num_G + 1
									else:
										if (this_char == 'H'):
											num_H = num_H + 1
										else:
											if (this_char == 'I'):
												num_I = num_I + 1
											else:
												if (this_char == 'S'):
													num_S = num_S + 1
												else:
													if (this_char == 'T'):
														num_T = num_T + 1
													else:
														if (this_char == '_'):
															num__ = num__ + 1
														else:
															print this_id, 'got wrong char', this_char, 'in pdb_struc2D', this_pdb_struc2D
					this_opm_TMsubunit_is_TMH = 0
					this_opm_TMsubunit_is_TMBB = 0
					#if ((num_H >= num_B) and (num_H >= num_E) and (num_H >= num_G) and (num_H >= num_I) and (num_H >= num_S) and (num_H >= num_T) and (num_H >= num__)):
					#	this_opm_TMsubunit_is_TMH = 1
					highest_num = 0
					highest_num_char = ''
					num_H = num_H + num_I + num_G
					num_B = num_B + num_E
					if (num_B >= highest_num):
						highest_num = num_B
						highest_num_char = 'B'
					if (num_H >= highest_num):
						highest_num = num_H
						highest_num_char = 'H'
					if (highest_num_char == 'B'):
						this_opm_TMsubunit_is_TMBB = 1
					if (highest_num_char == 'H'):
						this_opm_TMsubunit_is_TMH = 1
					if (this_opm_TMsubunit_is_TMBB == 1):
						for i in range( this_from_minus_1, this_to ):
							output_opm_tmh_segments_array[i] = 'B'
					if (this_opm_TMsubunit_is_TMH == 1):
						for i in range( this_from_minus_1, this_to ):
							output_opm_tmh_segments_array[i] = 'H'

			output_opm_tmh_segments = ''
			for i in range( 0, len(this_pdb_seq) ):
				output_opm_tmh_segments = output_opm_tmh_segments + output_opm_tmh_segments_array[i]
			output_line = this_id + ':::::' + output_opm_tmh_segments + ':::::' + "\r\n"
			outfile.write( output_line )
			#print output_opm_tmh_segments
			#print
			#print
			#print
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

