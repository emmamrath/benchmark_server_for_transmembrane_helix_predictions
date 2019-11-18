#!/usr/bin/env python

# python ../create_opm_pdbtm_topology_files_all_types.py alpha_polytopic_bitopic.opm_tm_segments alpha_polytopic_bitopic.opm_adj_tm_segments alpha_polytopic_bitopic.pdbtm_segments alpha_polytopic_bitopic.opm_tm_segments_topology alpha_polytopic_bitopic.opm_adj_tm_segments_topology alpha_polytopic_bitopic.pdbtm_segments_topology

# python create_opm_pdbtm_topology_files_all_types.py [input-file-1-opm_tm_segments_topology] [input-file-2-opm_adj_tm_segments_topology] [input-file-3-pdbtm_segments_topology] [output-file-1-opm_tm_segments_topology] [output-file-2-opm_adj_tm_segments_topology] [output-file-3-pdbtm_segments_topology]


# alpha_polytopic_bitopic.opm_tm_segments
#1a11_A:::::iiiiMMMMMMMMMMMMMMMMMMooo:::::
#1a91_A:::::iiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMMMMMoooooooooooomoMMMMMMMMMMMMMMMMMMMMMMMMMiiii:::::
#1afo_A:::::_________________________________________________________________mmooomMMMMMMMMMMMMMMMMMMMMMMiiiiiiii:::::
#1ar1_A:::::________________mmmiimmmiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMoooooooooooommooooooooooooMMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiimmMMMMMMMMMMMMMMMMMMMMMMoooooooooooMMMMMMMMMMMMMMMMMMMMMMMMMMiiiiimmMMMMMMMMMMMMMMMMMMMMMMMMmoooooooooooMMMMMMMMMMMMMMMMMMMMMMMmiiiiiiiiiimmMMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii_____________:::::

# alpha_polytopic_bitopic.opm_adj_tm_segments
#1ijd_A:::::_iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMmoommommo_______:::::
#1ijd_B:::::iiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMmom__:::::
#1j4n_A:::::iiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooMMMMMMMMMMMMMMMMMMMMmmmmmmmmMMMMMMMMMMmmmmMMMMMMMMMMMMMMMMMMMMMMMomooooooooommmoooooooommmMMMMMMMMMMMMMMMMMMMmiiiiimmmmMMMMMMMMMMMMMMMMMoommmmmmMMMMMMMMMMooooooMMMMMMMMMMMMMMMMMMMMMMmmmiiiiiiiiiiimii______________________:::::
#1jb0_A:::::____________iiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiimiimmimmiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooooooooooooooooooooooooooooooomoooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooooooooooooooooo___ommmmoooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMmmmiiimmmiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooomommmmoomoooooooooooooooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooooooooooooooooomoooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMmmmimmiimmiimiiiiiiiiiiiiiiiiimiiMMMMMMMMMMMMMMMMMMMMMoooooooooo:::::

# alpha_polytopic_bitopic.pdbtm_segments
#1h68_A:::::ALPHA__side_1,,,,,28-35,92-95,147-159,216-219,;;;;;ALPHA__side_2,,,,,2-5,57-69,118-124,180-193,;;;;;ALPHA__alpha_helix,,,,,6-27,36-56,70-91,96-117,125-146,160-179,194-215,;;;;;ALPHA__unknown,,,,,1-1,220-239,;;;;;:::::
#1h6i_A:::::ALPHA__side_1,,,,,9-13,69-71,87-93,157-167,232-233,;;;;;ALPHA__side_2,,,,,33-49,115-138,185-187,203-209,;;;;;ALPHA__alpha_helix,,,,,14-32,50-68,94-114,139-156,168-184,210-231,;;;;;ALPHA__membrane_loop,,,,,72-86,188-202,;;;;;ALPHA__unknown,,,,,1-8,234-269,;;;;;:::::
#1hgz_A:::::-:::::

# 1j4n_A : Crystal Structure of the AQP1 water channel
#SEQUENCE :MASEFKKKLFWRAVVAEFLAMILFIFISIGSALGFHYPIKSNQTTGAVQDNVKVSLAFGLSIATLAQSVGHISGAHLNPAVTLGLLLSCQISVLRAIMYIIAQCVGAIVATAILSGITSSLPDNSLGLNALAPGVNSGQGLGIEIIGTLQLVLCVLATTDRRRRDLGGSGPLAIGFSVALGHLLAIDYTGCGINPARSFGSSVITHNFQDHWIFWVGPFIGAALAVLIYDFILAPRSSDLTDRVKVWTSGQVEEYDLDADDINSRVEMKPK 
#OPMADJ   :iiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooMMMMMMMMMMMMMMMMMMMMmmmmmmmmMMMMMMMMMMmmmmMMMMMMMMMMMMMMMMMMMMMMMomooooooooommmoooooooommmMMMMMMMMMMMMMMMMMMMmiiiiimmmmMMMMMMMMMMMMMMMMMoommmmmmMMMMMMMMMMooooooMMMMMMMMMMMMMMMMMMMMMMmmmiiiiiiiiiiimii______________________ 
#OPM      :iiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMooooooooooooooooooMMMMMMMMMMMMMMMMMMMMmmmmmmmmMMMMMMMMMMmmmmMMMMMMMMMMMMMMMMMMMMMMMomooooooooommmoooooooommmMMMMMMMMMMMMMMMMMMMmiiiiimmmmMMMMMMMMMMMMMMMMMoommmmmmMMMMMMMMMMooooooMMMMMMMMMMMMMMMMMMMMMMmmmiiiiiiiiiiimii______________________ 
#PDBTM    :11111111111HHHHHHHHHHHHHHHHHHHHH2222222222222222222HHHHHHHHHHHHHHHHHHH11LLLLLLLLLLLLLLLL11111HHHHHHHHHHHHHHHHHHHHH222222222222222222222222222HHHHHHHHHHHHHHHHHH111111111HHHHHHHHHHHHHHHHHH2222LLLLLLLLLLLLLL2222222HHHHHHHHHHHHHHHHHHHHHH1111111111111111UUUUUUUUUUUUUUUUUUUUUU 


import sys
import re


# global variables
save_id = []
opm_tm_segments_topology = []
opm_adj_tm_segments_topology = []
pdbtm_segments_topology = []


######################################################
def convert_pdbtm_segments_to_seq( this_segments, this_length ):

	if (this_segments == '-'):
		this_seq = '-'

	else:
		seq_array = []
		for i in range( 0, this_length ):
			seq_array.append( '_' )

		bits = this_segments.rsplit( ';;;;;' )
		for bit in bits:
			if (bit != ''):
				bits2 = bit.split( ',,,,,' )
				this_type = bits2[0]
				bits2b = this_type.split( '__' )
				this_type2 = bits2b[1]
				this_type_segments = bits2[1]
				this_char = '_'
				if (this_type2 == 'side_1'):
					this_char = '1'
				elif (this_type2 == 'side_2'):
					this_char = '2'
				elif (this_type2 == 'alpha_helix'):
					this_char = 'H'
				elif (this_type2 == 'beta_strand'):
					this_char = 'B'
				elif (this_type2 == 'coil'):
					this_char = 'C'
				elif (this_type2 == 'membrane_inside'):
					this_char = 'I'
				elif (this_type2 == 'membrane_loop'):
					this_char = 'L'
				elif (this_type2 == 'unknown'):
					this_char = 'U'
				bits3 = this_type_segments.rsplit( ',' )
				for bit3 in bits3:
					if (bit3 != ''):
						if ( bit3.find('--') == -1 ):
							this_from = -1
							if ( bit3[0:1] == '-' ):
								bit3 = bit3[1:]
							bits4 = bit3.split( '-' )
							if (this_from == -1):
								this_from = int(bits4[0])
							if (this_from < 1):
								this_from = 1
							this_to = int(bits4[1])
							#print str(this_from), str(this_to), this_char
							for i in range( (this_from - 1), this_to ):
								if (i < len(seq_array)):
									seq_array[i] = this_char

		this_seq = ''
		#print seq_array
		for i in range( 0, this_length ):
			this_seq = this_seq + seq_array[i]

	return this_seq

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	input_file_3 = sys.argv[3]
	output_file_1 = sys.argv[4]
	output_file_2 = sys.argv[5]
	output_file_3 = sys.argv[6]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# read in alpha_polytopic_bitopic.opm_tm_segments

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			save_id.append( this_id )
			opm_tm_segments_topology.append( this_seq )
			opm_adj_tm_segments_topology.append( '' )
			pdbtm_segments_topology.append( '' )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	# read in alpha_polytopic_bitopic.opm_adj_tm_segments

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			i = 0
			found_id = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					found_id = 1
					opm_adj_tm_segments_topology[i] = this_seq
				i = i + 1
			if (found_id == 0):
				print 'opm_adj_tm_segments_topology has', this_id, 'and opm_tm_segments_topology does not have it'

	infile3 = open( input_file_3, "r" )
	inlines3 = infile3.readlines()
	infile3.close()

	# read in alpha_polytopic_bitopic.pdbtm_segments

	for inline in inlines3:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_segments = bits[1]
			i = 0
			found_id = 0
			for this_save_id in save_id:
				if (this_save_id == this_id):
					found_id = 1
					this_length = len( opm_tm_segments_topology[i] )
					this_seq = convert_pdbtm_segments_to_seq( this_segments, this_length )
					pdbtm_segments_topology[i] = this_seq
				i = i + 1
			if (found_id == 0):
				print 'pdbtm_segments_topology has', this_id, 'and opm_tm_segments_topology does not have it'

	i = 0
	for this_save_id in save_id:
		if (opm_adj_tm_segments_topology[i] == ''):
			print 'opm_tm_segments_topology has', this_id, 'and opm_adj_tm_segments_topology does not have it'
		if (pdbtm_segments_topology[i] == ''):
			print 'opm_tm_segments_topology has', this_id, 'and pdbtm_segments_topology does not have it'
		i = i + 1

	i = 0
	for this_save_id in save_id:

		this_opm_tm_segments_topology = opm_tm_segments_topology[i]
		this_opm_adj_tm_segments_topology = opm_adj_tm_segments_topology[i]
		this_pdbtm_segments_topology = pdbtm_segments_topology[i]
		this_opm_tm_segments_topology_array = []
		this_opm_adj_tm_segments_topology_array = []
		this_pdbtm_segments_topology_array = []

		for j in range( 0, len(opm_tm_segments_topology[i]) ):
			k = j + 1
			this_opm_tm_segments_topology_array.append( this_opm_tm_segments_topology[j:k] )
			this_opm_adj_tm_segments_topology_array.append( this_opm_adj_tm_segments_topology[j:k] )
			this_pdbtm_segments_topology_array.append( this_pdbtm_segments_topology[j:k] )

		# first pass at filling opm_tm_segments_topology

		prev_char = ''
		for j in range( 0, len(this_opm_tm_segments_topology_array) ):
			this_char = this_opm_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_tm_segments_topology_array[j] = this_char
			prev_char = this_char
		prev_char = ''
		for j in range( (len(this_opm_tm_segments_topology_array) - 1), -1, -1 ):
			this_char = this_opm_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_tm_segments_topology_array[j] = this_char
			prev_char = this_char

		# first pass at filling opm_adj_tm_segments_topology

		prev_char = ''
		for j in range( 0, len(this_opm_adj_tm_segments_topology_array) ):
			this_char = this_opm_adj_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_adj_tm_segments_topology_array[j] = this_char
			prev_char = this_char
		prev_char = ''
		for j in range( (len(this_opm_adj_tm_segments_topology_array) - 1), -1, -1 ):
			this_char = this_opm_adj_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_adj_tm_segments_topology_array[j] = this_char
			prev_char = this_char

		# assign sides 1 and 2 of pdbtm_segments_topology to inside or outside by comparing to opm_adj_tm_segments

		num_segments_1_is_inside = 0
		num_segments_2_is_inside = 0
		for j in range( 0, len(this_pdbtm_segments_topology_array) ):
			this_char = this_pdbtm_segments_topology_array[j]
			opm_adj_char = this_opm_adj_tm_segments_topology_array[j]
			if (this_char == '1'):
				if (opm_adj_char == 'i'):
					num_segments_1_is_inside = num_segments_1_is_inside + 1
				elif (opm_adj_char == 'o'):
					num_segments_2_is_inside = num_segments_2_is_inside + 1
			elif (this_char == '2'):
				if (opm_adj_char == 'i'):
					num_segments_2_is_inside = num_segments_2_is_inside + 1
				elif (opm_adj_char == 'o'):
					num_segments_1_is_inside = num_segments_1_is_inside + 1
		pdbtm_1_is = ''
		if (num_segments_1_is_inside > num_segments_2_is_inside):
			pdbtm_1_is = 'inside'
		elif (num_segments_2_is_inside > num_segments_1_is_inside):
			pdbtm_1_is = 'outside'

		if (pdbtm_1_is == ''):
			for j in range( 0, len(this_pdbtm_segments_topology_array) ):
				this_pdbtm_segments_topology_array[j] = '_'
		else:
			for j in range( 0, len(this_pdbtm_segments_topology_array) ):
				this_char = this_pdbtm_segments_topology_array[j]
				new_char = this_char
				if (this_char == '1'):
					if (pdbtm_1_is == 'inside'):
						new_char = 'i'
					elif (pdbtm_1_is == 'outside'):
						new_char = 'o'
				elif (this_char == '2'):
					if (pdbtm_1_is == 'inside'):
						new_char = 'o'
					elif (pdbtm_1_is == 'outside'):
						new_char = 'i'
				this_pdbtm_segments_topology_array[j] = new_char

		# fill any remaining m in opm_tm_segments_topology by inspecting pdbtm

		prev_char = ''
		curr_from = 0
		curr_to = 0
		for j in range( 0, len(this_opm_tm_segments_topology_array) ):
			this_char = this_opm_tm_segments_topology_array[j]
			if ((this_char == 'm') and (prev_char == 'M')):
				curr_from = j
			if ((this_char == 'M') and (prev_char == 'm')):
				curr_to = j - 1
				num_pdbtm_i = 0
				num_pdbtm_o = 0
				for k in range( curr_from, (curr_to + 1) ):
					if (this_pdbtm_segments_topology_array[k] == 'i'):
						num_pdbtm_i = num_pdbtm_i + 1
					elif (this_pdbtm_segments_topology_array[k] == 'o'):
						num_pdbtm_o = num_pdbtm_o + 1
				if (num_pdbtm_i >  num_pdbtm_o):
					this_opm_tm_segments_topology_array[j] = 'i'
				elif (num_pdbtm_o >  num_pdbtm_i):
					this_opm_tm_segments_topology_array[j] = 'o'
			prev_char = this_char

		# fill any remaining m in opm_adj_tm_segments_topology by inspecting pdbtm

		prev_char = ''
		curr_from = 0
		curr_to = 0
		for j in range( 0, len(this_opm_adj_tm_segments_topology_array) ):
			this_char = this_opm_adj_tm_segments_topology_array[j]
			if ((this_char == 'm') and (prev_char == 'M')):
				curr_from = j
			if ((this_char == 'M') and (prev_char == 'm')):
				curr_to = j - 1
				num_pdbtm_i = 0
				num_pdbtm_o = 0
				for k in range( curr_from, (curr_to + 1) ):
					if (this_pdbtm_segments_topology_array[k] == 'i'):
						num_pdbtm_i = num_pdbtm_i + 1
					elif (this_pdbtm_segments_topology_array[k] == 'o'):
						num_pdbtm_o = num_pdbtm_o + 1
				if (num_pdbtm_i >  num_pdbtm_o):
					this_opm_adj_tm_segments_topology_array[j] = 'i'
				elif (num_pdbtm_o >  num_pdbtm_i):
					this_opm_adj_tm_segments_topology_array[j] = 'o'
			prev_char = this_char

		# third pass at filling opm_tm_segments_topology

		prev_char = ''
		for j in range( 0, len(this_opm_tm_segments_topology_array) ):
			this_char = this_opm_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_tm_segments_topology_array[j] = this_char
			prev_char = this_char
		prev_char = ''
		for j in range( (len(this_opm_tm_segments_topology_array) - 1), -1, -1 ):
			this_char = this_opm_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_tm_segments_topology_array[j] = this_char
			prev_char = this_char

		# third pass at filling opm_adj_tm_segments_topology

		prev_char = ''
		for j in range( 0, len(this_opm_adj_tm_segments_topology_array) ):
			this_char = this_opm_adj_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_adj_tm_segments_topology_array[j] = this_char
			prev_char = this_char
		prev_char = ''
		for j in range( (len(this_opm_adj_tm_segments_topology_array) - 1), -1, -1 ):
			this_char = this_opm_adj_tm_segments_topology_array[j]
			if ((this_char == 'm') and ((prev_char == 'i') or (prev_char == 'o'))):
				this_char = prev_char
				this_opm_adj_tm_segments_topology_array[j] = this_char
			prev_char = this_char

		# serialise the arrays to strings

		this_seq = ''
		for j in range( 0, len(this_opm_tm_segments_topology_array) ):
			this_seq = this_seq + this_opm_tm_segments_topology_array[j]
		opm_tm_segments_topology[i] = this_seq

		this_seq = ''
		for j in range( 0, len(this_opm_adj_tm_segments_topology_array) ):
			this_seq = this_seq + this_opm_adj_tm_segments_topology_array[j]
		opm_adj_tm_segments_topology[i] = this_seq

		this_seq = ''
		for j in range( 0, len(this_pdbtm_segments_topology_array) ):
			this_seq = this_seq + this_pdbtm_segments_topology_array[j]
		pdbtm_segments_topology[i] = this_seq

		i = i + 1

	# write output for alpha_polytopic_bitopic.opm_tm_segments

	outfile1 = open( output_file_1, "w" )
	i = 0
	for this_seq in opm_tm_segments_topology:
		this_id = save_id[i]
		output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
		outfile1.write( output_line )
		i = i + 1
	outfile1.close()

	# write output for alpha_polytopic_bitopic.opm_adj_tm_segments

	outfile2 = open( output_file_2, "w" )
	i = 0
	for this_seq in opm_adj_tm_segments_topology:
		this_id = save_id[i]
		output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
		outfile2.write( output_line )
		i = i + 1
	outfile2.close()

	# write output for alpha_polytopic_bitopic.pdbtm_segments

	outfile3 = open( output_file_3, "w" )
	i = 0
	for this_seq in pdbtm_segments_topology:
		this_id = save_id[i]
		output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
		outfile3.write( output_line )
		i = i + 1
	outfile3.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

