#!/usr/bin/env python

# python ../field2_ids_choose_exist_file2_then_file3_resolution.py opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.ids_for_dupl_seqs pdbtm_feb2012.restricted_by_pdbtm.pdbtm_segments opm_feb2012.sorted_list1_from_pdb_method_resolution.pdb_method_resolution opm_feb2012.opm_category_list8_alpha_polytopic_bitopic_peptide.best_id_for_dupl_seqs

# python list_ids_having_seq.py [input-file-1] [input-file-2] [input-file-3] [output-file]

import sys
import re

# 1fjk_A:::::MDKVQYLTRSAIRRASTIEMPQQARQNLQNLFINFCLILIFLLLICIIVMLL:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::
# 1grm_B:::::VGALAVVVWLWLWLWX:::::
# 1grm_A:::::VGALAVVVWLWLWLWX:::::


# global variables
save_file2_id = []
save_file3_id = []
save_file3_resolution = []


######################################################
def find_id_in_file2( this_id ):

	is_in_file2 = 0
	for file2_id in save_file2_id:
		if (file2_id == this_id):
			is_in_file2 = 1
	return

######################################################
def get_resolution( this_id ):

	this_resolution = 999
	i = 0
	for file3_id in save_file3_id:
		if (file3_id == this_id):
			this_resolution = save_file3_resolution[i]
		i = i + 1
	return this_resolution

######################################################
def get_best_id( this_array_of_ids ):

	i = 0
	best_id = ''
	best_resolution = 9999
	best_in_file2_id = ''
	best_in_file2_resolution = 9999
	ouitput_id = ''
	output_comment = '-'
	for this_id in this_array_of_ids:
		if (this_id != ''):
			this_id_is_in_list = find_id_in_file2( this_id )
			this_id_resolution = get_resolution( this_id )
			if (best_id == ''):
				best_id = this_id
				best_resolution = this_id_resolution
			if (this_id_resolution < best_resolution):
				best_id = this_id
				best_resolution = this_id_resolution
			if (this_id_resolution < best_in_file2_resolution):
				best_in_file2_id = this_id
				best_in_file2_resolution = this_id_resolution
	if (best_in_file2_resolution <= best_resolution):
		output_id = best_in_file2_id
		if (best_in_file2_id == ''):
			output_id = best_id
			output_comment = 'BEST RESOLUTION NOT IN PDBTM;;;;;'
	else:
		output_id = best_id
		output_comment = 'BEST RESOLUTION NOT IN PDBTM;;;;;'
		if (best_in_file2_id != ''):
			output_comment = output_comment + 'BEST PDBTM ID IS,,,,,' + best_in_file2_id + ';;;;;'
	output_fields = output_id + ';;;;;' + output_comment
	return output_fields

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	input_file_3 = sys.argv[3]
	output_file_name = sys.argv[4]

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			save_file2_id.append( this_id )

	infile3 = open( input_file_3, "r" )
	inlines3 = infile3.readlines()
	infile3.close()
	for inline in inlines3:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_resolution = ''
			this_resolution_info = bits[1]
			# METHOD,,,,,X-RAY DIFFRACTION;;;;;RESOLUTION,,,,,2.80;;;;;
			# METHOD,,,,,SOLUTION NMR;;;;;
			bits2 = this_resolution_info.rsplit( ';;;;;' )
			this_resolution_method_info = bits2[0]
			bits3 = this_resolution_method_info.rsplit( ',,,,,' )
			this_resolution_method = bits3[1]
			if (this_resolution_method == 'X-RAY DIFFRACTION'):
				this_resolution_number_info = bits2[1]
				bits4 = this_resolution_number_info.rsplit( ',,,,,' )
				this_resolution = bits4[1]
			else: # this_resolution_method == 'SOLUTION NMR','ELECTRON CRYSTALLOGRAPHY'
				if (this_resolution_method == 'SOLUTION NMR'):
					this_resolution = 999
				else: # EPR, X-RAY DIFFRACTION,SOLID-STATE NMR,ELECTRON CRYSTALLOGRAPHY,ELECTRON MICROSCOPY,FIBER DIFFRACTION,INFRARED SPECTROSCOPY
					this_resolution = 9999
			save_file3_id.append( this_id )
			save_file3_resolution.append( this_resolution )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile = open( output_file_name, "w" )
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_seq = bits[0]
			this_list_of_ids = bits[1]
			this_array_of_ids = this_list_of_ids.rsplit( ',' )
			output_fields = get_best_id( this_array_of_ids )
			output_line =  this_seq + ':::::' + output_fields + ':::::' + "\r\n"
			outfile.write( output_line )
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

