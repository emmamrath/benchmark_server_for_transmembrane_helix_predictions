#!/usr/bin/env python

# python ../remove_duplicate_seqs_by_resolution.py opm_feb2012.list4b.pdb_seq opm_feb2012.list2_from_pdb_method_resolution.pdb_method_resolution opm_feb2012.list5.pdb_seq

# python remove_duplicate_seqs_by_resolution.py [input-file-1] [input-file-2] [output-file]

import sys
import re


# global variables
save_id = []
save_method = []
save_resolution = []
save_seq = []
save_output = []


######################################################
def look_for_id( this_id ):

	found_id = -1
	i = 0
	for this_save_id in save_id:
		if (this_save_id == this_id):
			found_id = i
		i = i + 1
	return found_id

######################################################
def find_best_id( this_seq ):

	best_id = ''
	best_method = ''
	best_resolution = ''

	i = 0
	for this_save_seq in save_seq:
		if (this_save_seq == this_seq):
			this_save_id = save_id[i]
			this_save_method = save_method[i]
			this_save_resolution = save_resolution[i]
			if (best_id == ''):
				best_id = this_save_id
				best_method = this_save_method
				best_resolution = this_save_resolution
			else:
				if (this_save_resolution < best_resolution):
					best_id = this_save_id
					best_method = this_save_method
					best_resolution = this_save_resolution
		i = i + 1

	return best_id

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_name = sys.argv[3]

	outfile = open( output_file_name, "w" )

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	# save the method and resolution of every id (pdbid+chain)

	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_method = bits[1]
			this_resolution = 999
			this_seq = ''
			this_output = 0
			if ( len(bits) >= 3 ):
				if (bits[2] != ''):
					this_resolution = bits[2]
			save_id.append( this_id )
			save_method.append( this_method )
			save_resolution.append( this_resolution )
			save_seq.append( this_seq )
			save_output.append( this_output )

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	# save the sequence of every id (pdbid+chain)

	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			bits = inline.rsplit( ':::::' )
			this_id = bits[0]
			this_seq = bits[1]
			found_id = look_for_id( this_id )
			if (found_id != -1):
				save_seq[found_id] = this_seq

	# flag which seq is the best according to resolution criteria

	i = 0
	for this_id in save_id:
		this_seq = save_seq[i]
		best_id = find_best_id( this_seq )
		if (best_id == this_id):
			save_output[i] = 1
		i = i + 1

	# output the best seq

	i = 0	
	for this_id in save_id:
		this_seq = save_seq[i]
		if (this_seq != ''):
			this_output = save_output[i]
			if (this_output == 1):
				output_line = this_id + ':::::' + this_seq + ':::::' + "\r\n"
				outfile.write( output_line )
		i = i + 1
	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

