#!/usr/bin/env python

# python print_pymol_for_MMMMMM_line.py C __________immmMMMMMMMMMMMMMMMMMMMMMoooooMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMmooooooooooooooooooooooooooooooooooooooooMMMMMMMMMMMMMMMMMMMiiiiimmiMMMMMMMMMMMMMMMMMMMMMMooooooooooooooommoooooomoMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMoooooooooooooooooommmommoomooomooommMMMMMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMmooMMMMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMMMMMMooooooooooMMMMMMMMMMMMMMMMMMMMiiiiiiiii________

# python print_pymol_for_MMMMMM_line.py [input-line]

import sys
import re


# global variables
save_id = []
save_mpstruc_category = []
save_opm_tm_segments = []
save_stride_struc2D = []

######################################################
def main():

	input_chain = sys.argv[1]
	input_line = sys.argv[2]

	print 'bg white'
	print 'select dum, resn dum'
	print 'select ' + input_chain + input_chain + ', chain ' + input_chain
	print 'center chain ' + input_chain
	print 'hide all'
	print 'show nb_spheres, dum'
	print 'show cartoon, chain ' + input_chain

	prev_char = ''
	this_start = ''
	this_end = ''
	for i in range( 0, len(input_line) ):
		j = i + 1
		this_char = input_line[i:j]
		if ((this_char == 'M') and (prev_char != 'M')):
			this_start = j
		if ((this_char != 'M') and (prev_char == 'M')):
			this_end = i
			output_line = 'select ' + str(this_start) + '-' + str(this_end) + ', chain ' + input_chain + ' and resi ' + str(this_start) + '-' + str(this_end)
			print output_line
			this_start = ''
			this_end = ''
		prev_char = this_char



if __name__ == "__main__":
	# Someone is launching this directly
	main()

