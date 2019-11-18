#!/usr/bin/env python

# python read_one_opm_pdb.py 1r3j.pdb alpha_polytopic.noseq.opm_feb2012 opm_feb2012.alpha_polytopic.membrane_posn opm_feb2012.alpha_polytopic.seq opm_feb2012.alpha_polytopic.opm_membrane_near

# python read_list_opm_pdb.py [input-list-of-pdb-files] [output-file-for-segments] [output-file-for-sequences] [output-file-for-membrane-near]

# This program reads a list of PDB files from OPM where the membrane position has been defined by DUM records.
# It outputs the segments as being M (membrane), I (inside the cell or on the inside side of the membrane)
# or O (outside the cell or on the outside side of the membrane).
# This program does not judge the secondary structure of the segments,
# which could be alpha helices, beta sheets or coils.

import sys
import os


######################################################
# global variables

# ATOM      1  N   ALA A  23     -16.378   6.310 -19.289  1.00181.62           N  
# ATOM      2  CA  ALA A  23     -15.241   5.381 -19.555  1.00181.62           C  
# ATOM      3  C   ALA A  23     -15.333   4.149 -18.663  1.00181.62           C  
# ATOM      4  O   ALA A  23     -16.317   3.409 -18.705  1.00181.62           O  
# ATOM      5  CB  ALA A  23     -13.914   6.100 -19.317  1.00 74.09           C  
# ATOM      6  N   LEU A  24     -14.296   3.936 -17.858  1.00163.39           N  
# ATOM      7  CA  LEU A  24     -14.240   2.806 -16.940  1.00163.39           C  
# ATOM      8  C   LEU A  24     -13.111   3.009 -15.938  1.00163.39           C  
# ATOM      9  O   LEU A  24     -13.348   3.278 -14.765  1.00163.39           O  
# ATOM     10  CB  LEU A  24     -14.010   1.503 -17.709  1.00126.45           C  
# ATOM     11  CG  LEU A  24     -13.979   0.232 -16.856  1.00126.45           C  
# ATOM     12  CD1 LEU A  24     -15.243  -0.576 -17.100  1.00126.45           C  
# ATOM     13  CD2 LEU A  24     -12.744  -0.586 -17.192  1.00126.45           C  
# HETATM 2146  N   DUM  2146     -26.000   0.000 -16.200
# HETATM 2146  O   DUM  2146     -26.000   0.000  16.200
# HETATM 2222  N   DUM  2222     -24.000 -10.000 -16.200
# HETATM 2222  O   DUM  2222     -24.000 -10.000  16.200
# HETATM 2223  N   DUM  2223     -24.000  -8.000 -16.200
# HETATM 2223  O   DUM  2223     -24.000  -8.000  16.200
# HETATM 2224  N   DUM  2224     -24.000  -6.000 -16.200
# HETATM 2224  O   DUM  2224     -24.000  -6.000  16.200


######################################################
def main():

	global outfile
	input_list_file = sys.argv[1]
	output_file_name_segments = sys.argv[2]
	output_file_name_sequences = sys.argv[3]
	output_file_name_membrane_near = sys.argv[4]

	outfile_segments = open( output_file_name_segments, "w" )
	outfile_sequences = open( output_file_name_sequences, "w" )
	outfile_membrane_near = open( output_file_name_membrane_near, "w" )

	infile = open( input_list_file, "r" )
	inlines = infile.readlines()
	infile.close()

	for inline in inlines:
		inline = inline.strip()
		if (inline != ''):
			this_pdb_file = inline
			this_command = 'python ../../../read_one_opm_pdb.py ' + this_pdb_file + ' ' + output_file_name_segments + ' ' + output_file_name_sequences + ' ' + output_file_name_membrane_near
			os.system( this_command )

	outfile_segments.close()
	outfile_sequences.close()
	outfile_membrane_near.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

