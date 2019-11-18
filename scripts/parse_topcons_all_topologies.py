#!/usr/bin/env python

# python parse_topcons_all_topologies.py alpha_polytopic_bitopic_noblanks.1line.topcons alpha_polytopic_bitopic_noblanks.1line.scampi_seq_topcons alpha_polytopic_bitopic_noblanks.1line.scampi_msa_topcons alpha_polytopic_bitopic_noblanks.1line.prodiv_topcons alpha_polytopic_bitopic_noblanks.1line.pro_topcons alpha_polytopic_bitopic_noblanks.1line.octopus_topcons alpha_polytopic_bitopic_noblanks.1line.topcons 

# python parse_topcons_all_topologies.py [input-file] [output-file-scampi-seq] [output-file-scampi-msa] [output-file-prodiv] [output-file-pro] [output-file-octopus] [output-file-topcons]


import sys
import re


# input file :
# 2b6o_A:::::
# ##############################################################################
# TOPCONS result file
# Generated from http://topcons.cbr.su.se/ at 2012-03-13 13:18:40
# Total request time: 19.98 seconds.
# ##############################################################################
# 
# 
# Sequence name: query_sequence
# Sequence length: 263 aa.
# Sequence:
# MWELRSASFWRAIFAEFFATLFYVFFGLGASLRWAPGPLHVLQVALAFGLALATLVQAVG
# HISGAHVNPAVTFAFLVGSQMSLLRAICYVVAQLLGAVAGAAVLYSVTPPAVRGNLALNT
# LHPGVSVGQATIVEIFLTLQFVLCIFATYDERRNGRLGSVALAVGFSLTLGHLFGMYYTG
# AGMNPARSFAPAILTRNFTNHWVYWVGPVIGAGLGSLLYDFLLFPRLKSVSERLSILKGT
# RPSESNGQPEVTGEPVELKTQAL
# 
# SCAMPI-seq predicted topology:
# iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoMMMMMMMMMMMMMMMMMMMMMi
# iiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiMMMMMMMMMMMMMMMMMMMMMooo
# oooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# SCAMPI-msa predicted topology:
# iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooMMMMMMMMMMMMMMMMMMMMM
# iiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiMMMMMMMMMMMMMMMMMMMMMo
# oooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# PRODIV predicted topology:
# iiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMoooooooooooooMMMMMMMMMMMMMMM
# MMMMMMiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMM
# MooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# PRO predicted topology:
# iiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooMMMMMMMMMMMMMMM
# MMMMMMiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMM
# ooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# OCTOPUS predicted topology:
# iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooMMMMMMMMMMMMMMMMMMMMM
# iiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiMMMMMMMMMMMMMMMMMMMMMo
# oooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# TOPCONS predicted topology:
# iiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooMMMMMMMMMMMMMMMMMMMMM
# iiiiiiiiiiiiiiiiiiiiiiiiiiMMMMMMMMMMMMMMMMMMMMMooooooooooooo
# oooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiMMMMMMMMMMMMMMMMMMMMMo
# oooooooooooooooooooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiiiiiiiii
# iiiiiiiiiiiiiiiiiiiiiii
# 
# Predicted Z-coordinates (left column=sequence position; right column=Z-coordinate)
# 1 	 20.50
# 2 	 21.46
# 3 	 23.05
# 4 	 22.63
# ...
# 
# Predicted Delta-G-values (kcal/mol) (left column=sequence position; right column=Delta-G)
# 11	5.13
# 12	4.59
# 13	3.76
# 14	2.73
# ...
# 
# Predicted TOPCONS reliability (left column=sequence position; right column=reliability)
# 11	0.914
# 12	0.914
# 13	0.914
# 14	0.914
# ...
# 2zt9_G:::::
# ##############################################################################
# TOPCONS result file
# Generated from http://topcons.cbr.su.se/ at 2012-03-13 13:19:03
# Total request time: 4.01 seconds.
# ##############################################################################
# 
# 
# Sequence name: query_sequence
# Sequence length: 37 aa.
# Sequence:
# MVEPLLSGIVLGLIVVTLAGLFYAAYKQYKRPNELGG
# 
# SCAMPI-seq predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# SCAMPI-msa predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# PRODIV predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# PRO predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# OCTOPUS predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# TOPCONS predicted topology:
# ooooMMMMMMMMMMMMMMMMMMMMMiiiiiiiiiiii
# 
# Predicted Z-coordinates (left column=sequence position; right column=Z-coordinate)
# 1 	 13.82
# 2 	 12.60


# global variables

######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	output_file_name_scampi_seq = sys.argv[2]
	output_file_name_scampi_msa = sys.argv[3]
	output_file_name_prodiv = sys.argv[4]
	output_file_name_pro = sys.argv[5]
	output_file_name_octopus = sys.argv[6]
	output_file_name_topcons = sys.argv[7]

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile_scampi_seq = open( output_file_name_scampi_seq, "w" )
	outfile_scampi_msa = open( output_file_name_scampi_msa, "w" )
	outfile_prodiv = open( output_file_name_prodiv, "w" )
	outfile_pro = open( output_file_name_pro, "w" )
	outfile_octopus = open( output_file_name_octopus, "w" )
	outfile_topcons = open( output_file_name_topcons, "w" )

	curr_id = ''
	in_an_id = 0
	in_scampi_seq = 0
	in_scampi_msa = 0
	in_prodiv = 0
	in_pro = 0
	in_octopus = 0
	in_topcons = 0
	output_for_scampi_seq = ''
	output_for_scampi_msa = ''
	output_for_prodiv = ''
	output_for_pro = ''
	output_for_octopus = ''
	output_for_topcons = ''
	for inline in inlines1:
		inline = inline.strip()
		if (inline != ''):
			# 2bg9_B:::::
			inline_predicted_text_found = ''
			if (len(inline) >= 9):
				inline_predicted_text_found = inline[0:9]
			if (len(inline) == 11):
				if (inline[6:11] == ':::::'):
					if (in_an_id == 1):
						if (output_for_scampi_seq == ''):
							output_for_scampi_seq = '-'
						if (output_for_scampi_msa == ''):
							output_for_scampi_msa = '-'
						if (output_for_prodiv == ''):
							output_for_prodiv = '-'
						if (output_for_pro == ''):
							output_for_pro = '-'
						if (output_for_octopus == ''):
							output_for_octopus = '-'
						if (output_for_topcons == ''):
							output_for_topcons = '-'
						output_line_for_scampi_seq = curr_id + ':::::' + output_for_scampi_seq + ':::::' + "\r\n"
						output_line_for_scampi_msa = curr_id + ':::::' + output_for_scampi_msa + ':::::' + "\r\n"
						output_line_for_prodiv = curr_id + ':::::' + output_for_prodiv + ':::::' + "\r\n"
						output_line_for_pro = curr_id + ':::::' + output_for_pro + ':::::' + "\r\n"
						output_line_for_octopus = curr_id + ':::::' + output_for_octopus + ':::::' + "\r\n"
						output_line_for_topcons = curr_id + ':::::' + output_for_topcons + ':::::' + "\r\n"
						outfile_scampi_seq.write( output_line_for_scampi_seq )
						outfile_scampi_msa.write( output_line_for_scampi_msa )
						outfile_prodiv.write( output_line_for_prodiv )
						outfile_pro.write( output_line_for_pro )
						outfile_octopus.write( output_line_for_octopus )
						outfile_topcons.write( output_line_for_topcons )
					curr_id = inline[0:6]
					in_an_id = 1
					in_scampi_seq = 0
					in_scampi_msa = 0
					in_prodiv = 0
					in_pro = 0
					in_octopus = 0
					in_topcons = 0
					output_for_scampi_seq = ''
					output_for_scampi_msa = ''
					output_for_prodiv = ''
					output_for_pro = ''
					output_for_octopus = ''
					output_for_topcons = ''
			elif (inline == 'SCAMPI-seq predicted topology:'):
				in_scampi_seq = 1
				in_scampi_msa = 0
				in_prodiv = 0
				in_pro = 0
				in_octopus = 0
				in_topcons = 0
			elif (inline == 'SCAMPI-msa predicted topology:'):
				in_scampi_seq = 0
				in_scampi_msa = 1
				in_prodiv = 0
				in_pro = 0
				in_octopus = 0
				in_topcons = 0
			elif (inline == 'PRODIV predicted topology:'):
				in_scampi_seq = 0
				in_scampi_msa = 0
				in_prodiv = 1
				in_pro = 0
				in_octopus = 0
				in_topcons = 0
			elif (inline == 'PRO predicted topology:'):
				in_scampi_seq = 0
				in_scampi_msa = 0
				in_prodiv = 0
				in_pro = 1
				in_octopus = 0
				in_topcons = 0
			elif (inline == 'OCTOPUS predicted topology:'):
				in_scampi_seq = 0
				in_scampi_msa = 0
				in_prodiv = 0
				in_pro = 0
				in_octopus = 1
				in_topcons = 0
			elif (inline == 'TOPCONS predicted topology:'):
				in_scampi_seq = 0
				in_scampi_msa = 0
				in_prodiv = 0
				in_pro = 0
				in_octopus = 0
				in_topcons = 1
			elif (inline_predicted_text_found == 'Predicted'):
				in_scampi_seq = 0
				in_scampi_msa = 0
				in_prodiv = 0
				in_pro = 0
				in_octopus = 0
				in_topcons = 0
			elif (in_scampi_seq == 1):
				output_for_scampi_seq = output_for_scampi_seq + inline
			elif (in_scampi_msa == 1):
				output_for_scampi_msa = output_for_scampi_msa + inline
			elif (in_prodiv == 1):
				output_for_prodiv = output_for_prodiv + inline
			elif (in_pro == 1):
				output_for_pro = output_for_pro + inline
			elif (in_octopus == 1):
				output_for_octopus = output_for_octopus + inline
			elif (in_topcons == 1):
				output_for_topcons = output_for_topcons + inline

	if (in_an_id == 1):
		if (output_for_scampi_seq == ''):
			output_for_scampi_seq = '-'
		if (output_for_scampi_msa == ''):
			output_for_scampi_msa = '-'
		if (output_for_prodiv == ''):
			output_for_prodiv = '-'
		if (output_for_pro == ''):
			output_for_pro = '-'
		if (output_for_octopus == ''):
			output_for_octopus = '-'
		if (output_for_topcons == ''):
			output_for_topcons = '-'
		output_line_for_scampi_seq = curr_id + ':::::' + output_for_scampi_seq + ':::::' + "\r\n"
		output_line_for_scampi_msa = curr_id + ':::::' + output_for_scampi_msa + ':::::' + "\r\n"
		output_line_for_prodiv = curr_id + ':::::' + output_for_prodiv + ':::::' + "\r\n"
		output_line_for_pro = curr_id + ':::::' + output_for_pro + ':::::' + "\r\n"
		output_line_for_octopus = curr_id + ':::::' + output_for_octopus + ':::::' + "\r\n"
		output_line_for_topcons = curr_id + ':::::' + output_for_topcons + ':::::' + "\r\n"
		outfile_scampi_seq.write( output_line_for_scampi_seq )
		outfile_scampi_msa.write( output_line_for_scampi_msa )
		outfile_prodiv.write( output_line_for_prodiv )
		outfile_pro.write( output_line_for_pro )
		outfile_octopus.write( output_line_for_octopus )
		outfile_topcons.write( output_line_for_topcons )

	outfile_scampi_seq.close()
	outfile_scampi_msa.close()
	outfile_prodiv.close()
	outfile_pro.close()
	outfile_octopus.close()
	outfile_topcons.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()

