#!/usr/bin/env python

# python parse_memsatsvm_results.py alpha_polytopic_bitopic_noblanks.fasta.memsatsvm1.memsatsvm_link.memsatsvm_results alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.memsatsvm1.memsatsvm alpha_polytopic_bitopic_noblanks.memsatsvm1.memsat3

# python parse_memsatsvm_results.py [input-file] [input-file-2] [output-file-memsatsvm] [output-file-memsat3]

import sys
import re


# global variables
save_id = []
save_seq_length = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]
	input_file_2 = sys.argv[2]
	output_file_1 = sys.argv[3]
	output_file_2 = sys.argv[4]

	# read in the lengths of each sequence, contained in input file 2

	infile2 = open( input_file_2, "r" )
	inlines2 = infile2.readlines()
	infile2.close()

	in_a_seq = 0
	curr_id = ''
	curr_seq = ''
	for inline in inlines2:
		inline = inline.strip()
		if (inline != ''):
			if (inline[0:1] == '>'):
				if (in_a_seq == 1):
					save_id.append( curr_id )
					save_seq_length.append( len(curr_seq) )
				in_a_seq = 1
				curr_id = inline[1:7]
				curr_seq = ''
			else:
				curr_seq = curr_seq + inline
	if (in_a_seq == 1):
		save_id.append( curr_id )
		save_seq_length.append( len(curr_seq) )

	# process input file 1 that contains the results to parse

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

	outfile1 = open( output_file_1, "w" )
	outfile2 = open( output_file_2, "w" )

	curr_id = ''
	curr_string_array = []
	curr_seq_length = 0
	helix_count = 0
	curr_memsatsvm_string = ''
	curr_memsat3_string = ''
	i = 0
	for inline in inlines1:

		is_id_line = 0
		if ( len(inline) >= 11):
			if (inline[6:11] == ':::::'):
				is_id_line = 1

		# 2xzb_A:::::
		if (is_id_line == 1):

			# finish previous id
			if (curr_id != ''):
				for j in range( 0, curr_seq_length ):
					curr_memsatsvm_string = curr_memsatsvm_string + curr_string_array[j]
				output_line_1 = curr_id + ':::::' + curr_memsatsvm_string + ':::::' + "\r\n"
				outfile1.write( output_line_1 )
				output_line_2 = curr_id + ':::::' + curr_memsat3_string + ':::::' + "\r\n"
				outfile2.write( output_line_2 )

			# initialise for new id
			curr_id = inline[0:6]
			curr_string_array = []
			curr_seq_length = 0
			helix_count = 0
			curr_memsatsvm_string = ''
			curr_memsat3_string = ''
			found_id = 0
			j = 0
			for this_save_id in save_id:
				if (this_save_id == curr_id):
					found_id = 1
					curr_seq_length = save_seq_length[j]
				j = j + 1
			if (found_id == 1):
				for j in range( 0, curr_seq_length ):
					curr_string_array.append( '_' )
			else:
					print curr_id, 'did not find its sequence length'

		# <td style="width: 20%">Signal peptide</td>
		# <td style="width: 20%"><b>Not detected.</b></td>
		# <td style="width: 20%"><b>1-12</b></td>
		elif ( inline.find('<td style="width: 20%">Signal peptide</td>') != -1 ):
			next_line = inlines1[ i + 1 ]
			bits1 = next_line.split('<b>')
			bit1 = bits1[1]
			bits2 = bit1.split('</b>')
			this_text = bits2[0]

		# <td style="width: 20%">Signal score</td>
		# <td style="width: 20%"><b>0.189</b></td>

		# <td style="width: 20%">Topology</td>
		# <td style="width: 20%"><b>11-32,42-59,87-107,128-149,161-176,203-222</b></td>
		elif ( inline.find('<td style="width: 20%">Topology</td>') != -1 ):
			next_line = inlines1[ i + 1 ]
			bits1 = next_line.split('<b>')
			bit1 = bits1[1]
			bits2 = bit1.split('</b>')
			this_text = bits2[0]
			if (this_text != 'Not detected.'):
				bits3 = this_text.rsplit(',')
				for this_helix in bits3:
					helix_count = helix_count + 1
					bits4 = this_helix.split('-')
					this_from = bits4[0]
					this_from = int(this_from)
					this_to = bits4[1]
					this_to = int(this_to)
					for j in range( (this_from - 1), this_to ):
						if (j < len(curr_string_array)):
							curr_string_array[j] = 'H'

		# <td style="width: 20%">Re-entrant helices</td>
		# <td style="width: 20%"><b>63-81,180-194</b></td>
		elif ( inline.find('<td style="width: 20%">Re-entrant helices</td>') != -1 ):
			next_line = inlines1[ i + 1 ]
			bits1 = next_line.split('<b>')
			bit1 = bits1[1]
			bits2 = bit1.split('</b>')
			this_text = bits2[0]
			if (this_text != 'Not detected.'):
				bits3 = this_text.rsplit(',')
				for this_helix in bits3:
					bits4 = this_helix.split('-')
					this_from = bits4[0]
					this_from = int(this_from)
					this_to = bits4[1]
					this_to = int(this_to)
					for j in range( (this_from - 1), this_to ):
						curr_string_array[j] = 'R'

		# <td style="width: 20%">Helix count</td>
		# <td style="width: 20%"><b>6</b></td>
		elif ( inline.find('<td style="width: 20%">Helix count</td>') != -1 ):
			next_line = inlines1[ i + 1 ]
			bits1 = next_line.split('<b>')
			bit1 = bits1[1]
			bits2 = bit1.split('</b>')
			this_helix_count = bits2[0]
			this_helix_count = int(this_helix_count)
			if (this_helix_count != helix_count):
				print curr_id, 'the added helix count', helix_count, 'is not the same as the read helix helix count', this_helix_count

		# <td style="width: 20%">N-terminal</td>
		# <td style="width: 20%"><b>in</b></td>
		elif ( inline.find('<td style="width: 20%">N-terminal</td>') != -1 ):
			next_line = inlines1[ i + 1 ]
			bits1 = next_line.split('<b>')
			bit1 = bits1[1]
			bits2 = bit1.split('</b>')
			this_n_terminal = bits2[0]
			curr_char = '?'
			if (this_n_terminal == 'in'):
				curr_char = 'i'
			elif (this_n_terminal == 'out'):
				curr_char = 'o'
			else:
				print curr_id, 'unexpected N-terminal :', this_n_terminal
			for j in range( 0, curr_seq_length ):
				if (curr_string_array[j] == '_'):
					curr_string_array[j] = curr_char
				elif (curr_string_array[j] == 'H'):
					if (j > 0):
						if (curr_string_array[j - 1] == 'i'):
							curr_char = 'o'
						elif (curr_string_array[j - 1] == 'o'):
							curr_char = 'i'
						elif (curr_string_array[j - 1] == 'H'):
							do_nothing = 1
						else:
							print curr_id, 'unexpected residue position :', curr_string_array[j - 1]

		# <td class="strand">S</td>
		# <td class="inner_loop">+</td>
		# <td class="inner_cap">I</td>
		# <td class="helix">X</td>
		# <td class="outer_cap">O</td>
		# <td class="outer_loop">-</td>
		elif ( inline.find('<td class="strand">S</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + 'S'
		elif ( inline.find('<td class="inner_loop">+</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + '+'
		elif ( inline.find('<td class="inner_cap">I</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + 'I'
		elif ( inline.find('<td class="helix">X</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + 'X'
		elif ( inline.find('<td class="outer_cap">O</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + 'O'
		elif ( inline.find('<td class="outer_loop">-</td>') == 0 ):
			curr_memsat3_string = curr_memsat3_string + '-'

		i = i + 1

	# finish last id
	if (curr_id != ''):
		for j in range( 0, curr_seq_length ):
			curr_memsatsvm_string = curr_memsatsvm_string + curr_string_array[j]
		output_line_1 = curr_id + ':::::' + curr_memsatsvm_string + ':::::' + "\r\n"
		outfile1.write( output_line_1 )
		output_line_2 = curr_id + ':::::' + curr_memsat3_string + ':::::' + "\r\n"
		outfile2.write( output_line_2 )

	outfile1.close()
	outfile2.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()


# input_file :

# <tr><td><h1>MEMSAT-SVM Prediction</h1></td></tr>
# <tr class="form_subtitle"><td><h4>Summary of MEMSAT-SVM Topology Analysis</h4></td></tr>
# <tr><td>
# <table>
# <tr>
# <td style="width: 20%">Signal peptide</td>
# <td style="width: 20%"><b>Not detected.</b></td>
# </tr>
# <tr>
# <td style="width: 20%">Signal score</td>
# <td style="width: 20%"><b>0.189</b></td>
# </tr>
# <tr>
# <td style="width: 20%">Topology</td>
# <td style="width: 20%"><b>11-32,42-59,87-107,128-149,161-176,203-222</b></td>
# </tr>
# <tr>
# <td style="width: 20%">Re-entrant helices</td>
# <td style="width: 20%"><b>63-81,180-194</b></td>
# </tr>
# <tr>
# <td style="width: 20%">Helix count</td>
# <td style="width: 20%"><b>6</b></td>
# </tr>
# <tr>
# <td style="width: 20%">N-terminal</td>
# <td style="width: 20%"><b>in</b></td>
# </tr>
# <tr>
# <td style="width: 20%">Score</td>
# <td style="width: 20%"><b>9.74787</b></td>
# </tr>
# </table>
# </td></tr>
# <tr><td><h1>MEMSAT3 Prediction</h1></td></tr>
# <tr class="form_subtitle"><td><h4>Summary of MEMSAT3 Topology Analysis</h4></td></tr>
# <tr><td>
# <table style="width: 40%">
# <tr><td><h4>Number</h4></td><td><h4>Type</h4></td><td><h4>Direction</h4></td><td><h4>Score</h4></td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>1</td><td>helix</td><td>+</td><td>38.043</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>1</td><td>helix</td><td>-</td><td>50.201</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>2</td><td>helices</td><td>+</td><td>88.244</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>2</td><td>helices</td><td>-</td><td>58.849</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>3</td><td>helices</td><td>+</td><td>102.912</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>3</td><td>helices</td><td>-</td><td>109.05</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>4</td><td>helices</td><td>+</td><td>153.113</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>4</td><td>helices</td><td>-</td><td>126.517</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>5</td><td>helices</td><td>+</td><td>170.58</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>5</td><td>helices</td><td>-</td><td>177.004</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>6</td><td>helices</td><td>+</td><td>221.067</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>8</td><td>helices</td><td>+</td><td>195.707</td></tr>
# </table>
# </td></tr>
# <tr class="form_subtitle"><td><h4>Final MEMSAT3 Prediction</h4></td></tr>
# <tr><td>
# <table style="width: 40%">
# <tr><td><h4>Segment</h4></td><td><h4>Range</h4></td><td><h4>Score</h4></td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>1</td><td>(in) 11-35</td>
# <td>34.05</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>2</td><td>42-66</td>
# <td>33.83</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>3</td><td>84-108</td>
# <td>26.63</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>4</td><td>126-149</td>
# <td>31.39</td></tr>
# <tr style="padding: 3px; background: #ccc; text-align: center;"><td>5</td><td>159-178</td>
# <td>20.07</td></tr>
# <tr style="padding: 3px; background: #fff; text-align: center;"><td>6</td><td>202-221</td>
# <td>25.89</td></tr>
# </table>
# </td></tr>
# <tr class="form_subtitle"><td><h4>Diagram</h4></td></tr>
# <tr><td>
# <h3>Key:</h3>
# <table>
# <tr><td class="helix">X</td><td>Central transmembrane helix segment</td><td class="strand">S</td><td>Possible N-terminal signal peptide</td></tr>
# <tr><td class="inner_loop">+</td><td>Inside loop</td><td class="inner_cap">I</td><td>Inside helix cap</td></tr>
# <tr><td class="outer_loop">-</td><td>Outside loop</td><td class="outer_cap">O</td><td>Outside helix cap</td></tr>
# </table>
# </td></tr>
# <tr><td>
# <table>
# <tr>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">3</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">4</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">5</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">6</td>
# <td class="alignment">0</td>
# </tr>
# <tr>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# </tr>
# <tr>
# <td class="alignment">M</td>
# <td class="alignment">W</td>
# <td class="alignment">E</td>
# <td class="alignment">L</td>
# <td class="alignment">R</td>
# <td class="alignment">S</td>
# <td class="alignment">A</td>
# <td class="alignment">S</td>
# <td class="alignment">F</td>
# <td class="alignment">W</td>
# <td class="alignment">R</td>
# <td class="alignment">A</td>
# <td class="alignment">I</td>
# <td class="alignment">F</td>
# <td class="alignment">A</td>
# <td class="alignment">E</td>
# <td class="alignment">F</td>
# <td class="alignment">F</td>
# <td class="alignment">A</td>
# <td class="alignment">T</td>
# <td class="alignment">L</td>
# <td class="alignment">F</td>
# <td class="alignment">Y</td>
# <td class="alignment">V</td>
# <td class="alignment">F</td>
# <td class="alignment">F</td>
# <td class="alignment">G</td>
# <td class="alignment">L</td>
# <td class="alignment">G</td>
# <td class="alignment">A</td>
# <td class="alignment">S</td>
# <td class="alignment">L</td>
# <td class="alignment">R</td>
# <td class="alignment">W</td>
# <td class="alignment">A</td>
# <td class="alignment">P</td>
# <td class="alignment">G</td>
# <td class="alignment">P</td>
# <td class="alignment">L</td>
# <td class="alignment">H</td>
# <td class="alignment">V</td>
# <td class="alignment">L</td>
# <td class="alignment">Q</td>
# <td class="alignment">V</td>
# <td class="alignment">A</td>
# <td class="alignment">L</td>
# <td class="alignment">A</td>
# <td class="alignment">F</td>
# <td class="alignment">G</td>
# <td class="alignment">L</td>
# <td class="alignment">A</td>
# <td class="alignment">L</td>
# <td class="alignment">A</td>
# <td class="alignment">T</td>
# <td class="alignment">L</td>
# <td class="alignment">V</td>
# <td class="alignment">Q</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">G</td>
# </tr>
# <tr><td colspan="60"><hr class="row_separator" /></td></tr>
# <tr>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">7</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">8</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">9</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">0</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">1</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">2</td>
# <td class="alignment">0</td>
# </tr>
# <tr>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# </tr>
# <tr>
# <td class="alignment">H</td>
# <td class="alignment">I</td>
# <td class="alignment">S</td>
# <td class="alignment">G</td>
# <td class="alignment">A</td>
# <td class="alignment">H</td>
# <td class="alignment">V</td>
# <td class="alignment">N</td>
# <td class="alignment">P</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">T</td>
# <td class="alignment">F</td>
# <td class="alignment">A</td>
# <td class="alignment">F</td>
# <td class="alignment">L</td>
# <td class="alignment">V</td>
# <td class="alignment">G</td>
# <td class="alignment">S</td>
# <td class="alignment">Q</td>
# <td class="alignment">M</td>
# <td class="alignment">S</td>
# <td class="alignment">L</td>
# <td class="alignment">L</td>
# <td class="alignment">R</td>
# <td class="alignment">A</td>
# <td class="alignment">I</td>
# <td class="alignment">C</td>
# <td class="alignment">Y</td>
# <td class="alignment">V</td>
# <td class="alignment">V</td>
# <td class="alignment">A</td>
# <td class="alignment">Q</td>
# <td class="alignment">L</td>
# <td class="alignment">L</td>
# <td class="alignment">G</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">A</td>
# <td class="alignment">G</td>
# <td class="alignment">A</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">L</td>
# <td class="alignment">Y</td>
# <td class="alignment">S</td>
# <td class="alignment">V</td>
# <td class="alignment">T</td>
# <td class="alignment">P</td>
# <td class="alignment">P</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">R</td>
# <td class="alignment">G</td>
# <td class="alignment">N</td>
# <td class="alignment">L</td>
# <td class="alignment">A</td>
# <td class="alignment">L</td>
# <td class="alignment">N</td>
# <td class="alignment">T</td>
# </tr>
# <tr><td colspan="60"><hr class="row_separator" /></td></tr>
# <tr>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">3</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">4</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">5</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">6</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">7</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">8</td>
# <td class="alignment">0</td>
# </tr>
# <tr>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# </tr>
# <tr>
# <td class="alignment">L</td>
# <td class="alignment">H</td>
# <td class="alignment">P</td>
# <td class="alignment">G</td>
# <td class="alignment">V</td>
# <td class="alignment">S</td>
# <td class="alignment">V</td>
# <td class="alignment">G</td>
# <td class="alignment">Q</td>
# <td class="alignment">A</td>
# <td class="alignment">T</td>
# <td class="alignment">I</td>
# <td class="alignment">V</td>
# <td class="alignment">E</td>
# <td class="alignment">I</td>
# <td class="alignment">F</td>
# <td class="alignment">L</td>
# <td class="alignment">T</td>
# <td class="alignment">L</td>
# <td class="alignment">Q</td>
# <td class="alignment">F</td>
# <td class="alignment">V</td>
# <td class="alignment">L</td>
# <td class="alignment">C</td>
# <td class="alignment">I</td>
# <td class="alignment">F</td>
# <td class="alignment">A</td>
# <td class="alignment">T</td>
# <td class="alignment">Y</td>
# <td class="alignment">D</td>
# <td class="alignment">E</td>
# <td class="alignment">R</td>
# <td class="alignment">R</td>
# <td class="alignment">N</td>
# <td class="alignment">G</td>
# <td class="alignment">R</td>
# <td class="alignment">L</td>
# <td class="alignment">G</td>
# <td class="alignment">S</td>
# <td class="alignment">V</td>
# <td class="alignment">A</td>
# <td class="alignment">L</td>
# <td class="alignment">A</td>
# <td class="alignment">V</td>
# <td class="alignment">G</td>
# <td class="alignment">F</td>
# <td class="alignment">S</td>
# <td class="alignment">L</td>
# <td class="alignment">T</td>
# <td class="alignment">L</td>
# <td class="alignment">G</td>
# <td class="alignment">H</td>
# <td class="alignment">L</td>
# <td class="alignment">F</td>
# <td class="alignment">G</td>
# <td class="alignment">M</td>
# <td class="alignment">Y</td>
# <td class="alignment">Y</td>
# <td class="alignment">T</td>
# <td class="alignment">G</td>
# </tr>
# <tr><td colspan="60"><hr class="row_separator" /></td></tr>
# <tr>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">1</td>
# <td class="alignment">9</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">0</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">1</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">2</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">3</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">4</td>
# <td class="alignment">0</td>
# </tr>
# <tr>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_loop">-</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="outer_cap">O</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="helix">X</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_cap">I</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# </tr>
# <tr>
# <td class="alignment">A</td>
# <td class="alignment">G</td>
# <td class="alignment">M</td>
# <td class="alignment">N</td>
# <td class="alignment">P</td>
# <td class="alignment">A</td>
# <td class="alignment">R</td>
# <td class="alignment">S</td>
# <td class="alignment">F</td>
# <td class="alignment">A</td>
# <td class="alignment">P</td>
# <td class="alignment">A</td>
# <td class="alignment">I</td>
# <td class="alignment">L</td>
# <td class="alignment">T</td>
# <td class="alignment">R</td>
# <td class="alignment">N</td>
# <td class="alignment">F</td>
# <td class="alignment">T</td>
# <td class="alignment">N</td>
# <td class="alignment">H</td>
# <td class="alignment">W</td>
# <td class="alignment">V</td>
# <td class="alignment">Y</td>
# <td class="alignment">W</td>
# <td class="alignment">V</td>
# <td class="alignment">G</td>
# <td class="alignment">P</td>
# <td class="alignment">V</td>
# <td class="alignment">I</td>
# <td class="alignment">G</td>
# <td class="alignment">A</td>
# <td class="alignment">G</td>
# <td class="alignment">L</td>
# <td class="alignment">G</td>
# <td class="alignment">S</td>
# <td class="alignment">L</td>
# <td class="alignment">L</td>
# <td class="alignment">Y</td>
# <td class="alignment">D</td>
# <td class="alignment">F</td>
# <td class="alignment">L</td>
# <td class="alignment">L</td>
# <td class="alignment">F</td>
# <td class="alignment">P</td>
# <td class="alignment">R</td>
# <td class="alignment">L</td>
# <td class="alignment">K</td>
# <td class="alignment">S</td>
# <td class="alignment">V</td>
# <td class="alignment">S</td>
# <td class="alignment">E</td>
# <td class="alignment">R</td>
# <td class="alignment">L</td>
# <td class="alignment">S</td>
# <td class="alignment">I</td>
# <td class="alignment">L</td>
# <td class="alignment">K</td>
# <td class="alignment">G</td>
# <td class="alignment">T</td>
# </tr>
# <tr><td colspan="60"><hr class="row_separator" /></td></tr>
# <tr>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">5</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment"> </td>
# <td class="alignment">2</td>
# <td class="alignment">6</td>
# <td class="alignment">0</td>
# <td class="alignment"> </td>
# </tr>
# <tr>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# <td class="inner_loop">+</td>
# </tr>
# <tr>
# <td class="alignment">R</td>
# <td class="alignment">P</td>
# <td class="alignment">S</td>
# <td class="alignment">E</td>
# <td class="alignment">S</td>
# <td class="alignment">N</td>
# <td class="alignment">G</td>
# <td class="alignment">Q</td>
# <td class="alignment">P</td>
# <td class="alignment">E</td>
# <td class="alignment">V</td>
# <td class="alignment">T</td>
# <td class="alignment">G</td>
# <td class="alignment">E</td>
# <td class="alignment">P</td>
# <td class="alignment">V</td>
# <td class="alignment">E</td>
# <td class="alignment">L</td>
# <td class="alignment">K</td>
# <td class="alignment">T</td>
# <td class="alignment">Q</td>
# <td class="alignment">A</td>
# <td class="alignment">L</td>
# </tr>
# <tr><td colspan="21"><hr class="row_separator" /></td></tr>
# </table>
# </td></tr>
# <tr><td>
# <h3>Key:</h3>
# <table>
# <tr><td class="helix">X</td><td>Central transmembrane helix segment</td><td class="strand">S</td><td>Possible N-terminal signal peptide</td></tr>
# <tr><td class="inner_loop">+</td><td>Inside loop</td><td class="inner_cap">I</td><td>Inside helix cap</td></tr>
# <tr><td class="outer_loop">-</td><td>Outside loop</td><td class="outer_cap">O</td><td>Outside helix cap</td></tr>
# </table>
# </td></tr>

