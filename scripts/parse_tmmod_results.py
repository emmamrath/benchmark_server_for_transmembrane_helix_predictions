#!/usr/bin/env python

# python parse_tmmod_results.py alpha_polytopic_bitopic_noblanks.1line.tmmod_link.tmmod_results alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.tmmod_results

# python parse_tmmod_results.py [input-file] [input-file-2] [output-file-tmmod]

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
	output_file = sys.argv[3]

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

	outfile = open( output_file, "w" )

	curr_id = ''
	curr_string_array = []
	curr_seq_length = 0
	helix_count = 0
	curr_string = ''
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
					curr_string = curr_string + curr_string_array[j]
				output_line = curr_id + ':::::' + curr_string + ':::::' + "\r\n"
				outfile.write( output_line )
				#print output_line

			# initialise for new id
			curr_id = inline[0:6]
			curr_string_array = []
			curr_seq_length = 0
			curr_string = ''
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

		else:
			rows = inline.split( '<tr>' )
			for this_row in rows:

				#<tr>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 106 <font></td>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 126 <font></td></tr>

				is_TMhelix = 0
				is_inside = 0
				is_outside = 0

				if (this_row.find( '&nbsp;&nbsp;&nbsp;&nbsp; TMhelix' ) != -1):
					is_TMhelix = 1
				elif (this_row.find( '&nbsp;&nbsp;&nbsp;&nbsp; inside' ) != -1):
					is_inside = 1
				elif (this_row.find( '&nbsp;&nbsp;&nbsp;&nbsp; outside' ) != -1):
					is_outside = 1

				if ((is_TMhelix == 1) or (is_inside == 1) or (is_outside == 1)):

					this_char = '?'
					if (is_TMhelix == 1):
						this_char = 'M'
					elif (is_inside == 1):
						this_char = 'i'
					elif (is_outside == 1):
						this_char = 'o'

					bits = this_row.split( '<font size= 2>' )
					got_from = bits[2]
					got_to = bits[3]
					got_from = got_from.replace( '&nbsp;', '' )
					got_from = got_from.replace( ' ', '' )
					got_to = got_to.replace( ' ', '' )
					from_bits = got_from.split( '<' )
					this_from = int(from_bits[0])
					to_bits = got_to.split( '<' )
					this_to = int(to_bits[0])
					for j in range( (this_from - 1), this_to ):
						curr_string_array[j] = this_char

		i = i + 1

	# finish last id
	if (curr_id != ''):
		for j in range( 0, curr_seq_length ):
			curr_string = curr_string + curr_string_array[j]
		output_line = curr_id + ':::::' + curr_string + ':::::' + "\r\n"
		outfile.write( output_line )
		#print output_line

	outfile.close()


if __name__ == "__main__":
	# Someone is launching this directly
	main()


# input_file :

#2xzb_A:::::
#0"><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Name <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 2xzb_A <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Annotation <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TM PROTEIN <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Length <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1033 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Number of predicted TMHs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 8 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number of AAs in TMHs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 157.123306 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number, first 60 AAs <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 0.000000 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Total prob of N-in <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 0.843466 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 105 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 106 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 126 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 127 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 137 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 138 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 158 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 159 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 302 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 303 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 327 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 328 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 328 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 329 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 353 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 354 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 798 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 799 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 819 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 820 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 924 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 925 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 945 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 946 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 963 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 964 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 984 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 985 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 992 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 993 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 1011 <font></td></tr><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1012 <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 1033 <font></td></tr><tr> <td ALIGN=LEFT VALIGN=CENTER > <font size= 2> <a href="frame.php?p=posterior&jobid=1633333921&index=0".> show posterior probabilities </a>  <font></td> </tr></table><font></td> <td ALIGN=RIGHT VALIGN=TOP > <font size= 2>&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  <font></td></tr></table><br><hr><br><br><tr><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td><td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#2b6o_A:::::



#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Name <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 2xzb_A <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Annotation <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TM PROTEIN <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Length <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1033 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Number of predicted TMHs <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 8 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number of AAs in TMHs <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 157.123306 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Exp number, first 60 AAs <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 0.000000 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> # Total prob of N-in <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 0.843466 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 1 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 105 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 106 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 126 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 127 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 137 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 138 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 158 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; inside <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 159 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 302 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 303 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 327 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; outside <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 328 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 328 <font></td></tr>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 329 <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 353 <font></td></tr>
#<tr> 
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> <a href="frame.php?p=posterior&jobid=1633333921&index=0".> show posterior probabilities </a>  <font></td> </tr></table><font></td> 
#	<td ALIGN=RIGHT VALIGN=TOP > <font size= 2>&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp;&nbsp;&nbsp;  <font></td></tr></table><br><hr><br><br>
#<tr>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp;  <font></td>
#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2>  <font></td></tr>



