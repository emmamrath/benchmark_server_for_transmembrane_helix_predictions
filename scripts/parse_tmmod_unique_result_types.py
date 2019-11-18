#!/usr/bin/env python

# python parse_tmmod_unique_result_types.py alpha_polytopic_bitopic_noblanks.1line.tmmod_link.tmmod_results

# python parse_tmmod_unique_result_types.py [input-file]

import sys
import re


# global variables
save_type = []


######################################################
def main():

	global outfile
	input_file_1 = sys.argv[1]

	# process input file 1 that contains the results to parse

	infile1 = open( input_file_1, "r" )
	inlines1 = infile1.readlines()
	infile1.close()

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
		if (is_id_line != 1):

			rows = inline.split( '<tr>' )
			for this_row in rows:

				#<tr>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; TMhelix <font></td>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> &nbsp;&nbsp;&nbsp;&nbsp; 106 <font></td>
				#	<td ALIGN=LEFT VALIGN=CENTER > <font size= 2> 126 <font></td></tr>

				cols = this_row.split( '<td' )
				if (len(cols) > 1):
					this_col = cols[1]
					if (this_col.find( '<a href="frame.php?' ) == -1):
						this_col = this_col.replace( 'ALIGN=LEFT VALIGN=CENTER >', '' )
						this_col = this_col.replace( '<font size= 2>', '' )
						this_col = this_col.replace( '&nbsp;&nbsp;&nbsp;&nbsp;', '' )
						this_col = this_col.replace( '<font>', '' )
						this_col = this_col.replace( '</td>', '' )
						this_col = this_col.strip()

						found_type = 0
						for this_save_type in save_type:
							if (this_col == this_save_type):
								found_type = 1
						if (found_type == 0):
							print this_col
							save_type.append( this_col )


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



