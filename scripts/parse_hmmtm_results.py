#!/usr/bin/env python

# python parse_hmmtm_results.py alpha_polytopic_bitopic_noblanks.fasta.hmmtm_output alpha_polytopic_bitopic_noblanks.fasta alpha_polytopic_bitopic_noblanks.hmmtm_results

# python parse_hmmtm_results.py [input-file] [input-file-2] [output-file-hmmtm]

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

				#<table>
				#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;9</span></td></tr></table><table>
				#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>10&nbsp;&nbsp;&nbsp;27</span></td></tr></table><table>
				#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>28&nbsp;&nbsp;&nbsp;46</span></td></tr></table><table>

				is_TMhelix = 0
				is_inside = 0
				is_outside = 0

				if (this_row.find( '<font color="RED">tm' ) != -1):
					is_TMhelix = 1
				elif (this_row.find( '<font color="GREEN">in' ) != -1):
					is_inside = 1
				elif (this_row.find( '<font color="BLUE">out' ) != -1):
					is_outside = 1

				if ((is_TMhelix == 1) or (is_inside == 1) or (is_outside == 1)):

					this_char = '?'
					if (is_TMhelix == 1):
						this_char = 'M'
					elif (is_inside == 1):
						this_char = 'i'
					elif (is_outside == 1):
						this_char = 'o'

					bits1 = this_row.split( '</font>' )
					bit1 = bits1[1]
					bits2 = bit1.split( '</span>' )
					bit2 = bits2[0]
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bit2 = bit2.replace( '&nbsp;&nbsp;', '&nbsp;' )
					bits3 = bit2.split( '&nbsp;' )
					this_from = int(bits3[0])
					this_to = int(bits3[1])
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

#2b6o_A:::::
#<td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;7</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>8&nbsp;&nbsp;&nbsp;&nbsp;26</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>27&nbsp;&nbsp;&nbsp;38</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>39&nbsp;&nbsp;&nbsp;61</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>62&nbsp;&nbsp;&nbsp;80</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>81&nbsp;&nbsp;&nbsp;104</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>105&nbsp;&nbsp;124</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>125&nbsp;&nbsp;146</span></td></tr></table></td><td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>147&nbsp;&nbsp;156</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>157&nbsp;&nbsp;177</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>178&nbsp;&nbsp;202</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>203&nbsp;&nbsp;224</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>225&nbsp;&nbsp;263</span></td></tr></table></td>
#2zt9_G:::::
#<td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;4</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>5&nbsp;&nbsp;&nbsp;&nbsp;25</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>26&nbsp;&nbsp;&nbsp;37</span></td></tr></table></td>
#2cfp_A:::::
#<td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;9</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>10&nbsp;&nbsp;&nbsp;27</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>28&nbsp;&nbsp;&nbsp;46</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>47&nbsp;&nbsp;&nbsp;64</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>65&nbsp;&nbsp;&nbsp;74</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>75&nbsp;&nbsp;&nbsp;95</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>96&nbsp;&nbsp;&nbsp;102</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>103&nbsp;&nbsp;124</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>125&nbsp;&nbsp;144</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>145&nbsp;&nbsp;162</span></td></tr></table></td><td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>163&nbsp;&nbsp;167</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>168&nbsp;&nbsp;186</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>187&nbsp;&nbsp;221</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>222&nbsp;&nbsp;239</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>240&nbsp;&nbsp;259</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>260&nbsp;&nbsp;278</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>279&nbsp;&nbsp;290</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>291&nbsp;&nbsp;309</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>310&nbsp;&nbsp;314</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>315&nbsp;&nbsp;331</span></td></tr></table></td><td class=col5><span class="seqcla2">
#<table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>332&nbsp;&nbsp;346</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>347&nbsp;&nbsp;365</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>366&nbsp;&nbsp;380</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>381&nbsp;&nbsp;402</span></td></tr></table><table><tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>403&nbsp;&nbsp;417</span></td></tr></table></td>

#<table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;4</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>5&nbsp;&nbsp;&nbsp;&nbsp;25</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>26&nbsp;&nbsp;&nbsp;37</span></td></tr></table></td>

#<table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>1&nbsp;&nbsp;&nbsp;&nbsp;9</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>10&nbsp;&nbsp;&nbsp;27</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>28&nbsp;&nbsp;&nbsp;46</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>47&nbsp;&nbsp;&nbsp;64</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>65&nbsp;&nbsp;&nbsp;74</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>75&nbsp;&nbsp;&nbsp;95</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>96&nbsp;&nbsp;&nbsp;102</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>103&nbsp;&nbsp;124</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>125&nbsp;&nbsp;144</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>145&nbsp;&nbsp;162</span></td></tr></table></td><td class=col5><span class="seqcla2">
#<table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>163&nbsp;&nbsp;167</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>168&nbsp;&nbsp;186</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>187&nbsp;&nbsp;221</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>222&nbsp;&nbsp;239</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>240&nbsp;&nbsp;259</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>260&nbsp;&nbsp;278</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>279&nbsp;&nbsp;290</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>291&nbsp;&nbsp;309</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>310&nbsp;&nbsp;314</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>315&nbsp;&nbsp;331</span></td></tr></table></td><td class=col5><span class="seqcla2">
#<table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>332&nbsp;&nbsp;346</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>347&nbsp;&nbsp;365</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="BLUE">out&nbsp;&nbsp;</font>366&nbsp;&nbsp;380</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="RED">tm&nbsp;&nbsp;&nbsp;</font>381&nbsp;&nbsp;402</span></td></tr></table><table>
#<tr><td nowrap><span class="seqcla2"><font color="GREEN">in&nbsp;&nbsp;&nbsp;</font>403&nbsp;&nbsp;417</span></td></tr></table></td>

