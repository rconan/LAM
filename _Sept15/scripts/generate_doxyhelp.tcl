#! /usr/bin/tclsh

set file_name "./doxy_latex"
set file_contents [read [open $file_name r] ];

set file_name_output "./doxy_output"
set file_out_pointer [open $file_name_output w ];


### <-- Begin Procedure::: Count how many needles are in a haystack
proc countSubstrings {haystack needle} {
    regexp -all ***=$needle $haystack
}
### <-- END Procedure::: Count how many needles are in a haystack




### <-- Begin Procedure::: Replace latex math
proc doxy_replace_math { file_contents } {

	set new_output ""
	set lines [split $file_contents \n]

	foreach line $lines { 

	regsub -all {[.*]\x7E} $line "" line  ;# remove all tildas!

	regsub -all {\\cite\{(\w*)\}} $line "" line
	regsub -all {\\ref\{(\w*)\}} $line "" line


		set line [string map {\\begin\{equation\} "\$"} $line]
		set line [string map {\\end\{equation\} "\$"} $line]

		set line [string map {\\begin\{eqnarray\} "\$"} $line]
		set line [string map {\\end\{eqnarray\} "\$"} $line]

		set line [string map {"\$" "\\f\$"} $line] ;# replace $ with /f$ 
		set line [string map {"\\f\\f\$" "\\f\$"} $line] ;# replace two f with one

	set line [regsub -all {\\textbf\{(.*?)\}} $line {\1}]

		append new_output $line "\n"
	}


#	puts $new_output;
	return $new_output;

}  
### <-- END Procedure::: Replace latex math




### <-- Begin Procedure::: Replace the empty stings for %> doxygen comments
proc doxy_comment_add { file_contents } {

	set new_output ""
	set lines [split $file_contents \n]
	
	foreach line $lines { 
	### read the file line-by-line

		if { [countSubstrings $line "%>"] == 0 } {
			append new_output "%> " $line "\n"

		} else {
			append new_output $line "\n"
		}
	}

#	puts $new_output;
	return $new_output;
}  ;# end of doxy_comment_add Procedure
### <-- END Procedure::: Replace the empty stings for %> doxygen comments



set output_stuff [ doxy_replace_math $file_contents ]
#puts $output_stuff

set output_stuff [doxy_comment_add $output_stuff ]
#puts $output_stuff


puts $file_out_pointer $output_stuff 

close $file_out_pointer