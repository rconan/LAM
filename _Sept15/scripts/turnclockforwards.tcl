#! /usr/bin/tclsh

set years_back 6;

# date 082015062009
#root@dot:/home/think# date 082015062009
#Thu Aug 20 15:06:00 EST 2009

catch "exec date --date=today +%Y%%%m%%%d%%%H%%%M" result;  set temp_result_list [split $result "%"]
set kmv_thisyear   [lindex $temp_result_list 0 ]
set kmv_thismonth  [lindex $temp_result_list 1 ]
set kmv_thisday    [lindex $temp_result_list 2 ]
set kmv_thishour   [lindex $temp_result_list 3 ]
set kmv_thisminute [lindex $temp_result_list 4 ]

puts "Today is $kmv_thishour $kmv_thisminute at $kmv_thisyear $kmv_thismonth $kmv_thisday\n";

set kmv_pastyear [expr $kmv_thisyear+$years_back];

puts "Turning the clock back to::: $kmv_pastyear $kmv_thismonth $kmv_thisday .... \n";

[exec date $kmv_thismonth$kmv_thisday$kmv_thishour$kmv_thisminute$kmv_pastyear];
