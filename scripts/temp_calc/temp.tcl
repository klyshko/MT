set nf [molinfo top get numframes]
set sel [atomselect top "all"]
set outfile [open "T.dat" w]
$sel frame 0
set tubulins0 [$sel get {x y z}]
set tubulins0_avg [$sel get {x y z}]
set ts 200
set frequency 100000
set stride 10
set gam 860000
set k 0.002
for {set i 0} {$i < $nf} {incr i} { 
	puts "frame $i of $nf"
	set sum 0
	set sum_avg 0
	$sel frame $i
	set tubulins [$sel get {x y z}] 
 	for {set j 0} {$j < [$sel num]} {incr j} {
		set dr [vecsub [lindex $tubulins $j] [lindex $tubulins0 $j] ]
		set dr_avg [vecsub [lindex $tubulins $j] [lindex $tubulins0_avg $j] ]
		set dr2 [vecdot $dr $dr]
		set dr2_avg [vecdot $dr_avg $dr_avg]
		set part_sum [expr $sum + $dr2]
		set sum $part_sum
		set part_sum_avg [expr $sum_avg + $dr2_avg]
		set sum_avg $part_sum_avg
	}
	set T [expr $sum * $gam / ([$sel num] * 6 * $ts * $frequency * $k * $stride)]
	set T_avg [expr $sum_avg * $gam / ([$sel num] * 6 * $ts * $frequency * $k * $stride * ($i + 1))]
	puts $outfile "$i $T $T_avg"
	set tubulins0 $tubulins 
}
close $outfile
	
