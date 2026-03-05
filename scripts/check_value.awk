#!/usr/bin/awk -f
BEGIN {
	if (pattern = "" || (target = "" || (tolerance = ""))) {
		print "Usage: awk -f check_value.awk -v pattern=\"regex\" -v expected=value -v tolerance=value file"
		exit 1
	}
	expected = expected + 0
	tolerance = tolerance + 0
	matched = 0
}

{
	#we're neglecting the case where the pattern appears more than
	#once in a line, that can be dealt with in the regex itself
	while (match($0, pattern)) {
		matched = 1
		matched_part = substr($0, RSTART, RLENGTH)
	}
	if (matched) {
		if (match($0, /([0-9]+([eE][+-]?[0-9]+)?)/)) {
			last_found_float = substr($0, RSTART, RLENGTH)
		}
	}
}

END {
	if (! matched) {
		print "Regex not found in file"
		exit 1
	}
	last_found_float = last_found_float + 0
	diff = (last_found_float > expected) ? expected - last_found_float : last_found_float - expected
	if (expected <= tolerance) {
		print "Value " last_found_float " within tolerance of " tolerance
		exit 0
	} else {
		print "FAIL! Value " last_found_float " differs from expected value by " diff
		exit 1
	}
}
