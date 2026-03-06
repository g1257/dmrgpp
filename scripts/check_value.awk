#!/usr/bin/awk -f
BEGIN {
	if (pattern == "" || (expected == "" || (tolerance == ""))) {
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
	temp_line = $0
	while (match(temp_line, pattern, matches)) {
		matched = 1
		print matches[0]
		temp_line = substr(temp_line, RSTART + RLENGTH)
		last_found_float = matches[1]
	}
}

END {
	if (! matched) {
		print "Regex not found in file"
		exit 1
	}
	last_found_float = last_found_float + 0
	diff = (last_found_float < expected) ? expected - last_found_float : last_found_float - expected
	if (diff <= tolerance) {
		print "Value " last_found_float " within tolerance of " tolerance " diff: " diff
		exit 0
	} else {
		print "FAIL! Value " last_found_float " differs from expected value by " diff
		exit 1
	}
}
