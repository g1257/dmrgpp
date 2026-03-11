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
		printf "Value %e within tolerance of %e of expected value %e diff: %e\n", last_found_float, tolerance, expected, diff
		exit 0
	} else {
		printf "FAIL! Value %e differs from expected value %e by %e\n", last_found_float, expected, diff
		exit 1
	}
}
