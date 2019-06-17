
BEGIN {

	maxeff=0;
	maxcorr=0;
	minmark=100000;

}

{
	# efficiency automatically gets rid of the cases that explode
	if ($3/$4 > maxeff) {
		besteff = $0;
		maxeff = $3/$4
	}

	if ($3 > maxcorr) {
		bestcorr = $0;
		maxcorr = $3
	}

	if ($4 < minmark) {
		bestmark = $0;
		minmark = $4
	}
}

END {
	print besteff;
	print bestcorr;
	print bestmark;
}

