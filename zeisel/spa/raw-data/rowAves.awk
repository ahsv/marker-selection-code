#!/bin/awk -f

{
	T = 0;
	if ($1 !~ /#/) {
            for (i=1; i<=NF; i++) T += $i;
            T/=NF;
            print T
        } else print $0
}
