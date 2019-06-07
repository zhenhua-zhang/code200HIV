#!/bin/bash

zcat $1 \
	| awk 'BEGIN {
	FS = " "; OFS = "\t"; getline; printf $1;
	for(i=1; i<=NF; i++) {
		where=match($i, /_X/);
		if(where!=0) {
			sbjnm = substr($i, where+1, length($i));
			printf "\t"sbjnm;
			arr[length(arr)+1]=i;
		}
	}
	print "";
}
{
	printf $1;
	len=length(arr);
	for(j=1; j<=len; j++) {printf "\t"$arr[j];}
	print "";
}
END {}
' \
	| gzip \
	> $2
