#!/bin/sh

name=$(echo "$1" | sed -E 's/^(([^_]*_){2}[^_]*)_/\1 /' |  awk -F ' ' '{print $2".split.anno"}'); echo $name
awk -F '[_\t]' '{print $4,$0}' OFS='\t' collinear/processed/$name | sed 's/^g//' | sort -k1,1n  | awk -F '\t' '$5!="-" {print $0}' OFS='\t' | sort -k1,1 -k4,4g -k5,5r  | awk '!seen[$1]++'  
#awk '!seen[$1]++' $name | awk -F '[_\t]' '{print $4,$0}' | sed 's/^g//' | sort -k1,1n  | awk -F '\t' '{print $1,$3,$4,$5,$6}' OFS='\t'
