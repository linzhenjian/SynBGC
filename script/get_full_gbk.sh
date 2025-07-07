#!/bin/bash
#while read i; do id=$(echo $i | sed 's/_/ /3' | awk '{print $2}'); echo "$i $(grep PF00109 V1_output/collinear/processed/$id.txt.anno)"; done < list > 1
while read i; do
	strain=$(echo $i  | awk -F '_' '{print $1"_"$2}' )
	contig=$(echo $i | awk -F '_' '{print $3}')
	id=$(echo $i | sed 's/_/ /3' | awk '{print $NF}')
	awk -F _ '{print $NF}' ./collinear/processed/$id.txt > temp_id
	if [ ! -s ${strain}_${contig}.gbk ]; then
		python ~/../schmidt-group3/software/syntenic_BGC_v2/script/extract_contig_gbk.py  ${strain}_anno.GBK ${strain}_${contig}.gbk $contig
		python ~/../schmidt-group3/software/syntenic_BGC_v2/script/extract_records_gbk.py   ${strain}_${contig}.gbk  ${i}.gbk temp_id
                python ~/../schmidt-group3/software/syntenic_BGC_v2/script/cleanup_gbk.py ${i}.gbk temp2 --prefix "${strain}_"
                mv temp2 ${i}.gbk

	else
		python ~/../schmidt-group3/software/syntenic_BGC_v2/script/extract_records_gbk.py   ${strain}_${contig}.gbk  ${i}.gbk temp_id
		python ~/../schmidt-group3/software/syntenic_BGC_v2/script/cleanup_gbk.py ${i}.gbk temp2 --prefix "${strain}_"
		mv temp2 ${i}.gbk
	fi
		
done < hit
