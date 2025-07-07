#!/bin/bash
#SBATCH --account=schmidt-np
#SBATCH --partition=schmidt-np
export TERM=xterm-256color  # or another supported terminal type like 'vt100i'
export PATH="/uufs/chpc.utah.edu/common/home/schmidt-group3/software/miniconda/bin":$PATH

while read id; do
esearch -db assembly -query $id < /dev/null \
  | esummary \
  | xtract -pattern DocumentSummary -element FtpPath_RefSeq,FtpPath_GenBank \
  | awk '{print $1}' | while read -r url; do
     
     
        base=$(basename "$url")
        #wget "${url}/${base}_genomic.fna.gz"
	wget "${url}/${base}_genomic.gbff.gz"
      
      i
    done
done < $1


