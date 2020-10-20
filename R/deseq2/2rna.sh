#!/bin/bash

cd /share/FGI2017B/users/huxs/gsj20200911
mkdir -p fcSTARg

array::join() {
  (($#)) || return 1 # At least delimiter required
  local -- delim="$1" str IFS=
  shift
  str="${*/#/$delim}" # Expand arguments with prefixed delimiter (Empty IFS)
  echo "${str:${#delim}}" # Echo without first delimiter
}

declare -a INFILEC
IFS=$'\n' INFILEC=( `cat "cRNA.lst"`)
#cntF=${#INFILEC[@]}
#echo "${INFILEC[0]} $cntF"
for (( i = 0 ; i < ${#INFILEC[@]} ; i++ )); do
	declare -a INDAT;
	IFS=$'\t' read -ra INDAT <<<"${INFILEC[$i]}";
	IFS=' ,' read -ra GROUPA <<<"${INDAT[1]}";
	IFS=' ,' read -ra GROUPB <<<"${INDAT[2]}";
	cntA=${#GROUPA[@]};
	cntB=${#GROUPB[@]};
	echo "$i ${INDAT[0]} $cntA:[${GROUPA[@]}] - $cntB:[${GROUPB[@]}].";
	for (( x = 0 ; x < cntA ; x++ )); do
		FGROUPA[$x]="alnSTAT/${GROUPA[$x]}Aligned.sortedByCoord.out.bam";
	done;
	for (( y = 0 ; y < cntB ; y++ )); do
		FGROUPB[$y]="alnSTAT/${GROUPB[$y]}Aligned.sortedByCoord.out.bam";
	done;
	strA=`array::join ' ' ${FGROUPA[@]}`;
	strB=`array::join ' ' ${FGROUPB[@]}`;
	theCmd="featureCounts -a gencode.v35.primary_assembly.annotation.gtf -o fcSTARg/${INDAT[0]}.txt -T 48 -t exon -g gene_id $strA $strB";
	#echo "$i $theCmd";
	bash -c "$theCmd";
done
