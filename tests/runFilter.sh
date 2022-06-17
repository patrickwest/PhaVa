#!/bin/bash

for FILE in test_outputs/*ratio.tsv; do
filename=$(echo ${FILE});
echo ${filename};
while read line; do
rev_reads=$(echo "$line" | cut -f3);
rev_all_ratio=$(echo "$line" | cut -f4);
if [ ${rev_reads} -ge 4 ] && [ 1 -eq $(bc<<<"${rev_all_ratio} > 0.009") ]
then
echo "${line}" >> test_outputs/filtered_inverton_hits.tsv
fi
done < ${filename}
done
