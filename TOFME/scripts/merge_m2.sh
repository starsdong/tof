#!/bin/bash
# Merge the completed StTofMEMaker v2 (m2) jobs: logs ending normally with the maker summary.
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
OK=(); BAD=()
for log in logs/me_m2/me_*.log; do
  id=$(basename "$log" .log); id=${id#me_}
  if grep -q "This is the end" "$log" && grep -q "Finish events" "$log"; then OK+=("output/me_m2/me_${id}.root"); else BAD+=("$id"); fi
done
echo "completed jobs: ${#OK[@]}   excluded: ${BAD[*]}"
printf "%s\n" "${OK[@]}" > output/me_m2_merged_inputs.txt
rm -f output/me_m2.root
bash scripts/sl7.sh "hadd -f output/me_m2.root `cat output/me_m2_merged_inputs.txt | tr '\n' ' '` >&! logs/hadd_m2.log"
tail -1 logs/hadd_m2.log
for f in "${OK[@]}"; do id=$(basename $f .root); grep -h "Finish events" logs/me_m2/${id}.log; done \
  | awk '{a+=$5; s+=$7} END {printf "events read %d selected %d\n", a, s}'
echo END
