#!/bin/bash
# Merge the completed TOFME mixed-event jobs (logs ending normally) and solve the geometry from the merged sums.
cd /gpfs01/star/pwg/dongx/tof/TOFME || exit 1
OK=(); BAD=()
for log in logs/me_1M/me_*.log; do
  id=$(basename "$log" .log); id=${id#me_}
  if grep -q "This is the end" "$log" && grep -q "Finish events" "$log"; then OK+=("output/me_1M/me_${id}.root"); else BAD+=("$id"); fi
done
echo "completed jobs: ${#OK[@]}   excluded: ${BAD[*]}"
printf "%s\n" "${OK[@]}" > output/me_1M_merged_inputs.txt
rm -f output/me_1M.root
bash scripts/sl7.sh "hadd -f output/me_1M.root `cat output/me_1M_merged_inputs.txt | tr '\n' ' '` >&! logs/hadd_1M.log"
tail -2 logs/hadd_1M.log
grep -h "Finish events" $(printf "logs/me_1M/me_%s.log " $(sed 's|output/me_1M/me_||; s|.root||' output/me_1M_merged_inputs.txt)) \
  | awk '{a+=$5; s+=$7; p+=$10} END {printf "events read %d selected %d projected tracks %d\n", a, s, p}'
bash scripts/sl7.sh "root -b -q -l 'macros/makeGeom.C+(\"output/me_1M.root\",\"geom/btof_geom_1M.txt\")' >&! logs/makeGeom_1M.log"
grep -v -E "GARFIELD|^$" logs/makeGeom_1M.log | tail -6
echo END
