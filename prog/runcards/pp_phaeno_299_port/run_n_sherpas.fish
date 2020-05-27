#!/bin/fish

for i in (seq $argv[1])
    /home/hiro/src/sherpa_rel_2_2_9/build/install/bin/Sherpa -R $i $argv[2..-1] "ANALYSIS_OUTPUT:=analysis_$i" &
end
wait
