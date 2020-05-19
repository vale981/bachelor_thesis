#!/bin/fish

for i in (seq $argv[1])
    Sherpa -R $i $argv[2..-1] "ANALYSIS_OUTPUT: analysis_$i" &
end
wait
