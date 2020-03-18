#!/usr/bin/fish
while inotifywait -e modify ./**.tex ./figs/**
    make thesis
end
