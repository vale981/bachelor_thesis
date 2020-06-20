#!/usr/bin/fish
while inotifywait -e modify ./**.tex ./figs/**
    make slides
end
