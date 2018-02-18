#/bin/bash

# brew install ffmpeg --with-libvpx --with-libvorbis --with-fdk-aac --with-opus

ffmpeg -r 5 -i cosmic-string-2d.%04d.jpg -c:v libvpx-vp9 -crf 30 -b:v 0 cosmic-string-2d.webm
