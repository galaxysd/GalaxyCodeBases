#!/bin/bash

echo Input "[${1}].png"

convert "$1.png" -crop 3146x5674+0+0 -background navy -gravity south -extent 3146x6030 "$1.1.png"
convert "$1.png" -crop 3146x3180+0+5674 "$1.2.png"
convert "$1.png" -crop 3146x1590+0+8854 "$1.3.png"

convert "$1.1.png" -alpha set -virtual-pixel background -background navy \
-set option:distort:scale 4 -distort BilinearForward \
'0,0,1573,0  3146,0,1574,0  0,6030,681,2010  3146,6030,2467,2010' \
-scale 25% \
-crop 2400x2010+374+0 \
-background white -gravity south -extent 2400x2274 \
"$1.t1c.png"

convert "$1.2.png" -alpha set -virtual-pixel background -background navy \
-set option:distort:scale 4 -distort BilinearForward \
'0,0,680,0  3146,0,2466,0  0,3180,374,1060  3146,3180,2772,1060' \
-scale 25% \
-crop 2400x1060+374+0 "$1.t2c.png"

convert "$1.3.png" -alpha set -virtual-pixel transparent \
-set option:distort:scale 4 -distort BilinearForward \
'0,0,374,0  3146,0,2772,0  0,1590,373,530  3146,1591,2772,530' \
-scale 25% \
-crop 2400x530+374+0 "$1.t3c.png"

convert "$1.t1c.png" "$1.t2c.png" "$1.t3c.png" \
-background navy -append +repage \
-crop 2400x3600+0+0 "out.$1.png"

identify "out.$1.png"


