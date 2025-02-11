#!/bin/bash
rm -r t
bash ./build.sh
if [[ $? -ne 0 ]]; then
    echo "Error: build.sh failed."
    exit 1
fi
./GFM <<EOF
0.8
t
EOF

#rm -r /mnt/c/Users/olive/MPhilData/t 
#cp -r t /mnt/c/Users/olive/MPhilData.
./plot rho
