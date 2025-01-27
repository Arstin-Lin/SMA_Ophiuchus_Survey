#!/bin/bash

targets='DoAr_44 AS_206 DoAr_16 DoAr_24E DoAr_25 DoAr_33 GSS_26 GSS_39 HBC_266 IRS_37 IRS_39 IRS_41 IRS_51 VSSG_1 WSB_31 WSB_60 YLW_47 YLW_8'
sidebands='lsb usb'
rxs='rx240 rx345'

# Combine visibility data
for target in $targets
do
    for sideband in $sidebands
    do
        for rx in $rxs
        do
            vis1="${target}_track1.${rx}.${sideband}.cal.miriad.c"
            vis2="${target}_track2a.${rx}.${sideband}.cal.miriad"
            output="${target}_track.${rx}.${sideband}.cal.miriad"

            echo "Checking files: '$vis1' and '$vis2'"
            if [ -e "$vis1" ] && [ -e "$vis2" ]; then
                echo "Combining: $vis1 and $vis2 into $output"
                uvaver vis="${vis1},${vis2}" options=nocal,nopass,nopol, out="$output"
            else
                echo "One or both files are missing: $vis1 or $vis2"
            fi
            mv *_track.* /home/arstinlin/DATA/imaging_Oph/imaging/track1
        done
    done
done

