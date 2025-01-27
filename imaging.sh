#!/bin/bash
flagval="f"

targets='AS_206  DoAr_16  DoAr_24E  DoAr_25  DoAr_33  DoAr_44  GSS_26  GSS_39  HBC_266  IRS_37  IRS_39  IRS_41  IRS_51  VSSG_1  WSB_31  WSB_60  YLW_47  YLW_8'
tracks='track2'
sidebands='lsb usb'
rxs='rx230 rx240'

for track in $tracks
do
	cd ch0/
	rm -rf $track
  	mkdir $track
 	cd ../

	for rx in $rxs
	do
		for sideband in $sidebands
		do
			if [ $track = 'track1' ]
			then
				cellsize=0.125
				imsize=1024
			fi
			if [ $track = 'track2' ]
			then
				cellsize=0.125
				imsize=1024
			fi

			for target in $targets
			do
				datadir='../imaging/'$track'/'
 			        vis=$target'_'$track'.'$rx'.'$sideband'.cal.miriad'
          			cp -r $datadir$vis ./

				uvflag vis=$vis edge=64,64,0 flagval="f"
				#uvaver vis=$vis options=nocal,nopass,nopol, line=channel,1,1,1024,1, out=$vis'.c'

				# creating non-self-calibrated image
                                rm -rf $vis'.dirty'
                                rm -rf $vis'.beam'
                                invert vis=$vis \
					options=systemp,double,mfs robust=2.0 \
					map=$vis'.dirty' \
					beam=$vis'.beam' cell=$cellsize imsize=$imsize

				rm -rf $vis'.dirty.fits'
				rm -rf $vis'.beam.fits'
				fits in=$vis'.dirty' op=xyout out=$vis'.dirty.fits'
				fits in=$vis'.beam' op=xyout out=$vis'.beam.fits'

				rm -rf $vis'.10.model'
				rm -rf $vis'.10.model.fits'
				clean map=$vis'.dirty' \
					beam=$vis'.beam'\
					out=$vis'.10.model' cutoff=0.01 niters=10 \
					options=positive
				fits in=$vis'.10.model' op=xyout out=$vis'.10.model.fits'

				rm -rf $vis'.clean'
				rm -rf $vis'.clean.fits'
				restor map=$vis'.dirty' \
					beam=$vis'.beam' \
					model=$vis'.10.model' \
					mode=clean out=$vis'.clean'
				fits in=$vis'.clean' op=xyout out=$vis'.clean.fits'

				rm -rf $vis'.residual'
				rm -rf $vis'.residual.fits'
				restor map=$vis'.dirty' \
                                        beam=$vis'.beam' \
                                        model=$vis'.10.model' \
                                        mode=residual out=$vis'.residual'
				fits in=$vis'.residual' op=xyout out=$vis'.residual.fits'

				output=$(python get_rms.sel.py  $rx  $sideband  $target  $track)
				IFS='   ' read -r -a array <<< "$output"
				rms=${array[0]}
				cut=$(bc -l <<< "${array[0]}*1.5")
				echo "The obtained rms for $target is ${array[0]} Jy/beam"

				rm -rf $vis'.dirty'
				rm -rf $vis'.beam'
				invert vis=$vis \
					options=systemp,double,mfs robust=2.0 \
					map=$vis'.dirty' \
					beam=$vis'.beam' cell=$cellsize imsize=$imsize

				rm -rf $vis'.dirty.fits'
                                rm -rf $vis'.beam.fits'
                                fits in=$vis'.dirty' op=xyout out=$vis'.dirty.fits'
                                fits in=$vis'.beam' op=xyout out=$vis'.beam.fits'

				rm -rf $vis'.model'
				rm -rf $vis'.model.fits'
				clean map=$vis'.dirty' \
					beam=$vis'.beam' \
					out=$vis'.model' cutoff=$cut niters=1000
				fits in=$vis'.model' op=xyout out=$vis'.model.fits'

				rm -rf $vis'.clean'
				rm -rf $vis'.clean.fits'
				restor map=$vis'.dirty' \
					beam=$vis'.beam' \
					model=$vis'.model' \
					mode=clean out=$vis'.clean'
				fits in=$vis'.clean' op=xyout out=$vis'.clean.fits'

				rm -rf $vis'.residual'
				rm -rf $vis'.residual.fits'
				restor map=$vis'.dirty' \
					beam=$vis'.beam' \
					model=$vis'.model' \
					mode=residual out=$vis'.residual'
				fits in=$vis'.residual' op=xyout out=$vis'.residual.fits'

				rm -rf $vis'.clean.pbcor'
				rm -rf $vis'.clean.pbcor.fits'
				linmos in=$vis'.clean' out=$vis'.clean.pbcor'
				fits in=$vis'.clean.pbcor' op=xyout out=$vis'.clean.pbcor.fits'



				mv *.beam ./ch0/$track/
                                mv *.model ./ch0/$track/
                                mv *.clean ./ch0/$track/
                                mv *.pbcor ./ch0/$track/
                                mv *.dirty ./ch0/$track/
                                mv *.residual ./ch0/$track/
                                mv *.fits ./ch0/$track/
                                rm -rf *.miriad
				rm -rf *.c

			done
		done
	done
done
