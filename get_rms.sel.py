import sys
import os
import numpy as np
import math

from astropy.io.fits import getdata
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy import constants as con
from astropy.coordinates import SkyCoord

import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib


targ_list=['AS_206', 'DoAr_16', 'DoAr_24E', 'DoAr_25', 'DoAr33', 'DoAr44', 'GSS_26', 'GSS_39', 'HBC_266', 'IRS_37', 'IRS_39', 'IRS_41', 'IRS_51', 'VSSG_1', 'WSB_31', 'WSB_60', 'YLW_47', 'YLW_8']

ifband = str(sys.argv[1])
sideband = str(sys.argv[2])
field = str(sys.argv[3])
track = str(sys.argv[4])

def write_to_file(file, field, value):
    txt_file = file
    new_row = 1
    add_line = str(value)+'   '
    try:
        with open(txt_file, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if line.split()[0] == field:
                new_line=line[:-1]+add_line+'\n'
                lines[i]=new_line
                new_row = 0
                break
        if new_row == 1:
            with open(txt_file, 'a') as f:
                f.write('\n')
                f.write(''.join(field+'   '+add_line))
        else:
            with open(txt_file, 'w') as f:
                f.writelines(lines)
    except:
        with open(txt_file, 'w') as f:
            f.write(''.join(field+'   '+add_line))

#directory = ''
modelmap = field+'_'+track+'.'+ifband+'.'+sideband+'.cal.miriad.10.model.fits'
residualmap = field+'_'+track+'.'+ifband+'.'+sideband+'.cal.miriad.residual.fits'
dirtymap = field+'_'+track+'.'+ifband+'.'+sideband+'.cal.miriad.dirty.fits'
cleanmap = field+'_'+track+'.'+ifband+'.'+sideband+'.cal.miriad.clean.fits'

if_success = False
try:

    # importing FITS image to a HDU
    rhdu   = fits.open(residualmap)
    mhdu   = fits.open(modelmap)
    dhdu   = fits.open(dirtymap)
    chdu   = fits.open(cleanmap)

    # editing the FITS image by multiplying a scaling factor
    residual_img = rhdu[0].data[0][0]
    model_img = mhdu[0].data[0][0]
    dirty_img = dhdu[0].data[0][0]
    if_success = True

except:
    print('Unable to read the intensity FITS image. Please double check the image file.')

if ( if_success == True ):
      # Reading FITS header
    try:
        naxis1 = rhdu[0].header['naxis1']
        naxis2 = rhdu[0].header['naxis2']
        crval1 = rhdu[0].header['crval1']
        crpix1 = rhdu[0].header['crpix1']
        cdelt1 = rhdu[0].header['cdelt1']
        crval2 = rhdu[0].header['crval2']
        crpix2 = rhdu[0].header['crpix2']
        cdelt2 = rhdu[0].header['cdelt2']
        hduwcs = wcs.WCS( rhdu[0].header)
        hduwcs = hduwcs.dropaxis(dropax=2)
        hduwcs = hduwcs.dropaxis(dropax=2)
    except:
        print( 'Warning. No coordinate headers' )

    try:
        bmaj = chdu[0].header['bmaj']
        bmin = chdu[0].header['bmin']
        bpa  = chdu[0].header['bpa']
    except:
        print('Warnning. No header for synthesized beam size')

    mdl_s = list(zip(*np.where(model_img > 0)))
    # mdl_s = [(128,128)]
    rms_img = residual_img.copy()

    for cen in mdl_s:
        radius = bmaj*2/(cdelt2)
        y,x = np.ogrid[:naxis1, :naxis2]
        dist = np.sqrt((x-cen[0])**2 + (y-cen[1])**2)
        mask = dist <= radius
        rms_img[mask] = 0
    
    fig = plt.figure()
    plt.imshow(rms_img)

    rms = math.sqrt(rms_img[~np.isnan(rms_img)].std()**2 + rms_img[~np.isnan(rms_img)].mean()**2)
    sys.stdout.write(str(rms)+'   ')

else:
    rms = 0.0000000000000
    sys.stdout.write(str(rms)+'   ')

write_to_file('rms_'+track+'.sel.txt', field, rms)

#filename='center_track456.txt'
#file = open(filename, 'r')
#lines = file.readlines()
#for i, line in enumerate(lines):
#    if line.split()[0] == field:
#        for k in range(4):
#            if line.split()[k+1] != '(0.0,0.0)':
#                cen_cord = eval(line.split()[k+1])
#                cen_x = hduwcs.wcs_world2pix(cen_cord[0],cen_cord[1],0)[0]
#                cen_y = hduwcs.wcs_world2pix(cen_cord[0],cen_cord[1],0)[1]
#                break
#            else:
#               cen_x = naxis1/2
#                cen_y = naxis2/2

#box_trx = cen_x+(radius*2)
#box_try = cen_y+(radius*2)
#box_blx = cen_x-(radius*2)
#box_bly = cen_y-(radius*2)

# before obtaining track456
#r_blx = int(naxis1/2-radius*2)
#r_bly = int(naxis2/2-radius*2)
#r_trx = int(naxis1/2+radius*2)
#r_try = int(naxis2/2+radius*2)

#dirty_img = dirty_img[r_blx:r_trx,r_bly:r_try]
#peak_value = np.amax(dirty_img)
#peak_pos = np.where(dirty_img == peak_value)
#box_trx = r_blx + peak_pos[0][0]+(radius*2)
#box_try = r_bly + peak_pos[1][0]+(radius*2)
#box_blx = r_blx + peak_pos[0][0]-(radius*2)
#box_bly = r_bly + peak_pos[1][0]-(radius*2)

#sys.stdout.write(str(box_blx)+'   '+str(box_bly)+'   '+str(box_trx)+'   '+str(box_try))

