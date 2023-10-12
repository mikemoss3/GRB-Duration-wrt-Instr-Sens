import numpy as np
from astropy.io import fits

def find_grid_id(imx,imy):
	"""
	Method to find the Swift/BAT response matrix GridID based on the position of the source on the detector plane according to Lien et al 2012.
	"""

	# Load table of GridIDs and imx,imy positions
	gridnum_imx_imy = np.genfromtxt("./util_packages/files-det-ang-dependence/gridnum_imx_imy.txt",dtype=[("GRIDID","U3"),("imx",float),("imy",float),("theta",float)])
	# Based on imx and imy, determine which grid number to use
	imx_list_cut1 = np.argwhere(gridnum_imx_imy['imx']<=imx+0.25).T[0]
	imx_list_cut2 = np.argwhere(gridnum_imx_imy['imx'][imx_list_cut1]>=imx-0.25).T[0]
	imy_list_cut1 = np.argwhere(gridnum_imx_imy['imy'][imx_list_cut1][imx_list_cut2] <= imy+0.17).T[0]
	imy_list_cut2 = np.argwhere(gridnum_imx_imy['imy'][imx_list_cut1][imx_list_cut2][imy_list_cut1] >= imy-0.17).T[0]
	gridid = gridnum_imx_imy['GRIDID'][imx_list_cut1][imx_list_cut2][imy_list_cut1][imy_list_cut2][0]

	return gridid

def find_inc_ang(imx,imy):
	"""
	Method to calculate the incidence angle from a given position on the detector plane 
	"""

	theta = np.arctan( np.sqrt( imx**2 + imy**2 ) )
	return theta

def find_pcode(imx,imy):
	"""
	Method to calculate the partial coding fraction on the detector plane for a given position on the detector plane 
	"""

	# Load pcode map image 
	pcode_img = fits.getdata("./util_packages/files-det-ang-dependence/pcode-map.img",ext=0) # indexing as pcode_img[y-index, x-index]
	# Load header from file
	pcode_img_header = fits.getheader("./util_packages/files-det-ang-dependence/pcode-map.img",ext=0)

	# Make (imx, imy) grid based on the indices (i,j)
	# i and j are the indices of each pixel
	# From headher files (the T at the end of the field name indicates tangent position): 
	# crpix: Reference pixel position
	# cdelt: Pixel spacing in physical units
	# crval: Coordinate value at reference pixel position (seems to be zero most of the time)
	i = int( ( ( imx - pcode_img_header["CRVAL1T"]) / pcode_img_header["CDELT1T"]) + pcode_img_header["CRPIX1T"] )
	j = int( ( ( imy - pcode_img_header["CRVAL2T"]) / pcode_img_header["CDELT2T"]) + pcode_img_header["CRPIX2T"] )

	# The given imx,imy may be calculated to be in the center of a pixel. To make this compatible with calling an index, we force it to be a integer.
	# This has the result of rounding down.

	return pcode_img[j,i]