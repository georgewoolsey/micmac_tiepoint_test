# cd /Users/Shared/data/usfs/micmac_tiepoint_test
# cd "/Users/Shared/data/micmac/src/uti_phgrm/TiePHistorical/"
import os
from os.path import exists, join, basename, splitext
import numpy as np
import cv2
import matplotlib.pyplot as plt
import pip
# python3 -m pip install "gdown"
import gdown

# download data
os.getcwd()
os.chdir('/Users/Shared/data/usfs/micmac_tiepoint_test')
# os.mkdir('content')
os.chdir('./content/')

# download data
url = 'https://drive.google.com/uc?id=1W299CM_lZk2KIgLV0odjxanz74T6Humy' 
output = 'historical_data.tar.gz'
gdown.download(url, output, quiet=False)

# unpack
#!unzip historical_data.zip $YOUR_PATH
os.system('tar -xf historical_data.tar.gz -C /Users/Shared/data/usfs/micmac_tiepoint_test/content')

os.chdir('historical_data')

# utility functions to visualise tie-points
utils_url = 'https://drive.google.com/uc?id=1ATO1Nz_aXApxVnm6l7x1xappGXtcjuvp'
gdown.download(utils_url, 'mm3d_utils.py', quiet=False)

# Extract SIFT tie-points between images within the same epoch with command "Tapioca" (i.e., a version of SIFT). The input, output and parameter interpretation of the command "Tapioca" are listed below:
os.system('mm3d Tapioca MulScale OIS-Reech_IGNF_PVA_1-0__1981.*tif 500 -1 PostFix=_1981')

# Create a mask
# Since historical images contain fiducial marks and they would yield nonsensical tie-points, we create a mask to remove any points in the vicinity of the fiducial marks. To create the mask we use the "SaisieMasq" program:
os.system('mm3d SaisieMasq OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0083.tif Name=Fiducial_marks_masq-1981-99.tif')

##!!! ^ this opens a window to select the boundaries of the image : https://youtu.be/dMxhD7UBEGE
# https://micmac.ensg.eu/index.php/SaisieMasqQT


# Remove tie-points on the fiducial marks
os.system('mm3d HomolFilterMasq OIS-Reech_IGNF_PVA_1-0__1981.*tif GlobalMasq=Fiducial_marks_masq-1981-99.tif PostIn=_1981 PostOut=_1981-Masq')


# Binary to txt conversion
# The resulting tie-points are in binary format. To visualise them we will convert them to a txt format. We can use the same command "HomolFilterMasq" to perform this conversion.

os.system('mm3d HomolFilterMasq "OIS-Reech_IGNF_PVA_1-0__1981.*tif" PostIn=_1981 PostOut=_1981-TXT ANM=1 ExpTxt=0 ExpTxtOut=1')

# Visualize tie points?
import mm3d_utils 

ImgName1 = 'OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0064.tif'
ImgName2 = 'OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0065.tif'
aIm1 = cv2.imread(ImgName1,cv2.IMREAD_IGNORE_ORIENTATION)
aIm2 = cv2.imread(ImgName2,cv2.IMREAD_IGNORE_ORIENTATION) 
TPtsVec = mm3d_utils.ImportHom('Homol_1981-TXT/Pastis'+ImgName1+'/'+ImgName2+'.txt') 

#img 1
mm3d_utils.plot_images([np.asarray(aIm1),np.asarray(aIm2)])  # np = numpy
len(TPtsVec)
mm3d_utils.plot_tiepts2(np.asarray(TPtsVec,dtype=float))
plt.show()


# Tie-points reduction
# The intra-epoch tie-points are generally very abundant. We reduce their quantity for several reasons:

# quantity of tie points does not impact the quality of the bundle adjustement (only the distribution is crucial)
# we wish to speed-up the bundle adjustement
# we want to "balance" the number of inter- and intra-epoch tie-points

os.system('mm3d TestLib NO_AllOri2Im OIS-Reech_IGNF_PVA_1-0__1981.*tif SH=_1981-Masq')
os.system('mm3d Ratafia OIS-Reech_IGNF_PVA_1-0__1981.*tif SH=_1981-Masq Out=_1981-Ratafia')

# Binary to txt conversion
# The resulting tie-points are in binary format. To visualise them we will convert them to a txt format. We can use the same command "HomolFilterMasq" to perform this conversion.
os.system('mm3d HomolFilterMasq "OIS-Reech_IGNF_PVA_1-0__1981.*tif" PostIn=_1981-Ratafia PostOut=_1981-Ratafia-TXT ANM=1 ExpTxt=0 ExpTxtOut=1')

# Visualize tie points?
import mm3d_utils 

ImgName1 = 'OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0064.tif'
ImgName2 = 'OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0065.tif'
aIm1 = cv2.imread(ImgName1,cv2.IMREAD_IGNORE_ORIENTATION)
aIm2 = cv2.imread(ImgName2,cv2.IMREAD_IGNORE_ORIENTATION) 
TPtsVec = mm3d_utils.ImportHom('Homol_1981-Ratafia-TXT/Pastis'+ImgName1+'/'+ImgName2+'.txt') 

mm3d_utils.plot_images([np.asarray(aIm1),np.asarray(aIm2)]) 
mm3d_utils.plot_tiepts2(np.asarray(TPtsVec,dtype=float))
plt.show()

