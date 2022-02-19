# cd /Users/Shared/data/usfs/micmac_tiepoint_test
# cd "/Users/Shared/data/micmac/src/uti_phgrm/TiePHistorical/"
import os
from os.path import exists, join, basename, splitext
import numpy as np
import cv2
import matplotlib.pyplot as plt
import pip
# python3 -m pip install "gdown"
# python3 -m pip install "scikit-image"
import gdown


####################
# !!!!!!!!!! to resize all images using % method
# using ImageMagick
####################
# navigate to folder with images
# os.system('mogrify -resize 30%  *.tif')

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


os.chdir('/Users/Shared/data/usfs/micmac_tiepoint_test/content/historical_data')

##############################################################
# Epoch 1
##############################################################

# Extract SIFT tie-points between images within the same epoch with command "Tapioca" (i.e., a version of SIFT). The input, output and parameter interpretation of the command "Tapioca" are listed below:
os.system('mm3d Tapioca MulScale OIS-Reech_IGNF_PVA_1-0__1981.*tif 500 -1 PostFix=_1981')

# Create a mask
# Since historical images contain fiducial marks and they would yield nonsensical tie-points, we create a mask to remove any points in the vicinity of the fiducial marks. To create the mask we use the "SaisieMasq" program:
os.system('mm3d SaisieMasq OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0083.tif Name=Fiducial_marks_masq-1981-3.tif')

##!!! ^ this opens a window to select the boundaries of the image : https://youtu.be/dMxhD7UBEGE
# https://micmac.ensg.eu/index.php/SaisieMasqQT


# Remove tie-points on the fiducial marks
os.system('mm3d HomolFilterMasq OIS-Reech_IGNF_PVA_1-0__1981.*tif GlobalMasq=Fiducial_marks_masq-1981-3.tif PostIn=_1981 PostOut=_1981-Masq')


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

##############################################################
# Epoch 2
##############################################################

#Recover tie-points
os.system('mm3d Tapioca MulScale OIS-Reech_IGNF_PVA_1-0__1971.*tif 500 -1 PostFix=_1971')

# Create a mask
# Since historical images contain fiducial marks and they would yield nonsensical tie-points, we create a mask to remove any points in the vicinity of the fiducial marks. To create the mask we use the "SaisieMasq" program:
os.system('mm3d SaisieMasq OIS-Reech_IGNF_PVA_1-0__1971-06-21__C2844-0141_1971_FR2117_0972.tif Name=Fiducial_marks_masq-1971-3.tif')

##!!! ^ this opens a window to select the boundaries of the image : https://youtu.be/dMxhD7UBEGE
# https://micmac.ensg.eu/index.php/SaisieMasqQT


# Remove tie-points on the fiducial marks
os.system('mm3d HomolFilterMasq OIS-Reech_IGNF_PVA_1-0__1971.*tif GlobalMasq=Fiducial_marks_masq-1971-3.tif PostIn=_1971 PostOut=_1971-Masq')

# Binary to txt conversion
# The resulting tie-points are in binary format. To visualise them we will convert them to a txt format. We can use the same command "HomolFilterMasq" to perform this conversion.

os.system('mm3d HomolFilterMasq "OIS-Reech_IGNF_PVA_1-0__1971.*tif" PostIn=_1971 PostOut=_1971-TXT ANM=1 ExpTxt=0 ExpTxtOut=1')


# Tie-points reduction
# The intra-epoch tie-points are generally very abundant. We reduce their quantity for several reasons:

# quantity of tie points does not impact the quality of the bundle adjustement (only the distribution is crucial)
# we wish to speed-up the bundle adjustement
# we want to "balance" the number of inter- and intra-epoch tie-points

os.system('mm3d TestLib NO_AllOri2Im OIS-Reech_IGNF_PVA_1-0__1971.*tif SH=_1971-Masq')
os.system('mm3d Ratafia OIS-Reech_IGNF_PVA_1-0__1971.*tif SH=_1971-Masq Out=_1971-Ratafia')

# Binary to txt conversion
# The resulting tie-points are in binary format. To visualise them we will convert them to a txt format. We can use the same command "HomolFilterMasq" to perform this conversion.
os.system('mm3d HomolFilterMasq "OIS-Reech_IGNF_PVA_1-0__1971.*tif" PostIn=_1971-Ratafia PostOut=_1971-Ratafia-TXT ANM=1 ExpTxt=0 ExpTxtOut=1')

# Visualize tie points?
import mm3d_utils 

ImgName1 = 'OIS-Reech_IGNF_PVA_1-0__1971-06-21__C2844-0141_1971_FR2117_0972.tif'
ImgName2 = 'OIS-Reech_IGNF_PVA_1-0__1971-06-21__C2844-0141_1971_FR2117_0973.tif'
aIm1 = cv2.imread(ImgName1,cv2.IMREAD_IGNORE_ORIENTATION)
aIm2 = cv2.imread(ImgName2,cv2.IMREAD_IGNORE_ORIENTATION) 
TPtsVec = mm3d_utils.ImportHom('Homol_1971-Ratafia-TXT/Pastis'+ImgName1+'/'+ImgName2+'.txt') 

mm3d_utils.plot_images([np.asarray(aIm1),np.asarray(aIm2)]) 
mm3d_utils.plot_tiepts2(np.asarray(TPtsVec,dtype=float))
plt.show()


##############################################################
# Relative orientation
##############################################################

# Recover relative orientations of images within the same epoch using structure-from-motion implemented in MicMac under the command "Tapas".

############
# Epoch 1
############

os.system('mm3d Tapas FraserBasic OIS-Reech_IGNF_PVA_1-0__1981.*tif Out=1981 SH=_1981-Masq')

############
# Epoch 2
############

os.system('mm3d Tapas FraserBasic OIS-Reech_IGNF_PVA_1-0__1971.*tif Out=1971 SH=_1971-Masq')

##############################################################
# DSM generation
##############################################################

# Compute DSM of each epoch based on orientations using the command "Malt".

############
# Epoch 1
############

os.system('mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1981.*tif 1981 NbVI=2 MasqImGlob=Fiducial_marks_masq-1981-3.tif DirMEC=MEC-Malt_1981 EZA=1 ZoomF=2 DoOrtho=0')

# Visualize DSM
import mm3d_utils 
from skimage import io
aIm1 = io.imread('MEC-Malt_1981/Z_Num8_DeZoom2_STD-MALT.tif')
mm3d_utils.plot_images([np.asarray(aIm1)])
plt.show()


############
# Epoch 2
############

os.system('mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1971.*tif 1971 NbVI=2 MasqImGlob=Fiducial_marks_masq-1971-3.tif DirMEC=MEC-Malt_1971 EZA=1 ZoomF=2 DoOrtho=0')

##############################################################
# Inter-epoch processing
##############################################################

# Use the TiePHistoP program to launch the automated inter-epoch processing pipeline.
# The input, output and parameter interpretation of the command "TiePHistoP" are listed below:

############
# Option 1: SuperGlue
############
# This command will produce 2 kinds of results:
# (1) roughly co-registered orientations, which will be stored in the folder "./Ori-1981";
# (2) inter-epoch tie-points, which will be stored in the folder "./Homol-SuperGlue-3DRANSAC-CrossCorrelation".

os.system('mm3d TiePHistoP Ori-1971 Ori-1981 ImgList1971all.txt ImgList1981all.txt MEC-Malt_1971 MEC-Malt_1981 CoRegPatchLSz=[1280,960] CoRegPatchRSz=[1280,960] PrecisePatchSz=[1280,960] Feature=SuperGlue')

############
# Visualize all the inter-epoch tie-points
############

# Check the distribution of the tie-points across the entire image block by visualising them in a common reference frame with the "TestLib VisuTiePtIn3D".

# os.system('mm3d TestLib VisuTiePtIn3D ImgList1971all.txt ImgList1981all.txt Ori-1971 Ori-1981 DSMDirL=MEC-Malt_1971 DSMDirR=MEC-Malt_1981 InSH=-SuperGlue-3DRANSAC-CrossCorrelation OutFile=Visu-SuperGlue-3DRANSAC-CrossCorrelation-all.txt')

# import mm3d_utils
# from skimage import io

# ImgName1 = 'MEC-Malt_1971/Z_Num8_DeZoom2_STD-MALT_gray.tif_sfs.tif'
# ImgName2 = 'MEC-Malt_1981/Z_Num8_DeZoom2_STD-MALT_gray.tif_sfs.tif'
# aIm1 = io.imread(ImgName1)
# aIm2 = io.imread(ImgName2)

# def GetPtsInDSM(DSMTFWL, Pt3DL, DSMTFWR, Pt3DR):
#   dsmInfoL = np.loadtxt(DSMTFWL, dtype=float)
#   ptL = np.loadtxt(Pt3DL, dtype=float)
#   dsmInfoR = np.loadtxt(DSMTFWR, dtype=float)
#   ptR = np.loadtxt(Pt3DR, dtype=float)
#   Result = []
#   for i in range(len(ptL)):
#     x1_2DL = (ptL[i, 0] - dsmInfoL[4])/dsmInfoL[0]
#     y1_2DL = (ptL[i, 1] - dsmInfoL[5])/dsmInfoL[3]
#     x1_2DR = (ptR[i, 0] - dsmInfoR[4])/dsmInfoR[0]
#     y1_2DR = (ptR[i, 1] - dsmInfoR[5])/dsmInfoR[3]
#     Result.append([[x1_2DL, y1_2DL],[x1_2DR, y1_2DR]])
#   return Result

# DSMTFWL = 'MEC-Malt_1971/Z_Num8_DeZoom2_STD-MALT.tfw'
# Pt3DL = 'Visu-SuperGlue-3DRANSAC-CrossCorrelation-all_L.txt'
# DSMTFWR = 'MEC-Malt_1981/Z_Num8_DeZoom2_STD-MALT.tfw'
# Pt3DR = 'Visu-SuperGlue-3DRANSAC-CrossCorrelation-all_R.txt'
# TPtsVec = GetPtsInDSM(DSMTFWL, Pt3DL, DSMTFWR, Pt3DR)

# print('                     epoch 1971                                          epoch 1981')
# mm3d_utils.plot_images([np.asarray(aIm1),np.asarray(aIm2)])
# if len(TPtsVec)>0:
#   mm3d_utils.plot_tiepts2(np.asarray(TPtsVec,dtype=float))
# else:
#   print('tie-point number is 0')

############
# Option 2: SIFT
############

# Note: (1) We set Feature=SIFT to switch to option SIFT; (2) the rough co-registration has been performed in the previous step, so we can skip it here; (3) the PrecisePatchSz is the same with the previous step, therefore we can skip the step of getting patch pairs.

# os.system('mm3d TiePHistoP Ori-1971 Ori-1981 ImgList1971all.txt ImgList1981all.txt MEC-Malt_1971 MEC-Malt_1981 PrecisePatchSz=[1280,960] Feature=SIFT SkipCoReg=1  SkipGetPatchPair=1')

# # not sure what this call is for 
# # os.system('!mm3d TiePHistoP Ori-1971 Ori-1981 OIS-Reech_IGNF_PVA_1-0__1971-06-21__C2844-0141_1971_FR2117_1018.tif OIS-Reech_IGNF_PVA_1-0__1981-06-16__C2544-0021_1981_F2544-2644_0083.tif MEC-Malt_1971 MEC-Malt_1981 PrecisePatchSz=[1280,960] Feature=SIFT SkipCoReg=1  SkipGetPatchPair=0 CCTh=0.1')

##############################################################
# 3. Evaluation
##############################################################
# In this section we will compare the DoD (Difference of DSMs) for evaluation.

# We use 3 sets of image orientations to compute the DSMs in epoch 1971 and 1981 individually and display the difference between them. The 3 sets of orientations are:

# (1) Roughly co-registered orientations: which are the orientations we got at the end of rough co-registration.

# (2) SuperGlue refined orientations: based on the roughly co-registered orientations, we refine them in a boudle adjustment routine using the SuperGlue inter-epoch tie-points.

# (3) SIFT refined orientations: based on the roughly co-registered orientations, we refine them in a boudle adjustment routine using the SIFT inter-epoch tie-points.

############################
# 3.1. DoD of roughly co-registered result
############################
# Get DSM of epoch 1981
# As the co-registered orientations are based on the reference of epoch 1981, we can use directly the DSM of epoch 1981 resulted from section 1.3.1

# Get DSM of epoch 1971
!mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1971.*tif 1981 NbVI=2 DirMEC=MEC-Malt_1971_CoReg EZA=1 MasqImGlob=Fiducial_marks_masq-1971-3.tif ZoomF=4 DoOrtho=0

#### Calculate DoD
# We use the command "CmpIm" to generate the DoD.

!mm3d CmpIm MEC-Malt_1971_CoReg/Z_Num7_DeZoom4_STD-MALT.tif MEC-Malt_1981/Z_Num8_DeZoom2_STD-MALT.tif UseFOM=1 FileDiff=DoD-CoReg.tif 16Bit=1

#### Visualize DoD
# The resulted DoD is visulized below:

# import mm3d_utils 
# from skimage import io

# def Convert2Gray(data, minV, maxV):
#     data[data<minV] = minV
#     data[data>maxV] = maxV
#     basemaxV = np.max(data)
#     baseminV = np.min(data)
#     data = (data-baseminV)*255.0/(basemaxV-baseminV)
#     return data

# aIm1 = np.asarray(io.imread('DoD-CoReg.tif'))

# minV = -0.015
# maxV = 0.015
# aIm1 = Convert2Gray(aIm1, minV, maxV)

# aLegend = np.arange(0, 255, 5)
# aLegend = np.vstack((aLegend, aLegend))
# mm3d_utils.plot_DoD([np.asarray(aLegend)])
# print('Legend:')
# print(minV, '                                                                                    ', maxV)

# mm3d_utils.plot_DoD([np.asarray(aIm1)])
# print('                                             DoD-CoReg')

############################
# 3.2. DoD of SuperGlue refined result
############################
# We need to refine the roughly co-registered orientations in a bundle adjustment (BA) routine using the SuperGlue inter-epoch tie-points resulted from section 2.1, as well as the reduced intra-epoch tie-points resulted from section 1.1.

# Set weight of inter-epoch tie-points
# First of all, we use the command "TestLib TiePtAddWeight" to set the weight of the inter-epoch tie-points to be 10, so that they will play a more important role in BA. (Please notice that the weight of the intra-epoch tie-points is by default 1.)

!mm3d TestLib TiePtAddWeight 10 InSH=-SuperGlue-3DRANSAC-CrossCorrelation

# Txt to binary conversion
# The SuperGlue inter-epoch tie-points we got are in txt format, we should transform them into binary format with the help of "HomolFilterMasq", so that they can be recognized in the following process.

!mm3d HomolFilterMasq "O.*tif" PostIn=-SuperGlue-3DRANSAC-CrossCorrelation-W10 PostOut=-SuperGlue-3DRANSAC-CrossCorrelation-W10-dat ANM=1 ExpTxt=1 ExpTxtOut=0

# Merge intra- and inter-epoch tie-points
# Then we need to merge the intra- and inter-epoch tie-points from different folders together using the command "MergeHomol".

!mm3d MergeHomol "Homol_1971-Ratafia|Homol_1981-Ratafia|Homol-SuperGlue-3DRANSAC-CrossCorrelation-W10-dat" Homol_Merged-SuperGlue

# Run bundle adjustment
# Now it is time to run BA with the command "Campari".

!mm3d Campari "O.*tif" 1981 Campari_Refined-SuperGlue SH=_Merged-SuperGlue AllFree=1 NbIterEnd=20 SigmaTieP=0.25 

# Get DSM of epoch 1981

!mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1981.*tif Campari_Refined-SuperGlue NbVI=2 DirMEC=MEC-Malt_1981_Refined-SuperGlue EZA=1 MasqImGlob=Fiducial_marks_masq-1981-3.tif ZoomF=2 DoOrtho=0

# Get DSM of epoch 1971

!mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1971.*tif Campari_Refined-SuperGlue NbVI=2 DirMEC=MEC-Malt_1971_Refined-SuperGlue MasqImGlob=Fiducial_marks_masq-1971-3.tif EZA=1 ZoomF=4 DoOrtho=0

#### Calculate DoD
# Finally we use the command "CmpIm" to generate the DoD.

!mm3d CmpIm MEC-Malt_1971_Refined-SuperGlue/Z_Num7_DeZoom4_STD-MALT.tif MEC-Malt_1981_Refined-SuperGlue/Z_Num8_DeZoom2_STD-MALT.tif UseFOM=1 FileDiff=DoD-Refined-SuperGlue.tif 16Bit=1


#### Visualize DoD
# The resulted DoD is visulized and compared to the DoD of section 3.1 below:

# import mm3d_utils 
# from skimage import io

# def Convert2Gray(data, minV, maxV):
#     data[data<minV] = minV
#     data[data>maxV] = maxV
#     basemaxV = np.max(data)
#     baseminV = np.min(data)
#     data = (data-baseminV)*255.0/(basemaxV-baseminV)
#     return data

# aIm1 = np.asarray(io.imread('DoD-CoReg.tif'))
# aIm2 = np.asarray(io.imread('DoD-Refined-SuperGlue.tif'))

# minV = -0.015
# maxV = 0.015
# aIm1 = Convert2Gray(aIm1, minV, maxV)
# aIm2 = Convert2Gray(aIm2, minV, maxV)

# aLegend = np.arange(0, 255, 5)
# aLegend = np.vstack((aLegend, aLegend))
# mm3d_utils.plot_DoD([np.asarray(aLegend)])
# print('Legend:')
# print(minV, '                                                                                    ', maxV)

# mm3d_utils.plot_DoD([np.asarray(aIm1), np.asarray(aIm2)])
# print('                     DoD-CoReg                                    DoD-Refined-SuperGlue')

############################
# 3.3. DoD of GuidedSIFT refined result
############################

# Set weight of inter-epoch tie-points
# First of all, we use the command "TestLib TiePtAddWeight" to set the weight of the inter-epoch tie-points to be 10, so that they will play a more important role in BA. (Please notice that the weight of the intra-epoch tie-points is by default 1.)

!mm3d TestLib TiePtAddWeight 10 InSH=-GuidedSIFT-3DRANSAC-CrossCorrelation

# Txt to binary conversion
# The SIFT inter-epoch tie-points we got are in txt format, we should transform them into binary format with the help of "HomolFilterMasq", so that they can be recognized in the following process.

!mm3d HomolFilterMasq "O.*tif" PostIn=-GuidedSIFT-3DRANSAC-CrossCorrelation-W10 PostOut=-GuidedSIFT-3DRANSAC-CrossCorrelation-W10-dat ANM=1 ExpTxt=1 ExpTxtOut=0

# Merge intra- and inter-epoch tie-points
# Then we need to merge the intra- and inter-epoch tie-points from different folders together using the command "MergeHomol".

!mm3d MergeHomol "Homol_1971-Ratafia|Homol_1981-Ratafia|Homol-SuperGlue-3DRANSAC-CrossCorrelation-W10-dat" Homol_Merged-GuidedSIFT

# Run bundle adjustment
# Now it is time to run BA with the command "Campari".

!mm3d Campari "O.*tif" 1981 Campari_Refined-GuidedSIFT SH=_Merged-GuidedSIFT AllFree=1 NbIterEnd=20 SigmaTieP=0.25 

# Get DSM of epoch 1981
# Based on the GuidedSIFT refined orientations "Campari_Refined-GuidedSIFT", we compute the DSMs in epoch 1981 using the command "Malt":

!mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1981.*tif Campari_Refined-GuidedSIFT NbVI=2 DirMEC=MEC-Malt_1981_Refined-GuidedSIFT EZA=1 MasqImGlob=Fiducial_marks_masq-1981-3.tif ZoomF=2 DoOrtho=0

# Get DSM of epoch 1971
# Based on the GuidedSIFT refined orientations "Campari_Refined-GuidedSIFT", we compute the DSMs in epoch 1971 using the command "Malt":

!mm3d Malt Ortho OIS-Reech_IGNF_PVA_1-0__1971.*tif Campari_Refined-GuidedSIFT NbVI=2 DirMEC=MEC-Malt_1971_Refined-GuidedSIFT EZA=1 MasqImGlob=Fiducial_marks_masq-1971-3.tif ZoomF=4 DoOrtho=0

# Calculate DoD
# Finally we use the command "CmpIm" to generate the DoD.

!mm3d CmpIm MEC-Malt_1971_Refined-GuidedSIFT/Z_Num7_DeZoom4_STD-MALT.tif MEC-Malt_1981_Refined-GuidedSIFT/Z_Num8_DeZoom2_STD-MALT.tif UseFOM=1 FileDiff=DoD-Refined-GuidedSIFT.tif 16Bit=1

# Visualize DoD
# The resulted DoD is visulized and compared to the other DoDs below:

# import mm3d_utils 
# from skimage import io

# def Convert2Gray(data, minV, maxV):
#     data[data<minV] = minV
#     data[data>maxV] = maxV
#     basemaxV = np.max(data)
#     baseminV = np.min(data)
#     data = (data-baseminV)*255.0/(basemaxV-baseminV)
#     return data

# aIm1 = np.asarray(io.imread('DoD-CoReg.tif'))
# aIm2 = np.asarray(io.imread('DoD-Refined-SuperGlue.tif'))
# aIm3 = np.asarray(io.imread('DoD-Refined-GuidedSIFT.tif'))

# minV = -0.015
# maxV = 0.015
# aIm1 = Convert2Gray(aIm1, minV, maxV)
# aIm2 = Convert2Gray(aIm2, minV, maxV)
# aIm3 = Convert2Gray(aIm3, minV, maxV)

# aLegend = np.arange(0, 255, 5)
# aLegend = np.vstack((aLegend, aLegend))
# mm3d_utils.plot_DoD([np.asarray(aLegend)])
# print('Legend:')
# print(minV, '                                                                                    ', maxV)

# mm3d_utils.plot_DoD([np.asarray(aIm1), np.asarray(aIm2), np.asarray(aIm3)])
# print('            DoD-CoReg                DoD-Refined-SuperGlue                DoD-Refined-GuidedSIFT')