# from __future__ import print_function
# from __future__ import division # undo if in Python 2 
import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import time
import sys
import copy
#Quick fudge to make import from ../Scripts work
import sys
sys.path.append('../Scripts')
import string

import ClassCosmicIntegrator  as CI #Given settings and redshifts returns rates (2D arrays) Loads the data
from PostProcessingScripts import * 
from ClassFormationChannels_5mainchannels import * 

import pandas as pd
from astropy import units as u
from astropy import constants as const





MSSFRnameslist = []
MSSFRnameslist.append('000') # add phenomenological 

for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1
        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1

            MSSFRnameslist.append('%s%s%s'%(ind_x, ind_y, ind_z))


GSMFs = [1,2,3]
SFRs = [1,2,3]
MZs=[1,2,3]


MSSFRnameslistCSV = []
MSSFRnameslistCSV.append('.0.0.0') # add phenomenological 


for ind_GSMF, GSMF in enumerate(GSMFs):
    ind_y = ind_GSMF + 1
    for ind_MZ, MZ in enumerate(MZs):
        ind_z = ind_MZ +1

        for ind_SFR, SFR in enumerate(SFRs):
            ind_x = ind_SFR+1            




            MSSFRnameslistCSV.append('.%s.%s.%s'%(ind_x, ind_y, ind_z))



####################################
#path to the data







# def writeToRatesFile_FormationChannels(BPSmodelName='Z', DCOtype='BHNS'):






#     # DCOtype='BHNS'

#     if DCOtype=='BHNS':
#         DCOname='BHNS'
#     elif DCOtype=='BBH':
#         DCOname='BHBH'
#     elif DCOtype=='BNS':
#         DCOname='NSNS'



#     # path for files 
#     path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
#     path_ = path_dir
#     path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
#     path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
#     # read in data 
#     fdata = h5.File(path)

#     # set optimistic true if that is the variation (H) 
#     OPTIMISTIC=False
#     if BPSmodelName=='H':
#         OPTIMISTIC=True 



#     # get formation channel Seeds!
#     seedsPercentageClassic, seedsPercentageOnlyStableMT = returnSeedsPercentageClassicAndOnlyStableMT(pathCOMPASOutput=path_,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)
#     seedsClassic, percentageClassic = seedsPercentageClassic
#     seedsOnlyStableMT, percentageOnlyStableMT = seedsPercentageOnlyStableMT



#     seedsDoubleCE, percentageDoubleCE = returnSeedsPercentageDoubleCoreCEE(pathCOMPASOutput=path_,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)


#     seedsSingleCE, percentageSingleCE = returnSeedsPercentageSingleCoreCEE(pathCOMPASOutput=path_,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)



#     seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE]

#     seedsOther, percentageOther = returnSeedsPercentageOther(pathCOMPASOutput=path_,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1, channelsSeedsList=seedschannels)




#     dictChannelsBHNS = { 'classic':seedsClassic, \
#                         'immediate CE':seedsSingleCE,\
#                              'stable B no CEE':seedsOnlyStableMT, \
#                          r'double-core CE':seedsDoubleCE,  \
#                             'other':seedsOther\
#                            }


#     dictPercentages = { 'classic':percentageClassic, \
#                         'immediate CE':percentageSingleCE,\
#                              'stable B no CEE':percentageOnlyStableMT, \
#                          r'double-core CE':percentageDoubleCE,  \
#                             'other':percentageOther\
#                            } 




#     # # obtain BH and NS masses
#     # M1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
#     # M2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
#     # MBH, MNS = obtainM1BHandM2BHassymetric(M1, M2)
#     # del M1
#     # del M2


#     # get intrinsic weights

#     fparam_intrinsic = 'weights_intrinsic'
#     # get detected weights

#     fparam_detected = 'weights_detected'


#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################
#     intrinsicRates = np.zeros(len(MSSFRnameslist))
#     detectedRates = np.zeros(len(MSSFRnameslist))
#     namesEMlist = []


#     intrinsicRates_I = np.zeros(len(MSSFRnameslist))
#     detectedRates_I = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_II = np.zeros(len(MSSFRnameslist))
#     detectedRates_II = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_III = np.zeros(len(MSSFRnameslist))
#     detectedRates_III = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_IV = np.zeros(len(MSSFRnameslist))
#     detectedRates_IV = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_V = np.zeros(len(MSSFRnameslist))
#     detectedRates_V = np.zeros(len(MSSFRnameslist))

#     DCOSeeds = fdata['doubleCompactObjects']['seed'][...].squeeze()

#     for ind_mssfr, mssfr in enumerate(MSSFRnameslist):

#         # print('mssfr =',ind_mssfr, 'mssfr= ', mssfr)
#         weightheader = 'w_' + mssfr
#         # print(ind_mssfr, weightheader)
#         w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
#         w_det = fdata[fparam_detected][weightheader][...].squeeze()

#         intrinsicRates[ind_mssfr] = np.sum(w_int)
#         detectedRates[ind_mssfr] = np.sum(w_det)
        
#         for nrC, Channel in enumerate(dictChannelsBHNSList):
                        
# #             #Get the seeds that relate to sorted indices
#             seedsInterest = dictChannelsBHNS[Channel]
#             mask_C = np.in1d(DCOSeeds, np.array(seedsInterest))
#             if Channel=='classic':
#                 intrinsicRates_I[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_I[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='stable B no CEE':
#                 intrinsicRates_II[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_II[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='immediate CE':
#                 intrinsicRates_III[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_III[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='double-core CE':
#                 intrinsicRates_IV[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_IV[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='other':
#                 intrinsicRates_V[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_V[ind_mssfr] = np.sum(w_det[mask_C])



       

#     stringgg =  'AllDCOsimulation_formation_channels'

#     df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv', index_col=0)
#     namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'

#     namez0_I = BPSmodelName + 'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_I = BPSmodelName + 'channel I observed (design LVK) [yr^{-1}]'
#     namez0_II = BPSmodelName + 'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_II = BPSmodelName + 'channel II observed (design LVK) [yr^{-1}]'
#     namez0_III = BPSmodelName + 'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_III = BPSmodelName + 'channel III observed (design LVK) [yr^{-1}]'
#     namez0_IV = BPSmodelName + 'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_IV = BPSmodelName + 'channel IV observed (design LVK) [yr^{-1}]'
#     namez0_V = BPSmodelName + 'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_V = BPSmodelName + 'channel V observed (design LVK) [yr^{-1}]'



#     df[namez0] = intrinsicRates
#     df[nameObs] = detectedRates

#     df[namez0_I] = intrinsicRates_I
#     df[nameObs_I] = detectedRates_I
#     df[namez0_II] = intrinsicRates_II
#     df[nameObs_II] = detectedRates_II
#     df[namez0_III] = intrinsicRates_III
#     df[nameObs_III] = detectedRates_III 
#     df[namez0_IV] = intrinsicRates_IV
#     df[nameObs_IV] = detectedRates_IV 
#     df[namez0_V] = intrinsicRates_V
#     df[nameObs_V] = detectedRates_V  


#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg  + '.csv')


#     fdata.close() 

#     return









# #### FUNCTIONS TO INITIALIZE CSV FILES

# def initialize_CSV_files_general(DCOname='BHNS'):



#     namesEMlist=[]


#     iii=0
    

#     # CREATE PANDAS FILE 
#     nModels=26
#     BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#     NAMES = []
#     stringgg = 'AllDCOsimulation'

#     for ind_l, L in enumerate(BPSnameslist):
#         str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#         str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#         NAMES.append(str_z0)
#         NAMES.append(str_obs)
        
        


#     datas=[]

#     for i in range(len(BPSnameslist)):
#         datas.append(np.zeros_like(MSSFRnameslist))
#         datas.append(np.zeros_like(MSSFRnameslist))
        
        
#     df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#     df.columns =   df.columns.map(str)
#     df.index.names = ['xyz']
#     df.columns.names = ['m']

#     # print(df) 

#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')



# def writeToRatesFile_lightestFormsFirst(BPSmodelName='Z', DCOtype='BHNS'):
#     """writes NS-BH rate to CSV file for all models"""

#     if DCOtype=='BHNS':
#         DCOname='BHNS'
#     elif DCOtype=='BBH':
#         DCOname='BHBH'
#     elif DCOtype=='BNS':
#         DCOname='NSNS'



#     # path for files 
#     path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
#     path_ = path_dir
#     path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
#     path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
#     # read in data 
#     fdata = h5.File(path)

#     # set optimistic true if that is the variation (H) 
#     OPTIMISTIC=False
#     if BPSmodelName=='H':
#         OPTIMISTIC=True 

#     # obtain BH and NS masses
#     M1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
#     M2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
#     MBH, MNS = obtainM1BHandM2BHassymetric(M1, M2)




#     whichSN = fdata['supernovae']['whichStar'][...].squeeze() 
#     maskSNnot3 = (whichSN ==1) | (whichSN==2)


#     DCOseeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
#     seedsSN = fdata['supernovae']['randomSeed'][...].squeeze()
    
#     u, indices = np.unique(seedsSN, return_index=True)
#     # uniqueSN = seedsSN[indices]

#     # bools = np.in1d(uniqueSN, DCOseeds)



#     # maskSNspecial = (bools==1) & (maskSNnot3==1)

#     # print(len(uniqueSN), )
#     print(len(fdata['supernovae']['whichStar'][...].squeeze()))
#     print(seedsSN)



#     whichSN = fdata['supernovae']['whichStar'][...].squeeze()[indices] # get whichStar for first SN 

#     print(len(seedsSN), len(set(seedsSN)))
#     print(whichSN)
#     print(set(whichSN))
#     print(len(whichSN), len(M1), len(M2))

    
#     mask_temp = ((whichSN==2) & (M1>M2) )
#     mask_2 = ((whichSN==1) & (M1<M2))
#     print('nr of weird reversals = %s'%np.sum(mask_temp))
#     print('nr of normal reversals = %s'%np.sum(mask_2))





#     maskNSBH = ((whichSN==2) & (M1>M2) ) | ((whichSN==1) & (M1<M2) ) 



#     del M1
#     del M2




#     # get intrinsic weights

#     fparam_intrinsic = 'weights_intrinsic'
#     # get detected weights

#     fparam_detected = 'weights_detected'


#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################
#     intrinsicRates = np.zeros(len(MSSFRnameslist))
#     detectedRates = np.zeros(len(MSSFRnameslist))
#     namesEMlist = []

#     for ind_mssfr, mssfr in enumerate(MSSFRnameslist):
# #             print('mssfr =',ind_mssfr)
#         weightheader = 'w_' + mssfr
#         w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()[maskNSBH]
#         w_det = fdata[fparam_detected][weightheader][...].squeeze()[maskNSBH]

#         intrinsicRates[ind_mssfr] = np.sum(w_int)
#         detectedRates[ind_mssfr] = np.sum(w_det)



       


#     stringgg =  'lightestFormsFirst'

#     df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv', index_col=0)
#     namez0 = BPSmodelName + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs = BPSmodelName + ' observed (design LVK) [yr^{-1}]'

#     df[namez0] = intrinsicRates
#     df[nameObs] = detectedRates 


#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg  + '.csv')
#     fdata.close() 
#     return



# #### FUNCTIONS TO INITIALIZE CSV FILES

# def initialize_CSV_files_lightestBHfirst(DCOname='BHNS'):



#     namesEMlist=[]


#     iii=0
    

#     # CREATE PANDAS FILE 
#     nModels=26
#     BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#     NAMES = []
#     stringgg =  'lightestFormsFirst'

#     for ind_l, L in enumerate(BPSnameslist):
#         str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#         str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#         NAMES.append(str_z0)
#         NAMES.append(str_obs)
        
        


#     datas=[]

#     for i in range(len(BPSnameslist)):
#         datas.append(np.zeros_like(MSSFRnameslist))
#         datas.append(np.zeros_like(MSSFRnameslist))
        
        
#     df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#     df.columns =   df.columns.map(str)
#     df.index.names = ['xyz']
#     df.columns.names = ['m']

#     # print(df) 

#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')




# #####################




# def initialize_CSV_files_MRRformationChannels(DCOname='BHBH'):

#     # namesEMlist=[]

#     iii=0
    

#     # CREATE PANDAS FILE 
#     nModels=26
#     BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#     NAMES = []
#     stringgg =  'MRR_FormationChannels'

#     for ind_l, BPSmodelName in enumerate(BPSnameslist):
#         # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#         # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

#         namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'
#         namez0_MRR = BPSmodelName + 'All MRR intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_MRR = BPSmodelName + 'All MRR observed (design LVK) [yr^{-1}]'

#         namez0_I = BPSmodelName + 'MRR channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_I = BPSmodelName + 'MRR channel I observed (design LVK) [yr^{-1}]'
#         namez0_II = BPSmodelName + 'MRR channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_II = BPSmodelName + 'MRR channel II observed (design LVK) [yr^{-1}]'
#         namez0_III = BPSmodelName + 'MRR channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_III = BPSmodelName + 'MRR channel III observed (design LVK) [yr^{-1}]'
#         namez0_IV = BPSmodelName + 'MRR channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_IV = BPSmodelName + 'MRR channel IV observed (design LVK) [yr^{-1}]'
#         namez0_V = BPSmodelName + 'MRR channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#         nameObs_V = BPSmodelName + 'MRR channel V observed (design LVK) [yr^{-1}]'

#         NAMES.append(namez0)
#         NAMES.append(nameObs)

#         NAMES.append(namez0_MRR)
#         NAMES.append(nameObs_MRR)

#         NAMES.append(namez0_I)
#         NAMES.append(nameObs_I)
#         NAMES.append(namez0_II)
#         NAMES.append(nameObs_II)
#         NAMES.append(namez0_III)
#         NAMES.append(nameObs_III)
#         NAMES.append(namez0_IV)
#         NAMES.append(nameObs_IV)
#         NAMES.append(namez0_V)
#         NAMES.append(nameObs_V)


        
        


#     datas=[]

#     for i in range(len(BPSnameslist)):
#         for ii in range(7):
#             datas.append(np.zeros_like(MSSFRnameslist))
#             datas.append(np.zeros_like(MSSFRnameslist))
        
        
#     df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#     df.columns =   df.columns.map(str)
#     df.index.names = ['xyz']
#     df.columns.names = ['m']

#     # print(df) 

#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')







# def writeToRatesFile_MRR_FormationChannels(BPSmodelName='Z', DCOtype='BHNS'):
#     """writes NS-BH rate to CSV file for all models"""

#     if DCOtype=='BHNS':
#         DCOname='BHNS'
#     elif DCOtype=='BBH':
#         DCOname='BHBH'
#     elif DCOtype=='BNS':
#         DCOname='NSNS'


#     # path for files 
#     path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
#     path_ = path_dir
#     path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
#     path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
#     # read in data 
#     fdata = h5.File(path)


#     # set optimistic true if that is the variation (H) 
#     OPTIMISTIC=False
#     if (BPSmodelName=='F') | (BPSmodelName=='K') :
#         OPTIMISTIC=True 

#     # obtain DCO masses
#     M1 = fdata['doubleCompactObjects']['M1'][...].squeeze()
#     M2 = fdata['doubleCompactObjects']['M2'][...].squeeze()
#     MBH, MNS = obtainM1BHandM2BHassymetric(M1, M2)


#     whichSN = fdata['supernovae']['whichStar'][...].squeeze() 
#     maskSNnot3 = (whichSN ==1) | (whichSN==2)


#     DCOseeds = fdata['doubleCompactObjects']['seed'][...].squeeze()
#     seedsSN = fdata['supernovae']['randomSeed'][...].squeeze()
    
#     u, indices = np.unique(seedsSN, return_index=True)

#     print(len(fdata['supernovae']['whichStar'][...].squeeze()))
#     print(seedsSN)


#     whichSN = fdata['supernovae']['whichStar'][...].squeeze()[indices] # get whichStar for first SN 

#     print(len(seedsSN), len(set(seedsSN)))
#     print(whichSN)
#     print(set(whichSN))
#     print(len(whichSN), len(M1), len(M2))

    
#     mask_temp = ((whichSN==2) & (M1>M2) )
#     mask_2 = ((whichSN==1) & (M1<M2))
#     print('nr of weird reversals = %s'%np.sum(mask_temp))
#     print('nr of normal reversals = %s'%np.sum(mask_2))




#     maskMRR = ((whichSN==2) & (M1>M2) ) | ((whichSN==1) & (M1<M2) ) 



#     del M1
#     del M2


#     print('at formation channel calc.')
#     fdata.close()

#     # get formation channel Seeds!
#     seedsPercentageClassic, seedsPercentageOnlyStableMT = returnSeedsPercentageClassicAndOnlyStableMT(pathCOMPASOutput=path,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)
#     seedsClassic, percentageClassic = seedsPercentageClassic
#     seedsOnlyStableMT, percentageOnlyStableMT = seedsPercentageOnlyStableMT



#     seedsDoubleCE, percentageDoubleCE = returnSeedsPercentageDoubleCoreCEE(pathCOMPASOutput=path,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)


#     seedsSingleCE, percentageSingleCE = returnSeedsPercentageSingleCoreCEE(pathCOMPASOutput=path,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1)



#     seedschannels = [seedsClassic, seedsOnlyStableMT, seedsSingleCE, seedsDoubleCE]

#     seedsOther, percentageOther = returnSeedsPercentageOther(pathCOMPASOutput=path,\
#                                     types=DCOtype,  withinHubbleTime=True, optimistic=OPTIMISTIC, \
#                                     binaryFraction=1, channelsSeedsList=seedschannels)

#     print('00')


#     dictChannelsBHNS = { 'classic':seedsClassic, \
#                         'immediate CE':seedsSingleCE,\
#                              'stable B no CEE':seedsOnlyStableMT, \
#                          r'double-core CE':seedsDoubleCE,  \
#                             'other':seedsOther\
#                            }


#     dictPercentages = { 'classic':percentageClassic, \
#                         'immediate CE':percentageSingleCE,\
#                              'stable B no CEE':percentageOnlyStableMT, \
#                          r'double-core CE':percentageDoubleCE,  \
#                             'other':percentageOther\
#                            }




#     # get intrinsic weights

#     fparam_intrinsic = 'weights_intrinsic'
#     # get detected weights

#     fparam_detected = 'weights_detected'


#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################



#     ####################################################
#     ######### ITERATE  OVER  MSSFR  MODELS #############
#     ####################################################
#     intrinsicRates = np.zeros(len(MSSFRnameslist))
#     detectedRates = np.zeros(len(MSSFRnameslist))
#     # namesEMlist = []
#     intrinsicRates_MRR = np.zeros(len(MSSFRnameslist))
#     detectedRates_MRR = np.zeros(len(MSSFRnameslist))

#     intrinsicRates_I = np.zeros(len(MSSFRnameslist))
#     detectedRates_I = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_II = np.zeros(len(MSSFRnameslist))
#     detectedRates_II = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_III = np.zeros(len(MSSFRnameslist))
#     detectedRates_III = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_IV = np.zeros(len(MSSFRnameslist))
#     detectedRates_IV = np.zeros(len(MSSFRnameslist))
#     intrinsicRates_V = np.zeros(len(MSSFRnameslist))
#     detectedRates_V = np.zeros(len(MSSFRnameslist))

#     fdata = h5.File(path)

#     DCOSeeds = fdata['doubleCompactObjects']['seed'][...].squeeze()

#     for ind_mssfr, mssfr in enumerate(MSSFRnameslist):


#         weightheader = 'w_' + mssfr
#         # print(ind_mssfr, weightheader)
#         w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
#         w_det = fdata[fparam_detected][weightheader][...].squeeze()



#         # ALL BBH RATES 
#         intrinsicRates[ind_mssfr] = np.sum(w_int)
#         detectedRates[ind_mssfr] = np.sum(w_det)
#         # MASS RATIO REVERSAL RATE ALL CHANNELS 

#         intrinsicRates_MRR[ind_mssfr] = np.sum(w_int[maskMRR])
#         detectedRates_MRR[ind_mssfr] = np.sum(w_det[maskMRR])

        
#         for nrC, Channel in enumerate(dictChannelsBHNSList):
                        
#     #             #Get the seeds that relate to sorted indices
#             seedsInterest = dictChannelsBHNS[Channel]
#             mask_temp = np.in1d(DCOSeeds, np.array(seedsInterest))
#             mask_C = (mask_temp==1) & (maskMRR==1)
#             if Channel=='classic':
#                 intrinsicRates_I[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_I[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='stable B no CEE':
#                 intrinsicRates_II[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_II[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='immediate CE':
#                 intrinsicRates_III[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_III[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='double-core CE':
#                 intrinsicRates_IV[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_IV[ind_mssfr] = np.sum(w_det[mask_C])
#             elif Channel=='other':
#                 intrinsicRates_V[ind_mssfr] = np.sum(w_int[mask_C])
#                 detectedRates_V[ind_mssfr] = np.sum(w_det[mask_C])




#     fdata.close()  

#     stringgg =  'MRR_FormationChannels'

#     df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv', index_col=0)
#     namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'
#     namez0_MRR = BPSmodelName + 'All MRR intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_MRR = BPSmodelName + 'All MRR observed (design LVK) [yr^{-1}]'

#     namez0_I = BPSmodelName + 'MRR channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_I = BPSmodelName + 'MRR channel I observed (design LVK) [yr^{-1}]'
#     namez0_II = BPSmodelName + 'MRR channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_II = BPSmodelName + 'MRR channel II observed (design LVK) [yr^{-1}]'
#     namez0_III = BPSmodelName + 'MRR channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_III = BPSmodelName + 'MRR channel III observed (design LVK) [yr^{-1}]'
#     namez0_IV = BPSmodelName + 'MRR channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_IV = BPSmodelName + 'MRR channel IV observed (design LVK) [yr^{-1}]'
#     namez0_V = BPSmodelName + 'MRR channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#     nameObs_V = BPSmodelName + 'MRR channel V observed (design LVK) [yr^{-1}]'



#     df[namez0] = intrinsicRates
#     df[nameObs] = detectedRates
#     df[namez0_MRR] = intrinsicRates_MRR
#     df[nameObs_MRR] = detectedRates_MRR

#     df[namez0_I] = intrinsicRates_I
#     df[nameObs_I] = detectedRates_I
#     df[namez0_II] = intrinsicRates_II
#     df[nameObs_II] = detectedRates_II
#     df[namez0_III] = intrinsicRates_III
#     df[nameObs_III] = detectedRates_III 
#     df[namez0_IV] = intrinsicRates_IV
#     df[nameObs_IV] = detectedRates_IV 
#     df[namez0_V] = intrinsicRates_V
#     df[nameObs_V] = detectedRates_V  


#     df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg  + '.csv')


#     fdata.close()













def writeToRatesFile_MRR_Spins(BPSmodelName='Z', DCOtype='BHNS', spin_threshold=0.05):
    """writes NS-BH rate to CSV file for all models"""

    if DCOtype=='BHNS':
        DCOname='BHNS'
    elif DCOtype=='BBH':
        DCOname='BHBH'
    elif DCOtype=='BNS':
        DCOname='NSNS'


    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
    # read in data 
    fdata = h5.File(path, 'r')
    fDCO      = fdata['doubleCompactObjects'] # hdf5 file with the DCO information
    fSN       = fdata['supernovae']  # hdf5 file with the SN information
    #
    M1 = fDCO['M1'][...].squeeze()   # Compact object mass [Msun] of the initially more massive star
    M2 = fDCO['M2'][...].squeeze()  # Compact object mass [Msun] of the initially less massive star

    print('using indices')
    seedsDCO = fDCO['seed'][...].squeeze()  # get the seeds in the DCO file 
    seedsSN = fSN['randomSeed'][...].squeeze()    # get the seeds in the SN file 
    indices = np.sort(np.unique(seedsSN[1::2], return_index=True)[1])
    maskSNdco = np.in1d(seedsSN,  seedsDCO) # mask in the SNe files the SNe that correspond to our DCO
    whichSN = fSN['whichStar'][...].squeeze()[maskSNdco]  # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova     
    whichSN2 = whichSN[1::2][indices]

    # either SN2 = primary (1) and M1 is > M2, or SN2 = secondary & M1 < M2 
    # this takes into account (first term) rejuvenation 
    MRR_mask = ((whichSN2==1) & (M1>M2) ) | ((whichSN2==2) & (M1<M2)) 

    del M1
    del M2
    del whichSN2
    del whichSN
    del maskSNdco
    del indices
    del seedsSN
    del seedsDCO
    del fDCO
    del fSN

    fdata.close()



    spin = COspin(data_path=path, state='he_depletion')  # set class 
    spin.setCOMPASData() # reads in the COMPAS DCO parameters 
    spinMZAMS1, spinMZAMS2  = spin.BaveraSpin()


    spinLVKM1, spinLVKM2 = np.zeros_like(spinMZAMS1), np.zeros_like(spinMZAMS1)
    spinLVKM1[MRR_mask] = spinMZAMS2[MRR_mask]  # MRR so M1 comes from M2ZAMS, we assign it spin from M2ZAMS
    spinLVKM1[~MRR_mask] = spinMZAMS1[~MRR_mask]  # no MRR so M1 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[MRR_mask] = spinMZAMS1[MRR_mask]   # MRR so M2 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[~MRR_mask] = spinMZAMS2[~MRR_mask]   # no MRR so M2 comes from M2ZAMS, we assign it spin from M2ZAMS     

    # spin_threshold = 0.05 # definition of "spinning BH"
    mask_LVKM1_spinning = (spinLVKM1 > spin_threshold ) 
    mask_LVKM2_spinning = (spinLVKM2 > spin_threshold ) # definition of "spinning BH"
    mask_anySpin = (spinLVKM1 > spin_threshold ) | (spinLVKM2 > spin_threshold )


    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'




    ####################################################
    ######### ITERATE  OVER  MSSFR  MODELS #############
    ####################################################
    intrinsicRates = np.zeros(len(MSSFRnameslist))
    detectedRates = np.zeros(len(MSSFRnameslist))
    # namesEMlist = []
    intrinsicRates_spinning = np.zeros(len(MSSFRnameslist))
    detectedRates_spinning = np.zeros(len(MSSFRnameslist))

    intrinsicRates_spinning_LVKM1 = np.zeros(len(MSSFRnameslist))
    detectedRates_spinning_LVKM1 = np.zeros(len(MSSFRnameslist))
    intrinsicRates_spinning_LVKM2 = np.zeros(len(MSSFRnameslist))
    detectedRates_spinning_LVKM2 = np.zeros(len(MSSFRnameslist))


    fdata = h5.File(path)

    DCOSeeds = fdata['doubleCompactObjects']['seed'][...].squeeze()

    for ind_mssfr, mssfr in enumerate(MSSFRnameslist):

        weightheader = 'w_' + mssfr
        # print(ind_mssfr, weightheader)
        w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
        w_det = fdata[fparam_detected][weightheader][...].squeeze()

        # ALL BBH RATES 
        intrinsicRates[ind_mssfr] = np.sum(w_int)
        detectedRates[ind_mssfr] = np.sum(w_det)
        
        # MASS RATIO REVERSAL RATE ALL CHANNELS 
        intrinsicRates_spinning[ind_mssfr] = np.sum(w_int[mask_anySpin])
        detectedRates_spinning[ind_mssfr] = np.sum(w_det[mask_anySpin])

        intrinsicRates_spinning_LVKM1[ind_mssfr] = np.sum(w_int[mask_LVKM1_spinning])
        detectedRates_spinning_LVKM1[ind_mssfr] = np.sum(w_det[mask_LVKM1_spinning])       
        intrinsicRates_spinning_LVKM2[ind_mssfr] = np.sum(w_int[mask_LVKM2_spinning])
        detectedRates_spinning_LVKM2[ind_mssfr] = np.sum(w_det[mask_LVKM2_spinning])     



    fdata.close()  

    stringgg =  'MRR_Spins'
    if spin_threshold!=0.05:
        stringgg =  'MRR_Spins' +  '_threshold_%s_'%spin_threshold

    df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv', index_col=0)
    namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
    nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'
    namez0_spinAny = BPSmodelName + 'All: at least one spin >%s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
    nameObs_spinAny = BPSmodelName + 'All: at least one spin >%s observed (design LVK) [yr^{-1}]'%spin_threshold

    namez0_spinning_LVKM1 = BPSmodelName + 'LVKM1 spinning > %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
    nameObs_spinning_LVKM1 = BPSmodelName + 'LVKM1 spinning > %s observed (design LVK) [yr^{-1}]'%spin_threshold
    namez0_spinning_LVKM2= BPSmodelName + 'LVKM2 spinning > %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
    nameObs_spinning_LVKM2 = BPSmodelName + 'LVKM2 spinning > %s observed (design LVK) [yr^{-1}]'%spin_threshold



    df[namez0]                      = intrinsicRates
    df[nameObs]                     = detectedRates
    df[namez0_spinAny]              = intrinsicRates_spinning
    df[nameObs_spinAny]             = detectedRates_spinning

    df[namez0_spinning_LVKM1]       = intrinsicRates_spinning_LVKM1
    df[nameObs_spinning_LVKM1]      = detectedRates_spinning_LVKM1
    df[namez0_spinning_LVKM2]       = intrinsicRates_spinning_LVKM2
    df[nameObs_spinning_LVKM2]      = detectedRates_spinning_LVKM2
 


    df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg  + '.csv')


    fdata.close()







def initialize_CSV_files_MRRspins(DCOname='BHBH', spin_threshold=0.05):

    # namesEMlist=[]

    iii=0
    

    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'MRR_Spins'
    if spin_threshold!=0.05:
        stringgg =  'MRR_Spins' +  '_threshold_%s_'%spin_threshold

    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
        # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

        namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'
        namez0_spinAny = BPSmodelName + 'All: at least one spin >%s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_spinAny = BPSmodelName + 'All: at least one spin >%s observed (design LVK) [yr^{-1}]'%spin_threshold

        namez0_spinning_LVKM1 = BPSmodelName + 'LVKM1 spinning > %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_spinning_LVKM1 = BPSmodelName + 'LVKM1 spinning > %s observed (design LVK) [yr^{-1}]'%spin_threshold
        namez0_spinning_LVKM2= BPSmodelName + 'LVKM2 spinning > %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_spinning_LVKM2 = BPSmodelName + 'LVKM2 spinning > %s observed (design LVK) [yr^{-1}]'%spin_threshold

        NAMES.append(namez0)
        NAMES.append(nameObs)

        NAMES.append(namez0_spinAny)
        NAMES.append(nameObs_spinAny)

        NAMES.append(namez0_spinning_LVKM1)
        NAMES.append(nameObs_spinning_LVKM1)
        NAMES.append(namez0_spinning_LVKM2)
        NAMES.append(nameObs_spinning_LVKM2)




    datas=[]

    for i in range(len(BPSnameslist)):
        for ii in range(4):
            datas.append(np.zeros_like(MSSFRnameslist))
            datas.append(np.zeros_like(MSSFRnameslist))
        
        
    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']

    # print(df) 

    df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')







def initialize_CSV_files_MRR_nonMRR_ratio(DCOname='BHBH', spin_threshold=0.05, GWname='GW150629'):

    # namesEMlist=[]

    iii=0
    

    # CREATE PANDAS FILE 
    nModels=26
    BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    NAMES = []
    stringgg =  'MRR_nonMRR_ratio_' + GWname 


    for ind_l, BPSmodelName in enumerate(BPSnameslist):
        # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
        # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

        namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'
        
        namez0_MRR = BPSmodelName + 'MRR: matches GW contour intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_MRR = BPSmodelName + 'MRR: matches GW contour observed (design LVK) [yr^{-1}]'
        namez0_nonMRR = BPSmodelName + 'nonMRR: matches GW contour intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_nonMRR = BPSmodelName + 'nonMRR: matches GW contour observed (design LVK) [yr^{-1}]'
        
        namez0_MRR_spin= BPSmodelName + 'MRR: matches GW contour spin2 < %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_MRR_spin = BPSmodelName + 'MRR: matches GW contour spin2 < %s observed (design LVK) [yr^{-1}]'%spin_threshold
        namez0_nonMRR_spin= BPSmodelName + 'nonMRR: matches GW contour spin1 < %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_nonMRR_spin = BPSmodelName + 'nonMRR: matches GW contour spin1 < %s observed (design LVK) [yr^{-1}]'%spin_threshold



        NAMES.append(namez0)
        NAMES.append(nameObs)

        NAMES.append(namez0_MRR)
        NAMES.append(nameObs_MRR)
        NAMES.append(namez0_nonMRR)
        NAMES.append(nameObs_nonMRR)

        NAMES.append(namez0_MRR_spin)
        NAMES.append(nameObs_MRR_spin)
        NAMES.append(namez0_nonMRR_spin)
        NAMES.append(nameObs_nonMRR_spin)




    datas=[]

    for i in range(len(BPSnameslist)):
        for ii in range(5):
            datas.append(np.zeros_like(MSSFRnameslist))
            datas.append(np.zeros_like(MSSFRnameslist))
        
        
    df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    df.columns =   df.columns.map(str)
    df.index.names = ['xyz']
    df.columns.names = ['m']

    # print(df) 

    df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





def GW_credible_intervals(GW_name, mode):

    # GW_list = ['GW151226','GW170729', 'GW190517_055101', 'GW190412','GW191109_010717'  ,'GW191103_012549', 'GW191126_115259']

    
    dfCSVname= '/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/GWTC_posterior_samples/' 
    dfCSVname_ = dfCSVname + 'posteriorSamples_' + GW_name  + '.csv' 


    df = pd.read_csv(dfCSVname_, index_col=0, skiprows=[1])
    mass_ratio = df['M2']/df['M1']
    total_mass = df['M2'] + df['M1']
    spin1 = df['spin1'] 
    spin2 = df['spin2'] 
    
    chi_eff = df['chi_eff']
    chirp_mass = chirpmass(df['M2'], df['M1'])
    y_quantiles = [0.05, 0.5, 0.95]

    if mode=='normal':
    
        total_mass_CI = weighted_quantile(values=total_mass, quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass, quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio, quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1, quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2, quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff, quantiles=y_quantiles)    

    elif mode=='spin1_is_zero':
        mask_spin = (abs(spin1)<0.05) & (spin2>0.05) # non MRR, spin 2 is the spinning one, we want spin1 to be zero 
        total_mass_CI = weighted_quantile(values=total_mass[mask_spin], quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'][mask_spin], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'][mask_spin], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass[mask_spin], quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio[mask_spin], quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1[mask_spin], quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2[mask_spin], quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff[mask_spin], quantiles=y_quantiles)   

    elif  mode=='spin2_is_zero':
        mask_spin = (abs(spin2)<0.05) & (spin1>0.05)  # MRR spin1 is the spinning one, we want the other one to be zero
        total_mass_CI = weighted_quantile(values=total_mass[mask_spin], quantiles=y_quantiles)
        mass1_CI = weighted_quantile(values=df['M1'][mask_spin], quantiles=y_quantiles)
        mass2_CI = weighted_quantile(values=df['M2'][mask_spin], quantiles=y_quantiles)
        chirp_mass_CI = weighted_quantile(values=chirp_mass[mask_spin], quantiles=y_quantiles)
        mass_ratio_CI = weighted_quantile(values=mass_ratio[mask_spin], quantiles=y_quantiles)
        spin1_CI = weighted_quantile(values=spin1[mask_spin], quantiles=y_quantiles)    
        spin2_CI = weighted_quantile(values=spin2[mask_spin], quantiles=y_quantiles)    
        chi_eff_quantiles = weighted_quantile(values=chi_eff[mask_spin], quantiles=y_quantiles)   

    return total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_quantiles
    



def writeToRatesFile_MRR_nonMRR_ratio(BPSmodelName='Z', DCOtype='BHNS', spin_threshold=0.05, GWnameList=['GW150629']):
    """writes NS-BH rate to CSV file for all models"""


    if DCOtype=='BHNS':
        DCOname='BHNS'
    elif DCOtype=='BBH':
        DCOname='BHBH'
    elif DCOtype=='BNS':
        DCOname='NSNS'


    # path for files 
    path_dir = '/Volumes/Andromeda2/DATA/AllDCO_bugfix/'
    path_ = path_dir
    path_ = path_ + alphabetDirDict[BPSmodelName] +'/'
    path  = path_ + 'COMPASCompactOutput_'+ DCOtype +'_' + BPSmodelName + '.h5'
            
    # read in data 
    fdata = h5.File(path, 'r')
    fDCO      = fdata['doubleCompactObjects'] # hdf5 file with the DCO information
    fSN       = fdata['supernovae']  # hdf5 file with the SN information
    #
    M1 = fDCO['M1'][...].squeeze()   # Compact object mass [Msun] of the initially more massive star
    M2 = fDCO['M2'][...].squeeze()  # Compact object mass [Msun] of the initially less massive star

    print('using indices')
    seedsDCO = fDCO['seed'][...].squeeze()  # get the seeds in the DCO file 
    seedsSN = fSN['randomSeed'][...].squeeze()    # get the seeds in the SN file 
    indices = np.sort(np.unique(seedsSN[1::2], return_index=True)[1])
    maskSNdco = np.in1d(seedsSN,  seedsDCO) # mask in the SNe files the SNe that correspond to our DCO
    whichSN = fSN['whichStar'][...].squeeze()[maskSNdco]  # this is 1 if the initially primary star goes SN and 2 if the secondary goes supernova     
    whichSN2 = whichSN[1::2][indices]

    # either SN2 = primary (1) and M1 is > M2, or SN2 = secondary & M1 < M2 
    # this takes into account (first term) rejuvenation 
    MRR_mask = ((whichSN2==1) & (M1>M2) ) | ((whichSN2==2) & (M1<M2)) 

    M1LVK, M2LVK = obtainM1BHandM2BHassymetric(M1, M2)
    chirp_mass = chirpmass(M1LVK, M2LVK)
    mass_ratio_LVK =  M2LVK/M1LVK

    del M1
    del M2
    del whichSN2
    del whichSN
    del maskSNdco
    del indices
    del seedsSN
    del seedsDCO
    del fDCO
    del fSN

    fdata.close()



    spin = COspin(data_path=path, state='he_depletion')  # set class 
    spin.setCOMPASData() # reads in the COMPAS DCO parameters 
    spinMZAMS1, spinMZAMS2  = spin.BaveraSpin()


    spinLVKM1, spinLVKM2 = np.zeros_like(spinMZAMS1), np.zeros_like(spinMZAMS1)
    spinLVKM1[MRR_mask] = spinMZAMS2[MRR_mask]  # MRR so M1 comes from M2ZAMS, we assign it spin from M2ZAMS
    spinLVKM1[~MRR_mask] = spinMZAMS1[~MRR_mask]  # no MRR so M1 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[MRR_mask] = spinMZAMS1[MRR_mask]   # MRR so M2 comes from M1ZAMS, we assign it spin from M1ZAMS
    spinLVKM2[~MRR_mask] = spinMZAMS2[~MRR_mask]   # no MRR so M2 comes from M2ZAMS, we assign it spin from M2ZAMS     

    chi_eff = ((spinLVKM1*M1LVK) + (spinLVKM2*M2LVK)) / (M1LVK + M2LVK)


    # spin_threshold = 0.05 # definition of "spinning BH"





    # get intrinsic weights
    fparam_intrinsic = 'weights_intrinsic'
    # get detected weights
    fparam_detected = 'weights_detected'







    fdata = h5.File(path)

    DCOSeeds = fdata['doubleCompactObjects']['seed'][...].squeeze()


    for ind_GW, GWname in enumerate(GWnameList):  
        print(GWname)


        ####################################################
        ######### ITERATE  OVER  MSSFR  MODELS #############
        ####################################################
        intrinsicRates = np.zeros(len(MSSFRnameslist))
        detectedRates = np.zeros(len(MSSFRnameslist))
        # namesEMlist = []

        intrinsicRates_MRR = np.zeros(len(MSSFRnameslist))
        detectedRates_MRR = np.zeros(len(MSSFRnameslist))
        intrinsicRates_nonMRR = np.zeros(len(MSSFRnameslist))
        detectedRates_nonMRR = np.zeros(len(MSSFRnameslist))

        intrinsicRates_MRR_spin = np.zeros(len(MSSFRnameslist))
        detectedRates_MRR_spin = np.zeros(len(MSSFRnameslist))
        intrinsicRates_nonMRR_spin = np.zeros(len(MSSFRnameslist))
        detectedRates_nonMRR_spin = np.zeros(len(MSSFRnameslist))




        total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_CI = GW_credible_intervals(GWname, mode='normal')
        mask_GW =  (total_mass_CI[0]<=(M1LVK+M2LVK))  & ((M1LVK+M2LVK)<=total_mass_CI[2])  &  (chirp_mass_CI[0]<=chirp_mass) & (chirp_mass<=chirp_mass_CI[2])  & (mass_ratio_CI[0]<= mass_ratio_LVK) & (mass_ratio_LVK<=mass_ratio_CI[2]) & (chi_eff_CI[0]<=chi_eff) & (chi_eff<=chi_eff_CI[2]) & (mass1_CI[0]<= M1LVK) & (M1LVK<=mass1_CI[2]) & (mass2_CI[0]<= M2LVK) & (M2LVK<=mass2_CI[2])
        mask_GW_MRR = (mask_GW==1) & (MRR_mask==1)
        mask_GW_nonMRR = (mask_GW==1) & (MRR_mask==0)


        total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_CI = GW_credible_intervals(GWname, mode='spin1_is_zero')
        mask_spin1_zero =  (total_mass_CI[0]<=(M1LVK+M2LVK))  & ((M1LVK+M2LVK)<=total_mass_CI[2])  &  (chirp_mass_CI[0]<=chirp_mass) & (chirp_mass<=chirp_mass_CI[2]) & (mass_ratio_CI[0]<= mass_ratio_LVK) & (mass_ratio_LVK<=mass_ratio_CI[2]) & (chi_eff_CI[0]<=chi_eff) & (chi_eff<=chi_eff_CI[2]) & (mass1_CI[0]<= M1LVK) & (M1LVK<=mass1_CI[2]) & (mass2_CI[0]<= M2LVK) & (M2LVK<=mass2_CI[2])# &  (spin2_CI[0]<=spinLVKM2 ) & (spinLVKM2<=spin2_CI[2]) 
        total_mass_CI, mass1_CI, mass2_CI, chirp_mass_CI, mass_ratio_CI, spin1_CI, spin2_CI, chi_eff_CI = GW_credible_intervals(GWname, mode='spin2_is_zero')
        mask_spin2_zero =  (total_mass_CI[0]<=(M1LVK+M2LVK))  & ((M1LVK+M2LVK)<=total_mass_CI[2])  &  (chirp_mass_CI[0]<=chirp_mass) & (chirp_mass<=chirp_mass_CI[2]) & (mass_ratio_CI[0]<= mass_ratio_LVK) & (mass_ratio_LVK<=mass_ratio_CI[2]) & (chi_eff_CI[0]<=chi_eff) & (chi_eff<=chi_eff_CI[2]) & (mass1_CI[0]<= M1LVK) & (M1LVK<=mass1_CI[2]) & (mass2_CI[0]<= M2LVK) & (M2LVK<=mass2_CI[2]) # &  (spin1_CI[0]<=spinLVKM1 ) & (spinLVKM1<=spin1_CI[2]) 


        mask_spin2zero_MRR = (mask_spin2_zero==1) & (MRR_mask==1)
        mask_spin1zero_nonMRR = (mask_spin1_zero==1) & (MRR_mask==0)



        for ind_mssfr, mssfr in enumerate(MSSFRnameslist):

            weightheader = 'w_' + mssfr
            # print(ind_mssfr, weightheader)
            w_int = fdata[fparam_intrinsic][weightheader][...].squeeze()
            w_det = fdata[fparam_detected][weightheader][...].squeeze()

            # ALL BBH RATES 
            intrinsicRates[ind_mssfr] = np.sum(w_int)
            detectedRates[ind_mssfr] = np.sum(w_det)
            
            # MASS RATIO REVERSAL RATE ALL CHANNELS 
            intrinsicRates_MRR[ind_mssfr] = np.sum(w_int[mask_GW_MRR])
            detectedRates_MRR[ind_mssfr] = np.sum(w_det[mask_GW_MRR])
            intrinsicRates_nonMRR[ind_mssfr] = np.sum(w_int[mask_GW_nonMRR])
            detectedRates_nonMRR[ind_mssfr] = np.sum(w_det[mask_GW_nonMRR])

            intrinsicRates_MRR_spin[ind_mssfr] = np.sum(w_int[mask_spin2zero_MRR])
            detectedRates_MRR_spin[ind_mssfr] = np.sum(w_det[mask_spin2zero_MRR])       
            intrinsicRates_nonMRR_spin[ind_mssfr] = np.sum(w_int[mask_spin1zero_nonMRR])
            detectedRates_nonMRR_spin[ind_mssfr] = np.sum(w_det[mask_spin1zero_nonMRR])     



    

        stringgg =  'MRR_nonMRR_ratio_' + GWname 

        df = pd.read_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv', index_col=0)
        namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'

        namez0_MRR = BPSmodelName + 'MRR: matches GW contour intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_MRR = BPSmodelName + 'MRR: matches GW contour observed (design LVK) [yr^{-1}]'
        namez0_nonMRR = BPSmodelName + 'nonMRR: matches GW contour intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
        nameObs_nonMRR = BPSmodelName + 'nonMRR: matches GW contour observed (design LVK) [yr^{-1}]'
        
        namez0_MRR_spin= BPSmodelName + 'MRR: matches GW contour spin2 < %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_MRR_spin = BPSmodelName + 'MRR: matches GW contour spin2 < %s observed (design LVK) [yr^{-1}]'%spin_threshold
        namez0_nonMRR_spin= BPSmodelName + 'nonMRR: matches GW contour spin1 < %s intrinsic (z=0) [Gpc^{-3} yr^{-1}]'%spin_threshold
        nameObs_nonMRR_spin = BPSmodelName + 'nonMRR: matches GW contour spin1 < %s observed (design LVK) [yr^{-1}]'%spin_threshold



        df[namez0]                      = intrinsicRates
        df[nameObs]                     = detectedRates

        df[namez0_MRR]              = intrinsicRates_MRR
        df[nameObs_MRR]             = detectedRates_MRR
        df[namez0_nonMRR]              = intrinsicRates_nonMRR
        df[nameObs_nonMRR]             = detectedRates_nonMRR

        df[namez0_MRR_spin]          = intrinsicRates_MRR_spin
        df[nameObs_MRR_spin]         = detectedRates_MRR_spin
        df[namez0_nonMRR_spin]       = intrinsicRates_nonMRR_spin
        df[nameObs_nonMRR_spin]      = detectedRates_nonMRR_spin
     


        df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')







    fdata.close()










GW_MRR_list = ['GW151226','GW170729', 'GW190517_055101', 'GW190412','GW191109_010717'  ,'GW191103_012549', 'GW191126_115259']


GWTC1 = ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729',  'GW170809', 'GW170814', 'GW170818', 'GW170823'] # 'GW170817', 
#     GWTC2 = ['GW190408_181802','GW190412','GW190413_052954','GW190413_134308','GW190421_213856',\
#     'GW190424_180648','GW190503_185404','GW190512_180714',\
#     'GW190513_205428','GW190514_065416','GW190517_055101','GW190519_153544','GW190521_074359',\
#     'GW190521','GW190527_092055','GW190602_175927','GW190620_030421','GW190630_185205','GW190701_203306',\
#     'GW190706_222641','GW190707_093326','GW190708_232457','GW190720_000836',\
#     'GW190727_060333','GW190728_064510','GW190731_140936','GW190803_022701','GW190828_063405',\
#     'GW190828_065509','GW190910_112807','GW190915_235702','GW190924_021846','GW190929_012149',\
#      'GW190930_133541', 'GW190425', 'GW190814', 'GW190426_152155']
    
GWTC2 = [ 'GW190408_181802', 'GW190412', 'GW190413_052954', 'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', 'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', 'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514', 'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701', 'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', 'GW190924_021846', 'GW190929_012149', 'GW190930_133541'] #, 'GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092255', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']
#     GWTC3 = ['GW191103_012549', 'GW191126_115259','GW191109_010717' ]    
GWTC3 = ['GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092254', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']             
    


BBHlist_GWTC_123 = ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729', 'GW170809', \
                    'GW170814', 'GW170818', 'GW170823', 'GW190408_181802', 'GW190412', 'GW190413_052954',\
                    'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', \
                    'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', \
                    'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', \
                    'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514',\
                    'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701',\
                    'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', \
                    'GW190924_021846', 'GW190929_012149', 'GW190930_133541', 'GW191103_012549', 'GW191105_143521',\
                    'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', \
                    'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', \
                    'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', \
                    'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092254', 'GW200216_220804', \
                    'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421',\
                    'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']



# ['GW150914', 'GW151012', 'GW151226', 'GW170104', 'GW170608', 'GW170729', 'GW170809', 'GW170814', 'GW170818', 'GW170823', 'GW190408_181802', 'GW190412', 'GW190413_052954', 'GW190413_134308', 'GW190421_213856', 'GW190424_180648', 'GW190503_185404', 'GW190512_180714', 'GW190513_205428', 'GW190514_065416', 'GW190517_055101', 'GW190519_153544', 'GW190521', 'GW190521_074359', 'GW190527_092055', 'GW190602_175927', 'GW190620_030421', 'GW190630_185205', 'GW190701_203306', 'GW190706_222641', 'GW190707_093326', 'GW190708_232457', 'GW190719_215514', 'GW190720_000836', 'GW190727_060333', 'GW190728_064510', 'GW190731_140936', 'GW190803_022701', 'GW190814', 'GW190828_063405', 'GW190828_065509', 'GW190910_112807', 'GW190915_235702', 'GW190924_021846', 'GW190929_012149', 'GW190930_133541', 'GW191103_012549', 'GW191105_143521', 'GW191109_010717', 'GW191113_071753', 'GW191126_115259', 'GW191127_050227', 'GW191129_134029', 'GW191204_110529', 'GW191204_171526', 'GW191215_223052', 'GW191216_213338', 'GW191222_033537', 'GW191230_180458', 'GW200112_155838', 'GW200128_022011', 'GW200129_065458', 'GW200202_154313', 'GW200208_130117', 'GW200208_222617', 'GW200209_085452', 'GW200210_092255', 'GW200216_220804', 'GW200219_094415', 'GW200220_061928', 'GW200220_124850', 'GW200224_222234', 'GW200225_060421', 'GW200302_015811', 'GW200306_093714', 'GW200308_173609', 'GW200311_115853', 'GW200316_215756', 'GW200322_091133']







# INITIALIZE
INITIALIZE_GENERAL = False # True #False #True #False#True #False
INITIALIZE_lightestBHfirst = False #True
INITIALIZE_MRR_FormationChannels = False
INITIALIZE_MRR_Spins = False
INITIALIZE_runMRR_nonMRR_ratio = False
# 
spin_threshold=0.05


if INITIALIZE_GENERAL==True:
    # initialize_CSV_files_general(DCOname='BHNS')
    initialize_CSV_files_general(DCOname='BHBH')
    # initialize_CSV_files_general(DCOname='NSNS')


if INITIALIZE_runMRR_nonMRR_ratio==True:
    for GWname in BBHlist_GWTC_123:
        print('initializing ', GWname)
        # initialize_CSV_files_general(DCOname='BHNS')
        initialize_CSV_files_MRR_nonMRR_ratio(DCOname='BHBH', spin_threshold=spin_threshold, GWname=GWname)



if INITIALIZE_lightestBHfirst==True:
    # initialize_CSV_files_general(DCOname='BHNS')
    initialize_CSV_files_lightestBHfirst(DCOname='BHBH')
    # initialize_CSV_files_general(DCOname='NSNS')

if INITIALIZE_MRR_FormationChannels==True:
    # initialize MRR formation channel file 
    initialize_CSV_files_MRRformationChannels(DCOname='BHBH')


if INITIALIZE_MRR_Spins==True:
    # initialize MRR formation channel file 
    initialize_CSV_files_MRRspins(DCOname='BHBH', spin_threshold=spin_threshold)

#### RUN different simulation summaries : 
runMejecta = False 
runFormationChannels =False 
runNSBH = False
runGeneralBHNS = False
runGeneralBHBH = False
runGeneralNSNS = False


runLightestFormsFirst=False
runMRR_FormationChannels = False
runMRR_Spins = False
runMRR_nonMRR_ratio = True



if runMRR_nonMRR_ratio==True:
    # for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
    for BPS in ['C']:
    # for BPS in ['N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            # for GWname in GW_MRR_list:
            
            # print('at GW =', GWname)
            writeToRatesFile_MRR_nonMRR_ratio(BPSmodelName=BPS, DCOtype=DCOtype, spin_threshold=spin_threshold, GWnameList=BBHlist_GWTC_123)

            print('done with ', BPS)
            print('----------------')
            print()





if runMRR_Spins==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            writeToRatesFile_MRR_Spins(BPSmodelName=BPS, DCOtype=DCOtype, spin_threshold=spin_threshold)



if runMRR_FormationChannels==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_MRR_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)




# if runMRR_FormationChannels==True:
#   for BPS in ['F', 'H', 'K']:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_MRR_FormationChannels(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)



if runLightestFormsFirst==True:
    for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_lightestFormsFirst(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)



if runGeneralBHBH ==True:
    for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
        print(BPS)
        for DCOtype in ['BBH']:
            print('at DCOtype =', DCOtype)
            writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
            print('done with ', BPS)

    # # INITIALIZE FILE 
    # namesEMlist=[]

    # DCOname ='BHBH'
    # iii=0
    

    # # CREATE PANDAS FILE 
    # nModels=26
    # BPSnameslist = list(string.ascii_uppercase)[0:nModels]

    # NAMES = []
    # stringgg = 'lightestFormsFirst'

    # for ind_l, L in enumerate(BPSnameslist):
    #   str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
    #   str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
    #   NAMES.append(str_z0)
    #   NAMES.append(str_obs)
        
        


    # datas=[]

    # for i in range(len(BPSnameslist)):
    #   datas.append(np.zeros_like(MSSFRnameslist))
    #   datas.append(np.zeros_like(MSSFRnameslist))
        
        
    # df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
    # df.columns =   df.columns.map(str)
    # df.index.names = ['.x.y.z']
    # df.columns.names = ['m']

    # # print(df) 

    # df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/MRR_Project/dataFiles/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')

    # # run calculation 
    # for BPS in   ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
    #   print(BPS)
    #   for DCOtype in ['BBH']:
    #       print('at DCOtype =', DCOtype)
    #       writeToRatesFile_lightestFormsFirst(BPSmodelName=BPS, DCOtype=DCOtype)
    #       print('done with ', BPS)














# if runMejecta ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_Mejecta(BPSmodelName=BPS)
#           print('done with ', BPS)




# if runFormationChannels ==True:
#   # for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O' ]:
#   #   print(BPS)
#   #   for DCOtype in ['BNS']:
#   #       print('at DCOtype =', DCOtype)
#   #       writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype='BNS')
#   #       print('done with ', BPS)


#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_FormationChannels(BPSmodelName=BPS, DCOtype='BBH')
#           print('done with ', BPS)

# if runNSBH ==True:
#   for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' ]:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_NSBH(BPSmodelName=BPS)
#           print('done with ', BPS)




# if runGeneralBHNS ==True:
#   for BPS in ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q' , 'R', 'S', 'T']:
#       print(BPS)
#       for DCOtype in ['BHNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)


# if runGeneralNSNS ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
#       print(BPS)
#       for DCOtype in ['BBH']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)

# if runGeneralBHBH ==True:
#   for BPS in  ['A','B',  'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T' ]:
#       print(BPS)
#       for DCOtype in ['BNS']:
#           print('at DCOtype =', DCOtype)
#           writeToRatesFile_GENERAL(BPSmodelName=BPS, DCOtype=DCOtype)
#           print('done with ', BPS)














# Models to RUN 

# May 20: I am updating my data with the AllDCO focused runs :-) 

# this is an overwrite with better data (old ones are in BHNS copy)


# for DCOtype in ['BHNS', 'BBH', 'BNS']:
#   print('at DCOtype =', DCOtype)
#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
#   modelname = 'A'
#   writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
#   modelname = 'B'
#   writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)


#   pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
#   modelname = 'G'
#   writeToRatesFile(modelname=modelname, pa thCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# INITIALIZE_FormationChannels = False
# INITIALIZE_NSBH= False #False#True
# INITIALIZE=False #False #True 
# INITIALIZE_GW190814 = False
# INITIALIZE_EM =False

# # INITIALIZE_NSBH= True #False#True
# # INITIALIZE=True #False #True 
# # INITIALIZE_GENERAL = True #False#True #False


# # ['.0.0.0', '.1.1.1', '.2.1.1', '.3.1.1', '.1.1.2', '.2.1.2', '.3.1.2', '.1.1.3', '.2.1.3', '.3.1.3', '.1.2.1', '.2.2.1', '.3.2.1', '.1.2.2', '.2.2.2', '.3.2.2', '.1.2.3', '.2.2.3', '.3.2.3', '.1.3.1', '.2.3.1', '.3.3.1', '.1.3.2', '.2.3.2', '.3.3.2', '.1.3.3', '.2.3.3', '.3.3.3']


# if INITIALIZE_FormationChannels==True:

#   # namesEMlist=[]

#   DCOname ='NSNS' 
#   iii=0
    

#   # CREATE PANDAS FILE 
#   nModels=26
#   BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#   NAMES = []
#   stringgg =  'AllDCOsimulation_formation_channels'

#   for ind_l, BPSmodelName in enumerate(BPSnameslist):
#       # str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#       # str_obs = str(L + ' observed (design LVK) [yr^{-1}]')

#       namez0 = BPSmodelName + 'All intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs = BPSmodelName + 'All observed (design LVK) [yr^{-1}]'

#       namez0_I = BPSmodelName + 'channel I intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_I = BPSmodelName + 'channel I observed (design LVK) [yr^{-1}]'
#       namez0_II = BPSmodelName + 'channel II intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_II = BPSmodelName + 'channel II observed (design LVK) [yr^{-1}]'
#       namez0_III = BPSmodelName + 'channel III intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_III = BPSmodelName + 'channel III observed (design LVK) [yr^{-1}]'
#       namez0_IV = BPSmodelName + 'channel IV intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_IV = BPSmodelName + 'channel IV observed (design LVK) [yr^{-1}]'
#       namez0_V = BPSmodelName + 'channel V intrinsic (z=0) [Gpc^{-3} yr^{-1}]'
#       nameObs_V = BPSmodelName + 'channel V observed (design LVK) [yr^{-1}]'

#       NAMES.append(namez0)
#       NAMES.append(nameObs)

#       NAMES.append(namez0_I)
#       NAMES.append(nameObs_I)
#       NAMES.append(namez0_II)
#       NAMES.append(nameObs_II)
#       NAMES.append(namez0_III)
#       NAMES.append(nameObs_III)
#       NAMES.append(namez0_IV)
#       NAMES.append(nameObs_IV)
#       NAMES.append(namez0_V)
#       NAMES.append(nameObs_V)


        
        


#   datas=[]

#   for i in range(len(BPSnameslist)):
#       for ii in range(6):
#           datas.append(np.zeros_like(MSSFRnameslist))
#           datas.append(np.zeros_like(MSSFRnameslist))
        
        
#   df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#   df.columns =   df.columns.map(str)
#   df.index.names = ['xyz']
#   df.columns.names = ['m']

#   # print(df) 

#   df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





# if INITIALIZE_GW190814==True:

#   for dcotype in ['NSNS', 'BHBH', 'BHNS']:

#       namesEMlist=[]

#       DCOname=dcotype
#       iii=0
        

#       # CREATE PANDAS FILE 
#       nModels=26
#       BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#       NAMES = []
#       stringgg = 'GW190814rate'

#       for ind_l, L in enumerate(BPSnameslist):
#           str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#           str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#           NAMES.append(str_z0)
#           NAMES.append(str_obs)
            
            


#       datas=[]

#       for i in range(len(BPSnameslist)):
#           datas.append(np.zeros_like(MSSFRnameslist))
#           datas.append(np.zeros_like(MSSFRnameslist))
            
            
#       df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#       df.columns =   df.columns.map(str)
#       df.index.names = ['xyz']
#       df.columns.names = ['m']

#       # print(df) 

#       df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')





# if INITIALIZE_NSBH==True:


#   namesEMlist=[]

#   DCOname ='BHNS'
#   iii=0
    

#   # CREATE PANDAS FILE 
#   nModels=26
#   BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#   NAMES = []
#   stringgg = 'NSBH'

#   for ind_l, L in enumerate(BPSnameslist):
#       str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#       str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#       NAMES.append(str_z0)
#       NAMES.append(str_obs)
        
        


#   datas=[]

#   for i in range(len(BPSnameslist)):
#       datas.append(np.zeros_like(MSSFRnameslist))
#       datas.append(np.zeros_like(MSSFRnameslist))
        
        
#   df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#   df.columns =   df.columns.map(str)
#   df.index.names = ['.x.y.z']
#   df.columns.names = ['m']

#   # print(df) 

#   df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringgg + '.csv')









# # #### INITIALIZE::: 
# if INITIALIZE_EM==True:

#   namesEMlist=[]

#   DCOname ='BHNS'
#   iii=0
#   for ind_chi, chi in enumerate([0.0, .5, 'Qin']):
#       # print(chi)
#       iii+=1
#       BH_chi   = chi 
#       for ind_Rns, NSradii in enumerate([11.5,13.0]):
#           iii+=1
#           Rns = NSradii
#           # if ind_mssfr ==0:
#           #   # print(chi)
#           stringg = 'Rns_'+ str(NSradii) + 'km_' + 'spinBH_' + str(chi) 
#           namesEMlist.append(stringg)


#           # CREATE PANDAS FILE 
#           nModels=26
#           BPSnameslist = list(string.ascii_uppercase)[0:nModels]

#           NAMES = []

#           for ind_l, L in enumerate(BPSnameslist):
#               str_z0 = str(L + ' intrinsic (z=0) [Gpc^{-3} yr^{-1}]')
#               str_obs = str(L + ' observed (design LVK) [yr^{-1}]')
#               NAMES.append(str_z0)
#               NAMES.append(str_obs)
                


#           datas=[]

#           for i in range(len(BPSnameslist)):
#               datas.append(np.zeros_like(MSSFRnameslist))
#               datas.append(np.zeros_like(MSSFRnameslist))
                
#           print(MSSFRnameslist)
#           df = pd.DataFrame(data=datas, index=NAMES, columns=MSSFRnameslistCSV).T
#           df.columns =   df.columns.map(str)
#           df.index.names = ['xyz']
#           df.columns.names = ['m']

#           # print(df) 

#           df.to_csv('/Users/floorbroekgaarden/Projects/GitHub/Double-Compact-Object-Mergers/dataFiles/summary_data_Fig_2/rates_MSSFR_Models_'+DCOname+ '_' + stringg + '.csv')


# # print(namesEMlist)


# # for DCOtype in ['BHNS']:

# #     for Rns in enumerate()

# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
# #     modelname = 'G'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# #     print('at DCOtype =', DCOtype)
# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'A'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'B'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)




# # for DCOtype in ['BHNS']:

# #     for Rns in enumerate()

# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/zeroBHkick/'
# #     modelname = 'G'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)

# #     print('at DCOtype =', DCOtype)
# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'A'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)


# #     pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/fiducial/'
# #     modelname = 'B'
# #     print('modelname')
# #     writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=True)





# # pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/alpha0_5/'
# # modelname, DCOtype = 'M', 'BNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BHNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BBH'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



# # pathCOMPASOutput = '/Volumes/Andromeda2/DATA/AllDCO/alpha2_0/'
# # modelname, DCOtype = 'N', 'BNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BHNS'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)

# # DCOtype='BBH'
# # writeToRatesFile(modelname=modelname, pathCOMPASOutput=pathCOMPASOutput, DCOtype=DCOtype, Optmistic=False)



