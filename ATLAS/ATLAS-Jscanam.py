### H-ATLAS script which is a modified version of HIPE Jscanamorphos pipeline
##  edited by Matthew Smith (Matthew.Smith@astro.cf.ac.uk) - Sep 2015

# This script allows the user to pass more than one Herschel observation to Jscanamorphos, rather than 
# only taking just a pair of observations.

### Input Parameters ###

# select cameras to process (select blue, green or red)
masterFilter = {"blue":False, "green":True, "red":False}

# the script can run through a number of tiles, the example below will process tile "01" and "02"
# for each tile the nominal and othogonal obsids are listed
obsids = {"01":{"nominal":[1342211292, 1342222677],"orthogonal":[1342210931,1342210902]},\
          "02":{"nominal":[1342222626, 1342222676],"orthogonal":[1342210567,1342210903]}}

# set name of the map to produce
outName = "SGP"

# set pixel size of the image to produce
pixelsize = {"blue": 3.0, "green":3.0, "red": 4.0}

# set folder to output the final map (note no final '/')
outFolder = "/path/to/folder"

# The H-ATLAS team process their L1-data seperatly and save them in a fits file. If your level-1
# data is saved, set parameter below to True and adjust L1folder path 
# the files must contain the obsid in their file name!
# if you wish to use the HSA level-1 data set saved L1 to false (not well tested)
savedL1 = True
if savedL1:
    # select folder with level-1 files
    L1folder = "/home/spiredata/spxmws/over-drive/ATLAS/SGP/jscanam100/L1data"

################################################################################

#
#  This file is part of Herschel Common Science System (HCSS).
#  Copyright 2001-2013 Herschel Science Ground Segment Consortium
# 
#  HCSS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as
#  published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
# 
#  HCSS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
# 
#  You should have received a copy of the GNU Lesser General
#  Public License along with HCSS.
#  If not, see <http://www.gnu.org/licenses/>.
#

""" 
TIP: if you edit this script in hipe:
- all lines of code actually executed by this script will appear in black,
- all comments and explanations will appear in green,
- all line of codes that are not executed by default but that could be useful
  in some particular cases, will appear in red.

If you want to edit the script and remove comments and commented lines of code,
be careful in preserving the indentation. Otherwise if statements and loops
might not work.

This script invokes the complete Scanamorphos port for the PACS photometer 
pipeline. 

Description:
- This script takes your data from Level 1 to 2.5, starting from a Level 1 HSA 
  product.
- This script is written so that you can run it entirely in one go or 
  line-by-line.  
- Comments are included in this script, but for a full explanation see the PACS 
  Data Reduction Guide (PDRG): 
  Help Menu -> Help Contents -> PACS Data Reduction Guide: Photometry

Inputs:
- camera: Photometer camera ("red" or "blue"). Be careful with memory 
  allocation, blue requires about 4 times the memory used in red.
- obsids: List of observation ids with the scans and cross scans needed to 
  combine. The cross scan direction should be perpendicular to the scan 
  direction.
- solarSystemObject: Set it to True if the source is a Solar System object 
  (i.e. is a moving target). Note that if there is extended emission in the
  map, the main JScanam assumptions will fail, because the scan and cross 
  scan will be shifted to the object reference system.
- galactic: The galactic option should be set to True if the maps extended 
  emission is not confined to a small area. The safest is to set it always 
  True. 
- calculateRaDec: Calculate the coordinates of every pixel with time. The 
  tasks will run faster if this option is set to True, but they will require 
  a factor of 2-3 more memory.
- showMapsAfterTasks: If True, it will produce a map after each step in the 
  pipeline.
- debug: Debug mode. Use wisely, debug=True will clutter your desktop with 
  graphs and displays. It will also consume more memory.
- deglitch: Set it to True if you want to deglitch the frames with a lower
  threshold than the one used by default (nSigmaDeglitch = 5).
- nSigmaDeglitch: The new deglitch threshold to use in case the deglitch 
  parameter is set to True.
- makeFinalMap: Generate a final map.
- outputPixelsize: Pixel size to use for the final map. It will use 1.6" for 
  the blue camera and 3.4" for the red camera if it's set to -1. 
- pixfrac: This parameter is used for the final map. It fixes the ratio between 
  the input detector pixel size and the map pixel size.

########################### SECTION 0 ##########################################
################################################################################
########################### SETTINGS ###########################################
"""

import os
from os.path import join as pj

solarSystemObject = False
galactic = False
calculateRaDec = False
showMapsAfterTasks = False
debug = False
deglitch = True
nSigmaDeglitch = 3
makeFinalMap = True
pixfrac = 0.1

bands = []
for band in masterFilter.keys():
    if masterFilter[band] == True:
        bands.append(band)

for tile in obsids.keys():
    for band in bands:
        filter = {"blue":False, "green":False, "red":False}
        filter[band] = True
		

        """
        ########################### SECTION 1 ##########################################
        ################################################################################
        ########################### PROCESSING #########################################
        """
        print " -- Executing JScanam -- "
        print "Starting Tile: ", tile
        
        # check only one filter has been selected
        if SUM(filter.values()) != 1:
        	raise Exception("One band not selected")
        # set camera variable based on filter
        if filter["red"]:
        	camera = "red"
        else:
        	camera = "blue" 
        
        # set map pixel size
        if filter["red"]:
        	outputPixelSize = pixelsize["red"]
        elif filter["green"]:
        	outputPixelSize = pixelsize["green"]
        else:
        	outputPixelSize = pixelsize["blue"]
        
        print " Loading Level 1 data "
        
	if savedL1:
        	# list all files in folder
        	allfiles = os.listdir(L1folder)
        
        	# restrict to the desired obsid files
        	nomL1files = []
        	for obsid in obsids[tile]["nominal"]:
        		# see if file present
        		if filter["red"]:
        			posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("red") > 0]
        		elif filter["green"]:
        			posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("green") > 0]
        		else:
        			posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("blue") > 0]
        		if len(posFiles) != 1:
        			raise Exception("Found " + str(len(posFiles)) + " possible files for obsid: " + str(obsid))
        		nomL1files.append(posFiles[0])
       		# restrict to the desired obsid files
        	orthL1files = []
        	for obsid in obsids[tile]["orthogonal"]:
        		# see if file present
        		if filter["red"]:
        			posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("red") > 0]
	        	elif filter["green"]:
	        		posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("green") > 0]
	        	else:
	        		posFiles = [allfile for allfile in allfiles if allfile.count(str(obsid)) > 0 and allfile.count("blue") > 0]
	        	if len(posFiles) != 1:
	        		raise Exception("Found " + str(len(posFiles)) + " possible files for obsid: " + str(obsid))
	        	orthL1files.append(posFiles[0])
	        
	        
	        fits = FitsArchive()
	        scansList = [fits.load(pj(L1folder,nomL1files[0]))]
	        blueFilter1 = scansList[0].meta["blue"].value
	        for i in range(1,len(nomL1files)):
	        	# combine scans
	        	scansList.append(fits.load(pj(L1folder,nomL1files[i])))
	
	        cscansList = [fits.load(pj(L1folder,orthL1files[0]))]
	        blueFilter2 = cscansList[0].meta["blue"].value
	        for i in range(1,len(orthL1files)):
	        	# combine scans
	        	cscansList.append(fits.load(pj(L1folder,orthL1files[i])))
        else:
	        fits = FitsArchive()
		scansList = []
		for obsid in obsids[tile]["nominal"]:	
			scansObs = getObservation(obsid, useHsa=True, instrument="PACS")
		        level1 = PacsContext(scansObs.level1)
		        scans = level1.averaged.getCamera(camera).product.selectAll() 
		        blueFilter1 = scans.meta["blue"].value
	        	scansList.append(scans)
		del(scansObs, level1, scans)

		cscansList = []
		for obsid in obsids[tile]["orthogonal"]:	
			cscansObs = getObservation(obsid, useHsa=True, instrument="PACS")
		        level1 = PacsContext(cscansObs.level1)
		        cscans = level1.averaged.getCamera(camera).product.selectAll() 
		        blueFilter2 = cscans.meta["blue"].value
	        	cscansList.append(cscans)
		del(cscansObs, level1, cscans)


        #calTree = getCalTree(obs=cscansObs)
        calTree = getCalTree(time=scansList[0].startDate)
        
        if (camera == "blue" and blueFilter1 != blueFilter2):
          print ""
          print "  ------------------------------------------------------------------"
          print "  ALERT!!! You are trying to combine two blue observations obtained "  
          print "  with different filters, which is most probably wrong!             "
          print "  ------------------------------------------------------------------"
          print ""
        
        """
         Set the scans and cross scans coordinates to the object reference system if it's a SSO
        """
        if(solarSystemObject):
          print " Setting the scans and crossScans coordinates to the object reference system "
          pp = scansObs.auxiliary.pointing
          orbitEphem = scansObs.auxiliary.orbitEphemeris
          horizons = scansObs.auxiliary.horizons
          cal = getCalTree(obs=scansObs)
          scans = photAddInstantPointing(scans, pp, calTree=cal, orbitEphem=orbitEphem, horizonsProduct=horizons)
          timeOffset = scans.getStatus("FINETIME")[0]
          scans = correctRaDec4Sso(scans, timeOffset=timeOffset, orbitEphem=orbitEphem, horizonsProduct=horizons, linear=0)
          #
          pp = cscansObs.auxiliary.pointing
          orbitEphem = cscansObs.auxiliary.orbitEphemeris
          horizons = cscansObs.auxiliary.horizons
          cal = getCalTree(obs=cscansObs)
          cscans = photAddInstantPointing(cscans, pp, calTree=cal, orbitEphem=orbitEphem, horizonsProduct=horizons)
          cscans = correctRaDec4Sso(cscans, timeOffset=timeOffset, orbitEphem=orbitEphem, horizonsProduct=horizons, linear=0)
          del pp, orbitEphem, horizons, cal, timeOffset
        
        """
         Calculate the coordinates of every pixel with time
        """
        if(calculateRaDec):
           print " Calculating the coordinates of every pixel with time "
           for i in range(0,len(scansList)):
                scansList[i] = photAssignRaDec(scansList[i], calTree=calTree)
           for i in range(0,len(cscansList)):
                cscansList[i] = photAssignRaDec(cscansList[i], calTree=calTree)
        
        """
         Reduce the product size removing unnecessary status information
        """
        # print " Reducing the product size "
        # scans = scanamorphosReduceFramesSize(scans)
        # cscans = scanamorphosReduceFramesSize(cscans)
        
        """
         Remove turn arounds
        """
        print " Removing turn arounds "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosRemoveTurnarounds(scansList[i], limit=50.0, debug=debug)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosRemoveTurnarounds(cscansList[i], limit=50.0, debug=debug)
        
        """
         Mask long term glitches. 
         This task produces a mask called Scanam_LongTermGlitchMask. You should check 
         this mask and if the results are not optimal (some glitches are not detected 
         or some sources are flagged as glitches), you can try to modify the parameter 
         stepAfter in order to get better results.
        """
        print " Masking long term glitches "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosMaskLongTermGlitches(scansList[i], stepAfter=20, galactic=galactic, calTree=calTree, debug=debug)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosMaskLongTermGlitches(cscansList[i], stepAfter=20, galactic=galactic, calTree=calTree, debug=debug)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans: " +str(i) + " after masking long term glitches")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) + " after masking long term glitches")

        
        """
         Save the scans and cross scans for later use
        """
        print " Saving a copy of the scans and cross scans in a temporal pool "
        from herschel.pacs.share.util import PacsProductSinkWrapper
        scansRef = []
        cscansRef = []
        for i in range(0,len(scansList)):
             scansRef.append(PacsProductSinkWrapper.getInstance().saveAlways(scansList[i]))
        for i in range(0,len(cscansList)):
             cscansRef.append(PacsProductSinkWrapper.getInstance().saveAlways(cscansList[i]))
        
        """
         Subtract the baseline of every scanleg in every pixel
        """
        print " Initial baseline subtraction / per scanleg "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosScanlegBaselineFitPerPixel(scansList[i], nSigma=2)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosScanlegBaselineFitPerPixel(cscansList[i], nSigma=2)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans: " + str(i) + " after baseline subtraction / per scanleg")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) + " after baseline subtraction / per scanleg")
        
        """
         Create the source mask.
         Modify nSigma until you get an optimal mask. The mask should cover only a small
         fraction of the map (<~30%). It's not necessary that all the faint emission is
         masked, only the brightest regions.
        """
        print " Creating a source mask "
        scans = scansList[0].copy()
        scans.join(cscansList[0])
        for i in range(1,MAX([len(scansList),len(cscansList)])):
             if len(scansList) < i:
                  scans.join(scansList[i])
             if len(cscansList) < i:
                  cscans.join(cscansList[i])
        del(scansList, cscansList)
        sourceImage, scans = scanamorphosCreateSourceMask(scans, nSigma=4.0, createMask=False, galactic=galactic, calTree=calTree, debug=debug)
        
        if(showMapsAfterTasks):
           d = Display(sourceImage, title="Masked sources")
        
        """
         Replace the scans and cross scans by the saved ones
        """
        print " Loading the saved scans and cross scans from the temporal pool "
        del(scans)
        scansList = []
        cscansList = []
        for i in range(0,len(scansRef)):
             scansList.append(scansRef[i].product)
        for i in range(0,len(cscansRef)):
             cscansList.append(cscansRef[i].product)
        del(scansRef, cscansRef)
        System.gc()
        
        """
         Add the mask to the saved scans and cross scans
        """
        print " Adding the source mask to the scans and cross scans "
        for i in range(0,len(scansList)):
             maskImage, scansList[i] = scanamorphosCreateSourceMask(scansList[i], inputImage=sourceImage, createMask=True, calTree=calTree, debug=debug)
        for i in range(0,len(cscansList)):
             maskImage, cscansList[i] = scanamorphosCreateSourceMask(cscansList[i], inputImage=sourceImage, createMask=True, calTree=calTree, debug=debug)
        
        """
         Baseline subtraction. 
         Here we use galactic=True, because we only want to remove an offset
        """
        print " Baseline subtraction "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosBaselineSubtraction(scansList[i], galactic=True, debug=debug)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosBaselineSubtraction(cscansList[i], galactic=True, debug=debug)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans: " + str(i) +  " after baseline subtraction")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) +" after baseline subtraction")
        
        """
         Baseline pre-processing
        """
        print " Baseline pre-processing "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosBaselinePreprocessing(scansList[i], debug=debug)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosBaselinePreprocessing(cscansList[i], debug=debug)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans: " + str(i) + " after baseline pre-processing")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) + " after baseline pre-processing")
        
        """
         Baseline subtraction
        """
        print " Baseline subtraction "
        for i in range(0,len(scansList)):
             scansList[i] = scanamorphosBaselineSubtraction(scansList[i], galactic=galactic, debug=debug)
        for i in range(0,len(cscansList)):
             cscansList[i] = scanamorphosBaselineSubtraction(cscansList[i], galactic=galactic, debug=debug)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans: " + str(i) + " after baseline subtraction")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) + " after baseline subtraction")
        
        """
         Destriping
        """
        print " Destriping "
        scans = scansList[0].copy()
        for i in range(1,len(scansList)):
             scans.join(scansList[i])
        del(scansList)
        cscans = cscansList[0].copy()
        for i in range(1,len(cscansList)):
             cscans.join(cscansList[i])
        del(cscansList)
        scans, cscans = scanamorphosDestriping(scans, cscans, iterations=6, calTree=calTree, debug=debug)
        
        if(showMapsAfterTasks):
           for i in range(0,len(scansList)):
                map, mi = photProject(scansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Scans after destriping")
           for i in range(0,len(cscansList)):
                map, mi = photProject(cscansList[i], pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
                d = Display(map, title="Cross scans: " + str(i) + " after destriping")
        
        """
         Merge the scans and cross scans
        """
        print " Merging scans and cross scans in a frame product "
        System.gc()
        mergedScans = scans
        mergedScans.join(cscans)
        del(scans, cscans)
        System.gc()
        
        """
         Deglitch the merged scans
         
         We use nSigma=5, and not a lower value, to avoid masking signal with strong drifts 
         that will be corrected by the scanamorphosIndividualDrifts task
        """
        print " Deglitching the merged scans "
        scanamorphosDeglitch(mergedScans, nSigma=5, calTree=calTree, debug=debug)
        
        """
         Remove individual drifts
        """
        print " Removing individual drifts "
        mergedScans = scanamorphosIndividualDrifts(mergedScans, calTree=calTree, debug=debug)
        
        if(showMapsAfterTasks):
           map, mi = photProject(mergedScans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
           d = Display(map, title="Merged scans after individual drifts correction")
        
        """
         Deglitch the merged scans again with the user provided nSigmaDeglitch
        
         Be careful setting nSigmaDeglitch to very small values, because it can
         mask bright sources. You should always check the mask called Scanamorphos_GlitchMask
         to make sure that the task is masking only real glitches.
        """
        if(deglitch):
            print " Deglitching the merged scans "
            scanamorphosDeglitch(mergedScans, nSigma=nSigmaDeglitch, calTree=calTree, debug=debug)
        
        """
         Create a new source mask, and remove the background median offset 
        """
        # print " Obtaining a new source mask "
        # sourceImage, mergedScans = scanamorphosCreateSourceMask(mergedScans, nSigma=2.0, createMask=True, galactic=galactic, calTree=calTree, debug=debug)
        #
        # if(showMapsAfterTasks):
        #   d = Display(sourceImage, title="Masked sources")
        #
        # print " Removing the background median offset "
        # mergedScans = scanamorphosSubtractOffset(mergedScans)
        #
        # if(showMapsAfterTasks):
        #   map, mi = photProject(mergedScans, pixfrac=0, calTree=calTree, useMasterMaskForWcs=False)
        #   d = Display(map, title="Merged scans after the background median offset subtraction")
        
        
        """
        ########################### SECTION 2 ##########################################
        ################################################################################
        ########################### SEE THE RESULTS#####################################
        
         Save the frames
        """
        # outputFitsFile = "/your/Home/Dir/finalResult.fits"
        # FitsArchive().save(outputFitsFile, mergedScans)
        
        """
         Final projection with drizzle
        """
        if filter["red"]:
            outputPixelSize = pixelsize['red']
        elif filter["green"]:
            outputPixelSize = pixelsize['green']
        else:
            outputPixelSize = pixelsize['blue']
        
        if(makeFinalMap):
           print " Projecting the merged scans onto the final map "
           print " Projection with drizzle using: " + str(outputPixelSize) + " arcsec"
           finalMap, mi = photProject(mergedScans, outputPixelsize=outputPixelSize, pixfrac=pixfrac, calTree=calTree, useMasterMaskForWcs=False)
           #d = Display(finalMap, title="Final map")
           # outputMapFile = "/your/Home/Dir/finalMap_"+camera+".jpg"
           # d.saveAsJPG(outputMapFile)
        
           # save map fits file
           if filter["red"]:
           	outFile = outName +"-" + tile + "_pacs160.fits"
           elif filter["green"]:
           	outFile = outName +"-" + tile + "_pacs100.fits"
           else:
           	outFile = outName +"-" + tile + "_pacs70.fits"
           fits.save(pj(outFolder,outFile),finalMap)
        
           # delete variables
           del(finalMap)

print "Finished Successfully"
