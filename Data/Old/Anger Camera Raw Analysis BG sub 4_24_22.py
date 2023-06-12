# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 11:25:02 2019

@author: 8jd
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
import scipy.constants as sc
import struct
import scipy.optimize as opt
import glob
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import os




#img is the limits for filtering data outside of the image
imgpar = [[0,1500],[0,1500]]
#imgpar = [[200,350],[150,300]]


spatialscale = 512./116.
#returns raw data and BG data, filtered by the imgdim and BG dim
def loadData(openfile, imgdim):
    listdata = []
    BGdata = []
    maxt = 0
    dummy = np.zeros(shape=(3))
    with open(openfile, "rb") as f:
        bytesin = f.read(24)
        while bytesin:
            temp = struct.unpack('3d', bytesin)
            dummy[0] = temp[0]
            dummy[1] = temp[1] * spatialscale
            dummy[2] = temp[2] * spatialscale
            listdata.append([dummy[0], int(dummy[1]-imgpar[0][0]), int(dummy[2]-imgpar[1][0])])
            
            bytesin = f.read(24)
    return listdata, BGdata

#This framestack is full image (spatially integrated)
#Think the raw data format is time in ticks (100ns/tick), position scaled to some length
def framestack(inarray, period, Nperiods, nbins):
    t_hist = np.zeros(shape=(nbins))
    binsize = period / nbins
    dlen = len(inarray)
    curindex = 0
    for i in range(0, dlen):
        curindex = int(np.floor( (inarray[i][0] % period) / binsize ))
        t_hist[curindex] += 1
    plt.plot(t_hist)
    return

def plotTotal(inarray):
    t_hist = np.zeros(shape=(int(max(map(lambda x: x[0], inarray)))+1))
    dlen = len(inarray)
    curindex = 0
    for i in range(0, dlen):
        curindex = int(np.floor( (inarray[i][0])))
        t_hist[curindex] += 1
    plt.plot(t_hist)
    return
    
def BGsubtractIMG(inImg, BG):
    ydim = int(max(map(lambda x: x[1], BG)) +1)
    xdim = int(max(map(lambda x: x[2], BG)) + 1)
    BGave = len(BG)/ (ydim * xdim)
    result = np.array(inImg)
    xlen = len(inImg)
    ylen = len(inImg[0])
    for i in range(0, xlen):
        for j in range(0,ylen):
            result[i][j] -= BGave
    return result

def ShowImage(inImage):
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    im = ax.imshow(inImage, origin="lower", aspect="auto")
    fig.colorbar(im)
    fig.tight_layout()
    plt.show()
    
def MakeImage(inarray):
    ydim = int(max(map(lambda x: x[1], inarray)) + 1)
    xdim = int(max(map(lambda x: x[2], inarray)) + 1)
    img = np.zeros(shape=(xdim, ydim))
    dlen = len(inarray)
    for i in range(0, dlen):
        xind = int(inarray[i][1])
        yind = int(inarray[i][2])
        img[yind][xind] += 1
    
    return img

#This function combines bins in time, need to take care to change scale when using
def BinDataT(inData, n):
    bSize = n
    result = np.array(inData)
    nData = len(inData)
    for i in range(0, nData):
        result[i][0] = np.floor(result[i][0] / bSize)
    return result
    
#This function bins spatial dimensions, need to take care to change scale when using
def BinImg(inImg, xbinsize, ybinsize):
    xdim = len(inImg[0])
    ydim = len(inImg)
    ybins = int(np.floor(ydim/ybinsize))
    xbins = int(np.floor(xdim/xbinsize))
    result = np.zeros(shape=(ybins, xbins))
        
    for xb in range(0, xbins):
        curx = xb * xbinsize
        for yb in range(0, ybins):
              cury = yb * ybinsize
              for x in range(curx, curx + xbinsize):
                  for y in range(cury, cury+ybinsize):
                      result[yb][xb] += int(inImg[y][x])
    return result
    
def BinRaw(inData, xbinsize, ybinsize):
    result = []
    #lastXbin = np.floor(max(map(lambda x: x[1], inData)) / xbinsize)
    #lastYbin = np.floor(max(map(lambda x: x[2], inData)) / ybinsize)
    for datum in inData:
        newXbin = datum[1]/xbinsize
        newYbin = datum[2]/ybinsize
        #if (newXbin > lastXbin):
        #    continue
        #if (newYbin > lastYbin):
        #    continue
        result.append([datum[0], int(np.floor(newXbin)), int(np.floor(newYbin))])
    return np.array(result)

def fitfun(x, a, f, phase, c):
    return c + a * np.sin(2*sc.pi*f *x + phase)

def tfit(tData, par0):
    t = np.linspace(0, len(tData), len(tData))
    #frequency estimate is in kHz, converted to pixels
    par0[0] = np.std(tData)
    par0[3] = np.mean(tData)
    popt, pcov = opt.curve_fit(fitfun, t, tData, p0 = par0, maxfev=1000000,ftol=1e-8)
    tt = np.linspace(0, len(tData), 1000)
    plt.plot(tt, fitfun(tt, *popt), linestyle="-", color='b')
    plt.plot(tData, linestyle="None", marker="o", color='r')
    plt.show()
    return popt

def T_analysis_byPixel(inData, par0):
    maxt = int(max(map(lambda x: x[0], inData)))
    #using hash table to uniquely identify each pixel 
    maxx = int(max(map(lambda x: x[1], inData)))
    maxy = int(max(map(lambda x: x[2], inData)))
    hashmax = maxx*(maxy + 1) + maxy
    #hash: i_x*(n_y - 1) + i_y (or vice versa)
    #iterate through all data once, separating into different hash groups.
    #then analyze the hash lists separately
    hashlist = np.zeros(shape=(hashmax+1, maxt+1), dtype=np.int32)
    for datum in inData:
        curhash = int(datum[1]*(maxy + 1) + datum[2])
        hashlist[curhash][int(datum[0])] += 1
    imgOut = np.zeros(shape=(maxy+1, maxx+1))
    for i in range(0, hashmax):
        tlist = hashlist[i]
        pixeldata = tfit(tlist, par0)
        xpix = int(np.floor((i / (maxy + 1)) ))
        ypix = int(i % (maxy + 1))
        #currently outputs in degrees
        imgOut[ypix][xpix] = 360*pixeldata[2]/(2.*sc.pi) #select which item from the fit you want to view
    return imgOut

def totalCounts(inData):
    return len(inData)

def scanCounts(f_i, f_f, imgdim, BGdim, directory):
    result = []
    filelist = []
    filenum = ""
    filenumlist = []
    configlist = []
    for filename in glob.glob(directory+"MURR*.dat"):
        filenum = filename[filename.find("MURR")+4:filename.find(".dat")]
        if (int(filenum) >= f_i and int(filenum) <= f_f):
            filelist.append(filename)
            filenumlist.append(filenum)
            configlist.append(getConfig(filenum, directory))
    for filename in filelist:
        data, BG = loadData(filename, imgdim, BGdim)
        BGperPix = totalCounts(BG) / ((BGdim[0][1] - BGdim[0][0]) * (BGdim[1][1] - BGdim[1][0]))
        rawcount = totalCounts(data)
        finalcount = rawcount - BGperPix * ((imgdim[0][1] - imgdim[0][0]) * (imgdim[1][1] - imgdim[1][0]))
        result.append(finalcount)
    return filenumlist, np.array(configlist), np.array(result)
 
#easier if runID is input as a string
#Big run did 2D scan of B2 and G2, try 3D plot or heatmap type plot
def getConfig(runID, direct):
    result = []
    for line in open(direct + "logfile.txt", "r"):
        element = line.find("MURR" + runID)
        if element > 0:
            result.append(round(float(line[line.find("__B1")+5:line.find("_B2")]) + 1e-4, 2)) #append B1 value
            result.append(round(float(line[line.find("_B2")+4:line.find("_G1")])+ 1e-4, 2)) #append B2 value
            result.append(round(float(line[line.find("_G1")+4:line.find("_G2")])+ 1e-4, 2)) #append G1 value
            result.append(round(float(line[line.find("_G2")+4:line.find("_N1")])+ 1e-4, 2)) #append G2 value
            result.append(round(float(line[line.find("_N1")+4:line.find("_N2")])+ 1e-4, 2)) #append N1 value
            result.append(round(float(line[line.find("_N2")+4:line.find("_RF1Fre")])+ 1e-4, 2)) #append N2 value
            result.append(round(float(line[line.find("_RF1Fre")+8:line.find("_RF1Amp")])+ 1e-4, 2)) #append RF1Freq value
            result.append(round(float(line[line.find("_RF1Amp")+8:line.find("_RF1phase")])+ 1e-4, 2)) #append RF1Amp value
            result.append(round(float(line[line.find("_RF1phase")+10:line.find("_RF2Fre")])+ 1e-4, 2)) #append RF1phase value
            result.append(round(float(line[line.find("_RF2Fre")+8:line.find("_RF2Amp")])+ 1e-4, 2)) #append RF2Freq value
            result.append(round(float(line[line.find("_RF2Amp")+8:line.find("_RF2phase")])+ 1e-4, 2)) #append RF2Amp value
            result.append(round(float(line[line.find("_RF2phase")+10])+ 1e-4, 2)) #append RF2phase value
            return result
    return result
            
def plot2dparam(configlist, countlist, n1, n2):
    histx, histy = convertarraystohist(configlist[:,n1], configlist[:,n2], countlist)
    h, xe, ye, im = plt.hist2d(histx, histy, bins=[50,27])
    plt.clf()
    plt.close()
    ext = [xe[0], xe[-1], ye[0], ye[-1]]
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    plt.imshow(h.T, extent=ext, origin="lower", interpolation="None", aspect='auto')
    cb = fig.colorbar(im)
    cb.ax.tick_params(labelsize=16)
    plt.tick_params(which='both', labelsize=16)
    ax.set_ylabel("$B_{rf}$ [Vpp]", fontsize=18)
    #plt.yticks(np.linspace(0,len(ylab),4), np.round(np.array([ylab[0],ylab[int(len(ylab)/3)], ylab[int(len(ylab)*2/3)], ylab[len(ylab)-1]]),2), fontsize=16)
    ax.set_xlabel("Gradient field [A]", fontsize=18)
    fig.tight_layout()
    plt.show()
    
def convertarraystohist(x, y, z):
    outx = []
    outy = []
    ndata = len(z)
    for i in range(0, ndata):
        for j in range(0, int(z[i])):
            outx.append(x[i])
            outy.append(y[i])
    return outx, outy

def convertarraystoimg(x, y, z):
    sort1 = sorted(set(x))
    sort2 = sorted(set(y))
    stepx = round(abs(sort1[1] - sort1[0]),2)
    stepy = round(abs(sort2[1] - sort2[0]),2)
    minx = int(min(x)*100)
    miny = int(min(y)*100)
    result = np.zeros(shape=(len(sort2), len(sort1)))
    #flipping x and y d/t image being row/col opposite of standard x/y notation
    xlabel = np.round(np.copy(sort1),2)
    ylabel = np.round(np.copy(sort2),2)
    ndata = len(z)
    for i in range(0, ndata):
        result[int(np.floor((y[i]*100-miny)/(stepy*100))), int(np.floor((x[i]*100-minx)/(stepx*100)))] = z[i]
    return result, xlabel, ylabel
 
def show2DdataImage(img, xlab, ylab):
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    ext = [xlab[0], xlab[-1], ylab[0], ylab[-1]]
    im = ax.imshow(img, origin="lower", extent=ext, aspect="auto", interpolation="None")
    cb = fig.colorbar(im)
    cb.ax.tick_params(labelsize=16)
    plt.tick_params(which='both', labelsize=16)
    ax.set_ylabel("$B_{rf}$ [Vpp]", fontsize=18)
    #plt.yticks(np.linspace(0,len(ylab),4), np.round(np.array([ylab[0],ylab[int(len(ylab)/3)], ylab[int(len(ylab)*2/3)], ylab[len(ylab)-1]]),2), fontsize=16)
    ax.set_xlabel("Gradient field [A]", fontsize=18)
    fig.tight_layout()
    plt.show()    

first, last = [230, 235]    #will include this index
names = np.arange(first, last + 1, 1, dtype=int)
filenames = ["CG4B" + str(i) + ".dat" for i in names]
print("Files used:", filenames)

direct = "C:/Users/Larmor/Dropbox/AngerCamera/Wide_Angle_NSE_2022/"

binsize = 3    # number of pixels binned together, each pixel is 0.2 mm x 0.2 mm
ROIsumlist = np.zeros((len(names)))
BGsumlist = np.zeros((len(names)))
BGperpixellist = np.zeros((len(names)))

MissingFiles = []
for i in np.arange(len(names)):
    if os.path.exists(direct+filenames[i]) == True:
        data, BG = loadData(direct+filenames[i], imgpar)
        RawImg = MakeImage(data) 
        BinnedImg = BinImg(RawImg, binsize, binsize)    #returns 2d np array

        ShowImage(BinnedImg)
        
        zstart = 1 #read from vertical axis of images
        zend = 160
        xstart = 1 #read from horizontal axis image 
        xend = 160
    
        rowrange = np.arange(xstart, xend + 1, 1)
        columnrange = np.arange(zstart, zend + 1, 1)
        
        #Background BG subtraction
        #Find Sum of background counts in background area, divide by number of pixels, subtract from each pixel specified above
        BGzstart = 0 #read from vertical axis of images
        BGzend = 40
        BGxstart = 0 #read from horizontal axis image 
        BGxend = 40
    
        BGrowrange = np.arange(BGxstart, BGxend + 1, 1)
        BGcolumnrange = np.arange(BGzstart, BGzend + 1, 1)

        for c in BGcolumnrange:
            for r in BGrowrange:
                BGsumlist[i] += BinnedImg[c,r]     #total number of background counts
        numBGpixels = len(BGrowrange)*len(BGcolumnrange)

        BGperpixellist[i] =  BGsumlist[i]/numBGpixels
        print("Subscan name: ",filenames[i])
        print("BG per pixel:", BGperpixellist[i])
        
        BGcorrectionxrange = rowrange
        BGcorrectionzrange = columnrange
        for c in BGcorrectionzrange:
            for r in BGcorrectionxrange:
                BinnedImg[c,r] -= BGperpixellist[i]
        print(BinnedImg)
        """
        df.to_excel (r'C:/Users/Larmor/Dropbox/AngerCamera/Analysis/Wide_Angle_NSE_2022/IndividualScansinExcel/BGsubtracted%sbin%s.xlsx' %(filenames[i],binsize), index = False, header=False)
        ROIzstart = 2 #read from vertical axis of images
        ROIzend = 75
        #Used z 5,70 x 90 145 binsize 3 for 2022/ scan 231-
        ROIxstart = 60 #read from horizontal axis image 
        ROIxend = 145
    
        ROIrowrange = np.arange(xstart, xend + 1, 1)
        ROIcolumnrange = np.arange(zstart, zend + 1, 1)
        
        for c in ROIcolumnrange:
            for r in ROIrowrange:
                ROIsumlist[i] += df.iloc[c,r]
        print("ROI sum", ROIsumlist[i])
        
        db = pd.DataFrame(ROIsumlist)
        db.to_excel (r'C:/Users/Larmor/Dropbox/AngerCamera/Analysis/Wide_Angle_NSE_2022/IndividualScansinExcel/ROIsumlistinitial%sfinal%sbin%s.xlsx' %(first,last,binsize), index = False, header=False)
        """
    else:
        MissingFiles.append(filenames[i])
print("Missing Files",MissingFiles)
#makes plot of ROI vs scan    
'''
plt.figure()   
#plt.scatter(names,ROIsumlist)
#plt.title("bin %s rows %s - %s columns %s - %s"  %(binsize, xstart, xend, zstart, zend))
'''
#plt.savefig('/Users/stephenkuhn/Dropbox/CG4B correction coil/Analysis/DetectorIntensity_vs_Run/%sbin%s.png' %(filenames[0],binsize))