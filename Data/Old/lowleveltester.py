# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import socket
import struct
import time
import numpy
#import pylab
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import sys
import threading
import ftplib
from PyQt5 import  uic, QtWidgets, QtCore

global app
global debugon

qtCreatorFile = "C:/Users/Larmor/Desktop/AngerCamera/lowleveltest.ui" # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)


global myfile
global fileopen

global psum_histox  #only one of these
global psum_histoy
global psum_histo_xaxis
global ad_histo_xaxis
global ad_histox #declared in main as rawx_advalues[8,4,4096]
global ad_histoy
global tof_histo
global zero_histos
global zero_histos_toggle
global psum_scale
global ad_scale
global bm1counts
global gcounts
global grate
global global_acquireflag
global twoD_histo
global printonce
global psum3Dhisto


#for testing only
#global firstten
#global firstin



from scipy.optimize import leastsq
from scipy.special import erf
#below is for 1mm no grease
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.044*(x-p[1])**2)+0.2*numpy.exp(-0.003*(x-p[1])**2))
#below is for 2mm no grease
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.037*(x-p[1])**2)+0.2*numpy.exp(-0.0045*(x-p[1])**2))
#below is for 2mm? with 2.8mm spacer
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.023*(x-p[1])**2)+0.2*numpy.exp(-0.003*(x-p[1])**2))
#below is for 2mm? with 2.3 mm spacer
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.025*(x-p[1])**2)+0.2*numpy.exp(-0.0045*(x-p[1])**2))
#below is for 1mm? grease and cookie
gaussvar= lambda p,x: p[0]*(0.775*numpy.exp(-0.029*(x-p[1])**2)+0.225*numpy.exp(-0.005*(x-p[1])**2))
#below is for 1mm? just grease
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.028*(x-p[1])**2)+0.2*numpy.exp(-0.006*(x-p[1])**2))
#below is for 1mm? no grease + .068" additional spacer
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.045*(x-p[1])**2)+0.2*numpy.exp(-0.0015*(x-p[1])**2))
#below is for 1mm? with 2.3 standard  +  1mm additional spacer
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.019*(x-p[1])**2)+0.2*numpy.exp(-0.003*(x-p[1])**2))
#below is for 1mm? no grease + + 5.3mm spacer no grease or cookie
#gaussvar= lambda p,x: p[0]*(0.8*numpy.exp(-0.034*(x-p[1])**2)+0.2*numpy.exp(-0.005*(x-p[1])**2))


def fit_funcvar(func, p,x,y,ye):
    errfunc = lambda p,x,y,ye: numpy.abs(func(p,x) - y )
    par, cov = leastsq(errfunc,p,args=(x,y,ye),maxfev=5)
#    chi2 = numpy.sum(errfunc(par,x,y,y)**2)/(len(x) - len(par))
#    parerr=np.sqrt(np.abs(cov.diagonal()))
    return par #, chi2


def analysisvar(xdata,ydata):
#    ye=numpy.sqrt(ydata)
    ye=1.0
#    xf=numpy.linspace(xdata[0],xdata[-1],100)

    xmean=numpy.sum(xdata*ydata)/numpy.sum(ydata)
#    sig=np.sqrt(abs(np.sum((xdata-xmean)**2*ydata)/np.sum(ydata)))
    #sig is fixed.
    ampl=numpy.max(ydata)

    gpar = fit_funcvar(gaussvar,[ampl, xmean],xdata,ydata,ye)
    return gpar

def analysisvaredge(xdata,ydata,xmean):
#    ye=numpy.sqrt(ydata)
    ye=1.0
#    xf=numpy.linspace(xdata[0],xdata[-1],100)

#    xmean=numpy.sum(xdata*ydata)/numpy.sum(ydata)
#    sig=np.sqrt(abs(np.sum((xdata-xmean)**2*ydata)/np.sum(ydata)))
    #sig is fixed.
    ampl=numpy.max(ydata)

    gpar = fit_funcvar(gaussvar,[ampl, xmean],xdata,ydata,ye)
    return gpar



class streamthread (threading.Thread):
    def __init__(self, threadID, socketID):
        self.RAW_EVT_FMT='32I'
        self.RAW_EVT_SZ=struct.calcsize(self.RAW_EVT_FMT)
        self.NORMAL_EVT_FMT='2I'
        self.NORMAL_EVT_SZ=struct.calcsize(self.NORMAL_EVT_FMT)
        self.numprints=0
        self .killstream=0
        self.xs=numpy.zeros((8,4),dtype=numpy.uint32)
        self.ys=numpy.zeros((8,4),dtype=numpy.uint32)
        self.xcaloffsets=numpy.zeros(32,dtype=numpy.float32)
        self.ycaloffsets=numpy.zeros(32,dtype=numpy.float32)
        self.xcalnew=numpy.zeros(16)
        self.ycalnew=numpy.zeros(16)
        self.xnew=numpy.zeros(16)
        self.ynew=numpy.zeros(16)
        self.lengthscale=116.0/512.0

        threading.Thread.__init__(self)
        self.threadID=threadID
        self.socketID=socketID
        self.printonce=0
        self.savelc=0
        self.index=threadID-1
        self.tablesin=0
        self.saveraw=0
        self.flushme=0
        self._stop = threading.Event()
        self.positions=numpy.zeros(16)
#        self.positions[0]=3.9
#        self.positions[1]=11.1
#        self.positions[2]=18.3
#        self.positions[3]=25.5
#        self.positions[4]=32.7
#        self.positions[5]=39.9
#        self.positions[6]=47.1
#        self.positions[7]=54.3
#        self.positions[8]=61.8
#        self.positions[9]=69.0
#        self.positions[10]=76.2
#        self.positions[11]=83.4
#        self.positions[12]=90.6
#        self.positions[13]=97.8
#        self.positions[14]=105.0
#        self.positions[15]=112.2
        self.positions[0]=3.85
        self.positions[1]=11.05
        self.positions[2]=18.25
        self.positions[3]=25.45
        self.positions[4]=32.65
        self.positions[5]=39.85
        self.positions[6]=47.05
        self.positions[7]=54.25
        self.positions[8]=61.75
        self.positions[9]=68.95
        self.positions[10]=76.15
        self.positions[11]=83.35
        self.positions[12]=90.55
        self.positions[13]=97.75
        self.positions[14]=104.95
        self.positions[15]=112.15
#        for i in range(16):
#            self.positions[i]=3.9+i*7.2
        self.xdata=(self.positions[0:5]-3.85)
#        self.pixelscale=256.0/57.6
        self.pixelscale=512.0/116.0
        self.num=0
        self.pslimit_max=8000
        self.pslimit_min=500
    def killme(self):
        self.killstream=1
    def run(self):
        global myfile
        global fileopen
        global gcounts
        global global_acquireflag
        numin=0
        connected=0
        resyncs=0
        print( "starting stream thread %d\n ", self.socketID)
        streamsock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        streamsock.settimeout(2)
        try:
            streamsock.connect(self.socketID)
            connected=1
        except socket.error as e:
            print('connection error %s', e.args[0])
            print('you must restart program to obtain a connection')
        while connected==1:
            try:
                data=streamsock.recv(24)  #get header only
                numin=numin+1
                totalin=len(data)
                while totalin < 24:  #seems to occur when there is no data just rtdl header.
                    print('second read for header only %d bytes received' % totalin)
                    data2=streamsock.recv(24-totalin)
                    data=data+data2
                    totalin=totalin+len(data2)
            except socket.timeout as e:
                err=e.args[0]
                if err=='timed out':
                    if self.killstream==1:
                        print ('terminating due to killstream in socket time out')
                        break;
                    else:
                        continue
                elif self.killstream==1:
                    print ('terminating due to killstream')
                    break;
            except socket.error as e:
                print ("socket error %s on thread id =%d closed", e.args[0], self.threadID)
                break;
            if len(data)==0 or self.killstream==1:
                if len(data)==0:
                    print ('data length zero terminating stream thread')
                else:
                    print ('terminating due to killstream')
                print ("board on thread id = %d closed" , self.threadID)
                break;
            if fileopen==1:
                myfile.write(data)
                    #get the rest of the data...
            ival=numpy.zeros(6)
            ival=numpy.frombuffer(data,dtype=numpy.uint32)
            if ival[1] != 0x1515C2C2:
                #or there any 0x1515c2c2?
                if resyncs < 5:
                    print('attempting resync on thread %d' %self.threadID)
                    resyncs=resyncs+1
                found=0
                for i in range(6):
                    if ival[i]==0x1515C2C2:
                        found=1
                        break
                if found:
                    if i==0:
                        ival[5]=ival[4]
                        ival[4]=ival[3]
                        ival[3]=ival[2]
                        ival[2]=ival[2]
                        ival[1]=ival[0]
                        ival[0]=0
                    else:
                        k=i-1
                        newdata=streamsock.recv(4*k)
                        if k==1:
                            ival[0]=ival[1]
                            ival[1]=ival[2]
                            ival[2]=ival[3]
                            ival[3]=ival[4]
                            ival[4]=ival[5]
                            ival[5]=numpy.frombuffer(newdata,dtype=numpy.uint32)
                        elif k==2:
                            ival[0]=ival[2]
                            ival[1]=ival[3]
                            ival[2]=ival[4]
                            ival[3]=ival[5]
                            ival[4:]=numpy.frombuffer(newdata,dtype=numpy.uint32)
                        elif k==3:
                            ival[0]=ival[3]
                            ival[1]=ival[4]
                            ival[2]=ival[5]
                            ival[3:]=numpy.frombuffer(newdata,dtype=numpy.uint32)
                        else:
                            ival[0]=ival[4]
                            ival[1]=ival[5]
                            ival[2:]=numpy.frombuffer(newdata,dtype=numpy.uint32)
                                
                else:
                    print ('sync error in stream data')
                    continue
                #tbd call resync function
            resyncs=0
            if ival[3] != 0:
                try:
                    payload=streamsock.recv(ival[3])  #get the rest of the data,put in loop because may not recieve all
                    totalin=len(payload)
                    while  (totalin != ival[3]):  #tbd make while with check for timeout etc.
#                        print("short payload received from %s = %d, expected = %d" % (self.socketID[0],totalin,ival[3]))
                        payload2=streamsock.recv(ival[3]-totalin)
                        payload = payload + payload2
                        totalin=totalin+len(payload2)
                except socket.timeout as e:
                    err=e.args[0]
                    if err=='timed out':
                        if self.killstream==1:
                            print ('terminating stream thread at socket timeout')
                            break;
                        else:
                            continue
                    elif self.killstream==1:
                        print ('terminating due to killstream in payload request')
                        break;
                except socket.error as e:
                    print ("socket error %s on thread id =%d closed", e.args[0], self.threadID)
                    break;
                if len(payload)==0 or self.killstream==1:
                    if len(data)==0:
                        print ('data length zero terminating stream thread in payload request')
                    else:
                        print ('terminating due to killstream in payload request')
                    print ("board on thread id = %d closed" , self.threadID)
                    break;
    #            if numin==3600:
    #                now=datetime.datetime.now()
    #                numin=0
    #                print(str(now))
                if fileopen==1:
                    myfile.write(payload)
    #now deal with the payload.  
                if global_acquireflag==1:   
                    self.flushme=self.flushme+1
                    self.handledata(payload,ival) 
                else:
                    self.flushme=0
            
        streamsock.close()
        print ("stream thread id =%d closed\n"  %  self.threadID)
        if self.saveraw==1:
            self.rawfile.close()
        
    def setsavelightcone(self,writefilename):
        self.savelc=1
        self.myfile=open(writefilename,'w+t')
        self.myfile2=open('C:/Users/Larmor/Desktop/AngerCamera/altlightconedata_gauss.txt','w+t')
        
    def stopsavelightcone(self):
        self.savelc=0
        self.myfile.close()
        self.myfile2.close()
        
    def setpslimits(self, maxlim, minlim):
        self.pslimit_max=maxlim
        self.pslimit_min=minlim
        
    def setcaloffsets(self,xoffsets,yoffsets):
        self.xcaloffsets=xoffsets
        self.ycaloffsets=yoffsets
        self.xcalnew=self.xcaloffsets[0:16]+self.xcaloffsets[16:32]
        self.ycalnew[0:8]=self.ycaloffsets[0:8]+self.ycaloffsets[8:16]
        self.ycalnew[8:16]=self.ycaloffsets[16:24]+self.ycaloffsets[24:32]
        
    def loadtables(self):
        rf='C:/Users/Larmor/Desktop/AngerCamera/ptable_max.dat'        
        z=numpy.fromfile(rf,dtype=numpy.float32)
#        print('max10,100=',z[10*512+100])
        self.psummax=z.reshape(512,512)
#        print(self.psummax[10,100])
        rf='C:/Users/Larmor/Desktop/AngerCamera/ptable_min.dat'        
        z=numpy.fromfile(rf,dtype=numpy.float32)
#        print('min10,100=',z[10*512+100])
        self.psummin=z.reshape(512,512)
        self.tablesin=1
        print('psum tables read')
#        Zswap=numpy.swapaxes(self.psummin,0,1)
#        ax = pylab.subplot(111)
#        im = pylab.imshow(Zswap, cmap=pylab.cm.jet, origin='lower') #, vmin=0, vmax=400)
#        cbar=pylab.colorbar()
#        cs='flood'
#        pylab.title(cs)
#        pylab.show()
        
    def saverawpos(self,value):
        if value==1:
            self.saveraw=1
            self.rawfile=open('C:/Users/Larmor/Desktop/AngerCamera/cg4btest.dat','w+b')
            print('raw write file open')
        else:
            if self.saveraw==1:
                self.rawfile.close()
                print('raw write file closed')
            self.saveraw=0
        
        
    def handledata(self,payload,header_unpacked):
        #currently this only works with short_pid header and raw data.
        #tbd in normal mode.
        global ad_histox
        global ad_histoy
        global psum_histox
        global psum_histoy
        global tof_histo
        global tof_histo_fine
        global zero_histos
        global zero_histos_toggle
        global bm1counts
        global gcounts
        global psum_scale
        global ad_scale
        global twoD_histo
        global printonce
        global firstin
        global firstten
        global psum3Dhisto
        
        invscale=1.0/psum_scale
        invscale_ad=1.0/ad_scale
        
        if (zero_histos==1 and zero_histos_toggle==0):
            zero_histos_toggle=1
            psum_histox.fill(0)
            psum_histoy.fill(0)
            ad_histox.fill(0)
            ad_histoy.fill(0)
            twoD_histo.fill(0)
            psum3Dhisto.fill(0)
            tof_histo.fill(0)
            gcounts=0
        if (zero_histos==0):
            zero_histos_toggle=0
        
        bm1counts=header_unpacked[4]
        rtdl_size=struct.calcsize('6I')
        ts_nsec,ts_sec,pcharge,gen_info,tsyn_period,tsyn_delay=struct.unpack('6I',payload[0:rtdl_size])
#        if (printonce < 10):
#            print('%x,%d, %d' % (header_unpacked[2],header_unpacked[4],header_unpacked[5]))
#            printonce=printonce+1
        
        k=len(payload)
#tbd below is to get the tof and pid index number.
        if (header_unpacked[2] & 0xF000000)==0x2000000:   #this is normal mode
            for j in range(24,k,self.NORMAL_EVT_SZ):
                tof,pixelid=struct.unpack(self.NORMAL_EVT_FMT,payload[j:j+self.NORMAL_EVT_SZ])
#                if (printonce < 10):
#                    print('tof=%X' % tof)
#                    printonce = printonce + 1
                ii=int(tof*.01)
                if ii < 1700:
                    tof_histo[ii]=tof_histo[ii]+1
                ii=int(tof*0.1)
                if ii < 2000:
                    tof_histo_fine[ii]=tof_histo_fine[ii]+1
                iypos=pixelid & 0x1FF  # doesn't work if there is a pixeloffset..net to be able to read that.
                ixpos=(pixelid >> 9) & 0x1FF
                gcounts=gcounts+1
#                if (ixpos >= 0 and ixpos <= 511 and iypos >=0 and iypos <= 511):
                twoD_histo[ixpos,iypos]=twoD_histo[ixpos,iypos]+1
                if self.saveraw==1:
                    kpack=struct.pack('ddd',tof,self.lengthscale*ixpos,self.lengthscale*iypos)
                    self.rawfile.write(kpack)
                    if (self.flushme >= 60):
                        self.rawfile.flush()
                        self.flushme=0

        else:
            for j in range(24,k,self.RAW_EVT_SZ):
                iupacked=struct.unpack(self.RAW_EVT_FMT,payload[j:j+self.RAW_EVT_SZ])                
                tof=0
                for i in range(8):
                    self.xs[i,1]=iupacked[i] & 0xfff                
                    ii=(iupacked[i] & 0x0fff0000)>>16
                    tof = tof >> 4
                    tof=tof | (iupacked[i] & 0xF0000000)
                    self.xs[i,0]=ii
                for i in range(8,16,1):
                    self.xs[i-8,3]=iupacked[i] & 0xfff
                    ii=(iupacked[i] & 0x0fff0000)>>16
                    self.xs[i-8,2]=ii  
                for i in range(16,24,1):
                    self.ys[i-16,2]=iupacked[i] & 0xfff
                    ii=(iupacked[i] & 0x0fff0000)>>16
                    self.ys[i-16,0]=ii
                for i in range(24,32,1):
                    self.ys[i-24,3]=iupacked[i] & 0xfff
                    ii=(iupacked[i] & 0x0fff0000)>>16
                    self.ys[i-24,1]=ii
#                if (printonce < 50):
#                    printonce=printonce+1                        
#                    print('%X,' % tof)
            
            #i know only tube 1 (index 0)works for now so I only use that one
    #for now until summer board is fixed, x3,x4 and x5,x6 are swapped... 
#subtract cal offsets

                self.xnew[0:8]=(self.xs[:,0]+self.xs[:,2])
                self.xnew[8:]=(self.xs[:,1]+self.xs[:,3])
                self.ynew[0:8]=(self.ys[:,0]+self.ys[:,1])
                self.ynew[8:]=(self.ys[:,2]+self.ys[:,3])
                self.xnew=self.xnew-self.xcalnew
                self.ynew=self.ynew-self.ycalnew
#for testing only                firstten[firstin,:]=self.xnew
    ######################################################
#                t3=self.xnew[3]
#                t2=self.xnew[2]
#                self.xnew[3]=self.xnew[5]
#                self.xnew[2]=self.xnew[4]
#                self.xnew[5]=t3
#                self.xnew[4]=t2
#                t3=self.xnew[11]
#                t2=self.xnew[10]
#                self.xnew[11]=self.xnew[13]
#                self.xnew[10]=self.xnew[12]
#                self.xnew[13]=t3
#                self.xnew[12]=t2
    ######################################################
                ii=int(tof*.01)
                if ii < 1700:
                    tof_histo[ii]=tof_histo[ii]+1
                ii=int(tof*.1)
                if ii < 2000:
                    tof_histo_fine[ii]=tof_histo_fine[ii]+1
                maxo=numpy.argmax(self.xnew)
                mini=maxo-5
                maxi=maxo+5
                if mini < 2:
                    back=0.25*(numpy.sum(self.xnew[maxi:maxi+4]))
                elif maxi > 13:
                    back=0.25*(numpy.sum(self.xnew[mini-3:mini+1]))
                else:
                    back=0.25*(self.xnew[mini]+self.xnew[mini-1]+self.xnew[maxi]+self.xnew[maxi+1])
                
                mini=maxo-3
                maxi=maxo+3
                if maxi > 16:
                    maxi=16
                if mini < 0:
                    mini=0
                xback=self.xnew-back
                keepps=numpy.sum(self.xnew[mini:maxi])
                if keepps > self.pslimit_max or keepps < self.pslimit_min:
                    continue
                
                dpsum=int(keepps/20.0)
                
                ps=int(invscale*keepps)
                maxo=numpy.argmax(self.xnew)
                

                if (maxo==0): 
#                    xback[:3]=xback[:3]-back
                    gpar=analysisvar(self.positions[:3],xback[:3])
                elif (maxo==1):
                    gpar=analysisvar(self.positions[:4],xback[:4])
                elif (maxo==15):
#                    xback[13:]=xback[13:]-back
                    gpar=analysisvar(self.positions[13:],xback[13:])
                elif (maxo==14):
                    gpar=analysisvar(self.positions[12:],xback[12:])
                else:
                    gpar=analysisvar(self.positions[maxo-2:maxo+3],xback[maxo-2:maxo+3])


                xpos=gpar[1]
                
#                if (maxo==0 or maxo==1):
#                    ydata=xback[:5]
##                    mypos=self.positions[0]
#                    self.xdata=self.positions[:5]
#                    gpar=analysisvar(self.xdata,ydata)
#                    xpos=gpar[1] #+mypos
#                elif (maxo==14 or maxo==15):
#                    ydata=xback[11:]
#                    self.xdata=self.positions[11:]
##                    mypos=self.positions[11]
#                    gpar=analysisvar(self.xdata,ydata)
#                    xpos=gpar[1] #+mypos
#                else:
#                    ydata=xback[maxo-2:maxo+3]
#                    self.xdata=self.positions[maxo-2:maxo+3]
#                    gpar=analysisvar(self.xdata,ydata)
#                    xpos=gpar[1] #+mypos
                    
#for testing only                    
#                if (firstin < 10):
#                    print('first ten' ,xpos, gpar[0])
#                    firstin=firstin+1
                  
    #            gpar=analysisvar(self.xdata,ydata)
    #            xpos=gpar[1]+mypos
                
                    
    #            if (maxo==3):
    #                back=0.5*(self.xs[6,0]+self.xs[7,0])
    #                xback=self.xs[:,0]-back
    #                ydata=numpy.zeros(5)
    #                ydata=xback[1:6]
    #                gpar=analysisvar(self.xdata,ydata)
    #                if (gpar[1] >= 17.9 and gpar[1] <= 18.1 and self.savelc==1):
    #                    cs='%d,%d,%d,%d,%d,%d,%d,%d\n' %(self.xs[0,0],self.xs[1,0],self.xs[2,0],self.xs[3,0],self.xs[4,0],self.xs[5,0],self.xs[6,0],self.xs[7,0])
    #                    self.myfile.write(cs)
    #            elif (maxo==4):
    #                back=0.5*(self.xs[0,0]+self.xs[1,0])
    #                xback=self.xs[:,0]-back
    #                ydata=numpy.zeros(5)
    #                ydata=xback[2:7]
    #                gpar=analysisvar(self.xdata,ydata)
    #                if (gpar[1] >= 10.7 and gpar[1] <= 10.9 and self.savelc==1):
    #                    cs='%d,%d,%d,%d,%d,%d,%d,%d\n' %(self.xs[0,0],self.xs[1,0],self.xs[2,0],self.xs[3,0],self.xs[4,0],self.xs[5,0],self.xs[6,0],self.xs[7,0])
    #                    self.myfile.write(cs)
                    
    
                adps=int(invscale_ad*self.xnew[maxo])
                if adps >= 500:
                    adps=499
                if ps >= 1000:
                    ps=999
                if (ps<1000 and ps >=0):  #tbd make 1000 a numonic defined above
                    psum_histox[ps,maxo]=psum_histox[ps,maxo]+1 
                    
                maxox=maxo
                    
    #debug of x1
    #            adpx=int(invscale_ad*self.xs[0,0])
                
    #            ad_histox[adpx,0,0]=ad_histox[adpx,0,0]+1
                
                ps=int(invscale*numpy.sum(self.ynew))
                maxo=numpy.argmax(self.ynew)
                adpsy=int(invscale_ad*self.ynew[maxo])
                if (adpsy >= 500):  
                    adpsy=499
                if ps >= 1000:
                    ps=999
                if (ps<1000 and ps >=0):  #tbd make 1000 a numonic defined above
                    psum_histoy[ps,maxo]=psum_histoy[ps,maxo]+1 
                
                maxoy=maxo
                if maxoy<8 and maxox <8:                    
                    ad_histoy[adpsy,maxoy,0]=ad_histoy[adpsy,maxoy,0]+1
                    ad_histox[adps,maxox,0]=ad_histox[adps,maxox,0]+1
                elif maxoy<8 and maxox >= 8:                    
                    ad_histoy[adpsy,maxoy,1]=ad_histoy[adpsy,maxoy,1]+1
                    ad_histox[adps,maxox-8,1]=ad_histox[adps,maxox-8,1]+1
                elif maxoy >= 8 and maxox <8:                    
                    ad_histoy[adpsy,maxoy-8,2]=ad_histoy[adpsy,maxoy-8,2]+1
                    ad_histox[adps,maxox,2]=ad_histox[adps,maxox,2]+1
                else:                    
                    ad_histoy[adpsy,maxoy-8,3]=ad_histoy[adpsy,maxoy-8,3]+1
                    ad_histox[adps,maxox-8,3]=ad_histox[adps,maxox-8,3]+1
                
                
                
    #get background estimate

                maxo=numpy.argmax(self.ynew)
                mini=maxo-5
                maxi=maxo+5
                if mini < 2:
                    back=0.25*(numpy.sum(self.ynew[maxi:maxi+4]))
                elif maxi > 13:
                    back=0.25*(numpy.sum(self.ynew[mini-3:mini+1]))
                else:
                    back=0.25*(self.ynew[mini]+self.ynew[mini-1]+self.ynew[maxi]+self.ynew[maxi+1])
                
                mini=maxo-3
                maxi=maxo+4
                if maxi > 16:
                    maxi=16
                if mini < 0:
                    mini=0
                yback=self.ynew-back





                
                if (maxo==0): 
#                    yback[:3]=yback[:3]-back
                    gpar=analysisvar(self.positions[:3],yback[:3])
                elif (maxo==1):
                    gpar=analysisvar(self.positions[:4],yback[:4])
                elif (maxo==15):
#                    yback[13:]=yback[13:]-back
                    gpar=analysisvar(self.positions[13:],yback[13:])
                elif (maxo==14):
                    gpar=analysisvar(self.positions[12:],yback[12:])
                else:
                    gpar=analysisvar(self.positions[maxo-2:maxo+3],yback[maxo-2:maxo+3])
#                    ydata=yback[maxo-2:maxo+3]
#                    self.xdata=self.positions[maxo-2:maxo+3]
                    
#                gpar=analysisvar(self.xdata,ydata)
                ypos=gpar[1] #+mypos
                ixpos=int(xpos*self.pixelscale+0.5)
                iypos=int(ypos*self.pixelscale+0.5)
        
                if (ixpos >= 0 and ixpos <= 511 and iypos >=0 and iypos <= 511):
                    psum3Dhisto[ixpos,iypos,dpsum]=psum3Dhisto[ixpos,iypos,dpsum]+1
                    if (self.savelc==1 and ixpos==368 and iypos > 270  and iypos < 520):
                        cs=''
                        for kk in range(4):
                            for k in range(8):
                                cs2='%d\t' % self.xs[k,kk]  #note that 3456 are not swapped!
                                cs=cs+cs2
                        for kk in range(4):
                            for k in range(8):
                                cs2='%d\t' % self.ys[k,kk]
                                cs=cs+cs2
                        cs=cs + '\n'
                        self.myfile.write(cs)
                        cs=''
                        for k in range(16):
                            cs2='%f\t' % self.xnew[k]
                            cs=cs+cs2
                        cs=cs +'\n'
                        self.myfile2.write(cs)
                    if self.tablesin==1:
#                        if printonce < 20:
#                            print(keepps,self.psummax[ixpos,iypos],self.psummin[ixpos,iypos])
#                            printonce=printonce+1
                        if keepps <= self.psummax[ixpos,iypos] and keepps >= self.psummin[ixpos,iypos]:
                            gcounts=gcounts+1
                            if self.saveraw==1:
                                kpack=struct.pack('ddd',tof,xpos,ypos)
                                self.rawfile.write(kpack)
                                if (self.flushme >= 60):
                                    self.rawfile.flush()
                                    self.flushme=0
                            twoD_histo[ixpos,iypos]=twoD_histo[ixpos,iypos]+1
                    else:
                        gcounts=gcounts+1
                        twoD_histo[ixpos,iypos]=twoD_histo[ixpos,iypos]+1
                        if self.saveraw==1:
                            kpack=struct.pack('ddd',tof,xpos,ypos)
                            self.rawfile.write(kpack)

            
            


class MySiPMComm:
    def __init__(self,messagebox):
        cmdserver_addr=('192.168.10.19',9119)
        self.cmdsock1 = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        messagebox.setText('attempting to connect to siPM at 10.19')
        self.cmdsock1.connect(cmdserver_addr)
        self.cmdsock1.settimeout(1)
        self.numprints=0
        messagebox.setText('connected')
        while(1):
            try:
                rsp=self.cmdsock1.recv(24)
                if len(rsp)==24:
                    print('crap sitting in comm buffer')
            except socket.timeout:
                print('clean comm buffer')
                break
            
    def __del__(self):
        self.cmdsock1.close()
    def setzeros(self,dacnum,value):
        #note, value is a sixteen array or tuple.
        #x zeros are
        # 1=num 1, chan 1..note chan zero sets all the channels to the same.
        # 2=num 1, chan 4
        # 3=num 1, chan 3
        # 4=num 1, chan 2
        # 5=num 0, chan 1
        # 6=num 0, chan 4
        # 7=num 0, chan 3
        # 8=num 0, chan 2
        # y zeros same chan sequence except 1-4 = num 2
        # y5-8= num 3
        addr=[1,1,1,1,0,0,0,0,2,2,2,2,3,3,3,3]
        chan=[1,4,3,2,1,4,3,2,1,4,3,2,1,4,3,2]
        cmd =10 
        numbytes=16
        for k in range (16):
            ivalue=(int)(value[k])
            data=struct.pack('10I',0,0,cmd,numbytes,0,0,dacnum,addr[k],chan[k],ivalue)
            self.cmdsock1.send(data)
            try:
                rsp=self.cmdsock1.recv(24)
                ival=struct.unpack('6I',rsp)
#                if k==0:
#                    print (ival[4:5])
#                return ival[4:5]  need to check to see if ack....later....
            except socket.timeout:
                print('zero timeout k=%d,dacnum=%d' % (k,dacnum))
#            time.sleep(.02)
            
    def setoffsets(self,dacnum,value):
        #note, value is a sixteen array or tuple.
        #x zeros are
        # 1=num 1, chan 5..note chan zero sets all the channels to the same.
        # 2=num 1, chan 8
        # 3=num 1, chan 6
        # 4=num 1, chan 7
        # 5=num 0, chan 5
        # 6=num 0, chan 8
        # 7=num 0, chan 6
        # 8=num 0, chan 7
        # y zeros same chan sequence except 1-4 = num 2
        # y5-8= num 3
        addr=[1,1,1,1,0,0,0,0,2,2,2,2,3,3,3,3]
        chan=[5,8,6,7,5,8,6,7,5,8,6,7,5,8,6,7]
        cmd =10 
        numbytes=16
        print('setoffset')
        print(value)
        for k in range (16):
            ivalue=(int)(value[k])
            data=struct.pack('10I',0,0,cmd,numbytes,0,0,dacnum,addr[k],chan[k],ivalue)
            self.cmdsock1.send(data)
            try:
                rsp=self.cmdsock1.recv(24)
                ival=struct.unpack('6I',rsp)
#                return ival[4:5]
            except socket.timeout:
                print('offset timeout k=%d,dacnum=%d' % (k,dacnum))
 #           time.sleep(.02)
        
    def setdac(self,dacnum,addr,chan,value):
        ivalue=(int)(value)
        cmd=10
        numbytes=16
        data=struct.pack('10I',0,0,cmd,numbytes,0,0,dacnum,addr,chan,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
#            return ival[4:5]
        except socket.timeout:
            print('timeout in setdac')
            
        
    def readdisc(self):
        cmd=12
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(24)
        ival=struct.unpack('6I',rsp)
        return ival[4:6]
    
    def readcipi(self):
        cmd=13
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(24)
        ival=struct.unpack('6I',rsp)
        return ival[4:6]
    
    def readtemperature(self):
        cmd=35
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(24)
        ival=struct.unpack('6I',rsp)
#        return ival[4]
        if ival[4] & 0x2000:  # means negitive temp.
            tempvalue=((ival[4] & 0x1FFF)- 8192)/32.0
        else:
            tempvalue=ival[4]/32.0
        return tempvalue


    def readadcstatus(self):  #this also reads the counts status
        cmd=20
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(40)
        ival=struct.unpack('10I',rsp)
        return ival[6:10]
    
    def readmem(self,ci):
        cmd=14
        numbytes=4
        ivalue=(int)(ci)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(152)
        bytesin=len(rsp)
        if bytesin < 152:
            rsp2=self.cmdsock1.recv(152-bytesin)
            rsp=rsp+rsp2
        ival=struct.unpack('38I',rsp)
        return ival[6:38]
    
    def testadc(self):
        cmd=15
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        rsp=self.cmdsock1.recv(24)
        ival=struct.unpack('6I',rsp)
        return ival[4:5]
    
    def initadc(self):
        cmd=16
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('initadc error')
    
    def endisc(self,value):
        cmd=17
        numbytes=4
        ivalue = (int)(value)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            if (ival[2] != 17):
                return -1
            else:
                return 0
        except socket.timeout:  #make an aditional retry here
            print('endisc error')
        
    def setinttime(self,value):
        cmd=18
        numbytes=4
        ivalue = (int)(value)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('set inittimec error')
        
    def enfaketrig(self,value):
        cmd=19
        numbytes=4
        ivalue = (int)(value)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('en faketrig error')
            
    def clrstatuscounts(self):
        cmd=37
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('clr stauts counts error')
            
            
    def readparams(self):
        cmd=21
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            plbytes=int(ival[3])
#            print('read parms bytes' , plbytes)
#            print ('read params payload size=%d' % plbytes)
            rsp=self.cmdsock1.recv(plbytes)
            inbytes=len(rsp)
            while inbytes < plbytes:
                rsp2=self.cmdsock1.recv(plbytes-inbytes)
                inbytes=inbytes+len(rsp2)
                rsp=rsp+rsp2
            return rsp
        except socket.timeout:  #make an aditional retry here
            rsp=' '
            print('response timeout in readparams')
            
    def writeparams(self,payload):  #payload is expected to be prepaked, payload_size is number of bytes.
        cmd=22
        payload_size=len(payload)
        headerdata=struct.pack('6I',0,0,cmd,payload_size,0,0)
        data=headerdata+payload
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            time.sleep(0.2)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('no ack/nack on writeparams')
#to use remove comments and return -1  and restore global firstten....          
    def caltest(self,index):
        return -1
#        global firstten
#        cmd=23
#        numbytes=16*4
#        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
#        data2=firstten[index,:].tobytes()
#        data=data+data2
#        self.cmdsock1.send(data)
#        try:
#            rsp=self.cmdsock1.recv(36)
#            numbytes=len(rsp)
#            if (numbytes < 36):
#                rsp2=self.cmdsock1.recv(36-numbytes)
#                rsp=rsp+rsp2
#            ival=struct.unpack('6I3f',rsp)
#            if (ival[2] != 23):
#                return -1
#            else:
#                return ival[6:9]
#        except socket.timeout:  #make an aditional retry here
#            print('en ADC intr error')
#            return -2
        
        
        
        
    def enADCintr(self,value):
        cmd=25
        numbytes=4
        ivalue = (int)(value)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            if (ival[2] != 25):
                return -1
            else:
                return 0
        except socket.timeout:  #make an aditional retry here
            print('en ADC intr error')
            return -2
    
    def setrawmode(self,value):
        cmd=26
        numbytes=4
        ivalue = (int)(value)
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            if (ival[2] != 26):
                return -1
            else:
                return 0
        except socket.timeout:  #make an aditional retry here
            print('set raw mode error')
            return -2
            
    def writeRTDL(self,reg0,reg1):
        cmd=30
        numbytes=8
        data=struct.pack('8I',0,0,cmd,numbytes,0,0,reg0,reg1)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('response timeout in writertdl')
    
    def readRTDL(self):
        cmd=29
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            plbytes=int(ival[3])
            rsp=self.cmdsock1.recv(plbytes)
            return struct.unpack('8I',rsp)
        except socket.timeout:  #make an aditional retry here
            print('response timeout in readparams')
        
    def readRTDLmem(self):
        cmd=31
        numbytes=4
        ivalue=16
        data=struct.pack('7I',0,0,cmd,numbytes,0,0,ivalue)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            plbytes=int(ival[3])
            rsp=self.cmdsock1.recv(plbytes)
            return struct.unpack('16I',rsp)
        except socket.timeout:  #make an aditional retry here
            print('response timeout in readparams')
        

        
    def killthread(self):
        cmd=100
        numbytes=0
        data=struct.pack('6I',0,0,cmd,numbytes,0,0)
        self.cmdsock1.send(data)
        try:
            rsp=self.cmdsock1.recv(24)
            ival=struct.unpack('6I',rsp)
            return ival[4:5]
        except socket.timeout:  #make an aditional retry here
            print('en kill thread error')
            
        
#############read values for all tube at given int time...
    def startaquire(self,fake,inttime):
        self.clrstatuscounts()
        self.setinttime(inttime)
        if (fake):
            self.endisc(0)  #disable all regular disc....
            self.enfaketrig(1)        
        else:
#            self.endisc(0x70)  #for single tube in slot 0
            self.endisc(0x7F)
            
    def stopaquire(self):
        self.endisc(0)
        self.enfaketrig(0)

#############read values for all tube at given int time...            
    def set_and_read(self,ticks):
        global debugon
        global plt1
        if ticks <0:
            ticks=0
        self.setinttime(ticks)        
        resp=self.readcipi()
        ci=resp[0]
        pi=resp[1]
        self.enfaketrig(1)
        time.sleep(.5)
        self.enfaketrig(0)
        resp=self.readcipi()
        ci2=resp[0]
        pi2=resp[1]
        
#        print( 'pi=%d,%d,pi2=%d,%d' %(pi,ci,pi2,ci2))
        
        xmeans=numpy.zeros((8,4))
        xstds=numpy.zeros((8,4))
        ymeans=numpy.zeros((8,4))
        ystds=numpy.zeros((8,4))
    
        if pi2 < pi:
            tlen=4096-pi+pi2
        else:
            tlen=pi2-pi
            
        xs=numpy.zeros((8,4,tlen))
        ys=numpy.zeros((8,4,tlen))
    
        
        if pi2 < pi:
            myadder=4096-pi
            for k in range(pi,4096,1):
                resp=self.readmem(k)
                for i in range(8):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i,1,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i,0,k-pi]=j
        #            print 'tube 1 x%d=%d' % (i+1,j)
                for i in range(8,16,1):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i-8,3,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i-8,2,k-pi]=j
                    
                for i in range(16,24,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-16,2,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-16,0,k-pi]=j
        #            print 'tube 1 y%d=%d' % (i-15,j)
                for i in range(24,32,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-24,3,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-24,1,k-pi]=j
                    
            for k in range(0,pi2,1):
                resp=self.readmem(k)
                for i in range(8):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i,1,myadder+k]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i,0,myadder+k]=j
        #            print 'tube 1 x%d=%d' % (i+1,j)
                for i in range(8,16,1):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i-8,3,myadder+k]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i-8,2,myadder+k]=j
                    
                for i in range(16,24,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-16,2,myadder+k]= j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-16,0, myadder+k]=j
        #            print 'tube 1 y%d=%d' % (i-15,j)
        
                for i in range(24,32,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-24,3,myadder+k]= j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-24,1, myadder+k]=j
    
        else:
            for k in range(pi,pi2,1):
                resp=self.readmem(k)
                for i in range(8):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i,1,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i,0,k-pi]=j
        #            print 'tube 1 x%d=%d' % (i+1,j)
                for i in range(8,16,1):
        ##        print '%x' % (uint32)(resp[i])
                    j=(numpy.uint32)(resp[i])
                    xs[i-8,3,k-pi]=j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    xs[i-8,2,k-pi]=j
                    
                for i in range(16,24,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-16,2,k-pi]= j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-16,0,k-pi]=j
        #            print 'tube 1 y%d=%d' % (i-15,j)
                for i in range(24,32,1):
                    j=(numpy.uint32)(resp[i])
                    ys[i-24,3,k-pi]= j & 0xfff
                    j=(j & 0x0fff0000)>>16
                    ys[i-24,1,k-pi]=j
    #    print 'for int time=100'
    ##    plot(xs[0,:])
    ##    show()
    #    print 'for int time=100'
    ##    plot(xs[0,:])
    ##    show()
#        myfile=open('/home/rre/Documents/tempdata/xdata.txt','w+t')
#        for j in range(tlen):
#            cs='%d,%d,%d,%d,%d,%d,%d,%d\n' %(xs[0,0,j],xs[1,0,j],xs[2,0,j],xs[3,0,j],xs[4,0,j],xs[5,0,j],xs[6,0,j],xs[7,0,j])
#            myfile.write(cs)
#        myfile.close()
#        myfile=open('/home/rre/Documents/tempdata/ydata.txt','w+t')
#        for j in range(tlen):
#            cs='%d,%d,%d,%d,%d,%d,%d,%d\n' %(ys[0,0,j],ys[1,0,j],ys[2,0,j],ys[3,0,j],ys[4,0,j],ys[5,0,j],ys[6,0,j],ys[7,0,j])
#            myfile.write(cs)
#        myfile.close()
        for k in range (4):
            for i in range (8):
                xstds[i,k]=numpy.std(xs[i,k,:])
                xmeans[i,k]=numpy.average(xs[i,k,:])
#                print('x%dmean=%f,x%dstd=%f', i+1,mymean,i+1,mysd)
            for i in range (8):
                ystds[i,k]=numpy.std(ys[i,k,:])
                ymeans[i,k]=numpy.average(ys[i,k,:])
#                print('y%dmean=%f,y%dstd=%f',i+1,mymean,i+1,mysd)
        return xmeans,ymeans,xstds,ystds    

        

    def lowlevel(self,tube,orgint_time):        
        myoffset=0
        offsets=[myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset,myoffset]
        if tube < 1 or tube > 4:
            return offsets
        
        if tube==1:
            dacnum=1
        elif tube==2:
            dacnum=2
        elif tube==3:
            dacnum=3
        elif tube==4:
            dacnum=4
            
        self.endisc(0)
        zeroval=4095
        zerovalues=[zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval,zeroval]
        self.setzeros(dacnum,zerovalues)
        for i in range(16):
            offsets[i]=2048
        self.setoffsets(dacnum,offsets)
        xmean1,ymean1,xstd1,ystd1=self.set_and_read(5)
        xmean2,ymean2,xstd2,ystd2=self.set_and_read(40)
    ###full range of zeros is to add between 0 to 76 constant A/D value
        #change offsets and redo.
        for i in range(16):
            offsets[i]=1648
        self.setoffsets(dacnum,offsets)
        xmean3,ymean3,xstd3,ystd3=self.set_and_read(5)
        xmean4,ymean4,xstd4,ystd4=self.set_and_read(40)
        xslopes_1=numpy.zeros(8)
        yslopes_1=numpy.zeros(8)
        xintercepts_1=numpy.zeros(8)
        yintercepts_1=numpy.zeros(8)
        for i in range (8):
            xslopes_1[i]=(xmean2[i,tube-1]-xmean1[i,tube-1])/35.0
            xintercepts_1[i]=xmean2[i,tube-1]-xslopes_1[i]*40.0
        for i in range (8):
            yslopes_1[i]=(ymean2[i,tube-1]-ymean1[i,tube-1])/35.0
            yintercepts_1[i]=ymean2[i,tube-1]-yslopes_1[i]*40.0
            
        xslopes_2=numpy.zeros(8)
        yslopes_2=numpy.zeros(8)
        xintercepts_2=numpy.zeros(8)
        yintercepts_2=numpy.zeros(8)
        for i in range (8):
            xslopes_2[i]=(xmean4[i,tube-1]-xmean3[i,tube-1])/35.0
            xintercepts_2[i]=xmean4[i,tube-1]-xslopes_2[i]*40.0
        for i in range (8):
            yslopes_2[i]=(ymean4[i,tube-1]-ymean3[i,tube-1])/35.0
            yintercepts_2[i]=ymean4[i,tube-1]-yslopes_2[i]*40.0
    
    
    #more checks later...............
    #now attempt to get slopes =0;
        dx=numpy.zeros(8)
        dy=numpy.zeros(8)
        for i in range(8):
            dx[i]=(-400.0*xslopes_1[i])/(xslopes_2[i]-xslopes_1[i])
            dy[i]=(-400.0*yslopes_1[i])/(yslopes_2[i]-yslopes_1[i])
            
    
    
        for i in range(8):
            offsets[i]=(int)(2048.0-dx[i])
            offsets[i+8]=(int)(2048.0-dy[i])
            zerovalues[i]=4095
            zerovalues[i+8]=4095
        
        self.setoffsets(dacnum,offsets)
        self.setzeros(dacnum,zerovalues)
        self.setinttime(orgint_time)
    #    s=raw_input('offsets=8000 hit a key')
        
        return offsets

        
        
class MyApp(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        global twoD_histo
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.instripgainlogic=0
        self.acquireflag=0  #set to one if something is aquiring
        self.enablefake=0
        self.prevcounts=0
        self.runningtime=QtCore.QTime()
#bunch of tube varibles from the detector...local copies....
        self.countsarray=numpy.zeros(3,dtype=numpy.uint32)
        self.statusarray=numpy.zeros(5,dtype=numpy.int32)
        self.discvalues=numpy.zeros(8,dtype=numpy.int32)
        self.xoffsets=numpy.zeros(32,dtype=numpy.int32)
        self.yoffsets=numpy.zeros(32,dtype=numpy.int32)
        self.xzeros=numpy.zeros(32,dtype=numpy.int32)
        self.yzeros=numpy.zeros(32,dtype=numpy.int32)
        self.cal_xoffsets=numpy.zeros(32,dtype=numpy.float32)
        self.cal_yoffsets=numpy.zeros(32,dtype=numpy.float32)
        self.cal_xstripgains=numpy.zeros(32,dtype=numpy.float32)
        self.cal_ystripgains=numpy.zeros(32,dtype=numpy.float32)
        self.tof_offset=0
        self.int_time=0
        self.disc_en=0
        self.en_flatfield=0
        self.en_photosum=0
        self.event_output_type=0
        self.rtdl_header_type=0
        self.psum_max=0.0
        self.psum_min=0.0
        self.alpha1=0.0
        self.alpha2=0.0
        self.weight_alpha1=0.0
        self.weight_alpha2=0.0
        self.whatplot=0  #default is psum, 1= a/d
        self.xpltstate=[1,1,1,1,1,1,1,1]
        self.ypltstate=[1,1,1,1,1,1,1,1]
        self.setupUi(self)
        self.MySiPM=MySiPMComm(self.messagebox)
        self.streamserver_addr = ('192.168.10.19',9229)
        self.mystream1=streamthread(1,self.streamserver_addr)
        self.mystream1.start()
        self.saveraw=0
        k=self.ReadPSumTables('192.168.10.19')
        if k==0:
            self.mystream1.loadtables()
        time.sleep(1)
        r=self.MySiPM.endisc(0)
        if r < 0:
            print('endisc failed')
        r=self.MySiPM.setrawmode(0)
        if r<0:
            print('setrawmode failed')
        r=self.MySiPM.enADCintr(1)
        if r < 0:
            print('enADCntr failed')
        rtdlregs=self.MySiPM.readRTDL();
        print ('rtdlreg0,1=%x,%x' % (rtdlregs[0],rtdlregs[1]))
        rtdlreg0=int(rtdlregs[0])
        rtdlreg1=int(rtdlregs[1])
#to enable external pulse sync
#        rtdlreg0=rtdlreg0 | 0xC0
#disable
        rtdlreg0=rtdlreg0 & 0xFFFFFF3F
        self.MySiPM.writeRTDL(rtdlreg0,rtdlreg1)
#to enable bm1 set 18, clear with bit 16
#        rtdlreg0=rtdlreg0 | 0x50000
#        self.MySiPM.writeRTDL(rtdlreg0,rtdlreg1)
#        rtdlreg0=rtdlreg0 & 0xFFFEFFFF
#        self.MySiPM.writeRTDL(rtdlreg0,rtdlreg1)
#        rtdlregs=self.MySiPM.readRTDL()
#        print ('%x,%x' % (rtdlregs[0],rtdlregs[1]))
###############################3
        self.win=pg.GraphicsWindow()
        self.win.resize(500,300)
        self.win.setWindowTitle('PSum histo x columns 1-8')
        self.plt1=self.win.addPlot()
#        self.plt2=self.win.addPlot()
#        self.plt3=self.win.addPlot()
#        self.plt4=self.win.addPlot()
        self.win2=pg.GraphicsWindow()
        self.win2.resize(500,300)
        self.win2.setWindowTitle('PSum histo y rows 1-8')
        self.plt2=self.win2.addPlot()


        self.win3=pg.GraphicsWindow()
        self.win3.resize(400,400)
        self.win3.setWindowTitle('2D histo')
        self.viewb=self.win3.addViewBox()
        self.viewb.setAspectLocked()
        self.img=pg.ImageItem(border='w')
        self.viewb.addItem(self.img)
        

#        self.plt2.addLegend()
#        self.plt6=self.win2.addPlot()
#        self.plt7=self.win2.addPlot()
#        self.plt8=self.win2.addPlot()
#        self.threeDwidget=gl.GLViewWidget()
#        X=numpy.linspace(0,256,256)

#        self.surfplot=gl.GLSurfacePlotItem(x=X,y=X,z=twoD_histo)
#        self.threeDwidget.addItem(self.surfplot)
#        self.threeDwidget.show()

        
        
        self.mytimer=QtCore.QTimer()
        self.mytimer.timeout.connect(self.UpdatePlots)
        self.GetAD.clicked.connect(self.GetADValues)
        self.getcaloff.clicked.connect(self.GetCalOffsets)
        self.getstripgains.clicked.connect(self.GetStripGains)
        self.CheckStuff.clicked.connect(self.GetStuff)
        self.discset.clicked.connect(self.SetMyDisc)
        self.setinttime.clicked.connect(self.SetIntTime)
        self.startacq.clicked.connect(self.StreamAcquire)
        self.stopacq.clicked.connect(self.AcquireStop)
        self.enablefaketrig.clicked.connect(self.EnableFake)
        self.trimInput.clicked.connect(self.LowLevelCal)
        self.savedata.clicked.connect(self.SaveData)
        self.readtemperature.clicked.connect(self.ReadTemp)
        self.tube1_tableWidget.setHorizontalHeaderLabels(['x ave','x std','y ave', 'y std'])
        self.saverawposition.stateChanged.connect(self.checksaverawposition)
        self.pltx1.stateChanged.connect(self.checkxplt)
        self.pltx2.stateChanged.connect(self.checkxplt)
        self.pltx3.stateChanged.connect(self.checkxplt)
        self.pltx4.stateChanged.connect(self.checkxplt)
        self.pltx5.stateChanged.connect(self.checkxplt)
        self.pltx6.stateChanged.connect(self.checkxplt)
        self.pltx7.stateChanged.connect(self.checkxplt)
        self.pltx8.stateChanged.connect(self.checkxplt)
        self.plotpsum1to8_rb.toggled.connect(lambda:self.checkrb(self.plotpsum1to8_rb,0))
        self.plotpsum9to16_rb.toggled.connect(lambda:self.checkrb(self.plotpsum9to16_rb,1))
        self.plotad1_rb.toggled.connect(lambda:self.checkrb(self.plotad1_rb,2))
        self.plotad2_rb.toggled.connect(lambda:self.checkrb(self.plotad2_rb,3))
        self.plotad3_rb.toggled.connect(lambda:self.checkrb(self.plotad3_rb,4))
        self.plotad4_rb.toggled.connect(lambda:self.checkrb(self.plotad4_rb,5))
        self.plottof_rb.toggled.connect(lambda:self.checkrb(self.plottof_rb,6))
        self.tree_items=self.initconfig_tree(self.config_tree.invisibleRootItem())
        self.ReadTemp()
        parms=self.MySiPM.readparams()
        if len(parms) < 1152:
            parms=self.MySiPM.readparams()
        struct_parms=struct.unpack('3I 5i 8i 128i 128f 4I 6i 6f',parms)
#control starts at index. 277
        self.updatetreeData(struct_parms)
        self.mystream1.setcaloffsets(self.cal_xoffsets,self.cal_yoffsets)
        self.mystream1.setpslimits(self.psum_max,self.psum_min)
            
        
#add config info tree view.
    
    def initconfig_tree(self,parent):        
        config_item=self.addParent(parent,0,'config',True)
        cal_item=self.addParent(parent,0,'cal',True)
        i1=self.addChild(config_item,0,'inttime','value')
        i2=self.addChild(config_item,0,'flatfield','0,0')
        i11=self.addChild(config_item,0,'discrim','0')
        i12=self.addChild(config_item,0,'enables','0')
        i13=self.addChild(config_item,0,'ps_max','0')
        i14=self.addChild(config_item,0,'ps_min','0')
        i23=self.addChild(config_item,0,'output','raw')
        i24=self.addChild(config_item,0,'fitparms','0.0,0.0,0.0,0.0')
        i25=self.addChild(config_item,0,'ff,ps okay','0,0')
        xoff_item=self.addParent(cal_item,0,'xoffsets',False)
        i3=self.addChild(xoff_item,0,'tube1','0')
        i4=self.addChild(xoff_item,0,'tube2','0')
        i5=self.addChild(xoff_item,0,'tube3','0')
        i6=self.addChild(xoff_item,0,'tube4','0')
        yoff_item=self.addParent(cal_item,0,'yoffsets',False)
        i7=self.addChild(yoff_item,0,'tube1','0')
        i8=self.addChild(yoff_item,0,'tube2','0')
        i9=self.addChild(yoff_item,0,'tube3','0')
        i10=self.addChild(yoff_item,0,'tube4','0')
#####################
        calxoff_item=self.addParent(cal_item,0,'calxoffsets',False)
        i15=self.addChild(calxoff_item,0,'tube1','0')
        i16=self.addChild(calxoff_item,0,'tube2','0')
        i17=self.addChild(calxoff_item,0,'tube3','0')
        i18=self.addChild(calxoff_item,0,'tube4','0')
        calyoff_item=self.addParent(cal_item,0,'calyoffsets',False)
        i19=self.addChild(calyoff_item,0,'tube1','0')
        i20=self.addChild(calyoff_item,0,'tube2','0')
        i21=self.addChild(calyoff_item,0,'tube3','0')
        i22=self.addChild(calyoff_item,0,'tube4','0')
        
#####################        
        self.config_tree.setColumnWidth(0,110)
        self.config_tree.setColumnWidth(1,300)
        z={'inttime':i1,
           'flatfield':i2,
           'xoff_t1':i3,
           'xoff_t2':i4,
           'xoff_t3':i5,
           'xoff_t4':i6,
           'yoff_t1':i7,
           'yoff_t2':i8,
           'yoff_t3':i9,
           'yoff_t4':i10,
           'disc_levels':i11,
           'disc_en':i12,
           'ps_max':i13,
           'ps_min':i14,
           'calxoff_t1':i15,
           'calxoff_t2':i16,
           'calxoff_t3':i17,
           'calxoff_t4':i18,
           'calyoff_t1':i19,
           'calyoff_t2':i20,
           'calyoff_t3':i21,
           'calyoff_t4':i22,
           'output_mode':i23,
           'fit_parameters':i24,
           'ffpsok':i25
           }
        return z
    def addParent(self,parent,column,title,expand=False):
        item = QtWidgets.QTreeWidgetItem(parent, [title])
        item.setChildIndicatorPolicy(QtWidgets.QTreeWidgetItem.ShowIndicator)
        item.setExpanded(expand)
        return item
    def addChild(self,parent,column,title,data):
        item = QtWidgets.QTreeWidgetItem(parent, [title])
        cstext='='+data
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, data)
        return item
    def updatetreeData(self,struct_parms):
        self.countsarray=numpy.asarray(struct_parms[0:3],dtype=numpy.uint32)
        self.statusarray=numpy.asarray(struct_parms[3:8],dtype=numpy.int32)
        self.discvalues=numpy.asarray(struct_parms[8:16],dtype=numpy.int32)
        self.xoffsets=numpy.asarray(struct_parms[16:48],dtype=numpy.int32)
        self.yoffsets=numpy.asarray(struct_parms[48:80],dtype=numpy.int32)
        self.xzeros=numpy.asarray(struct_parms[80:112],dtype=numpy.int32)
        self.yzeros=numpy.asarray(struct_parms[112:144],dtype=numpy.int32)
        self.cal_xoffsets=numpy.asarray(struct_parms[144:176],dtype=numpy.float32)
        self.cal_yoffsets=numpy.asarray(struct_parms[176:208],dtype=numpy.float32)
        self.cal_xstripgains=numpy.asarray(struct_parms[208:240],dtype=numpy.float32)        
        self.cal_ystripgains=numpy.asarray(struct_parms[240:272],dtype=numpy.float32)        
        self.tof_offset=struct_parms[272]
        self.pixel_offset=struct_parms[273]
        self.frame=struct_parms[274]
        self.frame_rate=struct_parms[275]
        self.int_time=struct_parms[276]
        self.disc_en=struct_parms[277]
        self.en_flatfield=struct_parms[278]
        self.en_photosum=struct_parms[279]
        self.event_output_type=struct_parms[280]
        self.rtdl_header_type=struct_parms[281]
        self.psum_max=struct_parms[282]
        self.psum_min=struct_parms[283]
        self.alpha1=struct_parms[284]
        self.weight_alpha1=struct_parms[285]
        self.alpha2=struct_parms[286]
        self.weight_alpha2=struct_parms[287]
        
        self.disc_setpoint.setValue(self.discvalues[0])
        self.inttime_setpoint.setValue(self.int_time)
        
        cstext='Version=%.1f Date=%x' %(.1*self.statusarray[2],self.statusarray[1])
        self.versionlabel.setText(cstext)
        
        item=self.tree_items['xoff_t1']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.xoffsets[i])
            cstext=cstext+cs
        cs='%d' %(self.xoffsets[7])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.xoffsets[0:8])
        
        item=self.tree_items['xoff_t2']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.xoffsets[i+8])
            cstext=cstext+cs
        cs='%d' %(self.xoffsets[15])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.xoffsets[8:16])
        
        item=self.tree_items['xoff_t3']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.xoffsets[i+16])
            cstext=cstext+cs
        cs='%d' %(self.xoffsets[23])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.xoffsets[16:24])
        
        item=self.tree_items['xoff_t4']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.xoffsets[i+24])
            cstext=cstext+cs
        cs='%d' %(self.xoffsets[31])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.xoffsets[24:32])
        
        item=self.tree_items['yoff_t1']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.yoffsets[i])
            cstext=cstext+cs
        cs='%d' %(self.yoffsets[7])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.yoffsets[0:8])
        
        item=self.tree_items['yoff_t2']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.yoffsets[i+8])
            cstext=cstext+cs
        cs='%d' %(self.yoffsets[15])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.yoffsets[8:16])
        
        item=self.tree_items['yoff_t3']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.yoffsets[i+16])
            cstext=cstext+cs
        cs='%d' %(self.yoffsets[23])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.yoffsets[16:24])
        
        item=self.tree_items['yoff_t4']
        cstext='= '
        for i in range(7):
            cs='%d,' %(self.yoffsets[i+24])
            cstext=cstext+cs
        cs='%d' %(self.yoffsets[31])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.yoffsets[24:32])
        
        
        item=self.tree_items['calxoff_t1']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_xoffsets[i])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_xoffsets[7])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_xoffsets[0:8])
        
        item=self.tree_items['calxoff_t2']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_xoffsets[i+8])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_xoffsets[15])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_xoffsets[8:16])
        
        item=self.tree_items['calxoff_t3']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_xoffsets[i+16])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_xoffsets[23])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_xoffsets[16:24])
        
        item=self.tree_items['calxoff_t4']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_xoffsets[i+24])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_xoffsets[31])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_xoffsets[24:32])
        
        item=self.tree_items['calyoff_t1']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_yoffsets[i])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_yoffsets[7])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_yoffsets[0:8])
        
        item=self.tree_items['calyoff_t2']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_yoffsets[i+8])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_yoffsets[15])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_yoffsets[8:16])
        
        item=self.tree_items['calyoff_t3']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_yoffsets[i+16])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_yoffsets[23])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_yoffsets[16:24])
        
        item=self.tree_items['calyoff_t4']
        cstext='= '
        for i in range(7):
            cs='%.1f,' %(self.cal_yoffsets[i+24])
            cstext=cstext+cs
        cs='%.1f' %(self.cal_yoffsets[31])
        cstext=cstext+cs
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.cal_yoffsets[24:32])
        
        
        
        item=self.tree_items['disc_levels']
        cstext='= %d,%d,%d,%d,%d,%d,%d' %(self.discvalues[0],self.discvalues[1],self.discvalues[2],self.discvalues[3],self.discvalues[4],self.discvalues[5],self.discvalues[6])
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.discvalues[0:6])
        item=self.tree_items['inttime']
        cstext='= %d' %(self.int_time)
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.int_time)
        
        item=self.tree_items['ps_max']
        cstext='= %.1f' %(self.psum_max)
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.psum_max)
        
        item=self.tree_items['ps_min']
        cstext='= %.1f' %(self.psum_min)
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.psum_min)
        
        item=self.tree_items['fit_parameters']
        cstext='=%.3f,%.3f,%.3f,%.3f' %(self.alpha1,self.alpha2,self.weight_alpha1,self.weight_alpha2)
        item.setText(1,cstext)
        
        item=self.tree_items['ffpsok']
        cstext='= %d,%d' % (self.statusarray[3],self.statusarray[4])
        item.setText(1,cstext)
        
        item=self.tree_items['disc_en']
        cstext='=0x %X' %(self.disc_en)
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.disc_en)
        
        
        item=self.tree_items['flatfield']
        if self.en_flatfield==0:
            cstext='= False'
        else:
            cstext='=True'
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.en_flatfield)
        
        item=self.tree_items['output_mode']
        if self.event_output_type==1:
            cstext='= normal'
        else:
            cstext='= raw'
        item.setText(1,cstext)
        item.setData(1, QtCore.Qt.UserRole, self.event_output_type)
        
        
        
#    def updateconfig_tree(self,configdic):
            
    def checkrb(self,b,value):
        if b.isChecked():
            self.whatplot=value
            if value==0:
                self.win.setWindowTitle('PSum histo x columns 1-8')
                self.win2.setWindowTitle('PSum histo y rows 1-8')
            elif value==1:
                self.win.setWindowTitle('PSum histo x columns 9-16')
                self.win2.setWindowTitle('PSum histo y rows 9-16')
            elif value==2:
                self.win.setWindowTitle('A/D histo x tube 1')
                self.win2.setWindowTitle('A/D histo y tube 1')
            elif value==3:
                self.win.setWindowTitle('A/D histo x tube 2')
                self.win2.setWindowTitle('A/D histo y tube 2')
            elif value==4:
                self.win.setWindowTitle('A/D histo x tube 3')
                self.win2.setWindowTitle('A/D histo y tube 3')
            elif value==5 :
                self.win.setWindowTitle('A/D histo x tube 4')
                self.win2.setWindowTitle('A/D histo y tube 4')            
            else :
                self.win.setWindowTitle('TOF histo 10usec bins')
                self.win2.setWindowTitle('TOF histo 0.1usec bins')
                
            if self.acquireflag==0:
                self.UpdatePlots()
        
    

    def checksaverawposition(self,value):
        if self.saverawposition.isChecked():
            self.saveraw=1
        else:
            self.saveraw=0
            self.mystream1.saveraw(0)
    
    def checkxplt(self,value):
        if self.pltx1.isChecked():
            self.xpltstate[0]=1
        else:
            self.xpltstate[0]=0
            
        if self.pltx2.isChecked():
            self.xpltstate[1]=1
        else:
            self.xpltstate[1]=0
            
        if self.pltx3.isChecked():
            self.xpltstate[2]=1
        else:
            self.xpltstate[2]=0
            
        if self.pltx4.isChecked():
            self.xpltstate[3]=1
        else:
            self.xpltstate[3]=0
            
        if self.pltx5.isChecked():
            self.xpltstate[4]=1
        else:
            self.xpltstate[4]=0
            
        if self.pltx6.isChecked():
            self.xpltstate[5]=1
        else:
            self.xpltstate[5]=0
            
        if self.pltx7.isChecked():
            self.xpltstate[6]=1
        else:
            self.xpltstate[6]=0
            
        if self.pltx8.isChecked():
            self.xpltstate[7]=1
        else:
            self.xpltstate[7]=0
            
        if self.acquireflag==0:
            self.UpdatePlots()

            
        
        
    def closeEvent(self,event):
        global fileopen
        global myfile
        self.mystream1.killme()
        time.sleep(1)
        self.win.close()
        self.win2.close()
        self.win3.close()
        time.sleep(1)
        self.AcquireStop()
        time.sleep(1)
        
    def ReadTemp(self):
        self.temperature=self.MySiPM.readtemperature()
        cs='%.2f' % self.temperature
        self.temperature_sb.setText(cs)
        
    def StreamAcquire(self):
        global zero_histos
        global gcounts
        global global_acquireflag
        if self.acquireflag==1:
            global_acquireflag=1
            return
        self.acquireflag=1
        global_acquireflag=1
        zero_histos=1
        gcounts=0
#        self.mystream1.setsavelightcone('C:/Users/Larmor/Desktop/AngerCamera/lightconedata_gauss.txt')
        if self.saveraw==1:
            self.mystream1.saverawpos(1)

        self.prevcounts=0
        self.MySiPM.startaquire(self.enablefake,self.int_time)
        self.runningtime.start()
        self.mytimer.start(5000)  #update every 5 seconds
        
    def AcquireStop(self):
        global zero_histos
        global global_acquireflag
        global twoD_histo
        zero_histos=0
        self.mytimer.stop()
        time.sleep(0.5)
        global_acquireflag=0
    
        self.MySiPM.stopaquire()
        self.acquireflag=0
        self.mystream1.saverawpos(0)
        et=self.runningtime.elapsed()
        if self.mystream1.savelc==1:
            self.mystream1.stopsavelightcone()
        if self.instripgainlogic==1:
            self.instripgainlogic=0
            self.CalStripGains()
        etm=et*0.001
        cs='acumm. time = %.1f sec' % etm
        self.stufflabel.setText(cs)
        self.UpdatePlots()

    def skewed_gauss(self,p,x):
        y=p[0]*numpy.exp(-0.5*(x-p[1])**2/p[2]**2)*(1.0+erf(p[3]*((x-p[1])/(1.414*p[2]))))+p[4]
        return y

    def fit_func(self,func, p,x,y,ye):
        errfunc = lambda p,x,y,ye: abs(func(p,x) - y )
        par, cov, dum,dum,dum = leastsq(errfunc,p,args=(x,y,ye),full_output=True)
        chi2 = sum(errfunc(par,x,y,y)**2)/(len(x) - len(par))
#    parerr=np.sqrt(np.abs(cov.diagonal()))
        return par, chi2


    def analysisnew(self,xdata,ydata):
        ye=numpy.sqrt(ydata)
        back=ydata[0]
        imax=numpy.argmax(ydata)
        xmean=xdata[imax]
        temp=-4
#find a full width half max
        maxval=ydata[imax]
        for i in range(100):
            if ydata[imax+i] < .5*maxval:
                break
        hwfm=xdata[imax+i]-xdata[imax]            
        sig=hwfm/1.25

        gpar,gchi2 = self.fit_func(self.skewed_gauss,[maxval, xmean,sig,temp,back],xdata,ydata,ye)
        return gpar



    def CalStripGains(self):
        global ad_histox
        global ad_histoy
        global ad_histo_xaxis
        answer = QtWidgets.QMessageBox.question(self, "",
                "Do you wish to continue \n"
                "calculating Strip Gains?\n",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)

        if answer == QtWidgets.QMessageBox.Yes:
            xmean=numpy.zeros(8)
            ymean=numpy.zeros(8)
            for i in range(8):
                xpar=self.analysisnew(ad_histo_xaxis[150:350],ad_histox[150:350,i,0])
                xmean[i]=xpar[1]
            for i in range(8):
                ypar=self.analysisnew(ad_histo_xaxis[150:350],ad_histoy[150:350,i,0])
                ymean[i]=ypar[1]
            xgainave=numpy.average(xmean)
            ygainave=numpy.average(ymean)
            cs='xpeakave=%.2f, ypeakave=%.2f' %(xgainave,ygainave)
            self.stufflabel.setText(cs)
#            print(xgainave,ygainave)
            #just do tube 1
            for i in range(8):
                if xmean[i] > 500.0:
                    self.cal_xstripgains[i]=xgainave/xmean[i]
                else:
                    self.cal_xstripgains[i]=1.0
                    
                if ymean[i] > 500.0:
                    self.cal_ystripgains[i]=ygainave/ymean[i]
                else:
                    self.cal_ystripgains[i]=1.0                
            #I am doing tube 1 only for now.
            pv=self.PackValues()
            self.MySiPM.writeparams(pv)
            
            answer = QtWidgets.QMessageBox.question(self, "",
                    "detector.ini is stale. \n"
                    "Do you want to upgrade the file?",
                    QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
      
            if answer == QtWidgets.QMessageBox.Yes:
                print('you clicked yes')
                self.ReadConfigIni('192.168.10.19')
                wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
                myfile=open(wf,'r+t')
                allfile=myfile.read()
                myfile.close()
                k=allfile.find('[xcalstripgains]')
                j=allfile.find('[',k+10)
                changepart='[xcalstripgains]\n'
                for t in range(4):
                    cs='t%d=' %(t+1)
                    changepart=changepart+cs
                    for i in range(7):
                        cs='%.3f,' % self.cal_xstripgains[8*t+i]
                        changepart=changepart+cs
                    cs='%.3f\n' % self.cal_xstripgains[8*t+7]
                    changepart=changepart+cs
                if (j >= 0):
                    allfile=allfile[0:k]+changepart+allfile[j:]
                else:
                    allfile=allfile[0:k]+changepart
                k=allfile.find('[ycalstripgains]')
                j=allfile.find('[',k+10)
                changepart='[ycalstripgains]\n'
                for t in range(4):
                    cs='t%d=' %(t+1)
                    changepart=changepart+cs
                    for i in range(7):
                        cs='%.3f,' % self.cal_ystripgains[8*t+i]
                        changepart=changepart+cs
                    cs='%.3f\n' % self.cal_ystripgains[8*t+7]
                    changepart=changepart+cs
                if (j >= 0):
                    newallfile=allfile[0:k]+changepart+allfile[j:]
                else:
                    newallfile=allfile[0:k]+changepart
                myfile=open(wf,'w+t')
                myfile.write(newallfile)
                myfile.close()
                self.WriteConfigIni('192.168.10.19')
    


    
    def GetStripGains(self):
        answer = QtWidgets.QMessageBox.question(self, "",
                "After selecting to continue, you must start acquire\n"
                "When stop acquire is pressed a new dialog will open\n"
                "To continue with Stip Gains Click yes. \n",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
  
        if answer == QtWidgets.QMessageBox.Yes:
            self.instripgainlogic=1
        else:
            self.instripgainlogic=0
        #just to tube zero for now
        return
    
    def PackValues(self):  #this routine packs tube and config params
        b1=self.discvalues.tostring()
        b2=self.xoffsets.tostring()
        b3=self.yoffsets.tostring()
        b4=self.xzeros.tostring()
        b5=self.yzeros.tostring()
        b6=self.cal_xoffsets.tostring()
        b7=self.cal_yoffsets.tostring()
        b8=self.cal_xstripgains.tostring()
        b9=self.cal_ystripgains.tostring()
        b10=struct.pack('4I 6i 6f', self.tof_offset, self.pixel_offset,self.frame,self.frame_rate,self.int_time,self.disc_en,self.en_flatfield,
                                 self.en_photosum, self.event_output_type,self.rtdl_header_type,
                                 self.psum_max,self.psum_min,self.alpha1,self.weight_alpha1,self.alpha2,self.weight_alpha2)

        packedvalues=b1+b2+b3+b4+b5+b6+b7+b8+b9+b10
        return packedvalues
    
    def GetCalOffsets(self):
        mb = QtWidgets.QMessageBox()
        mb.setText('Turn on bias voltage\nand place away from neutron source')
        mb.exec()
#get current values in detector memory.
        xmean1,ymean1,xstd1,ystd1=self.MySiPM.set_and_read(self.int_time)
        self.cal_xoffsets[0:8]=xmean1[:,0]
        self.cal_xoffsets[8:16]=xmean1[:,1]
        self.cal_xoffsets[16:24]=xmean1[:,2]
        self.cal_xoffsets[24:32]=xmean1[:,3]
        self.cal_yoffsets[0:8]=ymean1[:,0]
        self.cal_yoffsets[8:16]=ymean1[:,1]
        self.cal_yoffsets[16:24]=ymean1[:,2]
        self.cal_yoffsets[24:32]=ymean1[:,3]
        self.mystream1.setcaloffsets(self.cal_xoffsets,self.cal_yoffsets)
#need an update tree data call.
        labellist=['A/D  at inttime', 'SD at inttime', 'old A/d @t=30','old SD @t=30']
        self.tube1_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube2_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube3_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube4_tableWidget.setHorizontalHeaderLabels(labellist)
        for i in range(8):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,0])
            self.tube1_tableWidget.setItem(i,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,0])
            self.tube1_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,0])
            self.tube1_tableWidget.setItem(i+8,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,0])
            self.tube1_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,1])
            self.tube2_tableWidget.setItem(i,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,1])
            self.tube2_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,1])
            self.tube2_tableWidget.setItem(i+8,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,1])
            self.tube2_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,2])
            self.tube3_tableWidget.setItem(i,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,2])
            self.tube3_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,2])
            self.tube3_tableWidget.setItem(i+8,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,2])
            self.tube3_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,3])
            self.tube4_tableWidget.setItem(i,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,3])
            self.tube4_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,3])
            self.tube4_tableWidget.setItem(i+8,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,3])
            self.tube4_tableWidget.setItem(i+8,1,entry)
        pv=self.PackValues()
        self.MySiPM.writeparams(pv)
        time.sleep(1)
        parms=self.MySiPM.readparams()  #this gets data from processor memory which has data from MysiPM.lowlevel call
        struct_parms=struct.unpack('3I 5i 8i 128i 128f 4I 6i 6f',parms)
        self.updatetreeData(struct_parms)
        
        answer = QtWidgets.QMessageBox.question(self, "",
                "detector.ini is stale. \n"
                "Do you want to upade the file?",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
  
        if answer == QtWidgets.QMessageBox.Yes:
            print('you clicked yes')
            self.ReadConfigIni('192.168.10.19')
            wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
            myfile=open(wf,'r+t')
            allfile=myfile.read()
            myfile.close()
            k=allfile.find('[xcaloffsets]')
            j=allfile.find('[',k+10)
            changepart='[xcaloffsets]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%.1f,' % self.cal_xoffsets[8*t+i]
                    changepart=changepart+cs
                cs='%.1f\n' % self.cal_xoffsets[8*t+7]
                changepart=changepart+cs
            allfile=allfile[0:k]+changepart+allfile[j:]
            k=allfile.find('[ycaloffsets]')
            j=allfile.find('[',k+10)
            changepart='[ycaloffsets]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%.1f,' % self.cal_yoffsets[8*t+i]
                    changepart=changepart+cs
                cs='%.1f\n' % self.cal_yoffsets[8*t+7]
                changepart=changepart+cs
            newallfile=allfile[0:k]+changepart+allfile[j:]
            myfile=open(wf,'w+t')
            myfile.write(newallfile)
            myfile.close()
            self.WriteConfigIni('192.168.10.19')
                
        return
    
    def LowLevelCal(self):
        orgint_time=self.int_time
        labellist=['DAC Offsets', 'A/D @t=5', 'A/D @t=30','SD @t=30']
        self.tube1_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube2_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube3_tableWidget.setHorizontalHeaderLabels(labellist)
        self.tube4_tableWidget.setHorizontalHeaderLabels(labellist)
        labellist=['x1','x2', 'x3', 'x4', 'x5', 'x6','x7','x8','y1','y2', 'y3', 'y4', 'y5', 'y6','y7','y8']
        self.tube1_tableWidget.setVerticalHeaderLabels(labellist)
        self.tube2_tableWidget.setVerticalHeaderLabels(labellist)
        self.tube3_tableWidget.setVerticalHeaderLabels(labellist)
        self.tube4_tableWidget.setVerticalHeaderLabels(labellist)
        offsets=self.MySiPM.lowlevel(1,orgint_time)       
        for i in range(16):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % offsets[i])
            self.tube1_tableWidget.setItem(i,0,entry)
        offsets=self.MySiPM.lowlevel(2,orgint_time)       
        for i in range(16):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % offsets[i])
            self.tube2_tableWidget.setItem(i,0,entry)
        offsets=self.MySiPM.lowlevel(3,orgint_time)       
        for i in range(16):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % offsets[i])
            self.tube3_tableWidget.setItem(i,0,entry)
        offsets=self.MySiPM.lowlevel(4,orgint_time)       
        for i in range(16):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % offsets[i])
            self.tube4_tableWidget.setItem(i,0,entry)
###########################################################################
        xmean1,ymean1,xstd1,ystd1=self.MySiPM.set_and_read(5)
        for i in range(8):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,0])
            self.tube1_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,0])
            self.tube1_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,1])
            self.tube2_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,1])
            self.tube2_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,2])
            self.tube3_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,2])
            self.tube3_tableWidget.setItem(i+8,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,3])
            self.tube4_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,3])
            self.tube4_tableWidget.setItem(i+8,1,entry)
        xmean1,ymean1,xstd1,ystd1=self.MySiPM.set_and_read(30)
        for i in range(8):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,0])
            self.tube1_tableWidget.setItem(i,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,0])
            self.tube1_tableWidget.setItem(i,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,0])
            self.tube1_tableWidget.setItem(i+8,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,0])
            self.tube1_tableWidget.setItem(i+8,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,1])
            self.tube2_tableWidget.setItem(i,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,1])
            self.tube2_tableWidget.setItem(i,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,1])
            self.tube2_tableWidget.setItem(i+8,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,1])
            self.tube2_tableWidget.setItem(i+8,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,2])
            self.tube3_tableWidget.setItem(i,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,2])
            self.tube3_tableWidget.setItem(i,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,2])
            self.tube3_tableWidget.setItem(i+8,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,2])
            self.tube3_tableWidget.setItem(i+8,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,3])
            self.tube4_tableWidget.setItem(i,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,3])
            self.tube4_tableWidget.setItem(i,3,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,3])
            self.tube4_tableWidget.setItem(i+8,2,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,3])
            self.tube4_tableWidget.setItem(i+8,3,entry)
        self.MySiPM.setinttime(orgint_time)
        parms=self.MySiPM.readparams()  #this gets data from processor memory which has data from MysiPM.lowlevel call
        struct_parms=struct.unpack('3I 5i 8i 128i 128f 4I 6i 6f',parms)
        self.updatetreeData(struct_parms)  #updata local data
        answer = QtWidgets.QMessageBox.question(self, "",
                "detector.ini is stale. \n"
                "Do you want to upade the file?",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
  
        if answer == QtWidgets.QMessageBox.Yes:
#            print('you clicked yes')
            self.ReadConfigIni('192.168.10.19')
            wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
            myfile=open(wf,'r+t')
            allfile=myfile.read()
            myfile.close()
            k=allfile.find('[xoffsets]')
            j=allfile.find('[',k+10)
            changepart='[xoffsets]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%d,' % self.xoffsets[8*t+i]
                    changepart=changepart+cs
                cs='%d\n' % self.xoffsets[8*t+7]
                changepart=changepart+cs
            allfile=allfile[0:k]+changepart+allfile[j:]
            k=allfile.find('[yoffsets]')
            j=allfile.find('[',k+10)
            changepart='[yoffsets]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%d,' % self.yoffsets[8*t+i]
                    changepart=changepart+cs
                cs='%d\n' % self.yoffsets[8*t+7]
                changepart=changepart+cs
            k=allfile.find('[xzeros]')
            j=allfile.find('[',k+10)
            changepart='[xzeros]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%d,' % self.xzeros[8*t+i]
                    changepart=changepart+cs
                cs='%d\n' % self.xzeros[8*t+7]
                changepart=changepart+cs
            allfile=allfile[0:k]+changepart+allfile[j:]
            k=allfile.find('[yzeros]')
            j=allfile.find('[',k+10)
            changepart='[yzeros]\n'
            for t in range(4):
                cs='t%d=' %(t+1)
                changepart=changepart+cs
                for i in range(7):
                    cs='%d,' % self.yzeros[8*t+i]
                    changepart=changepart+cs
                cs='%d\n' % self.yzeros[8*t+7]
                changepart=changepart+cs
            allfile=allfile[0:k]+changepart+allfile[j:]
            newallfile=allfile[0:k]+changepart+allfile[j:]
            myfile=open(wf,'w+t')
            myfile.write(newallfile)
            myfile.close()
            self.WriteConfigIni('192.168.10.19')
                
                    
    def ReadPSumTables(self,ip):        
        session=ftplib.FTP()
        rsp=session.connect(ip)
        wf='C:/Users/Larmor/Desktop/AngerCamera/ptable_max.dat'
        myfile=open(wf,'w+b')
        try:
            rsp=session.retrbinary('RETR /mnt/tables/ptable_max.dat',myfile.write)
        except:
            myfile.close()
            session.quit()
            return -1
        myfile.close()
        k=rsp.find('succ')
        if k < 0:
            session.quit()
            return k
        wf='C:/Users/Larmor/Desktop/AngerCamera/ptable_min.dat'
        myfile=open(wf,'w+b')
        rsp=session.retrbinary('RETR /mnt/tables/ptable_min.dat',myfile.write)
        myfile.close()
        k=rsp.find('succ')
        if k < 0:
            session.quit()
            return k
        session.quit()    
        return 0
        
    def ReadConfigIni(self,ip):
        session=ftplib.FTP()
        rsp=session.connect(ip)
        print(rsp)
        wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
        myfile=open(wf,'w+b')
        rsp=session.retrbinary('RETR /mnt/detector.ini',myfile.write)
        myfile.close()
        print(rsp)
        session.quit()    
        return 0

    def WriteConfigIni(self,ip):
        session=ftplib.FTP()
        rsp=session.connect(ip)
        print(rsp)
        wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
        myfile=open(wf,'r+b')
        rsp=session.storbinary('STOR /mnt/detector.ini',myfile)
        myfile.close()
        print(rsp)
        session.quit()    
        return 0

    def SaveData(self):
        global psum_histox
        global psum_histoy
        global psum_histo_xaxis 
        global ad_histo_xaxis
        global ad_histox
        global ad_histoy
        global twoD_histo
        global psum3Dhisto
        global tof_histo
        filename='C:/Users/Larmor/Dropbox/AngerCamera/Correction_Anger/psums.txt'
        writefile=open(filename,'w+t')
        str1='columns are psum scale followed by x 1-16 then y 1-16\n'
        writefile.write(str1)
        for i in range(1000):
            mystr='%d\t' % psum_histo_xaxis[i]
            for j in range(16):
                str1='%d\t' % psum_histox[i,j]
                mystr=mystr+str1
            for j in range(16):
                str1='%d\t' % psum_histoy[i,j]
                mystr=mystr+str1
            mystr=mystr+'\n'
            writefile.write(mystr)
        writefile.close()
        filename='C:/Users/Larmor/Dropbox/AngerCamera/Correction_Anger/adhistos.txt'
        writefile=open(filename,'w+t')
        str1='columns are ad scale followed by tube 1, x 1-8, tube 2, x1-8, tube 3 x1-8 tube 4 x1-8 then ys in same tube order\n'
        writefile.write(str1)
        for i in range(500):
            mystr='%d\t' % ad_histo_xaxis[i]
            for j in range(4):
                for k in range(8):
                    str1='%d\t' % ad_histox[i,k,j]
                    mystr=mystr+str1
            for j in range(4):
                for k in range(8):
                    str1='%d\t' % ad_histoy[i,k,j]
                    mystr=mystr+str1
            mystr=mystr+'\n'
            writefile.write(mystr)
        writefile.close()
        filename='C:/Users/Larmor/Dropbox/AngerCamera/Correction_Anger/sipm2d.dat'
        writefile=open(filename,'w+b')
        writefile.write(twoD_histo)
        writefile.close()
        filename='C:/Users/Larmor/Dropbox/AngerCamera/Correction_Anger/sipmtofhisto.dat'
        writefile=open(filename,'w+b')
        writefile.write(tof_histo)
        writefile.close()
        filename='C:/Users/Larmor/Dropbox/AngerCamera/Correction_Anger/psum3d.dat'
        writefile=open(filename,'w+b')
        writefile.write(psum3Dhisto)
        writefile.close()
        
    def UpdatePlots(self):
        global gcounts
        global grate
        global psum_histox
        global psum_histoy
        global psum_histo_xaxis 
        global ad_histo_xaxis
        global ad_histox
        global ad_histoy
        global twoD_histo
        global bm1counts
        global tof_histo
        global tof_histo_fine
        cs='%d' % bm1counts
        self.bm1countsvalue.setText(cs)
        cs='%d' % gcounts
        self.countsvalue.setText(cs)
        grate=(gcounts-self.prevcounts)/5.0
        self.prevcounts=gcounts
        cs='%.1f' % grate
        self.countsratevalue.setText(cs)
        self.img.setImage(twoD_histo)
#whatplot is set in checkrb.  0 means plots 1-16 photosums 1 means plots a/d values.  (up to two tubes and a time.)
        if self.whatplot==6:
            self.plt1.plot(tof_histo,clear=True)
            self.plt2.plot(tof_histo_fine,clear=True)
            return

        if self.whatplot==0 or self.whatplot==1:
            xcleared=0
            ycleared=0
            if (self.xpltstate[0]==1):
                self.plt1.plot(psum_histo_xaxis,psum_histox[:,0+8*self.whatplot],clear=True,pen=(0,8))
                xcleared=1
            for i in range(1,8,1):
                if (self.xpltstate[i]==1):
                    if xcleared==1:
                        self.plt1.plot(psum_histo_xaxis,psum_histox[:,i+8*self.whatplot],pen=(i,8))  
                    else:
                        xcleared=1
                        self.plt1.plot(psum_histo_xaxis,psum_histox[:,i+8*self.whatplot],clear=True,pen=(i,8))  
                        
            if (self.ypltstate[0]==1):    
                self.plt2.plot(psum_histo_xaxis,psum_histoy[:,0+8*self.whatplot],clear=True,pen=(0,8)) 
                ycleared=1
            for i in range(1,8,1):
                if (self.ypltstate[i]==1):
                    if ycleared==1:
                        self.plt2.plot(psum_histo_xaxis,psum_histoy[:,i+8*self.whatplot],pen=(i,8))  
                    else:
                        ycleared=1
                        self.plt2.plot(psum_histo_xaxis,psum_histoy[:,i*self.whatplot],clear=True,pen=(i,8))  
        else:
            xcleared=0
            ycleared=0
            if (self.xpltstate[0]==1):
                self.plt1.plot(ad_histo_xaxis,ad_histox[:,0,self.whatplot-2],clear=True,pen=(0,8))
                xcleared=1
            for i in range(1,8,1):
                if (self.xpltstate[i]==1):
                    if xcleared==1:
                        self.plt1.plot(ad_histo_xaxis,ad_histox[:,i,self.whatplot-2],pen=(i,8))  
                    else:
                        xcleared=1
                        self.plt1.plot(ad_histo_xaxis,ad_histox[:,i,self.whatplot-2],clear=True,pen=(i,8))  
                        
            if (self.ypltstate[0]==1):    
                self.plt2.plot(ad_histo_xaxis,ad_histoy[:,0,self.whatplot-2],clear=True,pen=(0,8)) 
                ycleared=1
            for i in range(1,8,1):
                if (self.ypltstate[i]==1):
                    if ycleared==1:
                        self.plt2.plot(ad_histo_xaxis,ad_histoy[:,i,self.whatplot-2],pen=(i,8))  
                    else:
                        ycleared=1
                        self.plt2.plot(ad_histo_xaxis,ad_histoy[:,i,self.whatplot-2],clear=True,pen=(i,8))  
                
        
        
    def GetADValues(self):
        xmean1,ymean1,xstd1,ystd1=self.MySiPM.set_and_read(self.int_time)
#need an update tree data call.
        labellist=['A/D  at inttime', 'SD at inttime', 'old data ','old data']
        self.tube1_tableWidget.setHorizontalHeaderLabels(labellist)
        for i in range(8):       
            entry=QtWidgets.QTableWidgetItem('%.2f' % xmean1[i,0])
            self.tube1_tableWidget.setItem(i,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % xstd1[i,0])
            self.tube1_tableWidget.setItem(i,1,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ymean1[i,0])
            self.tube1_tableWidget.setItem(i+8,0,entry)
            entry=QtWidgets.QTableWidgetItem('%.2f' % ystd1[i,0])
            self.tube1_tableWidget.setItem(i+8,1,entry)
        xave=numpy.average(xmean1[:,0])
        yave=numpy.average(ymean1[:,0])
        xstdave=numpy.average(xstd1[:,0])
        ystdave=numpy.average(ystd1[:,0])
        cs='xave,std=%.2f,  %.2f, yave,std=%.2f, %.2f' %(xave,xstdave,yave,ystdave)
        self.stufflabel.setText(cs)
        
    def CheckFirstTen(self):
        for i in range(10):
            fit=self.MySiPM.caltest(i)
            print('fit = ',fit)
             
    def GetStuff(self):
        i=self.MySiPM.readdisc()
        k=self.MySiPM.readadcstatus()
        myrtdl=self.MySiPM.readRTDL()
#        myrtdlmem=self.MySiPM.readRTDLmem()
        cs='first rtdl regs = x%X, x%X, x%X' %(myrtdl[0],myrtdl[1], myrtdl[4])
        print(cs)
        cs='dis=x%X, psum vetos=%d, pos vetos=%d' %(i[0],k[1],k[2])
        self.stufflabel.setText(cs)
        
    def EnableFake(self):
        if (self.enablefake==0):
            self.enablefake=1
            self.enablefaketrig.setText('Disable Fake Trigger')
        else:
            self.enablefake=0
            self.enablefaketrig.setText('Enable Fake Trigger')
#        self.MySiPM.enfaketrig(self.enablefake)
    def SetMyDisc(self):
        i=self.disc_setpoint.value()
        ii=int(i)
#        print(ii)
        for k in range(1,8,1):
            self.MySiPM.setdac(0,0,k,ii)
            
        parms=self.MySiPM.readparams()
        struct_parms=struct.unpack('3I 5i 8i 128i 128f 4I 6i 6f',parms)
#control starts at index. 277
        self.updatetreeData(struct_parms)
        
        answer = QtWidgets.QMessageBox.question(self, "",
                "detector.ini is stale. \n"
                "Do you want to upade the file?",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
  
        if answer == QtWidgets.QMessageBox.Yes:
#            print('you clicked yes')
            self.ReadConfigIni('192.168.10.19')
            wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
            myfile=open(wf,'r+t')
            allfile=myfile.read()
            myfile.close()
            k=allfile.find('[discriminators]')
            j=allfile.find('[',k+10)
            changepart='[discriminators]\nd='
            for i in range(6):
                cs='%d,' % self.discvalues[i]
                changepart=changepart+cs
            cs='%d\n' % self.discvalues[6]
            changepart=changepart+cs
            newallfile=allfile[0:k]+changepart+allfile[j:]
            myfile=open(wf,'w+t')
            myfile.write(newallfile)
            myfile.close()
            self.WriteConfigIni('192.168.10.19')
        
    def SetIntTime(self):
        i=self.inttime_setpoint.value()
        ii=int(i)
#        print(ii)
        self.MySiPM.setinttime(ii)
        parms=self.MySiPM.readparams()
        struct_parms=struct.unpack('3I 5i 8i 128i 128f 4I 6i 6f',parms)
#control starts at index. 277
        self.updatetreeData(struct_parms)
        
        answer = QtWidgets.QMessageBox.question(self, "",
                "detector.ini is stale. \n"
                "Do you want to upade the file?",
                QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
  
        if answer == QtWidgets.QMessageBox.Yes:
#            print('you clicked yes')
            self.ReadConfigIni('192.168.10.19')
            wf='C:/Users/Larmor/Desktop/AngerCamera/detector.ini'
            myfile=open(wf,'r+t')
            allfile=myfile.read()
            myfile.close()
            k=allfile.find('inttime')
            j=allfile.find('\n',k+7)+1
            changepart='inttime=%d\n' % self.int_time
            newallfile=allfile[0:k]+changepart+allfile[j:]
            myfile=open(wf,'w+t')
            myfile.write(newallfile)
            myfile.close()
            self.WriteConfigIni('192.168.10.19')


        
if __name__ == "__main__":

    debugon=0
    myfile=None
    #triggerHz=20500*2/16
    triggerHz=60
    fileopen=0
    printonce=0
    zero_histos=0
    psum_scale=10.0
    ad_scale=10.0
    zero_histos_toggle=0
    bm1counts=0
    gcounts=0
    grate=0.0
    global_acquireflag=0
    psum3Dhisto=numpy.zeros((512,512,400),dtype=numpy.uint32)
    psum_histo_xaxis=numpy.zeros(1000,dtype=numpy.uint32)
    twoD_histo=numpy.zeros((512,512),dtype=numpy.uint32)
    psum_histox=numpy.zeros((1000,16),dtype=numpy.uint32)  #density will be ten units.
    psum_histoy=numpy.zeros((1000,16),dtype=numpy.uint32)  #density will be ten units.
    ad_histo_xaxis=numpy.zeros(500,dtype=numpy.uint32)  #density will be ten units.
    for i in range(1000):
        psum_histo_xaxis[i]=psum_scale*i
    for i in range(500):
        ad_histo_xaxis[i]=ad_scale*i
    tof_histo=numpy.zeros(1700,dtype=numpy.uint32)  #density is 10 usec
    usec=1.0e6/triggerHz
    finebins=int(usec*20)
    tof_histo_fine=numpy.zeros(finebins,dtype=numpy.uint32) #density is 0.1 usec
    tof_histo_fine_xaxis=numpy.zeros(finebins)
    for i in range(finebins):
        tof_histo_fine_xaxis[i]=0.1*i
    ad_histox=numpy.zeros((500,8,4),dtype=numpy.uint32) 
    ad_histoy=numpy.zeros((500,8,4),dtype=numpy.uint32) 
#    firstten=numpy.zeros((11,16),dtype=numpy.float32)
#    firstin=0
    app = QtWidgets.QApplication(sys.argv)
    window = MyApp()
    window.show()
    tx=app.exec_()
    exit(tx)
    
