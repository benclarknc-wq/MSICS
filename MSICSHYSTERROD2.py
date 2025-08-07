# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 15:01:18 2023

@author: bencl
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
#%%
imagesize=100 #number of pixels per side of image
rowstart=0 #the starting row of the region of intrest
collumnstart=0 #the starting column region of intrest
ROIheight=100 #the number of rows in the region of intrest
ROIlength=100 #the number of collumns in the region of intrest
Sfactor=15 #the shift amount for images must be odd
K3=1.96 #hysteresis constant determined for particular image parameters
d=3.0 #hysteresis constant for associated shift
pixelsize=30.0e-9 #pixel side length
taup=1.7e-6 #pixel dwell time
ecc=3.333 #the aspect ratio of the confocal volume
w0=0.353e-6 #the lateral width of a 3D gaussian confocal volume
tauT=5.0e-6 #the triplet state decay time

#%%
"""
Hysteresis Correction

The ideal pixel array correction based on a shifted sinesoidal model of bidirectional hysteresis.
Every row is shifted slightly by a factor delta*K3*sin(pi*x/N) as well d number of pixels.
"""
x=np.linspace(0,imagesize,imagesize)
y=np.linspace(0,imagesize, imagesize)
x,y=np.meshgrid(x,y)
x[::2,:]=(x[::2,:]+K3*np.sin(np.pi*x[::2,:]/imagesize))*pixelsize
x[1::2,:]=(x[1::2,:]-K3*np.sin(np.pi*x[1::2,:]/imagesize)+d)*pixelsize
y=y*pixelsize
"""
Correlation Fitting Function

This function is a combination of two seporate parts. G(xi,psi) = S(xi,psi)*Gd(xi,psi)
Gd is a diffusion function which can take on a number of forms from pure diffusion,
triplet state diffusion, confirmational diffusion and protanation diffusion. In this 
algorithm the trilet state function is used. We have also transformed xi,psi into
delr and taud. delr is the displacement between two pixels described by xi,psi and 
any hysteresis correlation. taud is the time displacment between two individual
pixels, which is also well defined by xi and psi. The inputs for the function 'GDL'
are: GDL((delr,taud),Dif. Coeff., Mean Particle Num., Triplet State Pop. Prop.)
"""
def GDL(inp,D,N,T):
    delr,taud=inp
    return 0.3535*1.0/N*np.exp(-1.0*(delr**2)/(w0**2+4*D*(np.abs(taud))))/(1.0+4*D*(np.abs(taud))/w0**2)/np.sqrt(1.0+4*D*(np.abs(taud))/ecc**2/w0**2)*(1+T*(np.exp(-taud/tauT)-1.0))*(1.0-T)

#Flatten function which takes non normal nested arrays and flattens to a single axis array
def flat(lis):
    flatList = []
    # Iterate with outer list
    for element in lis:
        if type(element) is np.ndarray:
            # Check if type is list than iterate through the sublist
            for item in element:
                flatList.append(item)
        else:
            flatList.append(element)
    return flatList
#%%
"""
Data Preporation

Each image frame is opened and saved into a 3D array 'Gmstack' the files which 
denote each image frame should be csv files, where each file name is defined by 
'filepath'+3digit+'datafiletype'. We used .txt files. Files should be written 
into Gmstack.
"""
Gmstack=[]
for k in range(1,601):
    filename='C:\\Users\\bencl\\Documents\\LeBlanc Lab\\PythonFiles\\11142023\\MSICS_1.7us_'+f'{k:03}'+'.txt'
    intdata = pd.read_csv(filename,sep="\t",header=None,skiprows=3,skipfooter=imagesize+3,engine='python',encoding='unicode_escape')
    intdata=intdata.to_numpy()
    Odata=intdata
    Gmstack.append(Odata)
Gmstack=np.array(Gmstack)
#%%
"""
Immobile Removal

As part of the algorithm we can also remove immobile fluorescent features from 
the image stack. We have executed a box car method to do so which subtracts an 
average from each pixel based on a pixel average. 'immobilefactor' sets the number
of frames in the boxcar average, and 'immobilefactor2' sets the frame shift for the 
boxcar. These varibles must be even.
"""
shape=np.shape(Gmstack)
immobilefactor=20
immobilefactor2=10
AF1=np.array([])
ACFlen=len(Gmstack)-immobilefactor
for i in range(immobilefactor2,ACFlen+immobilefactor2):
    AF1=np.append(AF1,np.mean(Gmstack[i-immobilefactor2:i+immobilefactor2,:,:],axis=0))
AF1=np.reshape(AF1,(shape[0]-immobilefactor,shape[1],shape[2]))
braF=np.average(np.average(Gmstack[immobilefactor2:ACFlen+immobilefactor2,:,:],axis=1),axis=1)
#%%
"""
Image Correlation

For the set number of pixel shifts in the xi and psi direction we will execuate a 
correlation function. This function is defined by a product of the shifted and 
non shifted residuals of each pixel from the boxcar average divided by the frame 
average squared.

We have portioned out the correlation into three different sections which 
span only xi shifts, even only psi shifts, and odd only psi shifts. This could 
be modified to cover a larger set of correlation with both xi and psi shifts.
We found that taking these three portions of data into consideration yeilded
adaquate correlation results.
"""
AF2alist=[]
AF2blist=[]
AF2clist=[]
for i in range(0,Sfactor+1):
    if i%2==1:
        Sdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart:rowstart+ROIheight-i,collumnstart:collumnstart+ROIlength]-AF1[:ACFlen,rowstart:rowstart+ROIheight-i,collumnstart:collumnstart+ROIlength]
        NSdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart+i:rowstart+ROIheight,collumnstart:collumnstart+ROIlength]-AF1[:ACFlen,rowstart+i:rowstart+ROIheight,collumnstart:collumnstart+ROIlength]
        AF2a=np.average(np.multiply(Sdata,NSdata),axis=1)
        for k in range(len(braF)):
            AF2a[k,:]=AF2a[k,:]/braF[k]**2
        AF2a=np.average(AF2a,axis=0)
        AF2alist.append(AF2a)
    else:
        Sdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart:rowstart+ROIheight-i,collumnstart:collumnstart+ROIlength]-AF1[:ACFlen,rowstart:rowstart+ROIheight-i,collumnstart:collumnstart+ROIlength]
        NSdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart+i:rowstart+ROIheight,collumnstart:collumnstart+ROIlength]-AF1[:ACFlen,rowstart+i:rowstart+ROIheight,collumnstart:collumnstart+ROIlength]
        AF2c=np.average(np.average(np.multiply(Sdata,NSdata),axis=1),axis=1)
        AF2c=AF2c/(braF**2)
        AF2c=np.average(AF2c,axis=0)
        AF2clist.append(AF2c)
for j in range(1,Sfactor+1):
    Sdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart:rowstart+ROIheight,collumnstart:collumnstart+ROIlength-j]-AF1[:ACFlen,rowstart:rowstart+ROIheight,collumnstart:collumnstart+ROIlength-j]
    NSdata=Gmstack[immobilefactor2:ACFlen+immobilefactor2,rowstart:rowstart+ROIheight,collumnstart+j:collumnstart+ROIlength]-AF1[:ACFlen,rowstart:rowstart+ROIheight,collumnstart+j:collumnstart+ROIlength]
    Sdata[:,::2, :] = Sdata[:,::2, ::-1]
    NSdata[:,::2, :] = NSdata[:,::2, ::-1]
    AF2b=np.average(np.multiply(Sdata,NSdata),axis=1)
    AF2b=AF2b/(braF**2)[:,None]
    AF2b=np.average(AF2b,axis=0)
    AF2blist.append(AF2b)
#%%
"""
Displacement and Timelag Calculation

We Calculate the spatial displacement based on xi, psi shift as well as any hysteresis
correction. As well as the time lag between pixels based on xi,psi. We order
all these values into ordered arrays of matching size along with the correlation 
data.
"""
displacement=[]
for i in range(1,Sfactor+1):
    displacement.append(x[0,i+collumnstart:collumnstart+ROIlength]-x[0,collumnstart:collumnstart+ROIlength-i])
for i in range(0,Sfactor+1):
    if i%2==1:
        displacement.append(np.sqrt((x[rowstart+i,collumnstart:collumnstart+ROIlength]-x[rowstart,collumnstart:collumnstart+ROIlength])**2+(i*pixelsize)**2))
    else:
        displacement.append(i*pixelsize)
Timelag=[]
for i in range(1,Sfactor+1):
    Timelag.append(taup*i*np.ones(ROIlength-i))
for i in range(0,Sfactor+1):
    if i%2==1:
        Timelag.append(np.linspace(imagesize*(i+1)-1,imagesize*(i+1)-2*imagesize+1,imagesize)[collumnstart:collumnstart+ROIlength][::-1]*taup)
    else:
        Timelag.append(i*imagesize*taup)

CorrelationV=AF2blist
for i in range(0,Sfactor+1):
    if i%2==1:
        CorrelationV.append(AF2alist[int((i-1)/2)])
    else:
        CorrelationV.append(AF2clist[int(i/2)])
#%%
critnum=len(flat(displacement[0:Sfactor]))
timearray=flat(Timelag)
disarray=flat(displacement)
corarray=flat(CorrelationV)
plt.figure(dpi=500)
plt.scatter(timearray[critnum+1:],disarray[critnum+1:],c=corarray[critnum+1:],cmap='inferno')
plt.scatter(timearray[:critnum],disarray[:critnum],c=corarray[:critnum],cmap='inferno')
plt.xlabel("Timelag $\\tau$ (s)")
plt.ylabel("Displacement $\Delta$r (m)")
plt.xscale('log')
plt.colorbar()
plt.title('G($\Delta r,\\tau$) vs TimeLag $\\tau$ and Displacement $\Delta r$')
plt.show()
#%%
plt.figure(dpi=500)
plt.scatter(timearray[critnum+1:],corarray[critnum+1:],c=disarray[critnum+1:],cmap='inferno_r')
plt.scatter(timearray[:critnum],corarray[:critnum],c=disarray[:critnum],cmap='inferno_r')
plt.xlabel("Timelag $\\tau$ (s)")
plt.ylabel("Correlation G($\\tau$,$\Delta r$)")
plt.xscale('log')
plt.title('G($\Delta r,\\tau$) vs TimeLag $\\tau$')
plt.show()
#%%
"""
Data Fitting

We fit the spatial displacement, time lag, and correlation data using GDL. We 
can fit to different sections based on only xi shifts, psi shifts and both of these
sets. We fit using the nonlinear Levenburg Marquardt algorithm with numerical
Jacobian. Reported are the diffusion coefficent, calculated N value, calculated
T value and the calculated concentration for each fitting set.
"""
Na=6.022e23
Veff=1.62e-15
popt,pcov=curve_fit(GDL, (flat(displacement[0:16]),flat(Timelag[0:16])),flat(CorrelationV[0:16]),p0=[64.0e-12, 2.19193849,0.55],method='lm')
Con=1/popt[1]/Na/Veff
print('total fit para '+str(popt[0]/10.0**(-12)))
print('total fit para '+str(popt[1]))
print('total fit para '+str(popt[2]))
print('total fit cov '+str(pcov))
print(Con)
popt,pcov=curve_fit(GDL, (flat(displacement[15:]),flat(Timelag[15:])),flat(CorrelationV[15:]),p0=[64.0e-12, 1.19193849,0.55],method='lm')
Con=1/popt[1]/Na/Veff
print('total fit para '+str(popt[0]/10.0**(-12)))
print('total fit para '+str(popt[1]))
print('total fit para '+str(popt[2]))
print('total fit cov '+str(pcov))
print(Con)
popt,pcov=curve_fit(GDL, (disarray,timearray),corarray,p0=[64.0e-12, 2.19193849,0.55],method='lm')
Con=1/popt[1]/Na/Veff
print('total fit para '+str(popt[0]/10.0**(-12)))
print('total fit para '+str(popt[1]))
print('total fit para '+str(popt[2]))
print('total fit cov '+str(pcov))
print(Con)
#%%
plt.figure(dpi=500)
plt.scatter(timearray[critnum+1:],corarray[critnum+1:],c=disarray[critnum+1:],cmap='inferno_r')
plt.plot(np.array(timearray[critnum+1:]),GDL((np.array(disarray[critnum+1:]),np.array(timearray[critnum+1:])),popt[0],popt[1],popt[2]),c='b')
plt.scatter(timearray[:critnum],corarray[:critnum],c=disarray[:critnum],cmap='inferno_r')
plt.plot(np.array(timearray[:critnum]),GDL((np.array(disarray[:critnum]),np.array(timearray[:critnum])),popt[0],popt[1],popt[2]),c='b')
plt.xlabel("Timelag $\\tau$ (s)")
plt.ylabel("Correlation G($\\tau$,$\Delta r$)")
plt.xscale('log')
# plt.colorbar()
plt.title('G($\Delta r,\\tau$) vs TimeLag $\\tau$')
plt.show()
#%%
plt.figure(dpi=500)
plt.scatter(np.array(timearray),np.array(disarray),c=GDL((np.array(disarray),np.array(timearray)),popt[0],popt[1],popt[2]),cmap='inferno')
plt.xlabel("Timelag $\\tau$ (s)")
plt.ylabel("Displacement $\Delta$r (m)")
plt.xscale('log')
plt.colorbar()
plt.title('G($\Delta r,\\tau$) vs TimeLag $\\tau$ and Displacement $\Delta r$')
plt.show()