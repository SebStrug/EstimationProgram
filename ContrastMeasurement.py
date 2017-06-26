# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 12:19:15 2017
@author: Sebastian
"""

import warnings

warnings.filterwarnings("ignore",".*GUI is implemented.*")

from PIL import Image, ImageDraw
from pylab import ginput
import matplotlib.pyplot as plt
import numpy as np
import time
import sys, os
from Tkinter import Tk
from tkFileDialog import askopenfilename

#define how big the spacing is to be
spacing = (5000)
#define the number of points to average over in your contrast line
avgNo = (100)
materialName = 'Fluorophlogopite' #change the name for your saved file

##Convert the image to RGB if it is a .gif for example
##imag = imag.convert ('RGB')
plt.figure(1)

"""Using PIL library, messes up conversion to greyscale"""
Tk().withdraw()
impath = askopenfilename()
imag = Image.open(impath) #open the image normally
plt.figure(1)
imgplot = plt.imshow(imag)
plt.show()

print "\nPlease click three times on the substrate"
substratePts = ginput(3)
substratePts = [list(i) for i in substratePts]
for i in range(0,3):
    for j in range(0,2):
        substratePts[i][j] = np.round(substratePts[i][j])
"""Find the RGB data for each pixel along this line, for the standard image"""
subRGB = []
subR = np.zeros(3)
subG = np.zeros(3)
subB = np.zeros(3)
subLuminance = np.zeros(3)
subX = []
subY = []
for tuples in substratePts:
    subX.append(tuples[0])
    subY.append(tuples[1])
for i in range(0,3):
    subRGBget = imag.getpixel((subX[i],subY[i])) #get the RGB of each pixel along the line
    subRGB.append(subRGBget) #list of all RGB
    subR[i] = subRGB[i][0]
    subG[i] = subRGB[i][1]
    subB[i] = subRGB[i][2]
    subLuminance[i] = ((0.2126*subR[i]) + (0.7152*subG[i]) + (0.0722*subB[i])) #Luminance (standard for certain colour spaces): 
substrateLuminance = np.average(subLuminance)
 
print "Now either zoom in on your flake, or click anywhere to skip this step"
plt.waitforbuttonpress() 
print "Click two times to draw a line across your flake, starting from the silicon substrate\n"
pts = ginput(2) #take two points
pts = [list(i) for i in pts] #convert tuples to a list of lists
#round the values taken from ginput
for i in range(0,2):
    for j in range(0,2):
        pts[i][j] = np.round(pts[i][j])
x1 = pts[0][0] #assign the points to x and y values
x2 = pts[1][0]
y1 = pts[0][1]
y2 = pts[1][1]  

"""Calculate a line for the two points chosen by ginput"""
def calculateLine(x1, x2, y1, y2):
    m = (y2-y1)/(x2-x1) # gradient
    c = y2 - (m*x2) #intercept
    return m, c
def calculateThickLine(x,y,m): #to calculate the thick line
    return y-(m*x) #return the intercept of the new line
    
m, c = calculateLine(x1, x2, y1, y2)
x = np.linspace(x1,x2,spacing)
y = (m*x) + c

#need a 'thickness' number of (x,y) coords for each y coord.
for i in range(0,len(y)):
    foobar=4

thickness = 10 #no. of pixels the line is thick
thicknessArr = np.arange((-thickness/2),(thickness/2))
thickLine = np.zeros((2,len(y),thickness)) #2 because one for x, one for y
for i in range(0,thickness):
    for j in range(0,len(y)):
        for k in range(0,2):
            #let 0 by x, 1 be y
            thickLine[0][j][i] = x[j] + ((-1/m)*thicknessArr[i])
            thickLine[1][j][i] = y[j] + (-m*thicknessArr[i])
            
#these configures the points to draw the thicker line over
for i in range(0,thickness):
    globals()['x1Thick%s' % i] = x1+ ((-1/m)*(thicknessArr[i]))
    globals()['x2Thick%s' % i] = x2+ ((-1/m)*(thicknessArr[i]))
    globals()['y1Thick%s' % i] = y1+ ((-m)*(thicknessArr[i]))
    globals()['y2Thick%s' % i] = y2+ ((-m)*(thicknessArr[i]))

"""Find the RGB data for each pixel along this line, for the standard image"""
pixelRGB = []
R = np.zeros(spacing)
G = np.zeros(spacing)
B = np.zeros(spacing)
luminance = np.zeros(spacing)
brightness = np.zeros(spacing)
for i in range(0,spacing):
    pixelRGBget = imag.getpixel((x[i],y[i])) #get the RGB of each pixel along the line
    pixelRGB.append(pixelRGBget) #list of all RGB
    R[i] = pixelRGB[i][0]
    G[i] = pixelRGB[i][1]
    B[i] = pixelRGB[i][2]
    luminance[i] = ((0.2126*R[i]) + (0.7152*G[i]) + (0.0722*B[i])) #Luminance (standard for certain colour spaces): 
    brightness[i] = np.sum([R[i],G[i],B[i]])/3 #this is another, simpler version of brightness, test both
    
"""Other formulas for the luminance::
    Standard: LuminanceA = (0.2126*R) + (0.7152*G) + (0.0722*B)
    Percieved A: LuminanceB = (0.299*R + 0.587*G + 0.114*B)
    Perceived B, slower to calculate: LuminanceC = np.sqrt( 0.241*R^2 + 0.691*G^2 + 0.068*B^2 )
"""

"""Show the image with a line drawn"""
plt.close()
plt.figure(2)
draw = ImageDraw.Draw(imag) 
draw.line((x1, y1, x2, y2), fill=(255,255,255))
#The following command draws the thicker line
#for i in range(0,thickness):
#    draw.line((globals()['x1Thick%s' % i],globals()['y1Thick%s' % i],globals()['x2Thick%s' % i],globals()['y2Thick%s' % i]), (255,255,255))
imag.save("Temporary_Image.jpg") #temporarily save the file
img3plot = plt.imshow(imag)
os.remove("Temporary_image.jpg") #remove the saved file from the system

"""Now we have the luminance/brightness, we can calculate the contrast.
Take this as the difference between the luminance at some point and the luminance
    at the start of the line (which should be on the silicon substrate)."""
#set initial contrast as 0 and then get contrast from that point
contrast = np.zeros(len(luminance))
contrast_2 = np.zeros(len(brightness))
contrast_test = np.zeros(len(brightness))

#luminance[0] = substrateLuminance #change the zero factor to be an average of points on the substrate.

for i in range(1,len(luminance)): 
    contrast[i] = (substrateLuminance - luminance[i])/substrateLuminance   
    
    """Alternative contrast definition: the Michelson contrast,
        contrast[i] = (luminance[i]-luminance[0])/(luminance[i]+luminance[0]) """

contrast[0] = 0 #have to set initial contrast to 0, cannot do this in previous loop
     
#general number of points averaged
avgContrast = np.zeros(len(contrast)/avgNo)
for i in np.linspace(0,len(contrast)-avgNo,len(contrast)/avgNo):
    for j in np.linspace(avgNo,len(contrast),len(contrast)/avgNo):
        if j == i+avgNo:
            avgContrast[int(i)/avgNo] = float(np.sum(contrast[int(i):int(j)])/avgNo)

timestr = time.strftime("{}_%Y%m%d-%H%M".format(materialName)) #current date and time
"""Write these two contrast arrays to a file"""
#if we want to save to a different directory:
scriptDir = sys.path[0]
contrastAltDir = os.path.join(scriptDir, "ContrastFiles/contrast_{}.txt")
avgContrastAltDir = os.path.join(scriptDir, "ContrastFiles/avgContrast_{}")
np.savetxt(contrastAltDir.format(timestr),contrast, fmt='%1.4e')
np.savetxt(avgContrastAltDir.format(timestr),avgContrast, fmt='%1.4e')

fig, ax = plt.subplots(1,1)
l1 = ax.plot(range(0,spacing),contrast,label='No average')
l2 = ax.plot(np.linspace(0,spacing,spacing/avgNo), avgContrast, label = '%s point average' %(avgNo))
plt.legend(); plt.xlabel('Displacement along line drawn'); plt.ylabel('Average contrast')

"""Now we can use this contrast to calculate the thickness"""
plt.show()