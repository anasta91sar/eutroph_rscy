#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 12:28:31 2017

@author: anastasia
"""

import os 
import numpy as np
from osgeo import gdal
from osgeo import ogr

def imread(image):
#	img = gdal.Open(image)
#	return img.ReadAsArray(), img.GetProjection(), img.GetGeoTransform(), img.Rastercount 
    img = gdal.Open(image)
    return img.ReadAsArray()

def imwrite (fileName, frmt, projection, geotransform, data) :
    drv = gdal.GetDriverByName(frmt)
    rows = data.shape[1]
    cols = data.shape[0]
    out = drv.Create(fileName, rows, cols, 1, gdal.GDT_Float32)
    band = out.GetRasterBand(1)
    band.WriteArray(data)
    band = None
    out.SetProjection(projection)
    out.SetGeoTransform(geotransform)
    out = None 

#reading AOI's bbox coordinates
coords = open('bbox.txt')
lines = coords.read().replace('\n','').split(',')
bbox=np.zeros((4,1))

for i in range (0,bbox.shape[0]):
        bbox[i] = lines[i]
        
line =lines[4]
del i,lines
 
months = {'January':'01','February':'02','March':'03','April':'04','May':'05','June':'06','July':'07','August':'08','September':'09','October':'10','November':'11','December':'12'}
    
    
#ALL DATA:
#os.system("python2 /home/anastasia/motu-client-python-1.4.00-20170410143941999/src/python/motu-client.py -u asarelli -p AnastasiaCMEMS2017 -m http://puertos2.cesga.es/motu-web/Motu -s IBI_ANALYSIS_FORECAST_BIO_005_004-TDS -d dataset-ibi-analysis-forecast-bio-005-004-daily -x {0} -X {2} -y {1} -Y {3} -t '2017-'{4}'-01' -T '2017-'{4}'-31' -z 0.493 -Z 0.4942 -v phy -v fer -v no3 -v chl -v eup -v sil -v nh4 -v po4 -v prp -v oxi -o ./ -f output.nc".format(float(bbox[0]),float(bbox[1]),float(bbox[2]),float(bbox[3]),months[line]))

#os.system("python2 /home/anastasia/motu-client-python-1.4.00-20170410143941999/src/python/motu-client.py -u asarelli -p AnastasiaCMEMS2017 -m http://puertos2.cesga.es/motu-web/Motu -s IBI_ANALYSIS_FORECAST_BIO_005_004-TDS -d dataset-ibi-analysis-forecast-bio-005-004-daily -x {0} -X {2} -y {1} -Y {3} -t '2017-'{4}'-01' -T '2017-'{4}'-31' -z 0.493 -Z 0.4942 -v no3 -v chl -v eup -v sil -v po4 -o ./ -f output.nc".format(float(bbox[0]),float(bbox[1]),float(bbox[2]),float(bbox[3]),months[line]))

## extract data layers

os.system('gdal_translate -of GTiff NETCDF:output.nc:sil /tmp/sil.tif')
#os.system('gdalwarp -te -1.592 50.313 5.160 54.331 -r bilinear -of GTiff NETCDF:output.nc:sil sil.tif')
os.system('gdal_translate -of GTiff NETCDF:output.nc:no3 /tmp/no3.tif')
#os.system('gdalwarp -te -1.592 50.313 5.160 54.331 -r bilinear -of GTiff NETCDF:output.nc:no3 no3.tif')
os.system('gdal_translate -of GTiff NETCDF:output.nc:po4 /tmp/po4.tif')
#os.system('gdalwarp -te -1.592 50.313 5.160 54.331 -r bilinear -of GTiff NETCDF:output.nc:po4 po4.tif')
os.system('gdal_translate -of GTiff NETCDF:output.nc:chl /tmp/chlorophyl.tif')
#os.system('gdalwarp -te -1.592 50.313 5.160 54.331 -r bilinear -of GTiff NETCDF:output.nc:chl chlorophyl.tif')
os.system('gdal_translate -of GTiff NETCDF:output.nc:eup /tmp/eup_zone.tif')
#os.system('gdalwarp -te -1.592 50.313 5.160 54.331 -r bilinear -of GTiff NETCDF:output.nc:eup eup_zone.tif')

# open datasets

img = gdal.Open('/tmp/sil.tif')
sil =  img.ReadAsArray()
nitr = imread('/tmp/no3.tif')
pho = imread('/tmp/po4.tif')
chl = imread('/tmp/chlorophyl.tif')
eudepth = imread('/tmp/eup_zone.tif')


# calculate dataset extents
gt = img.GetGeoTransform()

pixelX = gt[1]
pixelY = -gt[5]

width = img.RasterXSize
height = img.RasterYSize

minx = gt[0]
miny = gt[3] + width*gt[4] + height*gt[5]
maxx = gt[0] + width*gt[1] + height*gt[2]
maxy = gt[3]


#read shapefiles:
##Exclusive Economic Zones (200 nautical miles (nm) buffer from coastline) 
##Territorial Zones (12 nm buffer from coastline)
##Contiguous Zones (24 nm buffer from coastline)

shp='/home/anastasia/Desktop/eutrophication-indicators/Maritime_Boundaries/contiguous_zones_iberia/shp/z24.shp'

driver = ogr.GetDriverByName('ESRI Shapefile')

source = driver.Open(shp)
layer = source.GetLayer()
featCount = layer.GetFeatureCount()

os.system('mkdir zone24')

path1 = '/home/anastasia/Desktop/eutrophication-indicators/data_{0}'.format(line)

fiLe = open('{0}/stats24_{1}.txt'.format(path1,line),'w+')
   
for feature in layer:
    print feature.GetField('GeoName')
    ID = feature.GetField('PolygonID')
    
    fiLe.write('{0}\nEEZ ID = {1}\n\nDay,No-problem,Tendency,Possibility,Problem areas \n\n'.format(feature.GetField('GeoName'),ID))
    
    os.system('mkdir {0}/zone24/{1}'.format(path1,ID))
    
    path2 = '/home/anastasia/Desktop/eutrophication-indicators/data_{1}/zone24/{0}'.format(ID,line)
        
    zone = imread('/home/anastasia/Desktop/eutrophication-indicators/Maritime_Boundaries/contiguous_zones_iberia/tif/contiguous_{0}.tif'.format(ID))
    zone = np.where(zone==0,9999,1)
        
    eutr_class = np.zeros(eudepth.shape)
    stats = np.zeros((img.RasterCount,4))
    
    for i in range (0, img.RasterCount):
        nitr_ = np.nan_to_num(nitr[i])
        nitr_ = nitr_*0.1+300
        pho_ = np.nan_to_num(pho[i])
        pho_ = pho_*0.0099999998 + 20
            
        eutr1 = np.zeros(sil[i].shape)
        eutr2 = np.zeros(sil[i].shape)
        eutr1 = np.where((nitr_>=sil[i]),1,eutr1)
        eutr2 = np.where((pho_>=sil[i]),1,eutr2)
        eutr = eutr1 + eutr2
        eutr = np.where(sil[i]<-1,-1,eutr)
        #imwrite('eutr_nutr.tif', 'GTiff', sil[1], sil[2], eutr) 
        
        
        ## classification of chlorophyll concentration
        ##according to a case study in baltic sea
        
        chl_ = np.nan_to_num(chl[i])
        chl_ = chl_*0.1
        chl_ = np.where((chl_ < 2.0 ),0,chl_)
        chl_ = np.where((chl_ >= 2.0) & (chl_< 3.2 ),1,chl_)
        chl_ = np.where((chl_ >= 3.2), 2, chl_)
        chl_ = np.where(sil[i]<-1,-1,chl_)
        #imwrite('eutr_chl.tif', 'GTiff', chl[1], chl[2], chl_) 
        
        
        zd_theor1 = eudepth[i]*0.1/2.81                       #according to French et al,1982
        zd_theor2 = (eudepth[i]*0.1-2.6626)/1.9322            #according to Lee et al, 2007
        
        zd_theor1 = np.where((zd_theor1 <= 0)|(zd_theor1 > 6),0,zd_theor1)
        zd_theor1 = np.where((zd_theor1 <= 6) & (zd_theor1 > 4),1,zd_theor1)
        zd_theor1 = np.where((zd_theor1 <= 4) & (zd_theor1 > 0),2,zd_theor1)
        zd_theor1 = np.where(sil[i]<-1,-1,zd_theor1)
        #imwrite('depth_class1.tif', 'GTiff', depth[1], depth[2], zd_theor) 
    
        zd_theor2 = np.where((zd_theor2 <= 0)|(zd_theor2 > 6),0,zd_theor2)
        zd_theor2 = np.where((zd_theor2 <= 6) & (zd_theor2 > 4),1,zd_theor2)
        zd_theor2 = np.where((zd_theor2 <= 4) & (zd_theor2 > 0),2,zd_theor2)
        zd_theor2 = np.where(sil[i]<-1,-1,zd_theor2)
        
        eutr_test = eutr+chl_+zd_theor2
        
        ## final classification:
        ## 0  -->   no-problem areas
        ## 1  -->   tendency for eutrophication events
        ## 2  -->   possibility of eutrophication events
        ## 3  -->   problem areas
        
        eutr_test = np.multiply(zone,eutr_test)
        
        eutr_class[i] = np.where((eutr_test<=2)&(zone==1),0,eutr_test)
        eutr_class[i] = np.where((eutr_class[i]==3),1,eutr_class[i])
        eutr_class[i] = np.where((eutr_class[i]>3)&(eutr_class[i]<=5),2,eutr_class[i])
        eutr_class[i] = np.where((eutr_class[i]==6),3,eutr_class[i])
        eutr_class[i] = np.where((eutr_class[i]<=0)&(zone==9999),-1,eutr_class[i])
        eutr_class[i] = np.where((eutr_class[i]>6),-1,eutr_class[i])
        
        stats[i] = [int(np.count_nonzero(eutr_class[i]==0)),int(np.count_nonzero(eutr_class[i]==1)),int(np.count_nonzero(eutr_class[i]==2)),int(np.count_nonzero(eutr_class[i]==3))]
        
        os.system('cd {0}'.format(path1))
        fiLe.write('{0},{1},{2},{3},{4}\n'.format(i+1,stats[i][0],stats[i][1],stats[i][2],stats[i][3]))
     
        imwrite('{0}/eutr_class{1}.tif'.format(path2,i+1), 'GTiff', img.GetProjection(), img.GetGeoTransform(), eutr_class[i]) 
    
    sum_pixels = np.sum(stats[1,:])
    class0 = np.mean(stats[:,0])/sum_pixels
    class1 = np.mean(stats[:,1])/sum_pixels
    class2 = np.mean(stats[:,2])/sum_pixels
    class3 = np.mean(stats[:,3])/sum_pixels

    icep = round((class1+class2+class3)*100)
    fiLe.write('\n ICEP = {0}% \n\n ----------- \n\n'.format(icep))
    
fiLe.close()
