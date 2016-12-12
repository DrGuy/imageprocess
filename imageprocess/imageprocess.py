#/usr/bin/python
import os, datetime, time, shutil, sys, glob, csv, ENVIfile, modistools, platform, numpy, numexpr
from xml.dom import minidom
from subprocess import Popen
from ENVIfile import *
computername = platform.node()
global fmaskdir, srdir, btdir, ingestdir, ndvidir, evidir, archdir, logdir, WRS1, WRS2, defaulterrorfile

if sys.version_info[0] == 2:
    import ConfigParser as configparser
else:
    import configparser

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'config', 'imageprocess.ini')
print(config_path)
# cfgfile = 'imageprocess.ini'
if os.path.isfile(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)
    fmaskdir = config['DEFAULT']['fmaskdir']
    srdir = config['DEFAULT']['srdir']
    btdir = config['DEFAULT']['btdir']
    ingestdir = config['DEFAULT']['ingestdir']
    ndvidir = config['DEFAULT']['ndvidir']
    evidir = config['DEFAULT']['evidir']
    archdir = config['DEFAULT']['archdir']
    logdir = config['DEFAULT']['logdir']
    WRS1 = config['DEFAULT']['WRS1'] # WRS-1, Landsats 1-3
    WRS2 = config['DEFAULT']['WRS2'] # WRS-2, Landsats 4-8
    defaulterrorfile = os.path.join(logdir, 'errors.csv')
else:
    print('ERROR: no configuration data found.')
    print(sys.argv)
    print(sys.argv[0])
    sys.exit()
    
if not 'linux' in sys.platform: # this way I can use the same library for processing on multiple systems
    if sys.version_info[0] !=3: # Attempt to load ArcPy and EnviPy libraries, if not, use GDAL.
        try:
            from arcenvipy import *
        except:
            print('There was an error loading either ArcPy or EnviPy. Functions requiring this library will not be available.')
   
# Import GDAL
    from osgeo import gdal, osr

else:
    
    import gdal
    sys.path.append('/usr/bin')
    sys.path.append('/usr/local/bin')
    

def logerror(f, message, *args, **kwargs):
    errorfile = kwargs.get('errorfile',defaulterrorfile) 
    dirname, basename = os.path.split(errorfile)
    if not os.path.isdir(dirname):
        errorfile = os.path.join(logdir, basename)
    if not os.path.exists(errorfile):
        with open(errorfile,'w') as output:
            output.write('Time, File, Error\n')
    now = datetime.datetime.now()
    with open(errorfile, 'a') as output:
        output.write('%s, %s, %s\n'%(now.strftime('%Y-%m-%d %H:%M:%S'), f, message))
    

def extract_xml(s):
    '''Take a string like
    "<tag>         something          </tag>"
    drop the opening and closing tags and leading or trailing whitespace.'''
    start = s.index(">") + 1
    end  = s.rindex("<")
    s = str (s[start:end])
    return (s.strip())

def get_landsat_fileparams(filename):
    import time
    filename=os.path.basename(filename)
    if '_' in filename: # Extract data from MTL filename
        fileparts=filename.split('_')
        if len(fileparts)!=2:  # This is for old format tar.gz MTL files from the USGS that included the acquisition date in YYYYMMDD format.
            b1, b2=fileparts[0], fileparts[1]
            landsat, path, row = b1[1:2], b1[-6:-3], b1[-3:]
            year, yearmonthday = b2[-8:-4], b2[-8:]
            datetuple=time.strptime(yearmonthday,"%Y%m%d")
            julian=time.strftime("%j",datetuple)
    else:    # This will work with a simple scene ID
        landsat, path, row, year, julian = filename[2:3], filename[3:6], filename[6:9], filename[9:13], filename[13:16]
        datetuple=time.strptime(year+julian,'%Y%j')
        yearmonthday=time.strftime('%Y%m%d',datetuple)
    return landsat, path, row, year, julian, yearmonthday



## Ireland specific functions    

def reproject_ITM(in_raster, out_raster, *args, **kwargs): # Converts raster to Irish Transverse Mercator
    rastertype = kwargs.get('rastertype', None)
    landsat = kwargs.get('landsat', None)
    sceneid = kwargs.get('sceneid', None)
    outdir = kwargs.get('outdir', None)
    rewriteheader = kwargs.get('rewriteheader', True)
    # out_raster_ds = out_raster
    if os.access(in_raster, os.F_OK):
        src_ds = gdal.Open(in_raster)
        gt = src_ds.GetGeoTransform()
        # src_srs = osr.SpatialReference()
        # src_srs.ImportFromWkt(src_ds.GetProjectionRef())
        # extent = [gt[0],gt[3],gt[0]+gt[1]*src_ds.RasterXSize,gt[3]+gt[5]*src_ds.RasterYSize]
        # 
        # # Reproject to EPSG: 2157
        # dest_srs = osr.SpatialReference()
        # dest_srs.ImportFromEPSG(2157)
        # 
        # warpoptions = gdal.WarpOptions(format = 'ENVI', xRes = abs(gt[1]), yRes = abs(gt[5]), srcSRS = src_srs, dstSRS =  dest_srs)
        # gdal.Warp(out_raster_ds, in_raster, options = warpoptions)
        
        # Close datasets
        src_ds = None
        # out_raster_ds = None
        # dest_ds = None
        p = Popen(['gdalwarp','-t_srs','EPSG:2157','-tr', str(abs(gt[1])), str(abs(gt[5])),'-of','ENVI',in_raster,out_raster])
        print(p.communicate())
        if rewriteheader:
            ENVIfile(out_raster, rastertype, SceneID = sceneid, outdir = outdir).WriteHeader()

def parsecsv(incsv, filelist=None, returnscenelist=False): # Only useful for Ireland
    urls=[]
    wrslist=[]
    scenelist=[]
    WRS2P=['205','205','205','206','206','206','206','207','207','207','207','208','208','208','208','209','209','209']
    WRS2R=['22','23','24','21','22','23','24','21','22','23','24','21','22','23','24','22','23','24']
    n=0
    for x in WRS2P:
        wrslist.append([WRS2P[n],WRS2R[n]])
        n+=1
    csvfile=open(incsv)
    reader=csv.DictReader(csvfile)
    print('Reading from '+incsv+' for scenes to download.')
    for line in reader:
        row=line['Frame'].strip('"')
        path=line['Track'].strip('"')
        starttime=line['Start'].strip('"')
        datetuple=time.strptime(starttime[0:10],'%Y-%m-%d')
        yearjulian=time.strftime('%Y%j',datetuple)
        if not ((path=='205' or path=='209') and row=='21'):
            if any(path and row in x for x in wrslist):
                if line['Mission']=='Landsat-5':
                    sensor='LT5'
                elif line['Mission']=='Landsat-4':
                    sensor='LT4'
                elif line['Mission']=='Landsat-7':
                    sensor='LE7'
                else:
                    sensor=''
                filebase=sensor+path+'0'+row+yearjulian
                if returnscenelist:
                    scenelist.append(filebase)
                if filelist==None:
                    urls.append(line['Download_URL'].strip('"'))
                elif not any(filebase in y for y in filelist):
                    urls.append(line['Download_URL'].strip('"'))
                    print('Found scene '+filebase+', writing to output.')
            
    csvfile.close()
    numurls=len(urls)
    print('Found '+str(numurls)+' scenes for download from EOLI.')
    if returnscenelist:
        return urls, numurls, scenelist
    else:
        return urls, numurls
    
def makeEOLIdownloadhtml(filename, urllist):
    print('Writing %s to %s'%(os.path.basename(filename),os.path.dirname(filename)))
    with open(filename,'w') as output:
        output.write('<html>\n<head><title>Download URLs</title></head>\n<body>\n<ol>\n<p>\n<a href="https://eo-sso-idp.eo.esa.int/idp/AuthnEngine">Login here first, then navigate back to this page</a></p>\n<p>')
        for url in urllist:
            output.write('<li><a href="%s">Link</a></li>\n'%url)
        output.write('</ol>\n</p>\n</body>\n</html>')
    print('HTML file generated.')

def checkscenelocation(scene,dst=50.0): # This function assessed geolocation accuracy of scenes warped to ITM
    misplaced=False # Boolean value for whether scene centre fits in acceptable tolerances  
    basename=os.path.basename(scene)
    if 'lndsr.' in basename:
        basename.replace('lndsr.','')
    if int(basename[2:3])<4: # Determine WRS type, Path, and Row
        WRS=1
        polygon=WRS1
    else:
        WRS=2
        polygon=WRS2
    path=int(basename[3:6])
    row=int(basename[7:9])
    # Get scene centre coordinates
    print('Checking scene centre location accuracy for %s centre to within %.1f km of WRS-%d Path %d Row %d standard footprint centre.'%(basename, dst, WRS, path, row))
    
    src_ds = gdal.Open(scene)
    geoTrans=src_ds.GetGeoTransform()
    latmax = geoTrans[3]
    latmin = latmax + geoTrans[5]*float(src_ds.RasterYSize)
    longmin = geoTrans[0]
    longmax= longmin + geotrans[1]*float(src_ds.RasterXSize)
    src_ds=None
    X=(float(longmax.getOutput(0))+float(longmin.getOutput(0)))/2.
    Y=(float(latmax.getOutput(0))+float(latmin.getOutput(0)))/2.
    print('ITM scene centre coordinates are %.1f E, %.1f N.'%(X,Y))
    
    driver = ogr.GetDriverByName("FileGDB")
    gdb,wrs=os.path.split(WRS)
    ds=driver.Open(gdb,0)
    layer=gs.Getlayer(wrs)
    found=False
    while not found:
        feature=layer.GetNextFeature()
        items=feature.items()
        if path==items['PATH'] and row==items['ROW']:
            geometry=feature.geometry()
            envelope=geometry.GetEnvelope()
            wX=(envelope[0]+envelope[1])/2.
            wY=(envelope[2]+envelope[3])/2.
            found=True
    ds=None
    
    print('ITM standard WRS-%d footprint centre coordinates are %.1f E, %.1f N.'%(WRS,wX,wY))
    offset=(((X-wX)**2+(Y-wY)**2)**0.5)/1000 # determine distance in km between scene and standard footprint centres
    print('Offset = %.1f km out of maximum distance of %.1f km.'%(offset,dst))
    if dst>=offset:
        print('Scene %s is appropriately placed, and is %.1f km from the standard WRS-%d scene footprint centre.'%(basename,offset,WRS))
    else:
        print('Scene %s is improperly located, and is %.1f km from the standard WRS-%d scene footprint centre.'%(basename,offset,WRS))
        misplaced=True
    return misplaced, offset
    
## File handling functions

def isfilepresent(url, filelist): # This checks to see if a file to be downloaded from an URL exists locally
    i=url.rfind('/')
    j=url.rfind('-')
    downloaded=False
    filedata=url[i+1:j]
    if any(filedata in x for x in filelist):
        downloaded=True
    return downloaded
    
def makeurllist(inurl, filelist): # This creates a list of URLs of scenes to be downloaded from USGS ESPA
    inpage=urllib.urlopen(inurl) # Reads input URL page
    urls=[]
    for line in inpage:
        if "a href='/status/" in line:
            i=line.find('/')
            j=line.find('>')-1
            searchurl='https://espa.cr.usgs.gov'+line[i:j]
            print('Found ordered data in '+searchurl+', looking for files that have not been downloaded.')
            subpage=urllib.urlopen(searchurl)
            for x in subpage:
                if 'tar.gz' in x:
                    lpos=x.find('http:')
                    rpos=x.find('tar.gz')+6
                    suburl=x[lpos:rpos]
                    downloaded=isfilepresent(suburl, filelist)
                    if not downloaded:
                        urls.append(suburl)
                        print('Found undownloaded data in '+suburl)
    numurls=len(urls)
    print('A total of '+str(numurls)+' are to be downloaded.')
    return urls, numurls

def makefilelist(scandir):
    scandirs=glob.glob(scandir+'/*')
    filelist=[]
    numfilesondisk=0
    for dir in scandirs:
        print('Searching '+dir+' for files on disk.')
        dirlist=glob.glob(dir+'/mask.dat')
        if len(dirlist)>0:
            for file in dirlist:
                filelist.append(file)
                numfilesondisk+=1
    print('Total files downloaded: '+str(numfilesondisk))
    return filelist    

## Landsat import and VI calculation functions

def envihdracqtime(hdr):
    acqtime = None
    with open(hdr, 'r') as lines:
        for line in lines:
            if line.startswith('acquisition time'):
                acqtime = line
    return acqtime
    

def calcvis(refitm, *args, **kwargs): # This should calculate a masked NDVI.
    # fmaskdir = kwargs.get('fmaskdir', fmaskdir)
    # ndvidir = kwargs.get('ndvidir', ndvidir)
    # evidir = kwargs.get('evidir', evidir)
    dirname, basename = os.path.split(refitm)
    sceneid = basename[:21]
    acqtime = envihdracqtime(refitm.replace('.dat','.hdr'))
    fmaskfile = os.path.join(fmaskdir,'%s_cfmask.dat'%sceneid)
    usefmask = True
    if not os.path.exists(fmaskfile):
        fmaskfile = fmaskfile.replace('_cfmask.dat','_fmask.dat')
        if not os.path.exists(fmaskfile):
            print('ERROR: Fmask file does not exist, returning.')
            logerror(fmaskfile, 'File not found.')
            usefmask = False
        
    print('Calculating NDVI for scene %s.'%sceneid)
    refobj = gdal.Open(refitm)
    
    # Get file geometry
    geoTrans = refobj.GetGeoTransform()
    ns = refobj.RasterXSize
    nl = refobj.RasterYSize    
    
    if usefmask:
        fmaskobj = gdal.Open(fmaskfile)
        fmask = fmaskobj.GetRasterBand(1).ReadAsArray()
        # print(isinstance(fmask, numpy.ndarray))
    else: 
        print('Warning: No Fmask file found for scene %s.'%sceneid)
        fmask = None
    if basename[2:3] == '8':
        NIR = refobj.GetRasterBand(5).ReadAsArray()
        red = refobj.GetRasterBand(4).ReadAsArray()
        blue = refobj.GetRasterBand(2).ReadAsArray()
    else:
        NIR = refobj.GetRasterBand(4).ReadAsArray()
        red = refobj.GetRasterBand(3).ReadAsArray()
        blue = refobj.GetRasterBand(1).ReadAsArray()
    
    # NDVI calculation 
    NDVI = NDindex(NIR, red, fmask = fmask)
    ENVIfile(NDVI, 'NDVI', outdir = ndvidir, geoTrans = geoTrans, SceneID = sceneid, acqtime = acqtime).Save()
    NDVI = None
    
    # EVI calculation
    evi = EVI(blue, red, NIR, fmask = fmask)
    ENVIfile(evi, 'EVI', outdir = evidir, geoTrans = geoTrans, SceneID = sceneid, acqtime = acqtime).Save()
    evi = None
    
    NIR = None
    red = None
    refobj = None
    fmask = None
    fmaskobj = None

def EVI(blue, red, NIR, *args, **kwargs):
    fmask = kwargs.get('fmask', None)
    # dims = NIR.shape
    if not isinstance(fmask, numpy.ndarray):
        mask = numexpr.evaluate('(NIR < 0) | (NIR > 10000) | (red < 0) | (red > 10000) | (blue < 0) | (blue > 10000)') # masks exclude invalid pixels
    else:
        mask = numexpr.evaluate('(fmask > 0) | ((NIR < 0) | (NIR > 10000) | (red < 0) | (red > 10000) | (blue < 0) | (blue > 10000))') # reevaluate mask for EVI
    C1 = 6
    C2 = 7.5
    G = 2.5
    L = 1
    # EVI = numpy.zeroes(dims, dtype = numpy.float32)
    EVI = 10000 * (G * ((NIR - red)/(NIR+ C1 * red - C2 * blue + L)))
    EVI[mask] = 0.0 # replace invalid pixels with zero
    mask = None
    return EVI.astype(numpy.int16)

def NDindex(A, B, *args, **kwargs):
    fmask = kwargs.get('fmask', None)
    # dims = A.shape
    if not isinstance(fmask, numpy.ndarray):   
        mask = numexpr.evaluate('(A < 0) | (A > 10000) | (B < 0) | (B > 10000)')  # masks exclude invalid pixels
    else:
        mask = numexpr.evaluate('(fmask > 0) | ((A < 0) | (A > 10000) | (B < 0) | (B > 10000))')
    # data = numpy.zeroes(dims, dtype = numpy.float32)
    data =  10000 * ((A - B)/(A + B))#numpy.divide(numpy.subtract(A, B),numpy.add(A, B))
    data[mask] = 0.0 # replace invalid pixels with zero
    mask = None
    return data.astype(numpy.int16)

def importespa(f, *args, **kwargs):
    # srdir = kwargs.get('srdir', srdir)
    # btdir = kwargs.get('btdir', btdir)
    # fmaskdir = kwargs.get('fmaskdir', fmaskdir)
    # evidir = kwargs.get('evidir', evidir)
    # archdir = kwargs.get('archdir', archdir)
    overwrite = kwargs.get('overwrite', False)
    remove = kwargs.get('remove', False)
    btimg = None
    basename = os.path.basename(f)
    dirname = os.path.dirname(f)
    landsat = basename[2:3]
    if landsat == '8': 
        bands=['1', '2', '3', '4', '5', '6', '7']
    elif basename[1:2] == 'M':
        print('Landsat MSS is not supported yet, returning.')
        return
    else:
        bands=['1', '2', '3', '4', '5', '7']
    if f.endswith('.tar.gz'):
        if '-' in basename:
            i = f.rfind('-')
        else:
            i = f.find('.tar.gz')
        outputdir = f[:i]
        filelist = untarfile(f, outputdir)
    else:
        filelist = glob.glob(os.path.join(dirname, '*'))
        outputdir = dirname
    if filelist == 0 or len(filelist) == 0:
        print('ERROR: there is a problem with the files, skipping.')
        if len(filelist) == 0:
            logerror(f, 'No files found.')
        return
    for (i, fname) in enumerate(filelist):
        if '.img' in fname:
            shutil.move(fname, fname.replace('.img', '.dat'))
            filelist[i] = fname.replace('.img', '.dat')
    if any(x.endswith('.tif') for x in filelist):
        ext = 'tif'
    else:
        ext = 'dat'
    xml = glob.glob(os.path.join(outputdir, '*.xml'))
    if len(xml)>0:
        sceneid = os.path.basename(xml[0]).replace('.xml', '')
    elif basename[:1] == 'L' and len(basename) > 21:
        sceneid = basename[:21]
    else: 
        print('No XML file found, returning.')
        logerror(f, 'No XML file found.')
        return
    
    if overwrite:
        print('Deleting existing output files.')
        for d in [fmaskdir, srdir, btdir, ndvidir, evidir]:
            dellist = glob.glob(os.path.join(d,'{}*.*'.format(sceneid)))
            if len(dellist) > 0:
                for entry in dellist:
                    os.remove(entry)
    
    # Fmask file
    in_raster = os.path.join(outputdir, '{}_cfmask.{}'.format(sceneid, ext))
    if not os.access(in_raster, os.F_OK):
        print('Error, CFmask file is missing. Returning.')
        logerror(f, 'CFmask file missing.')
        return
        # if os.access(in_raster.replace('.tif', '.img'),os.F_OK):
        #     in_raster = in_raster.replace('.tif', '.img')
        # else:
        #     in_raster = in_raster.replace('_cfmask.tif', '_cfmask.dat')
    out_raster = os.path.join(fmaskdir, '{}_cfmask.dat'.format(sceneid))
    # if out_raster.endswith('.tif'):
    #     out_raster = out_raster.replace('.tif', '.dat')
    # elif out_raster.endswith('.img'):
    #     out_raster = out_raster.replace('.img', '.dat')
    if not os.path.exists(out_raster): # overwrite or 
        print('Reprojecting %s Fmask to ITM.'%sceneid)
        # if overwrite and os.path.exists(out_raster):
        #     flist =  glob.glob(out_raster.replace('.dat', '.*'))
        #     for fl in flist:
        #         os.remove(fl)
        reproject_ITM(in_raster, out_raster, sceneid = sceneid, rastertype = 'Fmask')
    
    # Surface reflectance data
    out_itm = os.path.join(srdir,'{}_ref_ITM.dat'.format(sceneid))
    if overwrite or not os.path.exists(out_itm):
        # if overwrite:
        #     flist =  glob.glob(out_itm.replace('.dat', '.*'))
        #     for fl in flist:
        #         os.remove(fl)
        print('Compositing surface reflectance bands to single file.')
        srlist = []
        out_raster = xml[0].replace('.xml', '.vrt')  
        if not os.path.exists(out_raster): # overwrite or 
            # if overwrite and os.path.exists(out_raster):
            #     os.remove(out_raster)  
            for band in bands:
                srlist.append(os.path.join(outputdir, '{}_sr_band{}.{}'.format(sceneid, band, ext))
                # if os.access(os.path.join(outputdir,'%s_sr_band%s.tif'%(sceneid, band)), os.F_OK):
                #     srlist.append(xml[0].replace('.xml', '_sr_band%s.tif'%band))
                # else:
                #     srlist.append(xml[0].replace('.xml', '_sr_band%s.dat'%band))
            mergelist = ['gdalbuildvrt', '-separate', out_raster]
            for s in srlist:
                mergelist.append(s)
            p = Popen(mergelist)
            print(p.communicate())
        print('Reprojecting %s reflectance data to Irish Transverse Mercator.'%sceneid)
        reproject_ITM(out_raster, out_itm, rastertype = 'ref', sceneid = sceneid)
    
    # Thermal data
    if basename[2:3] != '8':
        outbtdir = btdir
        rastertype = 'Landsat Band6'
        btimg = os.path.join(outputdir,'{}_toa_band6.dat'.format(sceneid))
    else:
        outbtdir = os.path.join(btdir, 'Landsat8')
        rastertype = 'Landsat TIR'
        btimg = os.path.join(tempdir,'{}_BT.vrt'.format(sceneid))
        print('Stacking Landsat 8 TIR bands for scene {}.'.format(sceneID))
        mergelist = ['gdalbuildvrt', '-separate', in_raster]
        for band in [10, 11]:
            mergelist.append(os.path.join(outputdir,'{}_toa_band{}.{}'.format(sceneid, band, ext)))
        p = Popen(mergelist)
        print(p.communicate())
    # if basename[2:3] != '8' and os.access(os.path.join(outputdir,'%s_toa_band6.dat'%(sceneid)),os.F_OK):
    #     
    # elif basename[2:3] != '8' and os.access(os.path.join(outputdir,'%s_toa_band6.tif'%(sceneid)),os.F_OK):
    #     btimg = os.path.join(outputdir,'%s_toa_band6.tif'%(sceneid))
    if btimg:
        BT_ITM = os.path.join(btdir,'%s_BT_ITM.dat'%sceneid)
        if not os.path.exists(BT_ITM): # overwrite or 
            # if overwrite and os.path.exists(BT_ITM):
            #     flist =  glob.glob(BT_ITM.replace('.dat', '.*'))
            #     for fl in flist:
            #         os.remove(fl)
        print('Reprojecting %s brightness temperature data to Irish Transverse Mercator.'%sceneid)
        reproject_ITM(btimg, BT_ITM, rastertype = rastertype, sceneid = sceneid)
        
    
    # Calculate EVI and NDVI
    if not os.path.exists(os.path.join(evidir, '%s_EVI.dat'%sceneid)): # overwrite or 
        try:
            calcvis(out_itm)
        except Exception as e:
            print('An error has occurred calculating VIs for scene %s:'%sceneid)
            print(e)
            logerror(out_itm, e)
    
    # Clean up files.
    
    if basename.endswith('.tar.gz'):
        print('Moving %s to archive: %s'%(basename, archdir))
        if not os.access(os.path.join(archdir,os.path.basename(f)),os.F_OK):
            shutil.move(f,archdir)
    if remove:
        print('Cleaning up files in directory.')
        filelist = glob.glob(os.path.join(outputdir, '{}*.*'.format(sceneid)))
        try:
            for fname in filelist:
                if os.access(fname,os.F_OK):
                    os.remove(fname)
            os.rmdir(outputdir)
        except Exception as e:
            print('An error has occurred cleaning up files for scene %s:'%sceneid)
            print(e)
            logerror(f, e)
    
    print('Processing complete for scene %s.'%sceneid)
        

def ESPAreprocess(SceneID, listfile):
    print('Adding scene {} for ESPA reprocessing to: {}'.format(SceneID, listfile))
    with open(listfile, 'a') as output:
        output.write('{}\n'.format(SceneID))


## File compression/ decompression

def unzip(f, outdir):
    import zipfile
    basename=os.path.basename(f)
    print('Unzipping '+basename+' to '+outdir) 
    zip_ref=zipfile.ZipFile(f, 'r')
    zip_ref.extractall(outdir)
    zip_ref.close()
    print('Files extracted to '+outdir)

def maketarfile(f, archdir):
    import tarfile
    basename=os.path.basename(f)
    dirname=os.path.dirname(f)
    j=basename.rfind('_')
    dataset=basename[:j]
    outname=archdir+'/'+dataset+'.tar.gz'
    datalist=glob.glob(dirname+'/'+dataset+'*.*')
    tarlist=[d for d in datalist if (('.dat' not in d) and ('.hdr' not in d))]
    tarfiles=len(tarlist)
    tarnum=1
    print('Writing '+str(tarfiles)+' files to '+outname)
    tar=tarfile.open(outname, 'w:gz')
    boname=os.path.basename(outname)
    for d in tarlist:
        bname=os.path.basename(d)
        print('Writing '+bname+' to '+boname+', file '+str(tarnum)+' of '+str(tarfiles)+'.')
        tar.add(d, recursive=False)
        tarnum+=1
    tar.close()
    print('Now deleting files from disk.')
    for d in tarlist:
        os.remove(d)
    print('Processing for '+boname+' complete.')  

def untarfile(file, outdir):
    import tarfile
    basename = os.path.basename(file)
    outbasepath = os.path.basename(outdir)
    if outbasepath in basename:
        outpath = outdir
    else:
        outpath = os.path.join(outdir,basename.rstrip('.tar.gz'))
    if not os.access(outpath,os.F_OK):
        os.mkdir(outpath)
    os.chdir(outpath)
    print('Extracting %s to: %s'%(basename, outpath))
    try:
        with tarfile.open(file) as tar:
            tar.extractall()
        filelist = glob.glob(os.path.join(outpath,'*.*'))
        # tar.close()
        if os.access(os.path.join(outpath,'gap_mask'),os.F_OK):
            shutil.rmtree(os.path.join(outpath,'gap_mask'))
        print('Completed extracting '+file)
        return filelist
    except Exception as e:
        logerror(file, e)
        print(e)
        tar.close()
        os.remove(file) # delete bad tar.gz
        return 0
    
## LEDAPS-specific functions


def fixespahdr(f, *args, **kwargs):
    xml = kwargs.get('xml', None)
    rastertype = kwargs.get('rastertype', None)
    if not f.endswith('.hdr'):
        j=f.rfind('.')+1
        f=f.replace(f[j:],'hdr')
    now = datetime.datetime.now()
    bak = '%s.%s.bak'%(f, now.strftime('%Y%m%d_%H%M%S'))
    basename = os.path.basename(f)
    sceneid = basename[:21]
    landsat=sceneid[2:3]
    if xml:
        tags=['acquisition_date','scene_center_time']
        dom = minidom.parse(xml)
        acquisitiontime='acquisition time = '
        for tag in tags:
            try:
                value = extract_xml(
                dom.getElementsByTagName(tag)[0].toxml())
                if tag == 'acquisition_date':
                    acquisitiontime += value+'T'
                else:
                    acquisitiontime += value+'\n'
            except Exception as e:
                print(e)
                return
            print(tag, value)
    else: 
        acquisitiontime = None
    shutil.move(f,bak)
    ENVIfile(f, rastertype, SceneID = sceneid, acqtime = acquisitiontime).WriteHeader()
    '''
    output=open(f,'w')
    units='wavelength units = Micrometers\n'
    if 'fmask' in f:
        imagetype='fmask'
        description='description = { Fmask Cloud and Shadow Classification (%s)}\n'%sceneid
        bnames='band names = { Cloud and Shadow Classification}\n'
        classnum='classes = 6\n'
        classlookup='class lookup = {0, 255, 0, 0, 0, 255, 127, 127, 127, 0, 255, 255, 255, 255, 255, 0, 0, 0}\n'
        classnames='class names = {Clear land, Clear water, Cloud shadow, Snow/ ice, Cloud, No data}\n'
    elif '_qa' in f:
        imagetype='QA'
        print('QA files not supported at this time, returning.')
        shutil.move(bak,f)
        return
    elif 'band6' in f:
        imagetype='BT'
        description='description = { LEDAPS Brightness Temperature (%s)}\n'%sceneid
        if (basename[0:3]!='LT8') or (basename[0:3]!='LC8'):
            bnames='band names = {TIR}\n'
            wavelength='wavelength = { 11.450000}\n'
            fwhm='fwhm = { 2.100000}\n'
        else:
            bnames='band names = {TIR 1, TIR 2}\n'
            wavelength='wavelength = { 10.895000, 12.005000}\n'
            fwhm='fwhm = { 0.590000, 1.010000}\n'
    elif '_ref' in f:
        imagetype='ref' 
        description='description = { LEDAPS Surface Reflectance (%s)}\n'%sceneid
        bnames='band names = {Blue, Green, Red, NIR, SWIR 1, SWIR 2}\n'
        if basename[0:3]=='LE7':
            wavelength='wavelength = { 0.483000, 0.560000, 0.662000, 0.835000, 1.648000, 2.206000}\n'
            fwhm='fwhm = { 0.070000, 0.080000, 0.060000, 0.120000, 0.200000, 0.260000}\n'
        elif (basename[0:3]=='LT4') or (basename[0:3]=='LT5'):
            wavelength='wavelength = { 0.485000, 0.560000, 0.662000, 0.830000, 1.648000, 2.215000}\n'
            fwhm='fwhm = { 0.070000, 0.080000, 0.060000, 0.130000, 0.200000, 0.270000}\n'
        elif landsat=='8':
            bnames='band names = {Coastal aerosol, Blue, Green, Red, NIR, SWIR 1, SWIR 2}\n'
            wavelength='wavelength = { 0.443000, 0.482600, 0.561300, 0.654600, 0.864600, 1.609000, 2.201000}\n'
            fwhm='fwhm = { 0.016000, 0.060100, 0.057400, 0.037500, 0.028200, 0.084700, 0.186700}\n'
    with open(bak,'r') as bakfile:
        projcs=True
        ignorelines=False
        lines=bakfile.readlines()
        if not any('coordinate system string' in p for p in lines):
            projcs=False 
        for line in lines:
            if not ignorelines:
                if 'description' in line:
                    output.write(description)
                elif 'map info' in line:
                    if not projcs:
                        gcsstring=gcs(line)
                        output.write(line+gcsstring)
                    else:
                        output.write(line)
                elif 'band names' in line:
                    if '}' not in line:
                        ignorelines=True
                    else:
                        output.write(bnames)
                else:
                    output.write(line)
            elif ignorelines and '}' in line:
                output.write(bnames)
                ignorelines=False
    if imagetype== 'BT' or imagetype=='ref':
        output.write(wavelength+fwhm+units)
    if imagetype=='fmask':
        output.write(classnum+classlookup+classnames)
    output.write(acquisitiontime)
    output.close()
    '''
    print('Header file updated for %s'%basename)    
        

def makehdrfile(filename, b1, taglist): # This modifies GDAL generated an ENVI header into a more useful with improved metadata.
    hdrname=filename.replace('.dat','.hdr')
    dirname=os.path.dirname(filename)
    temphdr=os.path.join(dirname,'temp.hdr')
    b1dict=b1.GetMetadata_Dict()
    sensor=b1dict['Satellite']+' '+b1dict['Instrument']+'\n'
    linenum=1
    outlines=['ENVI\n']
    print('Processing '+hdrname)
    inhdr=open(hdrname, 'r')
    lines=inhdr.readlines()
    basename=os.path.basename(hdrname)
    i=basename.rfind('_')+1
    j=basename.rfind('.')
    dataset=basename[:i-1]
    if 'fmask' in basename:
        imagetype='fmask'
    else:
        imagetype=basename[i:j]
    if imagetype=='ref':
        description='description = { LEDAPS Surface Reflectance (%s)}\n'%dataset
        bnames='band names = {Band 1 (Blue), Band 2 (Green), Band 3 (Red), Band 4 (NIR), Band 5 (SWIR1), Band 7 (SWIR2)}\n'
        if basename[0:3]=='LE7':
            wavelength='wavelength = { 0.483000, 0.560000, 0.662000, 0.835000, 1.648000, 2.206000}\n'
            fwhm='fwhm = { 0.070000, 0.080000, 0.060000, 0.120000, 0.200000, 0.260000}\n'
        if (basename[0:3]=='LT4') or (basename[0:3]=='LT5'):
            wavelength='wavelength = { 0.485000, 0.560000, 0.662000, 0.830000, 1.648000, 2.215000}\n'
            fwhm='fwhm = { 0.070000, 0.080000, 0.060000, 0.130000, 0.200000, 0.270000}\n'
    elif imagetype=='BT':
        description='description = { LEDAPS Brightness Temperature (%s)}\n'%dataset
        if (basename[0:3]!='LT8') or (basename[0:3]!='LC8'):
            bnames='band names = {Band 6 (TIR)}\n'
            wavelength='wavelength = { 11.450000}\n'
            fwhm='fwhm = { 2.100000}\n'
        else:
            bnames='band names = {Band 10 (TIR 1), Band 11 (TIR 2)}\n'
            wavelength='wavelength = { 10.895000, 12.005000}\n'
            fwhm='fwhm = { 0.590000, 1.010000}\n'
    elif imagetype=='QA':
        description='description = { Quality Assurance Data (%s)}\n'%dataset
        bnames='band names = {'
        for t in taglist:
            bnames+=t+','
        bnames=bnames[:-1]+'}\n'
    elif imagetype=='fmask':
        description='description = { Fmask Cloud and Shadow Classification (%s)}\n'%dataset
        bnames='band names = { Cloud and Shadow Classification}\n'
        classnum='classes = 6\n'
        classlookup='class lookup = {0, 255, 0, 0, 0, 255, 127, 127, 127, 0, 255, 255, 255, 255, 255, 0, 0, 0}\n'
        classnames='class names = {Clear land, Clear water, Cloud shadow, Snow/ ice, Cloud, No data}\n'
    outlines.append(description)
    aggregator=False
    for line in lines:
        linenum+=1
        if 'samples' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'lines' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'bands' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'header offset' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'header offset' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'file type' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'data type' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'interleave' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'byte order' in line:
            outlines.append(line.rstrip()+'\n')
        elif 'map info' in line:
            outlines.append(line.rstrip()+'\n')
            projcs=b1.GetProjection()
            if 'Unknown datum based upon the WGS 84 ellipsoid' in projcs:
                projcs=projcs.replace('Unknown datum based upon the WGS 84 ellipsoid','WGS 84')
            if 'Not specified (based on WGS 84 spheroid)' in projcs:
                projcs=projcs.replace('Not specified (based on WGS 84 spheroid)','WGS_1984')
            outlines.append('coordinate system string = '+projcs+'}\n')
    if imagetype == 'ref':
        outlines.append('default bands = {4,3,2}\n')
    if (imagetype != 'QA') and (imagetype != 'fmask'):
        outlines.append(bnames)
        outlines.append(wavelength)
        outlines.append(fwhm)
        outlines.append('wavelength units = Micrometers\n')
        if imagetype=='ref':
            outlines.append('data ignore value = -9999\n')
            outlines.append('reflectance scale factor = 10000.0\n')
        else:
            outlines.append('calibration scale factor = 100.0\n')
            outlines.append('data units = Celsius\n')
    if imagetype=='QA':#This aggregates multiple lines of band names into one line
        outlines.append(bnames)
    if imagetype=='fmask':
        outlines.append(bnames)
        outlines.append(classnum)
        outlines.append(classlookup)
        outlines.append(classnames)
    outlines.append('sensor type = Landsat '+b1dict['Instrument']+'\n')
    outlines.append('acquisition time = '+b1dict['AcquisitionDate']+'\n')
    output=open(temphdr, 'w')    
    for x in outlines:
        output.write(x)    
    output.close()
    inhdr.close()
    shutil.copy(hdrname, hdrname+'.bak')
    os.remove(hdrname)
    shutil.copy(temphdr,hdrname)
    os.remove(temphdr)
    print('ENVI header updated for '+hdrname+'.')

def runledaps(filename, outdir):
    filename=os.path.abspath(filename)
    infiledir=os.path.dirname(os.path.abspath(filename))
    basename=os.path.basename(filename)
    outdir=os.path.abspath(outdir)
    rmfilelist=[]
    if '_' in os.path.basename(filename):
        i=filename.rfind('_')
    else:
        i=filename.find('.tar.gz')
    h=filename.rfind('/')
    sceneid=filename[h:i]
    if infiledir!=outdir:
        print('Copying files from '+infiledir+' to '+outdir)
        if filename.endswith('.tar.gz'): 
            shutil.copy(filename, outdir)
            rmfilelist=untarfile(outdir+'/'+basename,outdir)
            os.remove(outdir+'/'+basename)
            filename=outdir+'/'+basename.replace('.tar.gz','_MTL.txt')
        else:
            filelist=glob.glob(filename[:i]+'*')
            for f in filelist:
                shutil.copy(f,outdir)
            filename=outdir+'/'+basename
        print('Copying complete.')
    if not filename.endswith('_MTL.txt'):
        if os.access(filename[:i]+'_MTL.txt',os.F_OK):
            filename=filename[:i]+'_MTL.txt'
        elif filename.endswith('.tar.gz'):
            rmfilelist=untarfile(filename,outdir)
            filename=filename.replace('.tar.gz','_MTL.txt')
        else:
            print('Error: '+basename+' is not in the correct format, skipping.')
            return
    elif filename.endswith('_MTL.txt') and not os.access(filename,os.F_OK):
        print('Error: '+basename+' does not exist, skipping.')
        return
    print('Running LEDAPS for '+basename)
    p=Popen(['do_ledaps.py','-f',filename])
    print(p.communicate())
    print('Processing complete for '+sceneid+'.')   
    