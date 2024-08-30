import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, AiryDisk2DKernel, convolve_fft
from scipy.interpolate import interp1d
import os
import shutil
import re
from .utils import *

def get_throughputs(survey, config):
    
    filters = split(config['filters'])
    survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
    # need to check if throughputs exist
    
    throughputs = []
    for fil in survey_filters:
        strings = fil.split('.')
        survey = strings[0]
        band = strings[1]
        
        throughput_file = f'../Data/filters/{survey}/{band}.fil'
        
        wavelength_scale = get_wavelength_scale(throughput_file)
        
        throughput = np.loadtxt(throughput_file)
        
        throughput[:, 0] = throughput[:, 0] * wavelength_scale
        
        throughputs.append(throughput)
        
    return throughputs

def get_PSFs(survey, config):
    
    if config['PSFFromFile']:
        filters = split(config['filters'])
        survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
        
        PSFs = []
        for fil in survey_filters:
            strings = fil.split('.')
            survey = strings[0]
            band = strings[1]
            
            psf_file = f'../Data/PSFs/{survey}/{band}.npy'
            psf = np.load(psf_file)
            PSFs.append(psf)
    else:
        fwhms = split(config['PSFFWHM'], float)
        
        PSFs = []
        for fwhm in fwhms:
            kernel = AiryDisk2DKernel(fwhm) # should use Gaussian or AiryDisk?
            PSFs.append(kernel)
            
    return PSFs

def bandpass_images(dataCube, throughputs, PSFs=None, backgroundSigma=None):
    file = fits.open(dataCube)
    wave_sed = file[1].data['grid_points'] * 10**4
    bandpass_images = []
    for i, thr in enumerate(throughputs):
        wave_min = np.min(thr[:, 0])
        wave_max = np.max(thr[:, 0])
        interp = interp1d(thr[:, 0], thr[:, 1])
        idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
        wave_in = wave_sed[idx]
        img_in = file[0].data[idx]
        trans_in = interp(wave_in)
        image = np.trapezoid(img_in * trans_in.reshape(-1, 1, 1), wave_in, axis=0)\
            /np.trapezoid(trans_in, wave_in)
        bandpass_images.append(image)

    if PSFs is not None:
        images_with_psf = []
        for i, img in enumerate(bandpass_images):
            images_with_psf.append(convolve_fft(img, PSFs[i]))
        
        bandpass_images = images_with_psf
    
    if backgroundSigma is not None:
        images_with_bkg = []
        for i, img in enumerate(bandpass_images):
            noise = np.random.normal(loc=0, scale=backgroundSigma[i],
                                     size=img.shape)
            images_with_bkg.append(img + noise)
        
        bandpass_images = images_with_bkg
            
    bandpass_images = np.array(bandpass_images)
    
    return bandpass_images


def postprocess(subhaloID, property, config):
    '''
    subhaloID: subhaloID used to generate the galaxy
    property: properties obtain from function modify_ski_file
    config: configuration in config.ini
    '''

    filters = split(config['filters'])
    
    surveys = [name.split('.')[0] for name in filters]
    
    surveys, numfilters = np.unique(surveys, return_counts=True)
    
    for survey, numfilter in zip(surveys, numfilters):
        saveBaseDir = f'mock_{survey}/Subhalo_{subhaloID}'
        os.makedirs(saveBaseDir, exist_ok=True)
        
        workingDir = config['workingDir']
        
        shutil.copy(os.path.join(workingDir, 'skirt_parameters.xml'),
                   os.path.join(saveBaseDir, 'skirt_parameters.xml'))
        
        shutil.copy(os.path.join(workingDir, 'quenched_stars.txt'),
                   os.path.join(saveBaseDir, 'quenched_stars.txt'))
        
        shutil.copy(os.path.join(workingDir, 'starforming_stars.txt'),
                   os.path.join(saveBaseDir, 'starforming_stars.txt'))
        
        shutil.copy(os.path.join(workingDir, 'dusts.txt'),
                   os.path.join(saveBaseDir, 'dusts.txt'))
        
        numViews = np.int32(config['numViews'])
        dataCubeFilenames = [workingDir + f'/skirt_view_{i:02d}_total.fits' for i in range(numViews)]
        
        throughputs = get_throughputs(survey, config)
        
        PSFs = None
        if config['includePSF']:
            PSFs = get_PSFs(survey, config)
        
        backgroundSigma = None
        if config['includeBackground']:
            backgroundSigma = split(config['backgroundSigma'], float)
        
        images = []
        for i in range(numViews):
            images.append(bandpass_images(dataCubeFilenames[i], throughputs, PSFs, backgroundSigma))
        
        shape = images[0].shape
        
        survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
        filterNames = [name.split('.')[1] for name in survey_filters]
        
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        for i in range(numViews):
            header = fits.Header()
            header['NAXIS'] = (3, 'number of data axes')
            header['NAXIS1'] = (shape[0], 'length of data axis 1')
            header['NAXIS2'] = (shape[1], 'length of data axis 2')
            header['NAXIS3'] = (shape[2], 'length of data axis 3')
            header['SURVEY'] = (survey, 'Survey')
            header['NFILTERS'] = (numfilter, 'Number of filters')
            header['UNIT'] = ('MJy/sr', 'Physical units of the array values')
            header['INCLI'] = (property['inclinations'][i], 'Inclination angle, in deg')
            header['AZIMUTH'] = (property['azimuths'][i], 'Azimuth angle, in deg')
            header['REDSHIFT'] = (property['redshift'], 'Redshift')
            header['FoV'] = (property['FoV'], 'Field of view, in pc')
            header['lumiDis'] = (property['lumiDis'], 'luminosity distance, in Mpc')
            for count in range(numfilter):
                header[f'FIL_{count:02d}'] = (f'{survey}.{filterNames[count]}', f'filter_{count:02d}')
                
            hdu = fits.ImageHDU(data=images[i], header=header)
            hdulist.append(hdu)
            
        savedImageName = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_image.fits'
        hdulist.writeto(savedImageName, overwrite=True)
        
        sedFilenames = [workingDir + f'/skirt_view_{i:02d}_sed.dat' for i in range(numViews)]
        
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        for i in range(numViews):
            sed = np.loadtxt(sedFilenames[i])
            shape = sed.shape
            
            header = fits.Header()
            header['NAXIS'] = (2, 'number of data axes')
            header['NAXIS1'] = (shape[0], 'length of data axis 1')
            header['NAXIS2'] = (shape[1], 'length of data axis 2')
            header['WUNIT'] = ('micron', 'Units of wavelength')
            header['FUNIT'] = ('Jy', 'Units of flux in F_nu')
            header['INCLI'] = (property['inclinations'][i], 'Inclination angle, in deg')
            header['AZIMUTH'] = (property['azimuths'][i], 'Azimuth angle, in deg')
            header['REDSHIFT'] = (property['redshift'], 'Redshift')
            header['lumiDis'] = (property['lumiDis'], 'luminosity distance, in Mpc')
    
            hdu = fits.ImageHDU(data=sed, header=header)
            hdulist.append(hdu)
            
        savedSEDName = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_SED.fits'
        hdulist.writeto(savedSEDName, overwrite=True)
        
        
                            
        if config[f'imgDisplay_{survey}']:
                            
            if config[f'RGBImg_{survey}']:
                
                RGBFilters = split(config[f'RGBFilters_{survey}'])
                
                RGBidx = [filterNames.index(RGBFilters[i]) for i in range(3)]
                
                for i in range(numViews):
                    RGBImg = convert_to_rgb(images[i], RGBidx)
                    savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_view_{i:02d}.png'
                    plot_image(RGBImg, savedFilename)
                    
            else:
                    
                displayFilter = config[f'displayFilter_{survey}']
                filteridx = filterNames.index(displayFilter)
                
                for i in range(numViews):
                    img = images[i][filteridx]
                    savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_view_{i:02d}.png'
                    plot_image(img, savedFilename)
                    
            
            sed = np.loadtxt(sedFilenames[0])
            savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_SED.png'
            
            plot_sed(sed, config[f'displaySEDxlogscale_{survey}'], savedFilename)
            
            
    if config['saveDatacube'] == True:
        
        datacubeDir = f'dataCubes/Subhalo_{subhaloID}/'
        os.makedirs(datacubeDir, exist_ok=True)
        
        for i in range(numViews):
            shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_total.fits'),
                        os.path.join(datacubeDir, f'skirt_view_{i:02d}_total.fits'))
            
            shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_sed.dat'),
                        os.path.join(datacubeDir, f'skirt_view_{i:02d}_sed.dat'))
            
    shutil.rmtree(workingDir)