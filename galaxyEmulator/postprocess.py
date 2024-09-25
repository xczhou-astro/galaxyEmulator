import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve_fft
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interp1d
import os
import re
from termcolor import colored
from skimage.transform import rescale
from joblib import Parallel, delayed
from scipy.integrate import trapezoid
from .utils import *

def get_throughputs(survey, config):
    
    filters = split(config['filters'])
    survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
    
    filterDir = f'../Data/filters/{survey}'
    filterLs = os.listdir(filterDir)
    filterNames = [name.split('.')[0] for name in filterLs]
    
    throughputs = []
    for fil in survey_filters:
        strings = fil.split('.')
        survey = strings[0]
        band = strings[1]
        
        filename = filterLs[filterNames.index(band)]
        throughput_file = os.path.join(filterDir, filename)
        
        wavelength_scale = get_wavelength_scale(throughput_file)
        
        throughput = np.loadtxt(throughput_file)
        
        throughput[:, 0] = throughput[:, 0] * wavelength_scale
        
        throughputs.append(throughput)
        
    return throughputs

def get_PSFs(survey, config):
    
    if config['PSFFromFile']:
        filters = split(config['filters'])
        survey_filters = [filter for filter in filters if bool(re.match(survey, filter, re.I))]
        
        psfDir = f'../Data/PSFs/{survey}'
        psfLs = os.listdir(psfDir)
        psfNames = [name.split('.')[0] for name in psfLs]
        
        PSFs = []
        for fil in survey_filters:
            strings = fil.split('.')
            survey = strings[0]
            band = strings[1]
            
            filename = psfLs[psfNames.index(band)]
            psf_file = os.path.join(psfDir, filename)
            
            try:
                psf = np.load(psf_file)
            except:
                psf = np.loadtxt(psf_file)

            psf = psf / np.sum(psf)
                
            PSFs.append(psf)
    else:
        fwhms_in_arcsec = split(config['PSFFWHM'], float)
        pixelScales = split(config['pixelScales'], float)
        x_stds = [fwhm/ps for fwhm, ps in zip(fwhms_in_arcsec, pixelScales)]
        
        PSFs = []
        for xstd in x_stds:
            kernel = Gaussian2DKernel(xstd)
            PSFs.append(kernel)
            
    return PSFs

def calculate_bandpass(img, tran, wave, factor):
    bandpass_img = factor * trapezoid(img * tran.reshape(-1, 1, 1) \
                                        * wave.reshape(-1, 1, 1), wave, axis=0)
    bandpass_img = bandpass_img.to(u.sr ** -1)
    bandpass_img = np.array(bandpass_img)
    return bandpass_img

def bandpass_images(dataCube, throughputs, infos,
                    PSFs=None, backgroundSigma=None):
    file = fits.open(dataCube)
    wave_sed = file[1].data['grid_points'] * 10**4 # to angstrom
    # image array in MJr/sr, need to convert to electron counts, since CCDs are used.

    numExp = infos['numExposure']
    exposureTime = infos['exposureTime']
    areaMirror = infos['areaMirror']

    ps_in_sr = infos['ps_in_sr']
    interpolate_ratios = infos['ratios']

    numfils = len(throughputs)

    image_arrs = []
    trans = []
    waves = []
    factors = []
    for i, thr in enumerate(throughputs):
        wave_min = np.min(thr[:, 0])
        wave_max = np.max(thr[:, 0])
        interp = interp1d(thr[:, 0], thr[:, 1])
        idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
        wave_in = wave_sed[idx] * u.angstrom
        trans_in = interp(wave_in)
        image_arr = file[0].data[idx] * (u.MJy/u.sr) * const.c \
                    / wave_in.reshape(-1, 1, 1)**2 # in flambda/u.sr
        factor = numExp[i] * (exposureTime * u.s) * (areaMirror * u.m**2) / (const.c * const.h)
        image_arrs.append(image_arr)
        trans.append(trans_in)
        waves.append(wave_in)
        factors.append(factor)

    bandpass_images = Parallel(n_jobs=numfils)(delayed(calculate_bandpass)(img, tran, wave, factor)
                                               for img, tran, wave, factor in zip(image_arrs, trans, waves, factors))
        
    
    # obtain bandpass images, in e^- * sr^-1
    # bandpass_images = []
    # for i, thr in enumerate(throughputs):
    #     wave_min = np.min(thr[:, 0])
    #     wave_max = np.max(thr[:, 0])
    #     interp = interp1d(thr[:, 0], thr[:, 1])
    #     idx = np.where((wave_sed > wave_min) & (wave_sed < wave_max))[0]
    #     wave_in = wave_sed[idx] * u.angstrom
    #     img_in = file[0].data[idx] * (u.MJy/u.sr) * const.c / wave_in.reshape(-1, 1, 1)**2
    #     # img_in = img_in.to(u.erg * u.cm**-2 * u.s**-1 * u.angstrom**-1 / u.sr) 
    #     trans_in = interp(wave_in) # unit 1
    #     factor = numExp[i] * (exposureTime * u.s) * (areaMirror * u.m**2) / (const.c * const.h)
    #     image_in_e_sr = factor * np.trapz(img_in * trans_in.reshape(-1, 1, 1) \
    #                                       * wave_in.reshape(-1, 1, 1), wave_in, axis=0)
    #     image_in_e_sr = image_in_e_sr.to(u.sr ** -1) # in e^- / sr
    #     image_in_e_sr = np.array(image_in_e_sr)
    #     bandpass_images.append(image_in_e_sr)

    # resize images by ratios derived from pixelscales
    resized_imgs = []
    for i, ratio in enumerate(interpolate_ratios):
        resized_imgs.append(rescale(bandpass_images[i], ratio))

    bandpass_images = resized_imgs

    # convert unit of bandpass images to e^-
    converted_imgs = []
    for i, ps_sr in enumerate(ps_in_sr):
        converted_imgs.append(bandpass_images[i] * ps_sr)
    
    bandpass_images = converted_imgs

    # apply PSF effects
    if PSFs is not None:
        images_with_psf = []
        for i, img in enumerate(bandpass_images):
            images_with_psf.append(convolve_fft(img, PSFs[i]))
        
        bandpass_images = images_with_psf

    # add instrumental noise 
    if backgroundSigma is not None:
        images_with_bkg = []
        for i, img in enumerate(bandpass_images):
            noise = np.random.normal(loc=0, scale=backgroundSigma[i],
                                     size=img.shape)
            images_with_bkg.append(img + noise)
        
        bandpass_images = images_with_bkg
            
    bandpass_images = np.array(bandpass_images)
    
    return bandpass_images

def postprocess(subhaloID, properties, config):
    '''
    subhaloID: subhaloID used to generate the galaxy
    properties: properties obtain from function modify_ski_file
    config: configuration in config.ini
    '''
    

    if config['postProcessing']:

        filters = split(config['filters'])
        
        surveys = [name.split('.')[0] for name in filters]
        
        surveys_unique, numfilters = np.unique(surveys, return_counts=True)

        resolutions = properties['resolution']
        baseRes = properties['baseRes']

        ratios = [baseRes / res for res in resolutions] # down sample ratios from the finest resolution
        pixelscales_in_sr = properties['pixelscales_in_sr']
        numExp = split(config['numExposure'], int)

        print('Considered survey(s)', surveys_unique)
        
        for survey, numfilter in zip(surveys_unique, numfilters):

            idx_survey = [i for i, sur in enumerate(surveys) if sur == survey]
            interpolate_ratios = [ratios[i] for i in idx_survey]
            ps_in_sr = [pixelscales_in_sr[i] for i in idx_survey]
            numExposure = [numExp[i] for i in idx_survey]
            resol_in_pc = [resolutions[i] for i in idx_survey]

            exposureTime = np.float32(config[f'exposureTime_{survey}'])
            areaMirror = np.pi * np.float32(config[f'mirrorSize_{survey}'])**2
            

            infos = {}
            infos['ratios'] = interpolate_ratios
            infos['ps_in_sr'] = ps_in_sr
            infos['numExposure'] = numExposure
            infos['exposureTime'] = exposureTime
            infos['areaMirror'] = areaMirror
            
            saveBaseDir = f'mock_{survey}/Subhalo_{subhaloID}'
            os.makedirs(saveBaseDir, exist_ok=True)
            
            workingDir = config['workingDir']

            copyfile(os.path.join(workingDir, 'skirt_parameters.xml'),
                    os.path.join(saveBaseDir, 'skirt_parameters.xml'))
            copyfile(os.path.join(workingDir, 'quenched_stars.txt'),
                    os.path.join(saveBaseDir, 'quenched_stars.txt'))
            copyfile(os.path.join(workingDir, 'starforming_stars.txt'),
                    os.path.join(saveBaseDir, 'starforming_stars.txt'))
            copyfile(os.path.join(workingDir, 'dusts.txt'),
                    os.path.join(saveBaseDir, 'dusts.txt'))
            
            # shutil.copy(os.path.join(workingDir, 'skirt_parameters.xml'),
            #         os.path.join(saveBaseDir, 'skirt_parameters.xml'))
            
            # shutil.copy(os.path.join(workingDir, 'quenched_stars.txt'),
            #         os.path.join(saveBaseDir, 'quenched_stars.txt'))
            
            # shutil.copy(os.path.join(workingDir, 'starforming_stars.txt'),
            #         os.path.join(saveBaseDir, 'starforming_stars.txt'))
            
            # shutil.copy(os.path.join(workingDir, 'dusts.txt'),
            #         os.path.join(saveBaseDir, 'dusts.txt'))
            
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
                images.append(bandpass_images(dataCubeFilenames[i], throughputs, 
                                            infos, PSFs, backgroundSigma))

            plot_infos = {}
            
            shape = images[0].shape
            plot_infos['pixels'] = shape[1]

            # plot_infos['unit'] = 'electrons'

            plot_infos['z'] = np.around(properties['redshift'], 2)
            plot_infos['ID'] = subhaloID
            
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
                header['snapNum'] = (config['snapNum'], 'snapshot ID of IllustrisTNG')
                header['SURVEY'] = (survey, 'Survey')
                header['NFILTERS'] = (numfilter, 'Number of filters')
                header['UNIT'] = ('e', 'Unit of image array, in electron counts')
                header['INCLI'] = (properties['inclinations'][i], 'Inclination angle, in deg')
                header['AZIMUTH'] = (properties['azimuths'][i], 'Azimuth angle, in deg')
                header['REDSHIFT'] = (properties['redshift'], 'Redshift')
                header['FoV'] = (properties['FoV'], 'Field of view, in pc')
                header['lumiDis'] = (properties['lumiDis'], 'Luminosity distance, in Mpc')
                for count in range(numfilter):
                    header[f'FILT_{count:02d}'] = (f'{survey}.{filterNames[count]}',
                                                   f'filter_{count:02d} for index {count}')
                for count in range(numfilter):
                    header[f'RES_{count:02d}'] = (properties['resolution'][count], 
                                                  f'Pixel scale, in pc for index {count}')
                for count in range(numfilter):
                    header[f'PS_{count:02d}'] = (properties['angleRes'][count],
                                                 f'Pixel scale, in arcsec for index {count}')
                    
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
                header['INCLI'] = (properties['inclinations'][i], 'Inclination angle, in deg')
                header['AZIMUTH'] = (properties['azimuths'][i], 'Azimuth angle, in deg')
                header['REDSHIFT'] = (properties['redshift'], 'Redshift')
                header['lumiDis'] = (properties['lumiDis'], 'luminosity distance, in Mpc')
        
                hdu = fits.ImageHDU(data=sed, header=header)
                hdulist.append(hdu)
                
            savedSEDName = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_SED.fits'
            hdulist.writeto(savedSEDName, overwrite=True)
            
                                
            if config[f'imgDisplay_{survey}']:
                                
                if config[f'RGBImg_{survey}']:
                    
                    RGBFilters = split(config[f'RGBFilters_{survey}'])
                    
                    RGBidx = [filterNames.index(RGBFilters[i]) for i in range(3)]
                    plot_infos['resol_in_pc'] = resol_in_pc[RGBidx[0]]
                    
                    for i in range(numViews):
                        RGBImg = convert_to_rgb(images[i], RGBidx)
                        savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_view_{i:02d}.png'
                        plot_image(RGBImg, plot_infos , savedFilename)
                        
                else:
                        
                    displayFilter = config[f'displayFilter_{survey}']
                    filteridx = filterNames.index(displayFilter)
                    
                    for i in range(numViews):
                        img = images[i][filteridx]
                        savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_view_{i:02d}.png'
                        plot_image(img, plot_infos, savedFilename)
                
                sed = np.loadtxt(sedFilenames[0])
                savedFilename = f'mock_{survey}/Subhalo_{subhaloID}/galaxy_SED.png'
                
                plot_sed(sed, config[f'displaySEDxlogscale_{survey}'], savedFilename)
                
                
        if config['saveDatacube'] == True:

            saveDatacube(subhaloID, config)
            
            # datacubeDir = f'dataCubes/Subhalo_{subhaloID}/'
            # os.makedirs(datacubeDir, exist_ok=True)
            
            # for i in range(numViews):
            #     shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_total.fits'),
            #                 os.path.join(datacubeDir, f'skirt_view_{i:02d}_total.fits'))
                
            #     shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_sed.dat'),
            #                 os.path.join(datacubeDir, f'skirt_view_{i:02d}_sed.dat'))
    else:
        saveDatacube(subhaloID, config)
            
    shutil.rmtree(workingDir)

def saveDatacube(subhaloID, config):
    numViews = np.int32(config['numViews'])
    workingDir = config['workingDir']

    datacubeDir = f'dataCubes/Subhalo_{subhaloID}/'
    os.makedirs(datacubeDir, exist_ok=True)
    for i in range(numViews):
        shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_total.fits'),
                    os.path.join(datacubeDir, f'skirt_view_{i:02d}_total.fits'))
        
        shutil.move(os.path.join(workingDir, f'skirt_view_{i:02d}_sed.dat'),
                    os.path.join(datacubeDir, f'skirt_view_{i:02d}_sed.dat'))
        
    