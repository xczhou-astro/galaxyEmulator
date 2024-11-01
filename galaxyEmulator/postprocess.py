import numpy as np
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, Moffat2DKernel, convolve_fft
import astropy.units as u
import astropy.constants as const
from scipy.interpolate import interp1d
import os
from skimage.transform import rescale
from joblib import Parallel, delayed
from scipy.integrate import trapezoid
from .utils import *
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
from matplotlib_scalebar.dimension import _Dimension
from mpl_toolkits.axes_grid1 import make_axes_locatable

class ParsecDimension(_Dimension):
    def __init__(self):
        super().__init__('pc')
        self.add_units('kpc', 1000)

class PostProcessing:

    def __init__(self, properties, config):

        self.config = config
        self.dataDir = self.config['dataDir']
        self.properties = properties
        self.workingDir = self.config['workingDir']
        self.subhaloID = self.properties['subhaloID']

    def __load_method(self, format):
        if format == '.npy':
            method = np.load
        else:
            method = np.loadtxt

        return method

    def __get_throughputs(self, survey):

        filters = split(self.config[f'filters_{survey}'])
        filterDir = os.path.join(self.dataDir, f'filters/{survey}')
        filterLs = os.listdir(filterDir)
        filterNames = [name.split('.')[0] for name in filterLs]

        throughputs = []
        for fil in filters:
            filename = filterLs[filterNames.index(fil)]
            throughput_file = os.path.join(filterDir, filename)
            wavelength_scale = get_wavelength_scale(throughput_file)
            throughput = np.loadtxt(throughput_file)
            throughput[:, 0] = throughput[:, 0] * wavelength_scale
            throughputs.append(throughput)

        return throughputs

    def __get_PSFs(self, survey):

        if self.config[f'PSFFromFile_{survey}']:
            filters = split(self.config[f'filters_{survey}'])
            psfDir = os.path.join(self.dataDir, f'PSFs/{survey}')
            psfLs = os.listdir(psfDir)

            psfNames = [name.split('.')[0] for name in psfLs]

            PSFs = []
            for fil in filters:
                filename = psfLs[psfNames.index(fil)]
                psf_file = os.path.join(psfDir, filename)
                
                psf = self.__load_method(os.path.splitext(psf_file)[1])(psf_file)
                if psf.shape[0] % 2 == 0:
                    psf = np.pad(psf, ((0, 1), (0, 1)), 'constant')
                psf = psf / np.sum(psf)
                PSFs.append(psf)
        
        else:
            numfilters = len(split(self.config[f'filters_{survey}']))
            fwhm_in_arcsec = split(self.config[f'PSFFWHM_{survey}'], float)
            fwhm_in_arcsec = extend(fwhm_in_arcsec, numfilters)
            
            pixelScales = self.properties[f'angleRes_{survey}']

            x_stds = [fwhm / ps for fwhm, ps in zip(fwhm_in_arcsec, pixelScales)]
            psfmodel = {'Moffat': Moffat2DKernel, 'Gaussian': Gaussian2DKernel}
            psfmodel = psfmodel[self.config[f'PSFModel_{survey}']]
            PSFs = []
            for xstd in x_stds:
                kernel = psfmodel(xstd)
                PSFs.append(kernel)
            
        return PSFs

    def __calculate_bandpass(self, img, tran, wave, factor):
        bandpass_img = factor * trapezoid(img * tran.reshape(-1, 1, 1) \
                                          * wave.reshape(-1, 1, 1), wave, axis=0)
        return bandpass_img
    
    def __bandpass_images(self, dataCube, survey, throughputs,
                    PSFs=None, bkgNoise=None):
        
        file = fits.open(dataCube)
        wave_sed = file[1].data['grid_points'] * 10**4
        numExp = self.properties[f'numExposure_{survey}']
        exposureTime = self.properties[f'exposureTime_{survey}']
        areaMirror = np.pi * (self.properties[f'aperture_{survey}'] / 2)**2

        ps_in_sr = self.properties[f'ps_in_sr_{survey}']
        interpolate_ratios = self.properties[f'ratios_{survey}']

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
            wave_in = wave_sed[idx]# * u.angstrom
            trans_in = interp(wave_in)
            image_arr = file[0].data[idx] / wave_in.reshape(-1, 1, 1)**2 # in flambda/u.sr
            converter = (u.MJy/u.sr) * const.c / u.angstrom**2 \
                        * numExp[i] * (exposureTime[i] * u.s) * (areaMirror * u.m**2) / (const.c * const.h) * u.angstrom**2
            converter = converter.to(u.sr ** -1)
            converter = converter.value
            image_arrs.append(image_arr)
            trans.append(trans_in)
            waves.append(wave_in)
            factors.append(converter)

        # bandpass_images = []
        # for img, tran, wave, factor in zip(image_arrs, trans, waves, factors):
        #     bandpass_images.append(self.__calculate_bandpass(img, tran, wave, factor))

        
        bandpass_images = Parallel(n_jobs=numfils, prefer='threads')(delayed(self.__calculate_bandpass)(img, tran, wave, factor)
                                                   for img, tran, wave, factor in zip(image_arrs, trans, waves, factors))
        

        resized_imgs = []
        for i, ratio in enumerate(interpolate_ratios):
            resized_imgs.append(rescale(bandpass_images[i], ratio))

        bandpass_images = resized_imgs

        converted_imgs = []
        for i, ps_sr in enumerate(ps_in_sr):
            converted_imgs.append(bandpass_images[i] * ps_sr)
        
        bandpass_images = converted_imgs

        # apply PSF effects
        if PSFs is not None:
            images_with_psf = []
            for i, img in enumerate(bandpass_images):
                # print(img.shape, PSFs[i].shape)
                images_with_psf.append(convolve_fft(img, PSFs[i]))
            
            bandpass_images = images_with_psf

        # add instrumental noise 
        if bkgNoise is not None:
            images_with_bkg = []
            for i, img in enumerate(bandpass_images):
                noise = np.random.normal(loc=0, scale=bkgNoise[i],
                                        size=img.shape)
                images_with_bkg.append(img + noise)
            
            bandpass_images = images_with_bkg
                
        bandpass_images = np.array(bandpass_images)
        
        return bandpass_images
    
    def __saveDataCube(self):
        numViews = np.int32(self.config['numViews'])
        
        datacubeDir = f'dataCubes/Subhalo_{self.subhaloID}'
        os.makedirs(datacubeDir, exist_ok=True)
        for i in range(numViews):
            shutil.move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_total.fits'),
                        os.path.join(datacubeDir, f'skirt_view_{i:02d}_total.fits'))
        
            shutil.move(os.path.join(self.workingDir, f'skirt_view_{i:02d}_sed.dat'),
                        os.path.join(datacubeDir, f'skirt_view_{i:02d}_sed.dat'))
            
    def __saveBandpassImages(self, images, survey):

        shape = images[0].shape

        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        numfilters = self.properties[f'numfilters_{survey}']
        
        
        for i in range(self.properties['numViews']):
            header = fits.Header()
            header['NAXIS'] = (3, 'number of data axes')
            header['NAXIS1'] = (shape[0], 'length of data axis 1')
            header['NAXIS2'] = (shape[1], 'length of data axis 2')
            header['NAXIS3'] = (shape[2], 'length of data axis 3')
            header['snapNum'] = (self.config['snapNum'], 'snapshot ID of IllustrisTNG')
            header['SURVEY'] = (survey, 'Survey')
            header['NFILTERS'] = (numfilters, 'Number of filters')
            header['UNIT'] = ('e', 'Unit of image array, in electron counts')
            header['INCLI'] = (self.properties['inclinations'][i], 'Inclination angle, in deg')
            header['AZIMUTH'] = (self.properties['azimuths'][i], 'Azimuth angle, in deg')
            header['REDSHIFT'] = (self.properties['redshift'], 'Redshift')
            header['FoV'] = (self.properties['FoV'], 'Field of view, in pc')
            header['lumiDis'] = (self.properties['lumiDis'], 'Luminosity distance, in Mpc')
            for count in range(numfilters):
                header[f'FILT_{count:02d}'] = (self.properties[f'filters_{survey}'][count],
                                                f'filter_{count:02d} for index {count}')
            for count in range(numfilters):
                header[f'RES_{count:02d}'] = (self.properties[f'resolution_{survey}'][count], 
                                                f'Pixel scale, in pc for index {count}')
            for count in range(numfilters):
                header[f'PS_{count:02d}'] = (self.properties[f'angleRes_{survey}'][count],
                                                f'Pixel scale, in arcsec for index {count}')
                
            hdu = fits.ImageHDU(data=images[i], header=header)
            hdulist.append(hdu)
            
        savedImageName = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_image.fits'
        hdulist.writeto(savedImageName, overwrite=True)

    def __saveSEDs(self, sedFilenames, survey):
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU()
        hdulist.append(primary_hdu)
        
        for i in range(self.properties['numViews']):
            sed = np.loadtxt(sedFilenames[i])
            shape = sed.shape
            
            header = fits.Header()
            header['NAXIS'] = (2, 'number of data axes')
            header['NAXIS1'] = (shape[0], 'length of data axis 1')
            header['NAXIS2'] = (shape[1], 'length of data axis 2')
            header['WUNIT'] = ('micron', 'Units of wavelength')
            header['FUNIT'] = ('Jy', 'Units of flux in F_nu')
            header['INCLI'] = (self.properties['inclinations'][i], 'Inclination angle, in deg')
            header['AZIMUTH'] = (self.properties['azimuths'][i], 'Azimuth angle, in deg')
            header['REDSHIFT'] = (self.properties['redshift'], 'Redshift')
            header['lumiDis'] = (self.properties['lumiDis'], 'luminosity distance, in Mpc')

            hdu = fits.ImageHDU(data=sed, header=header)
            hdulist.append(hdu)
            
        savedSEDName = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED.fits'
        hdulist.writeto(savedSEDName, overwrite=True)

    def __plot_image(self, image, res, savedFilename=None):

        if isinstance(image, Image.Image):
            pixels = np.array(image).shape[0]
        else:
            pixels = image.shape[0]
        
        z = str(np.around(self.properties['redshift'], 2))
        
        fig, ax = plt.subplots()
        ax.axis('off')
        im = ax.imshow(image)
        scalebarSize = 0.25 * pixels * res # in pc
        if scalebarSize > 1000:
            scalebarUnit = 'kpc'
            scalebarSize = np.around(scalebarSize / 1000, 1)
        else:
            scalebarUnit = 'pc'
            scalebarSize = np.around(scalebarSize, 1)
        
        pc_dim = ParsecDimension()
        
        scalebar = ScaleBar(res, 'pc', dimension=pc_dim, 
                            fixed_value=scalebarSize, fixed_units=scalebarUnit, frameon=False,
                            location='lower right', scale_loc='top',
                            color='white', font_properties={'size': 12})
        ax.add_artist(scalebar)
        ax.text(x=0.05, y=0.1, s=fr'$z$={z}', fontsize=12,
                transform=ax.transAxes, color='white')
        ax.text(x=0.05, y=0.05, s=f'ID:{self.subhaloID}', fontsize=12,
                transform=ax.transAxes, color='white')
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes('right', size='5%', pad=0.05)
        # fig.colorbar(im, cax=cax, label=unit)
        if savedFilename is not None:
            plt.savefig(savedFilename)
            plt.close()
        else:
            plt.show()

    def __plot_sed(self, sed, logscale=True, savedFilename=None):
        plt.figure()
        plt.plot(sed[:, 0] * 10**4, sed[:, 1], label='Total')
        # plt.plot(sed[:, 0] * 10**4, sed[:, 2], label='Transparent')
        if logscale:
            plt.xscale('log')
        # plt.legend(frameon=False)
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.ylabel(r'$F_{\nu}\ [Jy]$')
        if savedFilename is not None:
            plt.savefig(savedFilename)
            plt.close()
        else:
            plt.show()

    def runPostprocess(self):
        
        print('Run Postprocessing')

        if not self.config['postProcessing']:
            self.__saveDataCube()
        else:
            surveys = split(self.config['surveys'])
            
            for survey in surveys:
                filters = self.properties[f'filters_{survey}']
                saveBaseDir = f'mock_{survey}/Subhalo_{self.subhaloID}'
                os.makedirs(saveBaseDir, exist_ok=True)
                copyfile(os.path.join(self.workingDir, 'skirt_parameters.xml'),
                        os.path.join(saveBaseDir, 'skirt_parameters.xml'))
                copyfile(os.path.join(self.workingDir, 'quenched_stars.txt'),
                        os.path.join(saveBaseDir, 'quenched_stars.txt'))
                copyfile(os.path.join(self.workingDir, 'starforming_stars.txt'),
                        os.path.join(saveBaseDir, 'starforming_stars.txt'))
                copyfile(os.path.join(self.workingDir, 'dusts.txt'),
                        os.path.join(saveBaseDir, 'dusts.txt'))

                dataCubeFilenames = [os.path.join(self.workingDir, f'skirt_view_{i:02d}_total.fits')
                                     for i in range(self.properties['numViews'])]
                
                throughputs = self.__get_throughputs(survey)
                PSFs = None
                if self.config[f'includePSF_{survey}']:
                    PSFs = self.__get_PSFs(survey)

                bkgNoise = None
                if self.config[f'includeBkg_{survey}']:
                    bkgNoise = split(self.config[f'bkgNoise_{survey}'], float)

                images = []
                for i in range(self.properties['numViews']):
                    images.append(self.__bandpass_images(dataCubeFilenames[i], survey,
                                                       throughputs, PSFs, bkgNoise))

                self.__saveBandpassImages(images, survey)

                sedFilenames = [self.workingDir + f'/skirt_view_{i:02d}_sed.dat' 
                                for i in range(self.properties['numViews'])]
                self.__saveSEDs(sedFilenames, survey)

                if self.config[f'imgDisplay_{survey}']:

                    if self.config[f'RGBImg_{survey}']:
                        RGBFilters = split(self.config[f'RGBFilters_{survey}'])
                        RGBidx = [filters.index(RGBFilters[i]) for i in range(3)]
                        res = self.properties[f'resolution_{survey}'][RGBidx[0]]

                        for i in range(self.properties['numViews']):
                            RGBImg = convert_to_rgb(images[i], RGBidx)
                            savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                            self.__plot_image(RGBImg, res, savedFilename=savedFilename)
                    
                    else:
                        displayFilter = self.config[f'displayFilters_{survey}']
                        filteridx = filters.index(displayFilter)
                        res = self.properties[f'resolution_{survey}'][filteridx]
                    
                        for i in range(self.properties['numViews']):
                            img = images[i][filteridx]
                            savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_view_{i:02d}.png'
                            self.__plot_image(img, res, savedFilename=savedFilename)

                sed = np.loadtxt(sedFilenames[0])
                savedFilename = f'mock_{survey}/Subhalo_{self.subhaloID}/galaxy_SED.png'

                logscale = self.config['displaySEDxlogscale']
                # print(logscale)
                self.__plot_sed(sed, logscale=logscale,
                              savedFilename=savedFilename)
                
            if self.config['saveDataCube']:
                self.__saveDataCube()

        print('Postprocessing finished! Clearing workingDir!')
        shutil.rmtree(self.workingDir)