import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15
# import matplotlib.pyplot as plt
from PIL import Image
import sys
import os
# from matplotlib_scalebar.scalebar import ScaleBar
# from matplotlib_scalebar.dimension import _Dimension
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import shutil

def u2temp(u_energy, x_e):
    '''
    u_energy: InternelEnergy
    x_e: ElectronAbundance
    
    return:
    T: temperatrure in K
    '''
    
    X_H = 0.76
    u_energy = u_energy * (u.km/u.s)**2
    mu = 4 / (1 + 3 * X_H + 4 * X_H * x_e) * const.m_p
    T = ((5/3 - 1) * u_energy/const.k_B * mu).to(u.K)
    T = T.value
    return T

def fage(z):
    
    '''
    z: redshift
    
    return:
    time: age in Myr
    '''
    
    time = Planck15.age(z).to(u.Myr)
    time = time.value
    return time

class channel_RGB(object):
    def __init__(self, RED=None, GREEN=None, BLUE=None, backsub=False):
        self.red   = RED
        self.green = GREEN
        self.blue  = BLUE
        self.check_image_shapes()
        if backsub:
            self.subtract_background()

    def check_image_shapes(self):
        if (np.shape(self.red) != np.shape(self.green)) or \
            (np.shape(self.red) != np.shape(self.blue)):
            raise "Image arrays are of different shapes, exiting"
        else:
            self.NX, self.NY = np.shape(self.red)

    def subtract_background(self):
        self.red   -= np.median(self.red)
        self.green -= np.median(self.green)
        self.blue  -= np.median(self.blue)

    def apply_scale(self, scales=(1.0,1.0,1.0)):
        assert len(scales) == 3
        s1,s2,s3 = scales
        mean = (s1 + s2 + s3)/3.0
        self.red   *= (s1/mean)
        self.green *= (s2/mean)
        self.blue  *= (s3/mean)

    def pjm_offset(self, offset=0.0):
        if offset==None:
            pass
        else:
            self.red   += offset
            self.green += offset
            self.blue  += offset

    def pjm_mask(self,masklevel=None):
        if masklevel==None:
            pass
        else:
            tiny = 1e-12
            mask = self.red*0.0 + 1.0
            for image in (self.red, self.green, self.blue):
                image[np.isnan(image)] = 0.0
                image[np.isinf(image)] = 0.0
                mask[image < masklevel] = 0.0
                mask[(image > -tiny) & (image < tiny)] = 0.0
            self.red   *= mask
            self.green *= mask
            self.blue  *= mask

    def lupton_stretch(self, Q=1.0, alpha=1.0, itype='sum'):
        if itype == 'sum':
            I = (self.red+self.green+self.blue) + 1e-10
        elif itype == 'rms':
            I = np.sqrt(self.red**2.0+self.green**2.0+self.blue**2.0) + 1e-10
        stretch = np.arcsinh(alpha*Q*I) / (Q*I)
        self.red   *= stretch
        self.green *= stretch
        self.blue  *= stretch

    def lupton_saturate(self,threshold=1.0, saturation='white', unsat=0.995):
        if saturation=="white":
            pass
        elif saturation=="color":
            x = np.dstack((self.red, self.green,self.blue))
            maxpix = np.max(x, axis=-1)
            maxpix[maxpix<threshold] = 1.0
            self.red   /= maxpix
            self.green /= maxpix
            self.blue  /= maxpix
        else:
            print("Not a recognized type of saturation!!!")

        all_tmp = np.hstack([self.red.ravel(), self.green.ravel(), self.blue.ravel()])
        self.red    /= (all_tmp[all_tmp.argsort()[int(np.round(len(all_tmp)*unsat))]])
        self.green  /= (all_tmp[all_tmp.argsort()[int(np.round(len(all_tmp)*unsat))]])
        self.blue   /= (all_tmp[all_tmp.argsort()[int(np.round(len(all_tmp)*unsat))]])

    def pack_up(self, unsat=0.995):
        x = np.zeros([self.NX,self.NY,3])
        x[:,:,0] = np.flipud(self.red)
        x[:,:,1] = np.flipud(self.green)
        x[:,:,2] = np.flipud(self.blue)
        # x = x/(x.ravel()[x.ravel().argsort()[int(np.round(len(x.ravel())*unsat))]])
        x = np.clip(x,0.0,1.0)
        x = x*255
        self.imgRGB = Image.fromarray(x.astype(np.uint8))

def convert_to_rgb(bandpassImage, idx=[2, 3, 5]):
    
    '''
    bandpassImage: image in several bands
    idx: proper channels used to convert to RGB
    
    return:
    object_RGB.imgRGB: image in RGB
    '''
    img_g = bandpassImage[idx[0]]
    img_r = bandpassImage[idx[1]]
    img_z = bandpassImage[idx[2]]

    max_g = img_g.max()
    max_r = img_r.max()
    max_z = img_z.max()
    img_rscl = max([max_g, max_r, max_z])
    img_thds = 1e-6
    img_g_rscl = img_g/img_rscl
    img_r_rscl = img_r/img_rscl
    img_z_rscl = img_z/img_rscl
    img_g_rscl[np.where(img_g_rscl<=img_thds)] = img_thds
    img_r_rscl[np.where(img_r_rscl<=img_thds)] = img_thds
    img_z_rscl[np.where(img_z_rscl<=img_thds)] = img_thds
    
    scales, offset, Q, alpha, masklevel, saturation, itype = (0.5,1.05,1.5), 0.0, 3000, 0.3, -1, 'color', 'sum'
    
    #(0.5,1.05,1.5), 0.0, 200, 0.2, -1.0
    
    object_RGB = channel_RGB(RED=img_z_rscl, GREEN=img_r_rscl, BLUE=img_g_rscl)
    object_RGB.apply_scale(scales=scales)      
    object_RGB.lupton_stretch(Q=Q, alpha=alpha, itype=itype)
    object_RGB.pjm_mask(masklevel=masklevel)     
    object_RGB.pjm_offset(offset=offset)       
    object_RGB.lupton_saturate(saturation=saturation)
    object_RGB.pack_up()

    return object_RGB.imgRGB

# class ParsecDimension(_Dimension):
#     def __init__(self):
#         super().__init__('pc')
#         self.add_units('kpc', 1000)

# def plot_image(image, plot_infos, savedFilename=None):
#     '''
#     image: RGB image
#     savedFilename: directory to save image, if not None
#     '''
#     res = plot_infos['resol_in_pc']
#     pixels = plot_infos['pixels']
#     # unit = plot_infos['unit']
#     z = str(np.around(plot_infos['z'], 2))
#     ID = plot_infos['ID']
    
#     fig, ax = plt.subplots()
#     ax.axis('off')
#     im = ax.imshow(image)
#     scalebarSize = 0.25 * pixels * res # in pc
#     if scalebarSize > 1000:
#         scalebarUnit = 'kpc'
#         scalebarSize = np.around(scalebarSize / 1000, 1)
#     else:
#         scalebarUnit = 'pc'
#         scalebarSize = np.around(scalebarSize, 1)
    
#     pc_dim = ParsecDimension()
    
#     scalebar = ScaleBar(res, 'pc', dimension=pc_dim, 
#                         fixed_value=scalebarSize, fixed_units=scalebarUnit, frameon=False,
#                         location='lower right', scale_loc='top',
#                         color='white', font_properties={'size': 12})
#     ax.add_artist(scalebar)
#     ax.text(x=0.05, y=0.1, s=fr'$z$={z}', fontsize=12,
#             transform=ax.transAxes, color='white')
#     ax.text(x=0.05, y=0.05, s=f'ID:{ID}', fontsize=12,
#             transform=ax.transAxes, color='white')
#     # divider = make_axes_locatable(ax)
#     # cax = divider.append_axes('right', size='5%', pad=0.05)
#     # fig.colorbar(im, cax=cax, label=unit)
#     if savedFilename is not None:
#         plt.savefig(savedFilename)
#         plt.close()
#     else:
#         plt.show()


# def plot_sed(sed, logscale=True, savedFilename=None):
#     '''
#     sed: sed, only plot total (with dust) and transparent (without dust)
#     savedFilename: directory to save image, if not None
#     '''
#     plt.figure()
#     plt.plot(sed[:, 0] * 10**4, sed[:, 1], label='Total')
#     # plt.plot(sed[:, 0] * 10**4, sed[:, 2], label='Transparent')
#     if logscale:
#         plt.xscale('log')
#     # plt.legend(frameon=False)
#     plt.xlabel(r'Wavelength $[\AA]$')
#     plt.ylabel(r'$F_{\nu}\ [Jy]$')
#     if savedFilename is not None:
#         plt.savefig(savedFilename)
#         plt.close()
#     else:
#         plt.show()
        
        
def split(string, castType=None):
    splits = [inc for inc in "".join(string.split()).split(',')]
    if castType is not None:
        splits = [castType(sp) for sp in splits]
        
    return splits

def get_wavelength_scale(filename):
    with open(filename) as file:
        header = file.readline()
        
    if 'angstrom' in header or 'AA' in header:
        wavelength_scale = 1
    elif 'nm' in header:
        wavelength_scale = 10
    elif 'um' in header or 'micron' in header:
        wavelength_scale = 10 * 10**3
    else:
        print('unrecognized unit in wavelength!')
        sys.exit(0)
    return wavelength_scale
        

def calc_pivot(survey, filter):
    
    filterDir = f'../Data/filters/{survey}'
    filterLs = os.listdir(filterDir)
    filterNames = [name.split('.')[0] for name in filterLs]
    filename = filterLs[filterNames.index(filter)]
    filterName = os.path.join(filterDir, filename)
        
    wavelength_scale = get_wavelength_scale(filterName)
    
    try:
        transmission = np.loadtxt(filterName)
    except:
        transmission = np.load(filterName)
    
    transmission[:, 0] = transmission[:, 0] * wavelength_scale
    
    numerator = np.trapz(transmission[:, 1], transmission[:, 0])
    denomerator = np.trapz(transmission[:, 1] * transmission[:, 0]**-2,
                               transmission[:, 0])
    pivot = np.sqrt(numerator/denomerator)
    
    return pivot

def copyfile(src, tar):
    if os.path.exists(src) and not os.path.exists(tar):
        shutil.copyfile(src, tar)

def extend(values, nums):
    if len(values) == 1:
        values = [values][0] * nums
    else:
        values
    
    return values

# def calculate_spatialResolution(distance, properties, config):

#     if config['postProcessing']:
#         surveys = split(config['surveys'])
#         spatialResols = []
#         for survey in surveys:
#             filters = split(config[f'filters_{survey}'])
#             numfilters = len(filters)
#             # numfilters = len(split(config[f'filters_{survey}']))
#             if config[f'resolFromPix_{survey}']:
#                 pixelScales = split(config[f'pixelScales_{survey}'], float)
#                 pixelScales = extend(pixelScales, numfilters)
#                 resol = [distance * 10**6 * np.deg2rad(ang / 3600) for ang in pixelScales]
#                 pixelscales_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
#                 properties[f'resolution_{survey}'] = resol
#                 properties[f'angleRes_{survey}'] = pixelScales
#                 properties[f'ps_in_sr_{survey}'] = pixelscales_in_sr
                

#             else:
#                 resol = split(config[f'resolution_{survey}'], float)
#                 resol = extend(resol, numfilters)
#                 pixelScales = [np.rad2deg(res / (distance * 10**6)) * 3600 for res in resol]
#                 pixelscales_in_sr = [(ps**2 * u.arcsec**2).to(u.sr).value for ps in pixelScales]
#                 properties[f'resolution_{survey}'] = resol
#                 properties[f'angleRes_{survey}'] = pixelScales
#                 properties[f'ps_in_sr_{survey}'] = pixelscales_in_sr
            
#             spatialResols += resol
            
#             numExposure = split(config[f'numExposure_{survey}'], int)
#             numExposure = extend(numExposure, numfilters)
#             properties[f'numExposure_{survey}'] = numExposure
            
#             exposureTime = split(config[f'exposureTime_{survey}'], float)
#             exposureTime = extend(exposureTime, numfilters)
#             properties[f'exposureTime_{survey}'] = exposureTime
            
#             aperture = np.float32(config[f'aperture_{survey}'])
#             properties[f'aperture_{survey}'] = aperture
            
#             properties[f'numfilters_{survey}'] = numfilters
#             properties[f'filters_{survey}'] = filters
        
#         spatialResol = np.min(spatialResols)
        
#         for survey in surveys:
#             ratios = [spatialResol / res for res in properties[f'resolution_{survey}']]
#             properties[f'ratios_{survey}'] = ratios

        
#     else:
#         spatialResol = np.float32(config['spatialResol'])
    
#     return properties, spatialResol

# def get_properties(z, boxLength, config):
#     properties = {}
#     properties['redshift'] = z
#     distance = Planck15.luminosity_distance(z).value
#     properties['lumiDis'] = distance
    
#     properties, spatialResol = calculate_spatialResolution(distance, properties, config)
#     properties['baseRes'] = spatialResol
    
#     numViews = np.int32(config['numViews'])
#     properties['numViews'] = numViews
#     randomViews = config['randomViews']
#     if randomViews:
#         inclinations = np.random.uniform(0, 180, numViews).tolist()
#         azimuths = np.random.uniform(-360, 360, numViews).tolist()
#     else:
#         inclinations = split(config['inclinations'], float)
#         azimuths = split(config['azimuths'], float)
    
#     properties['inclinations'] = inclinations
#     properties['azimuths'] = azimuths
#     fovscale = np.float32(config['FoVboxLengthRatio'])
#     fov = fovscale * boxLength * 10**3
#     properties['FoV'] = fov
    
#     return properties