import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.cosmology import Planck15
from PIL import Image
import sys
import os
import shutil
import logging

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

