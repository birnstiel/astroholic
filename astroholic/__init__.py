from importlib import metadata as _md

from ._lic import lic

from .lic import contrast_enhance
from .lic import calc_2D_streamline
from .lic import LIC_twostage
from .lic import hsv_mix
from .lic import pcolormesh_rgb

__version__ = _md.version('astroholic')

__all__ = [
    'lic',
    'contrast_enhance',
    'calc_2D_streamline',
    'LIC_twostage',
    'hsv_mix',
    'pcolormesh_rgb',
]
