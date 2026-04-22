import numpy as _np
from scipy.ndimage import laplace as _laplace
import matplotlib.pyplot as _plt
from matplotlib.colors import Normalize as _Normalize, rgb_to_hsv as _rgb_to_hsv, hsv_to_rgb as _hsv_to_rgb

from ._lic import lic as _fortran
LIC = _fortran.flic
gen_noise_fast = _fortran.gen_noise_fast

__all__ = ['LIC', 'gen_noise_fast', 'contrast_enhance', 'calc_2D_streamline',
           'LIC_twostage', 'mix', 'pcolormesh_rgb', 'check_openmp', 'normalize_velocity']


def check_openmp():
    """Check if OpenMP is enabled in the Fortran code and how many threads are available."""
    enabled, nthreads = _fortran.check_openmp()
    if enabled:
        print(f"OpenMP is enabled with {nthreads} threads.")
    else:
        print("OpenMP is not enabled.")


def contrast_enhance(data, sig=2.0):
    minval, maxval = data.mean() + sig * _np.array([-1, 1]) * data.std()
    minval = max(0.0, minval)
    maxval = min(1.0, maxval)
    output = _np.clip((data - minval) / (maxval - minval), 0.0, 1.0)
    return output


def calc_2D_streamline(p, x, y, vel, data, length=1.0, n_steps=10, direction='forward'):
    """calculate streamline in 2D velocity field

    Parameters
    ----------
    p : array
        initial position as array or list of length 2, (x, y)
    x : array
        regular x grid, shape (nx)
    y : array
        regular y grid, shape (ny)
    vel : array
        velocity in x and y, shape (nx, ny, 2), where one point is (vx, vy)
    data : array
        scalar field of shape (nx, ny), values along the path are returned for LIC
    length : float, optional
        how long the path of the stream line should be, by default 1.0
    n_steps : int, optional
        number of steps to take along the path, by default 10
    direction : str, optional
        direction, by default 'forward', can be
        - 'forward': move in the direction of the velocity
        - 'backward': move in the opposite direction of the velocity
        - 'both': move half a length in both direction


    Returns
    -------
    _type_
        _description_

    Raises
    ------
    ValueError
        _description_
    """
    x0, y0 = p
    if direction == 'forward':
        path, ipath, values = _fortran.calc_2d_streamline_forward(
            x0, y0, x, y, vel, data, length=length, n_steps=n_steps)
    elif direction == 'backward':
        path, ipath, values = _fortran.calc_2d_streamline_forward(
            x0, y0, x, y, vel, data, length=-length, n_steps=n_steps)
    elif direction == 'both':
        path, ipath, values = _fortran.calc_2d_streamline_bothways(
            x0, y0, x, y, vel, data, length=length, n_steps=n_steps)
    else:
        raise ValueError("direction needs to be 'forward' or 'both'")
    return path, ipath, values


def LIC_twostage(x, y, vel, generate_plot=False, **kwargs):
    """computes a 2-stage LIC with some contrast enhancement

    Parameters
    ----------
    x : array
        regular x-grid of shape (nx)
    y : array
        regular y-grid of shape (ny)
    vel : array
        2D velocity (vx, vy) of shape (nx, ny, 2)
    generate_plot : bool, optional
        if true, plot the different stages, by default False

    kwargs : are passed to LIC, e.g. length of the streamline

    Returns
    -------
    array
        2D LIC pattern
    """

    nx = len(x)
    ny = len(y)
    length = kwargs.pop('length', min(x[-1] - x[0], y[-1] - y[0]) / 5.)
    length = abs(length / 2)

    noise = gen_noise_fast(nx, ny)
    noise_L = _Normalize()(LIC(noise, x, y, vel, length=length, **kwargs))
    noise_Ll = _laplace(noise_L)
    noise_LlL = _Normalize()(LIC(noise_Ll, x, y, vel, length=length, **kwargs))
    noise_LlLC = contrast_enhance(noise_LlL)

    if generate_plot:
        imgs = [
            'noise',
            'noise_L',
            'noise_Ll',
            'noise_LlL',
            'noise_LlLC',
        ]
        n = len(imgs)
        f, ax = _plt.subplots(n, 2, figsize=(8, 2 * n), dpi=100)

        _loc = locals()

        for i, name in enumerate(imgs):
            cc = ax[i, 0].imshow(_loc[name].T, cmap='gray', norm=_Normalize())
            _plt.colorbar(cc, ax=ax[i, 0])
            ax[i, 0].set_title(name)
            ax[i, 0].set_axis_off()

            ax[i, 1].hist(_loc[name].ravel(), bins=20, fc='k')

    return noise_LlLC


def mix(scalar, noise, cmap='magma', mode='hsv', norm=None, alpha=None):
    """Mixes a color-mapped scalar and a LIC pattern in HSV color space

    Parameters
    ----------
    scalar : array
        2d array of a scalar quantity used for defining the color
    noise : array
        2d array of a LIC pattern, same shape as scalar
    cmap : str, optional
        color map to be used, by default 'magma'
    norm : norm, optional
        norm used to normalize the scalar before color mapping, by default a linear norm is used

    Returns
    -------
    array
        RGB values for the shape of scalar
    """
    if norm is None:
        norm = _Normalize()

    assert scalar.shape == noise.shape, 'scalar and noise need to have the same shape'

    # get the RGB colors for the two images
    img_col = _plt.colormaps[cmap](norm(scalar))[..., :3]
    img_lic = _plt.colormaps['gray'](_Normalize()(noise))[..., :3]

    if mode == 'hsv':
        # convert to HSV (float64, no quantization)
        hsv_col = _rgb_to_hsv(img_col)
        hsv_lic = _rgb_to_hsv(img_lic)

        # blend value channel
        hsv = hsv_col.copy()
        if alpha is None:
            hsv[..., 2] = (hsv_col[..., 2] + hsv_lic[..., 2]) / 2.0
        else:
            V = hsv_lic[..., 2]
            fac = 1.0 + alpha * (V - 0.5)
            fac = _np.clip(fac, 0.85, 1.15)   # optional safeguard
            hsv[..., 2] = _np.clip(hsv_col[..., 2] * fac, 0.0, 1.0)

        # convert back to RGB as uint8
        result = (_hsv_to_rgb(hsv) * 255).astype(_np.uint8)

    elif mode == 'rgb':

        if alpha is not None:
            beta = alpha
        else:
            beta = 0.3

        lic = (noise - noise.min()) / (noise.max() - noise.min())
        lic_tex = 1.0 + beta * (lic - 0.5)

        result = img_col * lic_tex[..., None]
        result = _np.clip(result, 0, 1)

    return result


def pcolormesh_rgb(x, y, rgb, ax=None, **kwargs):
    """Makes a pcolormesh plot given RGB data

    Parameters
    ----------
    x : array
        x-array
    y : array
        y-array
    rgb : array
        color data, shape (len(x), len(y), 3 or 4)
    ax : axes, optional
        axes into which to plot, by default None

    kwargs : are passed to pcolormesh

    Returns
    -------
    f, ax
        figure and axes object
    """
    if ax is None:
        f, ax = _plt.subplots()
    else:
        f = ax.figure
    col_len = rgb.shape[-1]
    cc = ax.pcolormesh(x, y, rgb[:, :, 0].T, facecolors=rgb.transpose(
        1, 0, 2).reshape(-1, col_len) / 255, **kwargs)
    cc.set_array(None)
    return f, ax


def normalize_velocity(VX, VY, p=0.5, eps=1e-8, perc=99):
    """Normalize the velocity field for better visualization in LIC

    Parameters
    ----------
    VX : array
        velocity in x direction
    VY : array
        velocity in y direction
    p : float, optional
        power for normalization, by default 0.5
    eps : float, optional
        small value to avoid division by zero, by default 1e-8
    perc : float, optional
        percentile for clipping the velocity, by default 99

    Returns
    -------
    VXn, VYn : array
        normalized velocity in x and y direction
    """

    speed = (VX**2 + VY**2)
    sclip = _np.percentile(speed, perc)
    speed = _np.clip(speed, 0, sclip)

    f_norm = 1 / (speed**p + eps)

    VXn = VX * f_norm
    VYn = VY * f_norm

    return VXn, VYn
