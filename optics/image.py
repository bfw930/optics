
''' Image Translation

    image translation through path tracing engine


    functions:

        test_image = import_image(path, edge_len)

        patt = gen_image(edge_len, style = 'target')

        rays = gen_img_rays(edge_len, height, ss = 1.0)

        back_rays = gen_rev_rays(paths, opt_params)

        grid = translate_image(test_image, ss, paths, rev_paths, height, edge_len)

'''



''' Imports '''

# nd array manipulation
import numpy as np

# image import, processing, export
import PIL

# image manipulation
from scipy import ndimage




'''

D is dioptre at 1/m
optical power of lens, reciprocal of focal length


Myopia / Hyperopia

    long / short eyeball (cornea to retnia?), range -8.0 D to +8.0 D

    target ~ +- 3.0 D at 0.25 D step


Astigmatism

    cornea surface asymmetry, axis for rotation, range -5.0 D

    target ~ 3.0 D at 0.25 D step, axis over +-90 degrees at 1 degrees step


Presbyopia

    limit in lens focusing ability (close range), range +3.5 D

    target ~ 3.5 D in 0.25 D step

'''



''' test image preparation '''

def import_image(path, edge_len):

    ''' Import Image from File

    Args:
        path (str): path of image to import
        edge_len (int): image square length (pixels ~ mm)

    Returns:
        (np.array): 2d array of imported image data
    '''

    # import image data from file into np array
    test_image = np.array( PIL.Image.open(path) ).astype(np.float32)

## temporary pad 100 pxl img to 101
    test_image = np.pad(test_image, (1,0), 'edge')


    # scale image by longest edge to edge length
    test_image = ndimage.zoom(test_image, edge_len / np.max(test_image.shape), mode = 'nearest', order = 1)

    # pad minor edge to edge length
    test_image = np.pad(test_image, ((0, edge_len - test_image.shape[0]),
        (0, edge_len - test_image.shape[1])), 'edge')


    # return imported image
    return test_image



def gen_image(edge_len, style = 'target'):

    ''' Generate Target Image

    Args:
        edge_len (int): image square length (pixels ~ mm)
        style (str): image style keyword (not yet implimented)

    Returns:
        (np.array): 2d array of generated image data
    '''

    # generate target range
    div = [ (2.0 / (edge_len - 1)) for _ in range(2) ]
    rng = [ [-1., 1.] for _ in range(2) ]

    # prepare blank grid
    grid = [ np.arange(rng[i][0], rng[i][1] + div[i], div[i]) for i in range(2) ]
    zero = np.array([ [ 0. for xi in grid[0] ] for yi in grid[1] ])


    # set target origin
    orig = np.array([0., 0.])

    # calculate distance map
    dist = np.array([ [ np.sqrt(sum((np.array([xi, yi])-orig)**2)) for xi in grid[0] ] for yi in grid[1] ])


    # generate pattern from distance map
    patt = (np.cos(dist * 5 * 2 * np.pi) + 1.) / 2.0


    # trim edges
    patt[np.where(dist >= 0.9)] = 0.
    patt[np.where(dist <= 0.1)] = 0.


    # tophat floor/ciel filter for pattern
    patt[np.where(patt >= .5)] = 1.
    patt[np.where(patt < .5)] = 0.


    print(patt.shape[0]**2)

    # return generated pattern image
    return patt



''' image fransformation functions '''

def gen_img_rays(edge_len, height, ss = 1.0):

    ''' Generate Rays from Image

        Given image edge length and virtual height, generate rays over image plane for refraction

    Args:
        edge_len (int): image square length (pixels ~ mm)
        height (float): height of virtual image [mm]
        ss (float): supersampling factor for image subsample per pixel

    Returns:
        (list): list of initial rays (list np.array[px, py, pz, vx, vy, vz] )
    '''

    rays = []

    # set height of target image (mm)
    #height = 30.
    height = 5.

    # calculate ray initial positions
    rng = height / 2
    stp = 2 * rng / ((edge_len - 1) * ss)

    for Lz in np.arange(-rng, rng + stp, stp)[:]:
        for Ly in np.arange(-rng, rng + stp, stp)[:]:

            # origin ray position
            L0 = np.array([0., Ly, Lz])

            # normalised direction vector
            v = np.array([1., 0., 0.])
            v = v / np.linalg.norm(v)

            # store ray
            ray = np.concatenate([L0, v])
            rays.append(ray)


    # print number of rays and confirm edge length correct
    print( int(np.sqrt(len(rays))) )
    print( int(len(rays)) )


    # return list of generated rays
    return rays



def gen_rev_rays(paths, opt_params):

    ''' Generate Reverse Rays

        Generate rays for reverse path calculation, provided with calculated paths including inital and final ray
        positions and direction vectors

    Args:
        paths (list): initial ray paths through optics chain to retina
        opt_params (dict): optics chain parameters for given paths

    Returns:
        (list): list of reverse rays (list np.array[px, py, pz, vx, vy, vz] )
    '''


    # get index of only rays that reach retina
    j = [ i for i in range(len(paths)) if len(paths[i]) == 7 ]


    # calculate prime axis shift
    retina_dist = opt_params['eye_front'] + opt_params['cornea_thick'] + \
        opt_params['aqueous_thick'] + opt_params['lens_thick'] + \
        opt_params['retina_thick']# + 2*opt_params['retina_rad']

    #print(retina_dist)

    # reverse ray retinal position
    rev_rays_p = [ paths[i][-1][0] for i in range(len(paths)) if i in j ] # positions

    rev_rays_p = [ rev_rays_p[i] - [(retina_dist), 0, 0] for i in range(len(rev_rays_p)) ] # x axis to zero
    rev_rays_p = [ rev_rays_p[i] * np.array([-1, 1, 1]) for i in range(len(rev_rays_p)) ] # invert x axis


    # reverse ray direction vectors
    rev_rays_v = [ paths[i][-2][1] * np.array([1, -1, -1]) for i in range(len(paths)) if i in j ]


    # store reverse rays
    back_rays = [ np.concatenate([rev_rays_p[i], rev_rays_v[i]]) for i in range(len(rev_rays_p)) ]


    # return generated reverse rays
    return back_rays



def translate_image(test_image, ss, paths, rev_paths, height, edge_len):

    ''' Resample Translated Rays as Image

        Supersample image, get map of pixel to ray, average ray colour per pixel bin to generate image

    Args:
        paths (list): initial ray paths through optics chain to retina

    Returns:
        (np.array): 2d image array from translation through optics chain
    '''

    # get index of only rays that reach retina
    j = [ i for i in range(len(paths)) if len(paths[i]) == 7 ]

    # supersample image, get map of image pixel colour to each ray
    img_map = ndimage.zoom(test_image, (((test_image.shape[0] - 1) * ss) + 1) / test_image.shape[0]).flatten()[j]


    # get retinal image
    pts = np.stack([ path[-1][0] for path in rev_paths ])[:,1:]


    # calculate ray initial positions
    rng = height / 2
    # calculate step with subsample
    stp = 2 * rng / ((edge_len - 1))

    # define image ranges (pixel bins)
    xs = np.arange(-rng, rng + stp, stp)[:]
    ys = np.arange(-rng, rng + stp, stp)[:]


    # generate empty image map
    grid = np.zeros( (xs.shape[0],ys.shape[0]) )

    # iterate over pixel bins
    for i in range(len(xs)-1):
        for j in range(len(ys)-1):

            # get rays in pixel bin range
            k = np.where(
                (pts[:,0] >= xs[i]) &
                (pts[:,0] <= xs[i+1]) &
                (pts[:,1] >= ys[j]) &
                (pts[:,1] <= ys[j+1])
            )

            # average rays in pixel bin by image value
            if len(k[0]) > 0:
                #grid[j,i] = np.mean(img_map[k][:])
                grid[j,i] = np.median(img_map[k][:])


    # return generated image array
    return grid


