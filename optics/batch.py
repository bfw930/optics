
''' Batch Image Generation

    protocols for batch image generation protocols


    initialise test image
        define resolution, generate image

    generate initial rays
        define image virtual height, set image to ray sub-sampling

    define standard optics chain
        define standard optical component state

    calculate standard ray paths
        image rays path through optics chain to retina
        get initial rays for retina to image paths

    define stigmatism from standard
        select stigmatism, define adjust parameter ranges

    calculate return ray paths
        iterate each stigmatism parameter set
        calculate retina to virtual image paths

    resample rays to image
        iterate each virtual image
        resample return rays to generate image


    functions:

        test_image, rays, opt_params, opts, paths, back_rays = init_optics(edge_len, height, ss, opt_params)

        images = batch_translate(test_image, edge_len, height, ss, opt_params, paths, stig)

'''



''' Imports '''

# path tracing
from .engine import get_path

# optical model components
from .optics import std_opt_params, gen_optics, gen_optics_rev

# image translation
from .image import import_image, gen_image
from .image import gen_img_rays, gen_rev_rays, get_paths, translate_image


# nd array manipulation
import numpy as np

# plotting with matplotlib
import matplotlib.pyplot as plt


''' batch image generation protocols '''

def init_optics(edge_len, height, ss, opt_params):

    ''' Intialise Standard Optics

    Args:
        edge_len (int): image square length (pixels ~ mm)

    Returns:
        (np.array): 2d array of generated image data
    '''


    ''' Generate Target Image '''

    # generate pattern image (target)
    test_image = gen_image(edge_len)


    ''' Generate Initial Ray '''

    # generate rays for image translation (list np.array[px, py, pz, vx, vy, vz] )
    rays = gen_img_rays(edge_len, height, test_image, ss)


    ''' Define Standard Optics '''

    # overwrite standard optical parameters
    opt_params = { **std_opt_params(), **opt_params }


    # generate standard optics chain
    opts = gen_optics(opt_params)


    ''' Calculate Ray Paths through Optics '''

    # calculate ray paths through optics chain to retina
    paths = get_paths(rays, opts)


    ''' Generate Reverse Rays '''

    # generate reverse rays for back propagation through reverse optics chain
    back_rays = gen_rev_rays(paths, opt_params)


    # return initialised image, optics, and rays
    return test_image, rays, opt_params, opts, paths, back_rays



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


def batch_translate(test_image, edge_len, height, ss, opt_params, paths, stig, back_rays):

    ''' Batch Image Translation

    Args:
        edge_len (int): image square length (pixels ~ mm)

    Returns:
        (np.array): 2d array of generated image data
    '''

    # store image database
    images = []


    ''' Define Reverse Optics (Stigmatism) '''

    # generate range of stigmatism parameter values
    params = np.arange(stig['range'][0], stig['range'][1], stig['range'][2])

    print(len(params))

    # iterate over parameter delta range
    for i in range(len(params)):

        # get current stigmatism parameter value
        param = params[i]


        # define stigmatism optics chain by optical parameters
        rev_opt_params = {**opt_params,
            stig['param']: opt_params[stig['param']] + param,
        }

        # generate standard optics chain, overwrite existing params
        rev_opts = gen_optics_rev(rev_opt_params)


        ''' get ray paths through optics chain '''

        # calculate reverse ray paths from retina, set initial refractive index
        rev_paths = get_paths(back_rays, rev_opts, n0 = 1.337)


        ''' Resample Translated Rays as Image'''

        # build image by resample return rays over area
        grid = translate_image(test_image, ss, paths, rev_paths, height, edge_len)


        # collect image and parameters
        image = {
            'image': grid,
            'param': stig['param'],
            'delta': param
        }

        # store in database
        images.append(image)

        print('image {}/{}'.format(i+1, len(params)))


    # return generated pattern image
    return images



def store_images(images, out_path):

    ''' Store Image Batch

    Args:
        images (list): batch of generated images as list of dict (image, params)
        path (str): images output path

    Returns:
        (none): images written to file
    '''

    # iterate each image in batch
    for image in images:

        # get current image data
        _img = image['image']

        # initialise figure and axes, clean format
        _w = 4; _h = 4
        fig = plt.figure(figsize = (_w, _h))
        fig.canvas.layout.width = '{}in'.format(_w)
        fig.canvas.layout.height= '{}in'.format(_h)

        plt.subplot(111)
        plt.grid(''); plt.xticks([]); plt.yticks([])

        plt.imshow(_img, cmap = 'bone_r', vmin = 0., vmax = 1.)


        #plt.tight_layout()
        plt.subplots_adjust(left = .0, right = 1., top = 1., bottom = .0)

        plt.savefig('{}/output-{}-{:.3f}.png'.format(out_path, image['param'], image['delta']), dpi = 200)

        plt.close()


