
''' Batch Image Generation

    protocols for batch image generation protocols


    functions:

        images = batch_img_gen()

'''



''' Imports '''

# path tracing
from .engine import get_path


# nd array manipulation
import numpy as np



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



''' batch image generation protocols '''

def batch_image_gen():

    ''' Batch Image Generation

    Args:
        edge_len (int): image square length (pixels ~ mm)

    Returns:
        (np.array): 2d array of generated image data
    '''


    ''' Generate Target Image '''

    # set edge length; ensure odd
    edge_len = 151

    # generate pattern image (target)
    test_image = optics.gen_image(edge_len)


    ''' Generate Initial Ray '''

    # set height of target image (mm), and supersample factor
    height = 5.
    ss = 2.

    # generate rays for image translation (list np.array[px, py, pz, vx, vy, vz] )
    rays = optics.gen_img_rays(edge_len, height, test_image, ss)


    ''' Define Standard Optics '''

    # get standard optical parameters
    opt_params = optics.std_opt_params()

    # overwrite standard optical parameters
    opt_params = { **opt_params,

        #'eye_front': 300.,

        #'cornea_sph': 1.,
        #'cornea_axis': 0.,

        #'cornea_pow': np.sqrt(0.5),

        #'iris_dia': 4.,

        #'focus': 1.,

        #'lens_pow': np.sqrt(4.5),

        #'retina_thick': 17.2,
    }

    # generate standard optics chain
    opts = optics.gen_optics(opt_params)


    ''' Calculate Ray Paths through Optics '''

    # calculate ray paths through optics chain to retina
    paths = optics.get_paths(rays, opts)


    ''' Generate Reverse Rays '''

    # generate reverse rays for back propagation through reverse optics chain
    back_rays = optics.gen_rev_rays(paths, opt_params)


    ### begin iterate over parameter ranges

    ''' Define Reverse Optics (Stigmatism) '''

    # define stigmatism optics chain by optical parameters
    rev_opt_params = {**opt_params,
        'cornea_sph': opt_params['cornea_sph'] - 0.015,
    }

    # generate standard optics chain, overwrite existing params
    rev_opts = optics.gen_optics_rev(rev_opt_params)


    ''' get ray paths through optics chain '''

    # calculate reverse ray paths from retina, set initial refractive index
    rev_paths = optics.get_paths(back_rays, rev_opts, n0 = 1.337)



    ''' Resample Translated Rays as Image'''

    # build image by resample return rays over area
    grid = optics.translate_image(test_image, ss, paths, rev_paths, height, edge_len)




    # return generated pattern image
    return None



