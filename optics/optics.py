
''' Optical Model Components

    optical model definition and generation functions


### helper function to generate full optics chain (list of dicts with parameters) from minimal input parameter set
## have relative distance of optics

## include parameter variation aligned to optical aberations (prescription), eg retinal distance, cornea ashpericity


    functions:

        opts = gen_optics(params)

        opts = gen_optics_rev(params)

        opt_params = std_opt_params()

'''



''' Imports '''

# nd array manipulation
import numpy as np



''' Optics Chain Generation Functions '''

def std_opt_params():

    ''' Standard Optics Parameters

    Args:

    Returns:
        (dict): standard optical parameters
    '''

    # define standard eye parameters
    opt_params = {
        'eye_front': 300.,

        'cornea_sph': 1.,
        'cornea_axis': 0.,

        'cornea_pow': np.sqrt(0.5),
        'cornea_f_rad': 7.8,
        'cornea_r_rad': 6.4,
        'cornea_thick': 0.6,
        'aqueous_thick': 3.0,

        'iris_dia': 4.,

        'lens_thick': 4.,
        'lens_f_rad_max': 10.1,
        'lens_f_rad_min': 5.95,
        'lens_r_rad_max': 6.1,
        'lens_r_rad_min': 4.5,

        'focus': 1.,

        'lens_pow': np.sqrt(4.5),

        'retina_thick': 17.2,
        'retina_rad': 12.5,
    }

    # return standard optical parameters
    return opt_params


def gen_optics(params):

    ''' Generate Standard Optics

    Given set of variable parameters, generate and define parameters of each optics element in path

    Args:
        eye_dist (float): distance from image to cornea front surface
        iris (float): radius of iris opening

    Returns:
        (list of dict): list of dict for each optic in path
    '''

    eye_front = params['eye_front']


    # define aspericity by scale in direction normal to standard axis, rotation axis about centre
    cornea_sph = params['cornea_sph']
    cornea_axis = params['cornea_axis']

    cornea_pow = params['cornea_pow']
    cornea_f_rad = params['cornea_f_rad']
    cornea_r_rad = params['cornea_r_rad']

    cornea_thick = params['cornea_thick']

    cornea_f_dist = (eye_front) + cornea_f_rad * cornea_pow

    cornea_r_dist = (eye_front + cornea_thick) + cornea_r_rad


    aqueous_thick = params['aqueous_thick']

    iris_dia = params['iris_dia']
    iris_dist = (eye_front + cornea_thick + aqueous_thick) - 0.1


    lens_thick = params['lens_thick']
    focus = params['focus'] # 0.0 at inf., 1.0 at near limit

    lens_f_rad_max = params['lens_f_rad_max']
    lens_f_rad_min = params['lens_f_rad_min']

    lens_r_rad_max = params['lens_r_rad_max']
    lens_r_rad_min = params['lens_r_rad_min']

    lens_f_rad = lens_f_rad_max - (lens_f_rad_max - lens_f_rad_min) * focus
    lens_r_rad = lens_r_rad_max - (lens_r_rad_max - lens_r_rad_min) * focus

    lens_pow = params['lens_pow']

    lens_f_dist = (eye_front + cornea_thick + aqueous_thick) + lens_f_rad

    # subtract scaled radius for reverse lens
    lens_r_dist = (eye_front + cornea_thick + aqueous_thick + lens_thick) - lens_r_rad * lens_pow


    retina_thick = params['retina_thick']
    retina_rad = params['retina_rad']
    retina_dist = (eye_front + cornea_thick + aqueous_thick + lens_thick + retina_thick) - retina_rad


    opts = [
        { # cornea
            'centre': np.array([cornea_f_dist, 0., 0.]),
            'opt_den': 1.377,
            'scale': np.array([cornea_pow, cornea_sph, 1.]),
            'radius': cornea_f_rad,
            'rev': False,
        },
        { # aqueous
            'centre': np.array([cornea_r_dist, 0., 0.]),
            'opt_den': 1.336,
            'scale': np.array([1., 1., 1.]),
            'radius': cornea_r_rad,
            'rev': False,
        },
        { # iris
            'centre': np.array([iris_dist, 0., 0.]),
            'opt_den': 1.336,
            'scale': np.array([np.sqrt(1e-6), 1., 1.]),
            'radius': iris_dia / 2.,
            'rev': False,
        },
        { # lens front
            'centre': np.array([lens_f_dist, 0., 0.]),
            'opt_den': 1.411,
            'scale': np.array([1., 1., 1.]),
            'radius': lens_f_rad,
            'rev': False,
        },
        { # lens rear
            'centre': np.array([lens_r_dist, 0., 0.]),
            'opt_den': 1.337,
            'scale': np.array([lens_pow, 1., 1.]),
            'radius': lens_r_rad,
            'rev': True,
        },
        { # retina
            'centre': np.array([retina_dist, 0., 0.]),
            'opt_den': 1.337,
            'scale': np.array([1., 1., 1.]),
            'radius': retina_rad,
            'rev': True,
        },

    ]
    #print(retina_dist)


    # return optics chain
    return opts



def gen_optics_rev(params):

    ''' Generate Standard Optics

    Given set of variable parameters, generate and define parameters of each optics element in path0...........................................................................

    Args:
        eye_dist (float): distance from image to cornea front surface
        iris (float): radius of iris opening

    Returns:
        (list of dict): list of dict for each optic in path
    '''

    lens_thick = params['lens_thick']
    focus = params['focus'] # 0.0 at inf., 1.0 at near limit

    lens_f_rad_max = params['lens_f_rad_max']
    lens_f_rad_min = params['lens_f_rad_min']

    lens_r_rad_max = params['lens_r_rad_max']
    lens_r_rad_min = params['lens_r_rad_min']

    lens_f_rad = lens_f_rad_max - (lens_f_rad_max - lens_f_rad_min) * focus
    lens_r_rad = lens_r_rad_max - (lens_r_rad_max - lens_r_rad_min) * focus

    lens_pow = params['lens_pow']


    retina_thick = params['retina_thick']
    retina_rad = params['retina_rad']

    lens_r_dist = (retina_thick) + lens_r_rad * lens_pow


    lens_f_dist = (retina_thick + lens_thick) - lens_f_rad


    #iris_dia = params['iris_dia']
    #iris_dist = (retina_thick + lens_thick) + 0.1


    aqueous_thick = params['aqueous_thick']

    cornea_sph = params['cornea_sph']
    cornea_axis = params['cornea_axis']

    cornea_pow = params['cornea_pow']
    cornea_f_rad = params['cornea_f_rad']
    cornea_r_rad = params['cornea_r_rad']

    cornea_thick = params['cornea_thick']


    cornea_r_dist = (retina_thick + lens_thick + aqueous_thick) - cornea_r_rad

    cornea_f_dist = (retina_thick + lens_thick + aqueous_thick + cornea_thick) - cornea_f_rad * cornea_pow


    eye_front = params['eye_front']

    image_dist = (retina_thick + lens_thick + aqueous_thick + cornea_thick) + eye_front


    opts = [
        { # lens rear
            'centre': np.array([lens_r_dist, 0., 0.]),
            'opt_den': 1.411,
            'scale': np.array([lens_pow, 1., 1.]),
            'radius': lens_r_rad,
            'rev': False,
        },
        { # lens front
            'centre': np.array([lens_f_dist, 0., 0.]),
            'opt_den': 1.336,
            'scale': np.array([1., 1., 1.]),
            'radius': lens_f_rad,
            'rev': True,
        },
        { # aqueous
            'centre': np.array([cornea_r_dist, 0., 0.]),
            'opt_den': 1.377,
            'scale': np.array([1., 1., 1.]),
            'radius': cornea_r_rad,
            'rev': True,
        },
        { # cornea
            'centre': np.array([cornea_f_dist, 0., 0.]),
            'opt_den': 1.,
            'scale': np.array([cornea_pow, cornea_sph, 1.]),
            'radius': cornea_f_rad,
            'rev': True,
        },
        { # image
            'centre': np.array([image_dist, 0., 0.]),
            'opt_den': 1.,
            'scale': np.array([np.sqrt(1e-6), 1., 1.]),
            'radius': 50.,
            'rev': True,
        },

    ]


    # return optics chain
    return opts

