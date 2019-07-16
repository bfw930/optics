
''' Path Tracing Engine

    Core functions for path tracing engine


    chain of optics as functions in spatial dimensions
    rays cast from image source

    throw ray
        equation as line from point source

    ray intersection with optic
        calculate coordinates of intercept

    new ray path
        calculate tangent on optic function at intersection coordinates
        calculate refraction angle using snells law (use optical constants and tangent as interface)
        update ray path with new angle, convert to gradient of line with reference to dimensions axis

    throw new ray into next optical element
        propagate


        Need to extent initial ray positions for additional dimension

'''




''' Imports '''

## import functions from file
#from .file import func1, func2


# nd array manipulation
import numpy as np

# image import, processing, export
import PIL

# image manipulation
from scipy import ndimage




''' Core Path Tracing Calculation Functions '''

def get_intercept(C, r, e, L0, v):

    ''' Get Intercept of Line and Ellipsoid

        Given ellipsoid and line parameters, calculate point of intercept(s)
        ref: http://www.illusioncatalyst.com/notes_files/mathematics/line_nu_sphere_intersection.php
        test / verify: http://www.ambrsoft.com/TrigoCalc/Sphere/SpherLineIntersection_.htm

    Args:
        C (np.array[3,1]): ellipsoid centre point
        r (float): ellipsoid radius
        e (np.array[3,1]): ellipsoid dimensional scaling factors

        L0 (np.array[3,1]): line origin point
        v (np.array[3,1]): line direction unit vector

    Returns:
        (list of np.array[3,1]): list of line-ellipsoid intercept point(s)
    '''

    # ellipse centre, radius, and axis vectors
    C = np.array([C]).T

    # point on line (origin) and direction vector
    L0 = np.array([L0]).T
    v = np.array([v]).T


    # Translation Matrix
    T = np.identity(4)
    T[:-1,-1:] = C


## fix rotation if required in future

    # Rotation/Scale Matrix ?
    R = np.identity(4)
    R[:-1,:-1] = np.identity(3) * np.array(e).T

    a = np.array([[1., 0., 0.]]).T# * e[0]
    b = np.array([[0., 1., 0.]]).T# * e[1]
    c = np.array([[0., 0., 1.]]).T# * e[2]


    # Scale/Rotate Matrix ?
    #S = np.identity(4)
    #S[0,0] = np.linalg.norm(a)
    #S[1,1] = np.linalg.norm(b)
    #S[2,2] = np.linalg.norm(c)


    # Transformation Matrix
    M = T @ R #@ S


    # rescaled sphere centre
    C_ = np.array([[0, 0, 0, 1]]).T

    # transformed line origin
    L0_ = np.linalg.inv(M) @ np.row_stack([L0,[1]])
    # difference line vector
    w = L0_ - C_

    # transformed line vector
    v_ = np.linalg.inv(M) @ np.row_stack([v,[0]])


    # coefficients of quadratic intersection
    a_ = (v_.T @ v_)[0]
    b_ = 2*((v_.T @ w)[0])
    c_ = (w.T @ w)[0] - r**2

    # descriminant of quadratic
    d_ = b_**2 - 4*a_*c_


    # intersection points
    if d_ > 0:
        t1 = (-b_ + np.sqrt(d_)) / (2*a_)
        t2 = (-b_ - np.sqrt(d_)) / (2*a_)
        L_int1 = L0 + t1 * v
        L_int2 = L0 + t2 * v
        icept = [L_int1, L_int2]


    # tangent intercept
    elif d_ == 0:
        t = -b_ / (2*a_)
        L_int = L0 + t * v
        icept = [L_int]


    # no intercept
    else:
        icept = []


    # return intercepts
    return icept



def get_surface_normal(C, e, p):

    ''' Get Ellipsoid Tangent Plane Normal at Point

        Given ellipsoid parameters and point on surface, calculate tangent plane surface normal at point on ellipse

    Args:
        C (np.array[3,1]): ellipsoid centre point
        e (np.array[3,1]): ellipsoid dimensional scaling factors
        p (np.array[3,1]): point on ellipsoid surface

    Returns:
        (np.array[3,1]): surface normal unit vector
    '''

    # calculate normal vector using coefficients of equation of plane at intercept point
    #(x*x0/a**2) + (y*y0/b**2) + (z*z0/c**2) = r
    n = (p-C) / e**2

    # normalise to unit vector
    n = n / np.linalg.norm(n)


    # return surface normal
    return n



def get_refracted_vector(n1, n2, N, v):

    ''' Get Refraction Vector at Line-Ellipsoid Intersection Point

        Given incidies of refraction, surface normal vector, and incident line direction unit vector, calculate
        direction unit vector of refraction through surface using Snells Law
        ref: http://www.starkeffects.com/snells-law-vector.shtml

    Args:
        n1 (float): above surface index of refraction
        n2 (float): below surface index of refraction
        N (np.array[3,1]): refraction surface normal unit vector
        v (np.array[3,1]): incident line direction unit vector

    Returns:
        (np.array[3,1]): refraction direction unit vector
    '''

    # calculate vector of refraction using incident vector and surface normal
    V = (n1/n2) * np.cross(N, np.cross(-N, v)) - N * np.sqrt(1 - ((n1/n2)**2) * (np.cross(N, v) @ np.cross(N, v)))


    # return refraction direction unit vector
    return V



''' Path Tracing Orchestration Functions '''


## store points for initial and each intercept with direction unit vector

# build array of initial ray locations with direction vectors

# define optics in path

# calculate full path of each individual ray through optics chain, save path details

## finally plot optics and rays for display, use ray path details to build



def refraction(C, r, e, L0, v, n1, n2, rev):

    ''' Calculate Ellipsoid-Ray Intercept Point(s) and Refraction Vector(s)

        Given ray and optic parameters, calculate refraction points and direction unit vectors

    Args:
        C (np.array[3,1]): ellipsoid centre point
        r (float): ellipsoid radius
        e (np.array[3,1]): ellipsoid dimensional scaling factors
        rev (bool): reverse flag for ellipse intercept surface (front / rear) along primary axis

        L0 (np.array[3,1]): line origin point
        v (np.array[3,1]): line direction unit vector

        n1 (float): above surface index of refraction
        n2 (float): below surface index of refraction

    Returns:
        (list of list pair np.array[3,1]): list of intercept vector pairs (point and direction unit vector)
    '''

    # calculate ray ellipsoid intersection
    icepts = get_intercept(C, r, e, L0, v)


    # only calculate refraction for intersections (exclude tangents)
    if len(icepts) > 1:

        points = []

        # iterate intercept points
        for icept in icepts:

            # define point of intercept (on surface normal)
            p = icept[:,0]


            # calculate surface normal vector
            n = get_surface_normal(C, e, p)

            # invert surface normal if reversed optic
            if rev:
                n = -n


            # get unit vector of refracted ray using snells law
            V = get_refracted_vector(n1, n2, n, v)


            # store intercept point and refraction direction unit vector
            points.append([p, V])


        # return intercept point and refraction ray direction vector
        return points


    # if no intercepts found, return None
    return None



def get_paths(rays, opts, n0 = 1.0):

    ''' Calculate Ray Path through Optics

        Given rays and optics lists, calculate full ray path through optics, return all intercept points and vectors

    Args:
        rays (list of np.array[3,1]): list of ray starting points
        optics (list of dict(optics)): list of optics in path with parameters for each
        n0 (float): initial index of refraction at ray start location; default to air (1.0)

    Returns:
        (list of list pair np.array[3,1]): list of intercept vector pairs (point and direction unit vector)
    '''

    paths = []

    # iterate over rays
    for i in range(len(rays)):

        # set initial refractive index
        n1 = n0

        # get origin ray position and direction vector
        L0 = rays[i][:3]
        v = rays[i][3:]

        path = []

        # store initial ray path details
        path.append([L0, v])


        # iterate over optics in path
        for optic in opts:

            # get optic (ellipsoid): centre, radius, axes scale, refractive index, reverse flag
            C = optic['centre']
            r = optic['radius']
            e = optic['scale']
            n2 = optic['opt_den']
            rev = optic['rev']


            # for intercept ray with optic, get new refracted ray
            ray = refraction(C, r, e, L0, v, n1, n2, rev)


            # if refracted ray, unpack
            if ray is not None:

                # check which side intersection is wanted
                if rev:
                    p, V = ray[0]
                else:
                    p, V = ray[1]


                # store ray path details
                path.append([p, V])


                # update new ray start
                n1 = n2
                v = V
                L0 = p

            # if not refracted ray, terminate path
            else:
                break


        # store full ray path
        paths.append(path)


    # return list of ray paths
    return paths




''' Optics Chain Generation Functions '''


### helper function to generate full optics chain (list of dicts with parameters) from minimal input parameter set
## have relative distance of optics

## include parameter variation aligned to optical aberations (prescription), eg retinal distance, cornea ashpericity


## standard optics system


def gen_optics(params):

    ''' Generate Standard Optics

    Given set of variable parameters, generate and define parameters of each optics element in path0...........................................................................

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





''' Retinal Image Transformation Functions '''

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







def get_retinal_img(paths, r, d):

    ''' Get Ellipse Points for Plotting in 3D

        Given ellipsoid parameters, calculate points in 3d for plotting

    Args:
        paths (list of list pair np.array[3,1]): list of intercept vector pairs (point and direction unit vector)
        r (float): retinal ellipsoid radius
        d (float): primary axis distance of retinal optic centre

    Returns:
        x,y,z (np.array): arrays of ellipse points each dimension
    '''

    # get retinal image
    r_img = np.stack([ path[-1][0] for path in paths ])


## update filter based on input length of optics chain

    # get index of rays that reach retina, filter
    j = [ i for i in range(len(paths)) if len(paths[i]) == 7 ]
    r_img = r_img[j]


    # calculate stereographic projection onto 2d plane from 3d retinal image
    z = r_img[:,0] - d - r
    x = ( r_img[:,1] / (0.5 - z) ) * r
    y = ( r_img[:,2] / (0.5 - z) ) * r


    # return 2d retinal image array
    return np.stack([x, y]).T



'''

    import test image and supersample (1 px to 9 rays) initial positions

    define standard optical parameters

    calculate ray paths, final retinal image locations

    store retinal pixel map at standard optics


    adjust optical parameters for stigmatism

    set initial ray locations as standard retinal image

    calculate ray paths in reverse (retina to image location)

    save pixel map for optical configuration


    use pixel map to subsample over image grid, generate image pixels from average of standard pixels


    build database of pixel maps for given range and step of optical parameters each stigmatism

    generate image translation suit for given input image (greyscale) using pixel maps (parameter dependent)

'''







## import test image





## get mapping of test image to initial ray points


# define length of image dimensions (allow scale of image, bin average pixels per region)
# uniform divide pixels into spatial map of initial ray y,z coordinates







## calculate retinal projection of image, mapping of image to retina




## display retinal image projection: wide focus and fovea zoom




## generate inverted optics chain (ray start retinal image mapping)



## reverse path trace retinal image back to origin

#### confirm accurate mapping to original image




## vary optics parameters in forward / reverse optics chain, compare image to retinal image mapping, distortion



## develop methodology to recreate perceptive refraction image for clinical testing
# test image baseline retinal mapping
# inverse optic parameter variation through reverse optics chain
# generate perceptive distortion image for display


## generate test image sweep over parameter space for distorted image set




''' Display Helper Functions '''

def plot_3d_ellipsoid(C, r, e, rev, theta = 0., full = False):

    ''' Get Ellipse Points for Plotting in 3D

        Given ellipsoid parameters, calculate points in 3d for plotting

    Args:
        C (np.array[3,1]): ellipsoid centre point
        r (float): ellipsoid radius
        e (np.array[3,1]): ellipsoid dimensional scaling factors
        rev (bool): reverse flag for half ellipse along primary axis

    Returns:
        x,y,z (np.array): arrays of ellipse points each dimension
    '''

    # generate and plot lens (ellipsoid)
    rx, ry, rz = np.array(e) * r


    # set of all spherical angles:
    #u = np.linspace(0.5 * np.pi, 1.5 * np.pi, 50)
    if rev:
        u = np.linspace(1.5 * np.pi, 2.5 * np.pi, 50)
    else:
        u = np.linspace(0.5 * np.pi, 1.5 * np.pi, 50)

    if full:
        k = np.linspace(0., 2*np.pi, 50)
    else:
        k = np.linspace(0., np.pi, 50)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    x = rx * np.outer(np.cos(u), np.sin(k))# + C[0]
    y = ry * np.outer(np.sin(u), np.sin(k))# + C[1]
    z = rz * np.outer(np.ones_like(u), np.cos(k))# + C[2]


    # pure rotation about primary axis, input angle in degrees
    theta = np.deg2rad(theta)
    rotate = np.array([[1,0,0], [0, np.cos(theta), -np.sin(theta)], [0, np.sin(theta), np.cos(theta)]])

    for i in range(len(x)):
        for j in range(len(x)):
            [x[i,j],y[i,j],z[i,j]] = np.dot([x[i,j],y[i,j],z[i,j]], rotate) + C


    # return dimension points for ellipse to plot
    return x,y,z





def plot_3d_line(o, v, xl):

    ''' Get Line Points for Plotting in 3D

        Given line parameters, calculate points in 3d for plotting

    Args:
        c (np.array[3,1]): line origin point
        v (np.array[3,1]): line direction unit vector
        xl (float): length of line along primary axis

    Returns:
        x,y,z (np.array): arrays of line points each dimension
    '''

    # define ray line points and plot
    line = np.array([o]).T + np.arange(0., xl + .1, .1) * np.array([v]).T
    x = line[0,]
    y = line[1,]
    z = line[2,]


    # return dimension points for line to plot
    return x,y,z
