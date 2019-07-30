
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



### Path Tracing Orchestration Functions

## store points for initial and each intercept with direction unit vector

# build array of initial ray locations with direction vectors

# define optics in path

# calculate full path of each individual ray through optics chain, save path details

## finally plot optics and rays for display, use ray path details to build



    functions:

        icept = get_intercept(C, r, e, L0, v)

        n = get_surface_normal(C, e, p)

        V = get_refracted_vector(n1, n2, N, v)

        points = refraction(C, r, e, L0, v, n1, n2, rev)

        path = get_path(rays, opts, n0 = 1.0)

'''




''' Imports '''

## to import functions from file
#from .file import func1, func2


# nd array manipulation
import numpy as np



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



def get_path(ray, opts, n0):

    ''' Calculate Ray Path through Optics

        Given rays and optics lists, calculate full ray path through optics, return all intercept points and vectors

    Args:
        rays (list of np.array[3,1]): list of ray starting points
        optics (list of dict(optics)): list of optics in path with parameters for each
        n0 (float): initial index of refraction at ray start location; default to air (1.0)

    Returns:
        (list of list pair np.array[3,1]): list of intercept vector pairs (point and direction unit vector)
    '''

    # set initial refractive index
    n1 = n0

    # get origin ray position and direction vector
    L0 = ray[:3]
    v = ray[3:]

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
        _ray = refraction(C, r, e, L0, v, n1, n2, rev)


        # if refracted ray, unpack
        if _ray is not None:

            # check which side intersection is wanted
            if rev:
                p, V = _ray[0]
            else:
                p, V = _ray[1]


            # store ray path details
            path.append([p, V])


            # update new ray start
            n1 = n2
            v = V
            L0 = p

        # if not refracted ray, terminate path
        else:
            break


    # return full ray path
    return path

