
''' Helper Display Functions

    helper functions for optical model display and testing, framework matplotlib 3d


    functions:

        x,y,z = plot_3d_ellipsoid(C, r, e, rev, theta = 0., full = False)

        x,y,z = plot_3d_line(o, v, xl)

        img = get_retinal_img(paths, r, d)

'''



''' Imports '''

# nd array manipulation
import numpy as np



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

