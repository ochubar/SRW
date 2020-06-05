# -*- coding: utf-8 -*-
#############################################################################
# Generate a 2D random pattern of disks
# Authors/Contributors: Yugang Zhang, Rebecca Ann Coles
# March 01, 2019
# April 11, 2020
#############################################################################

import os
import sys
import numpy as np

from skimage.draw import polygon, circle, ellipse

# ********************** Get the distance between point jj and the other points
def get_r1j( px, py, jj):
    '''Used with Yugang random walk function.
    Get the distance between point jj and the other points

    :param px: mesgrid, coordination of x, shape will be [point number, point number]
    :param py: mesgrid, coordination of y, shape will be [point number, point number]
    :param jj: the index of the position
    :return: the distance between the j and other points, shape be point number -1
    '''

    pxf = np.ravel( px )
    pyf = np.ravel( py )
    pxn = np.delete( pxf, jj, 0)
    pyn = np.delete( pyf, jj, 0)
    dd = np.sqrt(  ( pxn - pxf[jj] )**2 + ( pyn - pyf[jj] )**2 )
    return dd


# ********************** Check the distance between point jj and other points 
def chk_r1j_rmin(px,py, jj, rmin):
    '''Used with Yugang random walk function.
    Check the distance between point j and other points if all the distances 
    larger than rmin, return True.

    :param px: mesgrid, coordination of x, shape will be [point number, point number]
    :param py: mesgrid, coordination of y, shape will be [point number, point number]
    :param jj: the index of the position
    :param rmin: the object minimum distance from other objects.
    :return: bool, True = distance is smaller than rmin, need to move position.
        False = distance larger than rmin, no need to move position.
    '''

    r1j = get_r1j( px, py, jj)

    if len(np.where( r1j < rmin)[0]) >= 1:  #distance smaller than rmin, need to move position 
        return True
    else:
        return False


# ********************** Move 2D points with random direction and with amplitue as amp for each step
def mv_2D( px, py, amp): 
    '''Used with Yugang random walk function.
    Move 2D points with random direction and with amplitue as amp for each step 
    
    :param px: mesgrid, coordination of x, shape will be [point number, point number]
    :param py: mesgrid, coordination of y, shape will be [point number, point number]
    :param amp: amplitue of each step.
    :return: px/py mesgrid, coordination of x/y, shape will be [point number, point number]
        or px/py: one float, in case of move one point
    '''

    NA=False
    if not isinstance(px, (list, np.ndarray) ):
        px = np.array( [px])
        py = np.array( [py] )
        NA = True

    theta = np.random.uniform(-0.5,0.5,(len(px),len(py))) * 2* np.pi
    dx = amp * np.cos( theta )
    dy = amp * np.sin( theta ) 

    if NA:
        return (px+dx)[0], (py+dy)[0]
    else:
        return px + dx, py + dy


# ********************** Simulate a 2D random walk 
def get_rnd_2D( px, py, amp, rmin=16, try_max=1000):
    '''Yugang's random walk function.
    Simulate a 2D random walk. Each of the 2D points, for each step, move with 
    random direction and with amplitue as amp.

    :param px: mesgrid, coordination of x, shape will be [point number, point number]
    :param py: mesgrid, coordination of y, shape will be [point number, point number]
    :param amp: amplitue of each step.
    :param rmin: the object minimum distance from other objects.
    :param try_max: try number to acheive rmin, if can't achevie print that point
    :return: px/py mesgrid, coordination of x/y, shape will be [point number, point number]
    '''

    px = np.ravel(px)
    py = np.ravel(py)
    for ii in range( px.size):
        trys = 0
        mov = True
        while mov:
            px[ii], py[ii] =  mv_2D( px[ii], py[ii], amp)
            mov = chk_r1j_rmin(px,py, ii, rmin)
            trys +=1
            if trys > try_max:
                mov = False
                print(ii, 'Not satify the rmin distance')
    return px,py


def uni_rnd_seed(num, rx_pixels, ry_pixels, obj_max_size, min_dist=16): #RAC25032020
    '''
    Uniform Random Seed Object Locations
    Evenly space "num" points over a (rx, ry) grid.

    :param num: particle number.
    :param rx_pixels: range of the horizontal coordinate [pixels] for which 
        the transmission is defined.
    :param ry_pixels: range of the vertical coordinate [pixels] for which 
        the transmission is defined.
    :param _obj_size_max: maximum allowed size of objects [pixels]
    :param min_dist = 16: minimum distance [pixels] between objects.    

    For num_x points in x-direction with spacing delta_x in a 2D array of width
    "w" and height "h":
        num_x = sqrt[(w/h)*num + (w-h)^2 / (4*h^2)] - [(w-h) / (2*h)]
        
    For num_y points in y-direction with spacing delta_y in a 2D array of width
    "w" and height "h":
        num_y = num / num_x

    To get (x,y) point pairs you need delta_x = delta_y. If delta_x != delta_y
    for num objects, num is reduced until delta_x = delta_y.
    '''
    # Get spacing (delta) between points. 
    #  If spacing can't be compleated with the requested
    #  number of points, the density of points will be reduced until 
    #  the points can all be uniformly spaced.
    final_num_points = num #initalize to maximum number of points
    delta_x = 0 #initalize to delta_x != delta_y
    delta_y = 1 #initalize to delta_x != delta_y
    dist_edge = obj_max_size*2 #minimum distance of an object from the edge of the surface
    rx_boundry = rx_pixels-dist_edge*2 #end edge boundry for objects
    ry_boundry = ry_pixels-dist_edge*2 #top edge boundry for objects

    while delta_x != delta_y: #delta_x=delta_y needed to evenly spread objects
        num_x = int(np.sqrt( (rx_boundry/ry_boundry)*final_num_points + 
                            np.power((rx_boundry-ry_boundry),2) / 
                            np.power((4*ry_boundry),2) ) - 
                            (rx_boundry-ry_boundry)/(2*ry_boundry))
        num_y = int(final_num_points/num_x)
        delta_x = int(rx_boundry/(num_x-1))
        delta_y = int(ry_boundry/(num_y-1))
        if delta_x != delta_y:
            # return a default set of objects if delta_x=delta_y fails
            if final_num_points <= 0:
                delta_x = delta_y = 1
                final_num_points = final_num_points/2
            # reduce num until delta_x=delta_y
            final_num_points = final_num_points-1
        else:
            # if delta_x=delta_y, set delta
            delta = delta_x

    # Construct px and py made up of uniformly spaced points in (rx, ry)
    px_edge = [] #x coord list of uniformly spaced points
    py_edge = [] #y coord list of uniformly spaced points
    rows = num_y #number of rows (y values)
    columns = num_x #number of columns (x values) 
    for ii in range(rows-1):
        #fill px and py with uniformly spaced points (row by row)
        row_value = ii*delta
        for jj in range(columns-1):
            px_edge.append(jj*delta)
            py_edge.append(row_value)

    # Enforce object distance from left an bottom surface boundry        
    px = [kk+dist_edge for kk in px_edge]
    py = [kk+dist_edge for kk in py_edge]

    # Create random noise to perturb points
    minimum_movement = obj_max_size + min_dist
    if delta <= minimum_movement:
        print("Warning: object_maximum_size + minimum_dist_between_obj is" 
              " greater than perfectly uniform object spaceing. To preserve"
              + " minimum_dist_between_obj you will need to reduce object_maximum_size"
              + " or decrease object density.")
    noise = np.random.uniform(low=minimum_movement,
                              high=delta,
                              size=(len(px), 2))

    #add noise to px and py
    px = np.asarray(px, dtype=np.float32)
    py = np.asarray(py, dtype=np.float32) 
    px += noise[:,0]
    py += noise[:,1]

    return py, px

# ********************** Place shape on meshgrid
def on_pxy(px, py, bx=100, by=None,
           _obj_type=1, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
           _ang_min = 0, _ang_max = 0, _ang_dist = 1,
           _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Place shapes on (px,py) points

    :param px: mesgrid, coordination of x, shape will be [point number, point number]
    :param py: mesgrid, coordination of y, shape will be [point number, point number]
    :param bx: box size (x), namely each frame size [array size (pixels)]
    :param by: box size (y), namely each frame size [array size (pixels)]
    :param _obj_type = 1: shape of objects. 
        Choices: 1=rectangle, 2=ellipse, 3=triangle, 4=polygon, 5=random_shapes
    :param _obj_size_min = 10: minimum allowed size of objects [pixels]
    :param _obj_size_max = 10: maximum allowed size of objects [pixels]
    :param _size_dist = 1: distribution of sizes. 
        Choices are: 1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of shape [degrees].
    :param _ang_max = 0: maximum rotation angle of shape [degrees].
    :param _ang_dist = 1: distribution of rotation angle. 
        Choices are: 1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: 
        rectangle: ratio of width to length (length is the longest dimension). 
            Leaving value as None will give 1.0 (square) as width to length 
            ratio unless _obj_par2 = True is selected to randomize the 
            width to length ratio.
        ellipse: ratio of minor to major semi-axes. 
            Leaving value as None will give 1.0 (circle) as minor to major 
            semi-axes ratio unless _obj_par2 = True is selected to 
            randomize the minor to major semi-axes ratio.
        triangle: ratio of height (y, or "opposite") to width (x, or adjunct). 
            Leaving value as None will give 1.0 (1/1) as height to width 
            ratio unless _obj_par2 = True is selected to randomize the 
            height to width ratio (default is an equilateral triangle).
        regular polygon: number of polygon vertices.
            Leaving value as None will give a hexagon (6 sides)
            unless _obj_par2 = True is selected to randomize the numer of 
            polygon vertices. Max vertices = 12.
        random shapes: Which shapes to randomly generate.
            Choices are:
                [1='rectangle',2='ellipse',3='triangle',4='polygon'].
            Leaving value as None will give all shape options as:[1,2,3,4]
    :param _obj_par2 = None:
        rectangle: value set to True will randomize the width to length ratio.
        ellipse: Value set to True will randomize the minor to major semi-axes ratio.
        triangle: value set to True will randomize the height (y, or "opposite")
            to width (x, or adjunct) ratio.
        regular polygon: Value set to True will randomize the numer of polygon 
            vertices.
    :return:array with shapes at (px,py) locations.
    '''

    if by is None:
        by = bx
    dd = np.zeros( [bx,by] , dtype = bool)   

    for i in range(px.size): 
        cx = np.ravel(px)[i]
        cy = np.ravel(py)[i]
        xx, yy = get_shape(cx, cy, _obj_type = _obj_type, 
                           _obj_size_min = _obj_size_min, _obj_size_max = _obj_size_max, _size_dist = _size_dist,
                           _ang_min = _ang_min, _ang_max = _ang_max, _ang_dist = _ang_dist,
                           _obj_par1 = _obj_par1, _obj_par2 = _obj_par2)
        ww = np.where( ( (xx>=0)  & (xx<bx) & (yy>=0) & (yy<by) ) )            
        dd[xx[ww],yy[ww]]=1 

    return dd


# ********************** Get shape in accordance with user specified parameters.
def get_shape(cx, cy, _obj_type = 1, 
              _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
              _ang_min = 0, _ang_max = 0, _ang_dist = 1,
              _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Get shape in accordance with user specified parameters.

    :param cx: x coordinate of object center.
    :param cy: y coordinate of object center.
    :param _obj_type = 1: shape of objects. Choices are: 
        1=rectangle, 2=ellipse, 3=triangle, 4=ploygon, 5=random_shapes
    :param _obj_size_min = 10: minimum allowed size of objects [pixels]
    :param _obj_size_max = 10: maximum allowed size of objects [pixels]
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of objects [degrees].
    :param _ang_max = 0: maximum rotation angle of objects [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: 
        rectangle: ratio of width to length (length is the longest dimension). 
            Leaving value as None will give 1.0 (square) as width to length 
            ratio unless _obj_par2 = True is selected to randomize the 
            width to length ratio.
        ellipse: ratio of minor to major semi-axes. 
            Leaving value as None will give 1.0 (circle) as minor to major 
            semi-axes ratio unless _obj_par2 = True is selected to 
            randomize the minor to major semi-axes ratio.
        triangle: ratio of height (y, or "opposite") to width (x, or adjunct). 
            Leaving value as None will give 1.0 (1/1) as height to width 
            ratio unless _obj_par2 = True is selected to randomize the 
            height to width ratio (default is an equilateral triangle).
        regular polygon: number of polygon vertices.
            Leaving value as None will give a hexagon (6 sides)
            unless _obj_par2 = True is selected to randomize the numer of 
            polygon vertices.
        random shapes: Which shapes to randomly generate.
            Choices are [1='rectangle',2='ellipse',3='triangle',4='polygon'].
            Leaving value as None will give all shape options as:
                [1,2,3,4]
    :param _obj_par2 = None:
        rectangle: value set to True will randomize the width 
            to length ratio.
        ellipse: Value set to True will randomize the minor to major 
            semi-axes ratio.
        triangle: value set to True will randomize the height (y, or "opposite")
            to width (x, or adjunct) ratio.
        regular polygon: Value set to True will randomize the numer of 
            polygon vertices. Max vertices = 12.
    '''

    if _obj_type == 1: #rectangle
        return get_rec(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist, 
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1, 
                       _obj_par2 = _obj_par2)   
    elif _obj_type == 2: #ellipse
        return get_elp(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist,
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1, 
                       _obj_par2 = _obj_par2)  
    elif _obj_type == 3: #triangle
        return get_tri(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist,
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1, 
                       _obj_par2 = _obj_par2)
    elif _obj_type == 4: #polygon
        return get_pol(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist,
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1, 
                       _obj_par2 = _obj_par2) 
    elif _obj_type == 5: #random
        return get_rnd(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist, 
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1)
    else: #ellipse
        return get_elp(cx, cy, 
                       _obj_size_min = _obj_size_min, 
                       _obj_size_max = _obj_size_max, 
                       _size_dist = _size_dist,
                       _ang_min = _ang_min, 
                       _ang_max = _ang_max, 
                       _ang_dist = _ang_dist, 
                       _obj_par1 = _obj_par1, 
                       _obj_par2 = _obj_par2)

# ********************** Create rectangle
def get_rec(cx, cy, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
            _ang_min = 0, _ang_max = 0, _ang_dist = 1,
            _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Rectangles

    :param cx: x coordinate of rectangle center.
    :param cy: y coordinate of rectangle center.
    :param _obj_size_min = 10: minimum allowed size of rectangle [pixels]
    :param _obj_size_max = 10: maximum allowed size of rectangle [pixels]
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of rectangle [degrees].
    :param _ang_max = 0: maximum rotation angle of rectangle [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: ratio of width to length (length is the 
        longest dimension). Leaving value as None will give 1.0 (square) as 
        width to length ratio unless _obj_par2 = True is selected to 
        randomize the width to length ratio.
    :param _obj_par2 = None: value set to True will randomize the width 
        to length ratio.
    '''

    #Use _obj_par1 (width to length ratio) or create random ratio 
    # if _obj_par2 = True
    default_side_ratio = 1 # 1/1 Square
    width, length = obj_opt_par(_obj_type = 1,
                                _obj_size_min = _obj_size_min,
                                _obj_size_max = _obj_size_max,
                                _size_dist = _size_dist,
                                _default_axis_ratio = default_side_ratio,
                                _obj_par1 = _obj_par1,
                                _obj_par2 = _obj_par2)
    
    #Create rectangle and get vertices
    _rows = np.array([cy-(width), cy-(width), 
                      cy+(width), cy+(width)])
    _columns = np.array([cx-(length), cx+(length), 
                         cx+(length), cx-(length)])

    ##Rotate the coordinates of vertices of the rectangle
    # Rotate a point counterclockwise by a given angle (degrees) around center
    rotated_rows, rotated_columns = rot_obj_vert(cx, cy, _rows, _columns,
                                                 _ang_min = _ang_min, _ang_max = _ang_max,
                                                 _ang_dist = _ang_dist)

    #Create rectangle (will use poloygon function so that we can create the rotated object)
    return polygon(rotated_rows, rotated_columns) 


# ********************** Create ellipse    
def get_elp(cx, cy, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
            _ang_min = 0, _ang_max = 0, _ang_dist = 1,
            _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Ellipses

    :param cx = x: coordinate of ellipse center.
    :param cy = y: coordinate of ellipse center.
    :param _obj_size_min = 10: minimum allowed size of ellipse [pixels]
    :param _obj_size_max = 10: maximum allowed size of ellipse [pixels]
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of ellipse [degrees].
    :param _ang_max = 0: maximum rotation angle of ellipse [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: ratio of minor to major semi-axes. Leaving 
        value as None will give 1.0 (circle) as minor to major semi-axes 
        ratio unless _obj_par2 = True is selected to randomize the ratio.
    :param _obj_par2 = None: value set to True will randomize the minor 
        to major semi-axes ratio.
    '''

    #Use _obj_par1 (minor to major semi-axes ratio) or create 
    #  random ratio if _obj_par2 = True
    default_axis_ratio = 1.0 #circle
    r_radius, c_radius = obj_opt_par(_obj_type = 2,
                                     _obj_size_min = _obj_size_min,
                                     _obj_size_max = _obj_size_max,
                                     _size_dist = _size_dist,
                                     _default_axis_ratio = default_axis_ratio,
                                     _obj_par1 = _obj_par1,
                                     _obj_par2 = _obj_par2) #minor axis, major axis

    try:
        _ang_max >= _ang_min
    except:
        print("Maximum angle of rotation is less than than minimum angle of rotation")

    #Select random rotation angle for ellipse
    _obj_ang = get_dist(_min_value = _ang_min, _max_value = _ang_max, _dist = _ang_dist)

    #Create ellipse
    return ellipse(cx, cy, r_radius, c_radius, rotation=np.deg2rad(-_obj_ang))


# ********************** Create triangle
def get_tri(cx, cy, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
            _ang_min = 0, _ang_max = 0, _ang_dist = 1, _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Triangles

    :param cx: x coordinate of triangle center.
    :param cy: y coordinate of triangle center.
    :param _obj_size_min = 10: minimum allowed size of triangle [pixels]
    :param _obj_size_max = 10: maximum allowed size of triangle [pixels]
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of triangle [degrees].
    :param _ang_max = 0: maximum rotation angle of triangle [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: ratio of height (y, or "opposite") to width 
        (x, or adjunct). Leaving value as None will give 1.0 (1/1) as 
        height to width ratio unless _obj_par2 = True is selected 
        to randomize the height to width ratio (default is an 
        equilateral triangle).
    :param _obj_par2 = None: value set to True will randomize the height 
        to width ratio.
    '''

    #Use _obj_par1 (height (y, or "opposite") to width (x, or adjunct) 
    #  ratio) or create random ratio if _obj_par2 = True
    default_side_ratio = 1.0 #1/1 equilateral triangle
    height, width = obj_opt_par(_obj_type = 3,
                                _obj_size_min = _obj_size_min,
                                _obj_size_max = _obj_size_max,
                                _size_dist = _size_dist,
                                _default_axis_ratio = default_side_ratio,
                                _obj_par1 = _obj_par1,
                                _obj_par2 = _obj_par2)

    #Create triangle and get vertices
    _rows = np.array([cy-(height), cy+(height), cy+(height)])
    _columns = np.array([cx-(width), cx+(width), cx-(width)])

    ##Rotate the coordinates of vertices of the triangle
    # Rotate a point counterclockwise by a given angle (degrees) around center
    rotated_rows, rotated_columns = rot_obj_vert(cx, cy, _rows, _columns,
                                                 _ang_min = _ang_min,
                                                 _ang_max = _ang_max,
                                                 _ang_dist = _ang_dist)

    #Create triangle (will use poloygon function so that we can create the rotated the object)
    return polygon(rotated_rows, rotated_columns)    


# ********************** Create polygon
def get_pol(cx, cy, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
            _ang_min = 0, _ang_max = 0, _ang_dist = 1, _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Regular Polygon

    :param cx: x coordinate of polygon center.
    :param cy: y coordinate of polygon center.
    :param _obj_size_min = 10: minimum allowed size of polygon [pixels].
    :param _obj_size_max = 10: maximum allowed size of polygon [pixels].
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of ploygon [degrees].
    :param _ang_max = 0: maximum rotation angle of polygon [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: number of polygon vertices.
        Leaving value as None will give a hexagon (6 sides) unless 
        _obj_par2 = True is selected to randomize the numer of polygon 
        vertices.
    :param _obj_par2 = None: Value set to True will randomize the numer 
        of polygon vertices. Max vertices = 12.
    '''

    #Use _obj_par1 (number of vertices) to get radius and num_polygon_vertices
    #  or create random number of vertices if _obj_par2 = True
    #  Function "obj_opt_par" sets the default number of regular polygon
    #  vertices to "num_polygon_vertices = 6" (hexagon).
    radius, num_polygon_vertices = obj_opt_par(_obj_type = 4,
                                               _obj_size_min = _obj_size_min,
                                               _obj_size_max = _obj_size_max,
                                               _size_dist = _size_dist,
                                               _obj_par1 = _obj_par1,
                                               _obj_par2 = _obj_par2)

    #Create regular polygon and get vertices
    # rr = radius, ii = vertex, nn = total vertices, aa = angle
    # vx = cx + rr * cos(2 * pi * ii / nn + aa)
    # vy = cy + rr * sin(2 * pi * ii / nn + aa)
    _columns = [] #x values
    _rows = [] #y values
    aa = 360/num_polygon_vertices
    for ii in range(num_polygon_vertices):
        _columns.append(int(cx + radius * np.cos(2 * np.pi * ii / num_polygon_vertices + aa)))
        _rows.append(int(cy + radius * np.sin(2 * np.pi * ii / num_polygon_vertices + aa)))

    ##Rotate the coordinates of vertices of the regular polygon
    # Rotate a point counterclockwise by a given angle (degrees) around center
    rotated_rows, rotated_columns = rot_obj_vert(cx, cy, 
                                                 _rows = np.array(_rows), 
                                                 _columns = np.array(_columns), 
                                                 _ang_min = _ang_min, 
                                                 _ang_max = _ang_max, 
                                                 _ang_dist = _ang_dist)

    #Create polygon
    return polygon(rotated_rows, rotated_columns)  


# ********************** Create random shape
def get_rnd(cx, cy, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
            _ang_min = 0, _ang_max = 0, _ang_dist = 1, _obj_par1 = None): #RAC25032020
    '''Random Shapes

    :param cx: x coordinate of shape center.
    :param cy: y coordinate of shape center.
    :param _obj_size_min = 10: minimum allowed size of shape [pixels].
    :param _obj_size_max = 10: maximum allowed size of shape [pixels].
    :param _size_dist = 1: distribution of sizes. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _ang_min = 0: minimum rotation angle of shape [degrees].
    :param _ang_max = 0: maximum rotation angle of shape [degrees].
    :param _ang_dist = 1: distribution of rotation angle. Choices are: 
        1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    :param _obj_par1 = None: Which shapes to randomly generate.
        Choices are [1='rectangle',2='ellipse',3='triangle',4='polygon'].
        Leaving value as None will give all shape options as: [1,2,3,4]
    '''

    #Use _obj_par1 to get random shape from list of shapes
    random_shape = obj_opt_par( _obj_type = 5, _obj_par1 = _obj_par1)

    #Generate shape
    if random_shape == 1: #rectangle
        return get_rec(cx, cy,
                       _obj_size_min = _obj_size_min,
                       _obj_size_max = _obj_size_max,
                       _size_dist = _size_dist,
                       _ang_min = _ang_min,
                       _ang_max = _ang_max,
                       _ang_dist = _ang_dist)
    elif random_shape == 2: #ellipse
        return get_elp(cx, cy,
                       _obj_size_min = _obj_size_min,
                       _obj_size_max = _obj_size_max,
                       _size_dist = _size_dist,
                       _ang_min = _ang_min,
                       _ang_max = _ang_max,
                       _ang_dist = _ang_dist)
    elif random_shape == 3: #triangle
        return get_tri(cx, cy,
                       _obj_size_min = _obj_size_min,
                       _obj_size_max = _obj_size_max,
                       _size_dist = _size_dist,
                       _ang_min = _ang_min,
                       _ang_max = _ang_max,
                       _ang_dist = _ang_dist)
    elif random_shape == 4: #polygon
        return get_pol(cx, cy,
                       _obj_size_min = _obj_size_min,
                       _obj_size_max = _obj_size_max,
                       _size_dist = _size_dist,
                       _ang_min = _ang_min,
                       _ang_max = _ang_max,
                       _ang_dist = _ang_dist)
    else: #ellipse
        return get_elp(cx, cy,
                       _obj_size_min = _obj_size_min,
                       _obj_size_max = _obj_size_max,
                       _size_dist = _size_dist,
                       _ang_min = _ang_min,
                       _ang_max = _ang_max,
                       _ang_dist = _ang_dist)


# ********************** Use optional object parameters to get optional shape parameters
def obj_opt_par(_obj_type = 1, _obj_size_min = 10, _obj_size_max = 10, _size_dist = 1,
                _default_axis_ratio = 1.0, _obj_par1 = None, _obj_par2 = None): #RAC25032020
    '''Object Optional Parameters

        :param _obj_type = 1: shape of objects. Choices are: 
            1=rectangle, 2=ellipse, 3=triangle, 4=polygon, 5=random_shapes
        :param _obj_size_min = 10: minimum allowed size of rectangle [pixels]
        :param _obj_size_max = 10: maximum allowed size of rectangle [pixels]
        :param _size_dist = 1: distribution of sizes. Choices are: 
            1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
        :param _obj_par1 = None: 
            rectangle: ratio of width to length (length is the 
                longest dimension). 
                Leaving value as None will give 1.0 (square) as width to length 
                ratio unless _obj_par2 = True is selected to randomize the 
                width to length ratio.
            ellipse: ratio of minor to major semi-axes. 
                Leaving value as None will give 1.0 (circle) as minor to major 
                semi-axes ratio unless _obj_par2 = True is selected to 
                randomize the minor to major semi-axes ratio.
            triangle: ratio of height (y, or "opposite") to width 
                (x, or adjunct). Leaving value as None will give 1.0 (1/1) as 
                height to width ratio unless _obj_par2 = True is selected 
                to randomize the height to width ratio (default is an 
                equilateral triangle).
            regular polygon: number of sides for polygon.
                Leaving value as None will give a hexagon (6 sides)
                unless _obj_par2 = True is selected to randomize the numer of 
                polygon vertices.
            random shapes: Which shapes to randomly generate.
                Choices are [1='rectangle',2='ellipse',3='triangle',4='polygon'].
                Leaving value as None will give all shape options as:
                    [1,2,3,4]
        :param _obj_par2 = None:
            rectangle: value set to True will randomize the width 
                to length ratio.
            ellipse: Value set to True will randomize the minor to major 
                semi-axes ratio.
            triangle: None: value set to True will randomize the height 
                to width ratio.
            regular polygon: Value set to True will randomize the numer of polygon 
                vertices. Max vertices = 12.
    '''

    #Default number of regular polygon vertices
    num_polygon_vertices = 6 #default: hexagon (all 6 sides)
    min_number_sides = 3
    max_number_sides = 12

    #RECTANGLE, ELLIPSE, or TRIANGLE
    #  Use _obj_par1 to get ratio, or create 
    #  random ratio if _obj_par2 = True
    if (_obj_type == 1 or _obj_type == 2 or 
        _obj_type == 3): #rectangle, ellipse, or triangle
        if ((_obj_par1 != None and _obj_par2 == True) or 
            (_obj_par1 == None and _obj_par2 == True)):
            #randomize width to length ratio
            larger = get_dist(_min_value = _obj_size_min, 
                              _max_value = _obj_size_max, 
                              _dist = _size_dist)
            smaller = get_dist(_min_value = _obj_size_min, 
                               _max_value = _obj_size_max, 
                               _dist = _size_dist)             
        elif ((_obj_par1 != None and _obj_par2 == None) or 
              (_obj_par1 != None and _obj_par2 != None) or 
              (_obj_par1 != None and _obj_par2 == False)):
            #use provided smaller to larger ratio
            #  get random int between min and max for larger
            #  use ratio to get associated int smaller
            if _obj_size_min == _obj_size_max:
                larger = _obj_size_max
            else:
                larger = get_dist(_min_value = _obj_size_min, 
                                  _max_value = _obj_size_max, 
                                  _dist = _size_dist)
            smaller = int(larger*_obj_par1)
        else: 
            #define ratio as default
            #  get random int between min and max for length
            #  use default ratio to get assoceated smaller
            if _obj_size_min == _obj_size_max:
                larger = _obj_size_max
            else:
                larger = get_dist(_min_value = _obj_size_min, 
                                  _max_value = _obj_size_max, 
                                  _dist = _size_dist)
            smaller = int(larger*_default_axis_ratio)
        return smaller, larger

    #POLYGON
    #  Use _obj_par1 to get number of vertices, or create 
    #  random number of vertices if _obj_par2 = True
    if _obj_type == 4: #polygon
        #Check if _obj_par2 == True
        if _obj_par2 == True:
            if _obj_size_min == _obj_size_max:
                radius = _obj_size_min
            else:
                radius = get_dist(_min_value = _obj_size_min, 
                                  _max_value = _obj_size_max, 
                                  _dist = _size_dist)
            # use random number of polygon vertices.
            num_polygon_vertices = np.random.randint(min_number_sides, 
                                                     max_number_sides)
        #else, check if _obj_par1 == an intiger
        elif isinstance(_obj_par1, int):
            if _obj_size_min == _obj_size_max:
                radius = _obj_size_min
            else:
                radius = get_dist(_min_value = _obj_size_min, 
                                  _max_value = _obj_size_max, 
                                  _dist = _size_dist)
            # use user provided number of polygon vertices.
            num_polygon_vertices = _obj_par1
        #else, set sides to default
        else:
            if _obj_size_min == _obj_size_max:
                radius = _obj_size_min
            else:
                radius = get_dist(_min_value = _obj_size_min, 
                                  _max_value = _obj_size_max, 
                                  _dist = _size_dist)
            # use hexagon (all 6 sides same size)
            num_polygon_vertices = num_polygon_vertices
        return radius, num_polygon_vertices

    #RANDOM
    #  Select random shape from list of _obj_par1 shapes or default 
    #  (all shapes).
    #  1  = rectangle, 2 = ellipse, 3 = triangle, 4 = polygon 
    if _obj_type == 5: #random shapes
        #Turn _obj_par1 into a list if user failed to enter it as a list
        if not isinstance(_obj_par1, list): 
            if isinstance(_obj_par1, str):
                split_string = _obj_par1.split(",")
                _obj_par1 = [int(i) for i in split_string]
            if isinstance(_obj_par1, tuple):
                _obj_par1 = list(_obj_par1)                

        #Create list of possible shapes (in case user accidentally added 
        #    non-viable options to their list). If "shape" is in _obj_par1, 
        #    append "shapes" with shape number
        shape_options = []
        if (_obj_par1 != None and ((1 in _obj_par1) or
                                   (2 in _obj_par1) or
                                   (3 in _obj_par1) or
                                   (4 in _obj_par1))):
            if 1 in _obj_par1:
                shape_options.append(1) #rectangle
            if 2 in _obj_par1:
                shape_options.append(2) #ellipse
            if 3 in _obj_par1:
                shape_options.append(3) #triangle
            if 4 in _obj_par1:
                shape_options.append(4) #polygon
        else:
            #[1='rectangle',2='ellipse',3='triangle',4='polygon']
            shape_options = [1,2,3,4]
 
        #Get random shape from list
        #[1='rectangle',2='ellipse',3='triangle',4='polygon']
        random_shape = np.random.choice(shape_options)
        if random_shape == 1:
            return(1) #rectangle
        elif random_shape == 2:
            return(2) #ellipse
        elif random_shape == 3:
            return(3) #triangle
        elif random_shape == 4:
            return(4) #polygon
        else:
            return(2) #ellipse


# ********************** Rotate vertices of an object counterclockwise by a given angle
def rot_obj_vert(cx, cy, _rows, _columns, _ang_min = 0, _ang_max = 0, _ang_dist = 1): #RAC25032020
    '''
    Rotate vertices counterclockwise by a given angle (degrees) around center.
    Used for shapes: rectangle, square, polygon

    Rotation Parameters
        :param cx = x: coordinate of object center.
        :param cy = y: coordinate of object center.
        :param _rows: x coordinates of object vertices (array).
        :param _columns: y coordinates of object vertices (array).
        :param _ang_min = 0: minimum rotation angle of shape [degrees].
        :param _ang_max = 0: maximum rotation angle of shape [degrees].
        :param _ang_dist = 1: distribution of rotation angle. Choices are: 
            1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    '''

    try:
        _rows.size == _columns.size
    except:
        print("Can not rotate object. Number of row and column vertex coordinates are not equal")
        
    try:
        _ang_max >= _ang_min
    except:
        print("Maximum angle of rotation is larger than minimum angle of rotation")

    #Rotate vertex coordinates
    rotated_columns = [] #x values
    rotated_rows = [] #y values
    if _ang_max == _ang_min:
        for ii in range(_rows.size):
            _angle = _ang_max
            rotated_columns.append(int(cx + np.cos(np.deg2rad(_angle)) * 
                                       (_columns[ii] - cx) - np.sin(np.deg2rad(_angle)) * 
                                       (_rows[ii] - cy)))
            rotated_rows.append(int(cy + np.sin(np.deg2rad(_angle)) * 
                                    (_columns[ii] - cx) + np.cos(np.deg2rad(_angle)) * 
                                    (_rows[ii] - cy)))
    else:
        _angle = get_dist(_min_value = _ang_min, _max_value = _ang_max, _dist = _ang_dist)
        for ii in range(_rows.size): #Uniform Distribution Parameters
            rotated_columns.append(int(cx + np.cos(np.deg2rad(_angle)) * 
                                       (_columns[ii] - cx) - np.sin(np.deg2rad(_angle)) * 
                                       (_rows[ii] - cy)))
            rotated_rows.append(int(cy + np.sin(np.deg2rad(_angle)) * 
                                    (_columns[ii] - cx) + np.cos(np.deg2rad(_angle)) * 
                                    (_rows[ii] - cy)))

    return rotated_rows, rotated_columns


# ********************** Get value using a given distribution
def get_dist(_min_value, _max_value, _dist = 1): #RAC30032020
    '''Get value using a given distribution.

    Distribution Parameters
        :param _min_value = 0: minimum value.
        :param _max_value = 0: maximum value.
        :param _dist = 1: distribution of value selection. Choices are: 
            1=uniform, 2=normal(Gaussian), 3=Flory–Schulz
    '''

    # load package srwl_uti_dist
    try:
        sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
        from uti_math import get_dist_uni, get_dist_norm, get_dist_schultz
    except:
        print('uti_math functions: get_uni_rand, get_norm_dist, get_schulz_dist failed to load.')

    # chose distribution method
    if _dist == 2: #Normal (Gaussian) Distribution Parameters
        value = get_dist_norm(_min_value, _max_value, _scale=((_max_value-_min_value)/3))
    elif _dist == 3: #Flory–Schulz Distribution Parameters
        value = get_dist_schultz(_min_value, _max_value)
    else: #Uniform Distribution Parameters
        value = get_dist_uni(_min_value, _max_value)  

    return value







