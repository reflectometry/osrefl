# Copyright (C) 2008 University of Maryland
# All rights reserved.
# See LICENSE.txt for details.
# Author: Christopher Metting

#Starting Date:6/5/2009

from numpy import  greater, less, greater_equal, less_equal, min, max
from numpy import array, size, shape, hstack, vstack, linalg, cross
from numpy import zeros, empty, sum, sort, searchsorted
from numpy import cos,sin, tan, arctan, pi, abs, Inf, degrees


def sphere_point_test(center,r,x,y,z):
    '''
    Overview:
        Determines whether a given point is in a sphere given the point being
    tested and a Sphere object.
        This module is much simpler than the calculation done for the k3d 
    module.
        
        
    Parameters:
    
    center:(float,[3]|angstroms) = The coordinates of the center of the sphere.
    This parameter is in the form (x center,y center, z center)
    
    r: (float|angstroms) = The radius of the sphere.
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    
    
    Note:
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    
    '''
    
    test_results = (((x - center[0])**2 + (y - center[1])**2 +
        (z - center[2])**2)) <=(r**2 )
    return test_results


def parallel_point_test(center,dim,x,y,z):
    '''
    Overview:
        Determines whether a given point is in a parallelapiped given the point
    being tested and the relevant parameters.
    
    
    Parameters:
    
    center:(float,[3]|angstroms) = The coordinates of the center of the 
    parallelapiped. This parameter is in the form (x center,y center, z center)
    
    dim:(float,[3]|angstroms) = The x, y and z dimensions of the parallelapiped
    object.
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    
    
    Note:
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    
    '''
    
    low_lim = (array(center) - (array(dim)/2.0))
    high_lim = (array(center) +(array(dim)/2.0))
    
    height_lim = greater_equal (z,low_lim[2])*less_equal (z,high_lim[2])
    length_lim = greater_equal (y,low_lim[1])*less_equal (y,high_lim[1])
    width_lim = greater_equal (x,low_lim[0])*less_equal (x,high_lim[0])
    
    test_results = height_lim * length_lim * width_lim

    return test_results
   
   

def ellipse_point_test(center,dim,x,y,z):
    '''
    Overview:
        Determines whether a given point is in an ellipse given the point being
    tested and the relevant parameters.
    
    
    Parameters:

    center:(float,[3]|angstroms) = The coordinates of the center of the 
    ellipse. This parameter is in the form (x center,y center, z center)
    
    dim:(float,[3]|angstroms) = The 'a' component, 'b' component and thickness
    of the Ellipse respectively. 'a' is the radius of the Ellipse in the x 
    direction and 'b' is the radius of the ellipsoid in the y direction. 
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    
    
    Notes:
    -To solve this equation more efficiently, the program takes in an array of
    x,y and z so that x[size(x),1,1], y[1,size(y),1], z[1,1,size(z)]. This
    module then solves each part of the test individually and takes the product.
    Only the points where all of the inquires are True will be left as true in
    the test_results array
    
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    
    '''

    low_height_lim = greater_equal (z,(center[2] - dim[2]/2 ))
    up_height_lim =  less_equal (z,(center[2] + dim[2]/2))
    
    xy_test = (((x-center[0])**2)/(dim[0]**2))+(((y-center[1])**2)
        /(dim[1]**2))
    
    in_plane_low_lim = less_equal (0.0,xy_test)
    in_plane_high_lim = greater_equal (1.0,xy_test)
        
    test_results = (low_height_lim * up_height_lim * in_plane_low_lim *
                     in_plane_high_lim)
        
    return test_results

def cone_point_test(center,dim,stub,x,y,z):
    '''
    Overview:
        Determines whether a given point is in an cone given the point being
    tested and the relevant parameters..
    
    
    Parameters:
    
    center:float,[3]|angstroms) = The x, y, and z component of the central
    point of the ellipsoid. In the case that the center is set to 
    [None,None,None] the shape will be put in the bottom corner of the unit cell 
    (the bounding box will start at (0,0,0).
    
    dim:(float,[3]|angstroms) = The x component, y component and thickness
    of the cone respectively. x is the radius of the cone base in the x 
    direction and b is the radius of the cone base in the y direction. 
    
    stub:(float|angstroms) = provides a hard cut-off for the thickness of the
    cone. this allows for the creation of a truncated cone object who side slope
    can be altered by using different z component values while keeping the stub
    parameter fixed.
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    
    
    Notes:
    -To solve this equation more efficiently, the program takes in an array of
    x,y and z so that x[size(x),1,1], y[1,size(y),1], z[1,1,size(z)]. This
    module then solves each part of the test individually and takes the product.
    Only the points where all of the inquires are True will be left as true in
    the test_results array
    
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    
    '''
    
    a_angle = arctan(dim[2]/dim[0])
    b_angle = arctan(dim[2]/dim[1])

    low_height_lim = greater_equal (z,(center[2] - dim[2]/2))
    
    if stub == None:
        up_height_lim =  less_equal (z,(center[2] + dim[2]/2))
    else:
        up_height_lim =  less_equal (z,(center[2] + stub/2))

    xy_test = ((((x-center[0])**2)/((((center[2] + 
           dim[2]/2)-z)/tan(a_angle))**2))+(((y-center[1])**2)/((((center[2] + 
           dim[2]/2)-z)/tan(b_angle))**2)))

    in_plane_low_lim = less_equal (0.0,xy_test)
    in_plane_high_lim = greater_equal (1.0,xy_test)
        
    test_results = (low_height_lim * up_height_lim * in_plane_low_lim * 
                    in_plane_high_lim)
        
    return test_results


def pyrimid_point_test(center,dim,stub,x,y,z):
    '''
    Overview:
        Determines whether a given point is in an pyramid given the point being
    tested and the relevant parameters..
    
    
    Parameters:
    
    center:float,[3]|angstroms) = The x, y, and z component of the central
    point of the ellipsoid. In the case that the center is set to 
    [None,None,None] the shape will be put in the bottom corner of the unit cell 
    (the bounding box will start at (0,0,0).
    
    dim:(float,[3]|angstroms) = The x component, y component and thickness
    of the cone respectively. x is the length of the Pyramid base and y is the 
    width of the Pyramid base. 
    
    stub:(float|angstroms) = provides a hard cut-off for the thickness of the
    Pyramid. this allows for the creation of a trapezoidal object who side slope
    can be altered by using different z component values while keeping the stub
    parameter fixed.
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    
    
    Notes:
    
    -To solve this equation more efficiently, the program takes in an array of
    x,y and z so that x[size(x),1,1], y[1,size(y),1], z[1,1,size(z)]. This
    module then solves each part of the test individually and takes the product.
    Only the points where all of the inquires are True will be left as true in
    the test_results array
    
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    '''
    a_angle = arctan(dim[2]/dim[0])
    b_angle = arctan(dim[2]/dim[1])
    
    low_height_lim = greater_equal(z,(center[2] - dim[2]/2))
    
    if stub == None:
        up_height_lim =  less_equal(z,(center[2] + dim[2]/2))
    else:
        up_height_lim =  less_equal(z,(center[2] + stub/2))
    
    
    
    test_results = (less_equal((center[0] - 
        ((center[2] + dim[2]/2)-z)/tan(a_angle)/2.0),x) *
        less_equal(x,(center[0] + ((center[2] + dim[2]/2)-z)/tan(a_angle)/2.0))*
        less_equal((center[0] - ((center[2] + dim[2]/2)-z)/tan(b_angle)/2.0),y)*
        less_equal(y,(center[0] + ((center[2] + dim[2]/2)-z)/tan(b_angle)/2.0)))

    return test_results

def layer_point_test(thickness,start_point,z):
    '''
    Overview:
        Creates an object that extends the length and width of the unit cell but
    is parameterized in the thickness direction. This is useful for over-layers
    or embed-layers.
    
    
    Parameters:
    
    thickness_value:(float|angstroms) = The thickness of the layer.
    
    start_point:(float|angstroms) = The starting position of layer in the z
    direction. This allows the layer to start at a height anywhere in the unit
    cell. This is useful for over-layers or passivation layers.
    
    z:(float|angstroms) = coordinate for the point being tested in  the z
    direction. Unlike the other tests, this only depends on the z coordinate.
    
    Notes:
    -The API is left intentionally independent of the class structures used in
    sample_prep.py to allow for code resuabilitiy.
    
    '''

    low_height_lim = greater_equal (z,start_point)
    up_height_lim =  less_equal (z,start_point + thickness)

    test_results = low_height_lim * up_height_lim
    return test_results


def ellipsoid_point_test(center, a, b, c, x, y, z):
    '''
    Overview:
        Uses the generic formula for a Ellipsoid feature to create a
    Ellipsoid object. This object can be used to create a sphere by setting
    dim[0] = dim[1] = dim[2]
    
    
    Parameters:
    
    a,b,c:(float|angstroms) = 'a' is the radius of the Ellipsoid 
    in the x direction, 'b' is the radius of the Ellipsoid in the y direction, 
    and 'c' is the radius of the Ellipsoid in the z direction. 
    
    center:float,[3]|angstroms) = The x, y, and z component of the central
    point of the ellipsoid. In the case that the center is set to 
    [None,None,None] the shape will be put in the bottom corner of the unit cell 
    (the bounding box will start at (0,0,0).
    
    x,y,z:(float|angstroms) = coordinates for the point being tested.
    '''
    test_one = less_equal(0,
        (((x - center[0])**2) / (a**2)) +
        (((y - center[1])**2) / (b**2)) +
        (((z - center[2])**2) / (c**2))
        )
    

    test_two = greater_equal(1,
        (((x - center[0])**2) / (a**2)) +
        (((y - center[1])**2) / (b**2)) +
        (((z - center[2])**2) / (c**2)))

    test = test_one * test_two
    return test


    
def K3D_point_test(self,x,y,z):
    '''
    This module makes a list of real space values that the matrix represents
    
    cell_count = the number of cells in a specific dimension
    step_size = the real space value of each voxile
    '''
    
    on_line_check = False

    #************Calculates the points used in the matrix representation*******
    unit_cell = zeros([size(x),size(y),size(z)])
    mag_cell = zeros([size(x),size(y),size(z)])
    x_neg = -2*z[-1]
    x_far = 2*z[-1]
    x_vector = array(x)

    for l in range(size(self.k3d_shapelist)):
        print 'VOXILIZING SHAPE.....',str(l+1)
        
        #**************Loop is done for each point in my data******************
        for k in range(size(z)):
            for j in range(size(y)):
                point_of_inter = zeros([3,self.k3d_shapelist[l].numpoly])
                
                l1,l2 = line_chooser([x[0],y[0],
                    z[0]],y[j],z[k],x[-1],
                    self.k3d_shapelist[l].edges,
                    self.k3d_shapelist[l].vertices,
                    self.k3d_shapelist[l].numpoly)
                
                y[j] = l1[1]
                z[k] = l1[2]
                inclusion = 0
                line_end_neg = [x_neg,y[j],z[k]]
                line_end_far = [x_far,y[j],z[k]]
                crosses_raw = zeros(self.k3d_shapelist[l].numpoly)
                
                
                for ii in range(self.k3d_shapelist[l].numpoly):
                    crosses_raw[ii],point_of_inter[:,ii] = (plane_test(
                                       self.k3d_shapelist[l].vertices,
                                       self.k3d_shapelist[l].edges[ii,:],
                                       line_end_neg,line_end_far))
                  

                point_of_inter_array = array(point_of_inter)
                crosses = point_of_inter_array[:,crosses_raw==True]
                
                if (sum(crosses))>0:
                    new_vector=point_inclusion_test(crosses,x_vector,
                                                    self.SLD_list[l])
                    current_vector = unit_cell[:,j,k]
                    unit_cell[current_vector==0,j,k]= (new_vector
                                                       [current_vector==0])
    #does not handle magnetic K3D feature. This mag_cell is a placeholder
    return unit_cell, mag_cell

    return
#******************Start of Line Test Functions*******************************

def plane_test(vertices,edges,line_end_neg,line_end_far):
    p1 = vertices[edges[0],:]
    p2 = vertices[edges[1],:]
    p3 = vertices[edges[2],:]
    p4 = vertices[edges[3],:]

    [crosses_raw,point_of_inter] = line_through_plane(p1,p2,p3,p4,
                                                      line_end_neg,line_end_far)
    return crosses_raw, point_of_inter



def approx(a,b,tol=1e-5):
    '''
    return true if a = b within tolerance
    '''
    return abs(a-b) < tol

def approx_between(a,b,c,tol=1e-5):

    '''
    returns true if b is approximately greater than or equal to a but
    approximately less than or equal to c
    '''

    if b >= (a - tol) and b <= (c + tol):
        return True
    else:
        return False
    
def vector_calc(point_one,point_two):
    ans = point_one - point_two
    return ans

def colinear(p1,p2,p3,p4):
    '''
    This module determines whether 4 points are colinear
    p1,p2,p3,p4 = points to be tested for linearity
        
    '''
    check = empty(2)
    p21 = p2 - p1
    p31 = p3 - p1
    p41 = p4 - p1
    
    num = sum(abs(cross((p21),(p31))))
    denom = sum(abs((p21)))
    
    if (approx(num,0)) and (approx(denom,0)):
        d = 0.0
    else:
        d = num/denom

    if (approx(d,0) == True):
        check[0] = True
    else:
        check[0] = False
    
    num = float(sum(abs(cross((p21),(p41)))))
    denom = float(sum(abs((p21))))
    
    if (approx(num,0)) and (approx(denom,0)):
        d = 0.0
    else:
        d = num/denom

    if (approx(d,0) == True):
        check[1] = True
    else:
        check[1] = False
    
    if (approx(check[0],1)) and (approx(check[1],1)):
        colinear = True
    else:
        colinear = False

    return colinear
    
def line_inter(p1,p2,p3,p4):
    '''
    This module determines whether a ray and a segment intersect
    
    p1,p2 defines a line segment from either a polyhedron face or a polygon side
    p3,p4 defines the ray from the test point(p3) to an outer sphere (4)
    '''
    A = empty((2,2))
    b = empty((2,1))
    x = empty ((2,1))
    
    p1_test = array([p1[0],p1[1]])
    p2_test = array([p2[0],p2[1]])
    p3_test = array([p3[0],p3[1]])
    p4_test = array([p4[0],p4[1]])

    #poorly written helps with the plane check
    
    if size(p1)==3:
        p1_test_2 = array([p1[1],p1[2]])
        p2_test_2 = array([p2[1],p2[2]])
        p3_test_2 = array([p3[1],p3[2]])
        p4_test_2 = array([p4[1],p4[2]])

    if (any(p1_test != p2_test)) and (any(p3_test != p4_test)):
        #This is the case if there are only two dimensions being test or if the
        #two lines are in the xy-plane
        for f in range (2):
            A[f,0] =  -(p4[f] - p3[f])
            A[f,1] = p2[f]- p1[f]
            b[f] = p3[f]-p1[f]

        try:
            x = linalg.solve(A,b)
            #x[0] is the position on the ray of the intersection
            #x[1] is the position on the line segment
            crossing = (0 <= x[0]) and (0 <= x[1] <= 1)
            if crossing and len(p1) == 3:
                crossing = approx((p3[2] - p1[2]),((p3[2]+p4[2])*x[0] - 
                                                   (p1[2]+p2[2])*x[1]))
        except linalg.LinAlgError:
            crossing = False
              
    elif(any(p1_test_2 != p2_test_2)) and (any(p3_test_2 != p4_test_2)):
        #This is the case if the two lines are in the yz plane
        
        for f in range (2):
            A[f,0] =  -(p4[f+1] - p3[f+1])
            A[f,1] = p2[f+1]- p1[f+1]
            b[f] = p3[f+1]-p1[f+1]
        try:
            x = linalg.solve(A,b)
            #x[0] is the position on the ray of the intersection
            #x[1] is the position on the line segment
            crossing = (0<= x[0]) and (0<= x[1] <= 1)
            if crossing and len(p1) == 3:
                crossing = approx((p3[0] - p1[0]),((p2[0]-p1[0])*x[1]-
                                                   (p4[0]-p3[0])*x[0]))
        except linalg.LinAlgError:
            crossing = False
    else:
        #This is the case if the two lines are in the xz plane
         
        for f in range (2):
            A[f,0] =  -(p4[f] - p3[f+1])
            A[f,1] = p2[f]- p1[f+1]
            b[f] = p3[f]-p1[f+1]
        try:
            x = linalg.solve(A,b)
            #x[0] is the position on the ray of the intersection
            #x[1] is the position on the line segment
            crossing = (-1e-14<= x[0]) and (-1e-14<= x[1] <= 1+1e-14)
            if crossing and len(p1) == 3:
                crossing = approx((p3[0] - p1[0]),((p3[0]+p4[0])*x[0] - 
                                                   (p1[0]+p2[0])*x[1]))
        except linalg.LinAlgError:
            crossing = False
    return crossing

def point_on_line(p1,p2,point):
    '''
        This module determines where a point falls on a line segment.
        
        p1, p2 = points which define the line being tested.
        point = the point which will be determined whether or not falls on the 
        line.
        
    '''

    high = empty(3)
    low = empty(3)
    p21 = p2- p1
    p1_point = p1- point
    
    for i in range (3):
        if p21[i] > 0:
            high[i] = p2[i]
            low [i] = p1[i]
        else:
            high[i] = p1[i]
            low [i] = p2[i]
    
    num = float(sum(abs(cross((p21),(p1_point)))))
    denom = float(sum(abs((p21))))

    if approx(num,0) and approx(denom,0):
        d = 0
    else:
        d = num/denom

    if ((approx(d,0) == True) and 
        (approx_between(low[0],point[0],high[0])) and 
        approx_between(low[1],point[1],high[1]) and 
        approx_between(low[2],point[2],high[2])):
        on_line = True
    else:
        on_line = False
    return on_line

def point_on_poly(p1,p2,p3,p4,point):
    '''
    This module determines if a point falls on a plane
    
    p1,p2,p3,p4 = four points which make up a plane
    V = a vector normal to the plane
    point = point being tested to determine if it lies on the plane
    '''
    v = (cross(vector_calc(p1,p2),vector_calc(p3,p2)))
    D = -(p4[0]*v[0]) - (p4[1]*v[1]) - (p4[2]*v[2])
    angle = array([pi,pi/2,pi/4])
    ang = 0
    check = empty(4)
    circle_point = empty(3)
    point_test = True
    
    plane_eq = hstack([v,D])
    distance = (plane_eq[0]*point[0] + plane_eq[1]*point[1] + 
                plane_eq[2]*point[2] + plane_eq[3])**2/(plane_eq[0]**2 + 
                                                        plane_eq[1]**2 + 
                                                        plane_eq[2]**2)
    
    if approx(distance,0):
        if not colinear(array([p1[0],p1[1],0]), array([p2[0],p2[1],0]), 
                        array([p3[0],p3[1],0]) , array([p4[0],p4[1],0])):
            proj_axis_one = 0
            proj_axis_two = 1
            zero_axis = 2
            
            
        elif not colinear(array([p1[0],0,p1[2]]), array([p2[0],0,p2[2]]), 
                          array([p3[0],0,p3[2]]) , array([p4[0],0,p4[2]])):
            proj_axis_one = 0
            proj_axis_two = 2
            zero_axis = 1
            
        elif not colinear( array([0,p1[1],p1[2]]), array([0,p2[1],p2[2]]), 
                           array([0,p3[1],p3[2]]) , array([0,p4[1],p4[2]])):
            proj_axis_one = 1
            proj_axis_two = 2
            zero_axis = 0
            
        
        poly = vstack([[p1],[p2],[p3],[p4]])
        
        if   (max(poly[:,proj_axis_one]) > point[proj_axis_one] > 
              min(poly[:,proj_axis_one])):
            
            if (max(poly[:,proj_axis_two]) > point[proj_axis_two] >
                 min(poly[:,proj_axis_two])):
                
                r = max(poly) + 1
                
                while point_test == True:

                    circle_point[proj_axis_one] =r*cos(angle[ang])
                    circle_point[proj_axis_two] = r*sin(angle[ang])
                    circle_point[zero_axis] = 0
                    
                    # Make a private copy for projection onto the zero_axis
                    point,p1,p2,p3,p4 = [array(m) for m in (point,p1,p2,p3,p4)]
                    point[zero_axis] = 0
                    p1[zero_axis] = 0
                    p2[zero_axis] = 0
                    p3[zero_axis] = 0
                    p4[zero_axis] = 0
                    
                    x = array([point[proj_axis_one],point[proj_axis_two],0])

                    check[0] = point_on_line(point,circle_point,p1)
                    check[1] = point_on_line(point,circle_point,p2)
                    check[2] = point_on_line(point,circle_point,p3)
                    check[3] = point_on_line(point,circle_point,p4)
                    if check[0]  == check[1] == check[2]== check[3] == False:
                        point_test = False
                    else:
                        ang =+ 1

                check[0] = line_inter([p1[proj_axis_one],p1[proj_axis_two]],
                                      [p2[proj_axis_one],p2[proj_axis_two]],
                                      [point[proj_axis_one],
                                       point[proj_axis_two]],
                                       [circle_point[proj_axis_one],
                                        circle_point[proj_axis_two]])
                
                check[1] = line_inter([p2[proj_axis_one],p2[proj_axis_two]],
                                      [p3[proj_axis_one],p3[proj_axis_two]],
                                      [point[proj_axis_one],
                                       point[proj_axis_two]],
                                       [circle_point[proj_axis_one],
                                        circle_point[proj_axis_two]])
                
                check[2] = line_inter([p3[proj_axis_one],p3[proj_axis_two]],
                                      [p4[proj_axis_one],p4[proj_axis_two]],
                                      [point[proj_axis_one],
                                       point[proj_axis_two]],
                                       [circle_point[proj_axis_one],
                                        circle_point[proj_axis_two]])
                
                check[3] = line_inter([p4[proj_axis_one],p4[proj_axis_two]],
                                      [p1[proj_axis_one],p1[proj_axis_two]],
                                      [point[proj_axis_one],
                                       point[proj_axis_two]],
                                       [circle_point[proj_axis_one],
                                        circle_point[proj_axis_two]])
                
                total_lines = sum(check)
                if total_lines == 0 or total_lines == 2:
                    on_plane = False
                else:
                    on_plane = True
            else:
                on_plane = False
        else:
            on_plane = False
    else:
        on_plane = False
    
    return on_plane

def line_through_plane(p1,p2,p3,p4,l1,l2): 
    
    '''
        This module if a line which is known not to cross a
        3D shape face at a line or a point and does not lie
        along a line on the face crosses the face
        
        p1, p2, p3, p4 = points which define the shape face in three dimensions
        l1,l2 = the line being tested for crossing
    '''
    f = 0
    point_inter = empty((3))
    A = empty((3,3))
    b = empty((3,1))
    x = empty ((3,1))
    for f in range (3):
        b[f] = l1[f] - p1[f]
        
        A[f,0] =  l1[f] - l2[f]
        A[f,1] =  p2[f] - p1[f]
        A[f,2] =  p3[f] - p1[f]
    
    
    try:
        x = linalg.solve(A,b)
        #x[0] is the position on the ray of the intersection

        ray_check = (-1e-5<= x[0])
        if ray_check:
            for f in range (3):
                point_inter[f] = l1[f]+(l2[f]-l1[f]) * x[0]
            crossing = point_on_poly(p1,p2,p3,p4,point_inter)
        else:
            crossing = False
    except linalg.LinAlgError:
        crossing = False
        
    return crossing,point_inter

def circle_test(x,y,z,sphere_r):
    '''
    This module is the final module which can determine whether a point falls
    inside a sphere which encompasses a feature
        
     x,y,z = the values of the point being tested
     sphere_r = a point which, from zero, creates a sphere that encompases
     the whole feature
    '''

    half_point = [(sphere_r[0]/2),(sphere_r[1]/2),(sphere_r[2]/2)]
    
    if (((x - half_point[0])**2 + 
         (y - half_point[1])**2 + 
         (z - half_point[2])**2) > 
         (half_point[0]**2 ) + (half_point[1]**2) + (half_point[2]**2)):
      
        inside = False
    else:
        inside = True
    return inside  
    
def line_chooser(step_list,y_value,z_value,Dx,poly_array,
                 point_array,num_polygons):
    '''
    This module shortens computation time by calculating a single line that can
    be utilized for a full row of voxels
        
    x_value,y_value,z_value = points which define the shape face in three
    dimensions
    '''
    fraction_tracker = 1
    y_start = y_value
    z_start = z_value
    x_neg = -2*Dx
    x_far = 2*Dx
    inner_point = array([x_neg,y_start,z_start])
    outer_point = array([x_far,y_start,z_start])
    point_determine = 0
    for ii in range(num_polygons): #LOOP FOR EACH PLANE    
        for iii in range (4): #LOOP FOR EACH LINE THAT MAKES UP EACH PLANE
            poly_point_one = point_array[poly_array[ii,iii],:]
            if ((iii+1) < (4)):
                poly_point_two = point_array[poly_array[ii,iii+1],:]
            else:
                poly_point_two = point_array[poly_array[ii,0]]                            
            #   test vector crosses a line or a node: chooses new vector
                
            intersect_check = line_inter(poly_point_one,poly_point_two,
                                         inner_point,outer_point)
            if (intersect_check  == True):
                fraction_tracker += 1
                y_start = y_value - (step_list[1]/fraction_tracker)
                z_start = z_value - (step_list[2]/fraction_tracker)
                inner_point = array([x_neg,y_start,z_start])
                outer_point = array([x_far,y_start,z_start])
                intersect_check  = False

                ii = 0
                iii = 4
    return inner_point,outer_point

def shape_builder(x_value,y_value,z_value,poly_array,point_array,
                  num_polygons,far_point):
    '''
    This module is the final module which can determine whether a point 
    falls inside a shape
        
    x_value,y_value,z_value = points which define the shape face in three 
    dimensions
    poly_array = array of numbers that indicate the four points in point_array 
    that make up a face
    point_array = array of real points 
    
    '''
    sphere_point = zeros(3)
    poly_point_one = zeros(3)
    poly_point_two = zeros(3)
    check = False
    on_line_check = False
    intersect_check = False
    test_result = False
    
    face_count = 0
    d = 0
    e = 0
    
#**Tests possible complications caused by the point falling on a specific feature
    test_point = [x_value,y_value,z_value]
    for ii in range(num_polygons): #LOOP FOR EACH PLANE
        #Test to see if the point falls on the polygon face
        if point_on_poly(point_array[poly_array[ii,0],:],
                         point_array[poly_array[ii,1],:],
                         point_array[poly_array[ii,2],:],
                         point_array[poly_array[ii,3],:],
                         test_point) == True:
            on_line_check = True
            break
        #Tests to see of the point falls on any of the lines that make up the polygon
        elif point_on_line(point_array[poly_array[ii,0],:],
                           point_array[poly_array[ii,1],:],
                           test_point):
            on_line_check = True
            break
        elif point_on_line(point_array[poly_array[ii,1],:],
                           point_array[poly_array[ii,2],:],
                           test_point):
            on_line_check = True
            break
        elif point_on_line(point_array[poly_array[ii,2],:],
                           point_array[poly_array[ii,3],:],
                           test_point):
            on_line_check = True
            break
        elif point_on_line(point_array[poly_array[ii,3],:],
                           point_array[poly_array[ii,0],:],
                           test_point):
            on_line_check = True
            break
        else:
            on_line_check = False
            
    if (on_line_check == False):
        for ii in range(num_polygons): #LOOP FOR EACH PLANE
            if line_through_plane(point_array[poly_array[ii,0],:],
                                  point_array[poly_array[ii,1],:],
                                  point_array[poly_array[ii,2],:],
                                  point_array[poly_array[ii,3],:],
                                  test_point,far_point) == True:
                face_count = face_count + 1 
        if (face_count == 1) or (face_count == 3) or (face_count == 5) or (face_count == 7):
            test_result = True
        else:
            test_result = False        
    else:
        test_result = True
        
    return test_result

def point_inclusion_test(poly_crossed,vector,SLD):
    '''
    This module determines if a point falls in a feature by a less 
    expensive method involving the calculation of the point of line/poly 
    intersections and determining of the point falls between these points inside
    or outside of the polyhedron.
        
    poly_crossed = is a list of which polygons were and were not crossed.
    point_of_inter = the point at which the line crosses the polygon for those
    that do pass through it.
    test_point = the point on which the inclusion test is being done.
    
    '''
    sorted_poly = sort(poly_crossed[0,:])
    
    
    locator = sorted_poly.searchsorted(vector)%2
    SLD_profile = zeros(size(locator))
    SLD_profile[locator==1] = SLD
    return SLD_profile

def test():
    '''
    this test contains an array of assertion statements to ensure the shapes
    are being properly treated.
    '''
    from sample_prep import Parallelapiped, Sphere, Ellipse, Cone
    
    first_cone = Cone(SLD = 9.4e-5,dim = [5.0,5.0,5.0],stub = None)
    assert cone_point_test(first_cone.center, first_cone.dim, first_cone.stub,2.5,2.5,2.0) == True, 'Cone calculation broke'
    
    test_sphere = Sphere(9.5e-4,10.0)
    assert sphere_point_test(test_sphere.center,test_sphere.r,10,10,10) == True, 'Sphere: point in sphere'
    assert sphere_point_test(test_sphere.center,test_sphere.r,1,1,1) == False, 'Sphere: point out of sphere'    
    
    test_parallel = Parallelapiped(9.87e-6,[10,10,10])
    assert parallel_point_test(test_parallel.center, test_parallel.dim,5,5,5) == True, 'parallel: Inside box'
    assert parallel_point_test(test_parallel.center, test_parallel.dim,11,5,5) == False, 'parallel: Outside Box'
    assert parallel_point_test(test_parallel.center, test_parallel.dim,11,11,11) == False, 'parallel: all axes Outside Box'
    
    test_ellipse = Ellipse(9.8e-6, [6.0,8.0,10.0])
    assert ellipse_point_test(test_ellipse.center,test_ellipse.dim,3,8,5) == True, 'ellipse: yaxis limit '
    assert ellipse_point_test(test_ellipse.center,test_ellipse.dim,3,8,11) == False, 'ellipse: z test '
    assert ellipse_point_test(test_ellipse.center,test_ellipse.dim,0,4,5)== True, 'ellipse: xaxis limit'
    assert ellipse_point_test(test_ellipse.center,test_ellipse.dim,5,4,5)== True, 'ellipse: random'

if __name__=="__main__":test()