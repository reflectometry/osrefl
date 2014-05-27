from numpy import *

EPSILON = 1e-10

def calculateB(kz, dz, rhoN, rhoM, mx, my, mz):
    """\
    Calculation of reflectivity in magnetic sample in framework that 
    also yields the wavefunction in each layer for DWBA.  
    Arguments:
    
     kz is the incoming momentum along the z direction
     dz is the thickness of each layer (top and bottom thickness ignored)
     rhoN is nuclear scattering length density array
     rhoM is magnetic scattering length density array
     mx, my, mz are components of unit vector along M for layer
    ###################################################################
    #  all of dz, rhoN, rhoM, mx, my, mz should have length equal to  #
    #  the number of layers in the sample including fronting and      #
    #  substrate.                                                     #
    ###################################################################
     AGUIDE is the angle of the z-axis of the sample with respect to the 
     z-axis (quantization) of the neutron in the lab frame.
     AGUIDE = 0 for perpendicular films (= 270 for field along y)

    """
    # sld is array([[sld, thickness, mu], [...], ...])
    # ordered from top (vacuum usually) to bottom (substrate)
    N = len(rhoN) 
    PI4 = complex(pi * 4.0)
    
    if abs(kz) < EPSILON:
        # we're at zero:
        YA = complex(-1.0);
        YB = complex(0.0);
        YC = complex(0.0);
        YD = complex(-1.0);
        return [YA, YB, YC, YD];
    
    
    global B 
    B = zeros((N-1,4,4), dtype='complex') # one matrix per layer
    
    newB = matrix([[1, 0, 0, 0], 
                  [0, 1, 0, 0], 
                  [0, 0, 1, 0], 
                  [0, 0, 0, 1]], dtype='complex')
    
    m_tot = sqrt(mx**2 + my**2 + mz**2)
    m_tot_nz = (m_tot != 0.0)
    
    m_ip = sqrt(mx**2 + my**2) # in-plane projection of m unit vector
    m_ip_nz = (m_ip != 0.0) # mask for nonzero
    
    #mx_ip = zeros_like(mx) 
    #my_ip = zeros_like(my) 
    #mx_ip[m_ip_nz] = mx[m_ip_nz] / m_ip[m_ip_nz]
    #my_ip[m_ip_nz] = my[m_ip_nz] / m_ip[m_ip_nz]
    
    thetaM = arctan2(my, mx)
    
    rhoM_ip = zeros_like(rhoM)
    rhoM_ip[m_tot_nz] = rhoM[m_tot_nz] * m_ip[m_tot_nz] / m_tot[m_tot_nz]
    
    KSQREL = kz**2 + PI4*rhoN[0] # fronting medium removed from effective kz
    
    S1 = sqrt(PI4*(rhoN + rhoM_ip)-KSQREL)
    S3 = sqrt(PI4*(rhoN - rhoM_ip)-KSQREL)
    mu = exp(1j * thetaM)
    
    l=0
    step = 1
    z = complex(0)
    
    if N>1:
        # chi in layer l
        lp = l+step
        chi_l = matrix(\
            [[  1,      1,      0,      0   ],
             [  0,      0,      1,      1   ],
             [S1[l], -S1[l],    0,      0   ],
             [  0,      0,    S3[l], -S3[l] ]])
        
        # inverse chi at layer l+1
        chi_inv_lp = matrix(\
            [[ 1,  1/mu[lp],  1/S1[lp],  1/(mu[lp]*S1[lp]) ],
             [ 1,  1/mu[lp], -1/S1[lp], -1/(mu[lp]*S1[lp]) ],
             [ 1, -1/mu[lp],  1/S3[lp], -1/(mu[lp]*S3[lp]) ],
             [ 1, -1/mu[lp], -1/S3[lp],  1/(mu[lp]*S3[lp]) ]])
        chi_inv_lp *= 0.25
        
        S_l_vec = array([exp( S1[l] * z),
                         exp(-S1[l] * z),
                         exp( S3[l] * z),
                         exp(-S3[l] * z)])
                          
        S_lp_inv_vec = array([exp(-S1[lp] * z),
                              exp( S1[lp] * z),
                              exp(-S3[lp] * z),
                              exp( S3[lp] * z)])
                                                       
        S_lp_inv = matrix(eye(4,4) * S_lp_inv_vec)
        S_l = matrix(eye(4,4) * S_l_vec)
        
        A = S_lp_inv * chi_inv_lp * chi_l * S_l
        #print A
        
        newB = A * newB
        #print newB
        B[l] = newB.copy()
        z += dz[l]
        l = lp
    
    for I in range(1, N-1):
        # chi in layer l
        lp = l+step
        chi_l = matrix(\
            [[1,              1,             1,              1            ],
             [mu[l],          mu[l],        -mu[l],         -mu[l]        ],
             [S1[l],         -S1[l],         S3[l],         -S3[l]        ],
             [mu[l]*S1[l],   -mu[l]*S1[l],  -mu[l]*S3[l],    mu[l]*S3[l]] ])
        
        # inverse chi at layer l+1
        chi_inv_lp = matrix(\
            [[ 1,  1/mu[lp],  1/S1[lp],  1/(mu[lp]*S1[lp]) ],
             [ 1,  1/mu[lp], -1/S1[lp], -1/(mu[lp]*S1[lp]) ],
             [ 1, -1/mu[lp],  1/S3[lp], -1/(mu[lp]*S3[lp]) ],
             [ 1, -1/mu[lp], -1/S3[lp],  1/(mu[lp]*S3[lp]) ]])
        chi_inv_lp *= 0.25
        
        S_l_vec = array([exp( S1[l] * z),
                         exp(-S1[l] * z),
                         exp( S3[l] * z),
                         exp(-S3[l] * z)])
                          
        S_lp_inv_vec = array([exp(-S1[lp] * z),
                              exp( S1[lp] * z),
                              exp(-S3[lp] * z),
                              exp( S3[lp] * z)])
                                                       
        S_lp_inv = matrix(eye(4,4) * S_lp_inv_vec)
        S_l = matrix(eye(4,4) * S_l_vec)
        
        A = S_lp_inv * chi_inv_lp * chi_l * S_l
        #print A
        
        newB = A * newB
        #print newB
        B[l] = newB.copy()
        z += dz[l]
        l = lp
    
    """
    denom = complex(1.0) / ((newB[3,3] * newB[1,1]) - (newB[1,3] * newB[3,1]))
    YA_sam = ((newB[1,3] * newB[3,0]) - (newB[1,0] * newB[3,3])) * denom # r++
    YB_sam = ((newB[1,0] * newB[3,1]) - (newB[3,0] * newB[1,1])) * denom # r+-
    YC_sam = ((newB[1,3] * newB[3,2]) - (newB[1,2] * newB[3,3])) * denom # r-+
    YD_sam = ((newB[1,2] * newB[3,1]) - (newB[3,2] * newB[1,1])) * denom # r--
    
    r_lab = unitary_LAB_SAM_LAB2(matrix([[YA_sam, YB_sam], [YC_sam, YD_sam]]), AGUIDE);
      
    YA_lab = r_lab[0,0];
    YB_lab = r_lab[0,1];
    YC_lab = r_lab[1,0];
    YD_lab = r_lab[1,1];
    
    return YA_lab, YB_lab, YC_lab, YD_lab, B
    """
    return B

def calculateR_sam(B):
    denom = complex(1.0) / ((B[3,3] * B[1,1]) - (B[1,3] * B[3,1]))
    YA_sam = ((B[1,3] * B[3,0]) - (B[1,0] * B[3,3])) * denom # r++
    YB_sam = ((B[1,0] * B[3,1]) - (B[3,0] * B[1,1])) * denom # r+-
    YC_sam = ((B[1,3] * B[3,2]) - (B[1,2] * B[3,3])) * denom # r-+
    YD_sam = ((B[1,2] * B[3,1]) - (B[3,2] * B[1,1])) * denom # r--
    return [YA_sam, YB_sam, YC_sam, YD_sam]

def calculateR_lab(R_sam, AGUIDE):
    r_lab = unitary_LAB_SAM_LAB2(matrix([[R_sam[0], R_sam[1]], [R_sam[2], R_sam[3]]]), AGUIDE);

    YA_lab = r_lab[0,0];
    YB_lab = r_lab[0,1];
    YC_lab = r_lab[1,0];
    YD_lab = r_lab[1,1];
    
    return YA_lab, YB_lab, YC_lab, YD_lab

def calculateRB(kz, dz, rhoN, rhoM, mx, my, mz, AGUIDE):
    B = calculateB(kz, dz, rhoN, rhoM, mx, my, mz)
    R_sam = calculateR_sam(B[-1])
    R_lab = calculateR_lab(R_sam, AGUIDE)
    return R_lab[0], R_lab[1], R_lab[2], R_lab[3], B

def calculateC_sam(C0_sam, B):
    """ take 4x1 matrix (row vector) of initial C
    and get C in each layer """
    layers = B.shape[0]
    C_all = zeros((layers, 4), dtype='complex')
    for l in range(layers):
        C_all[l] = B[l] * C0_sam
    return C_all

def calculateC0_lab(I_plus, I_minus, R_lab):
    # C0 is I+, r+, I-, r-
    # where I+ and I- are supplied by the user and 
    # r+ = (r++ * I+) + (r+- * I-)
    # r- = (r-+ * I+) + (r-- * I-)
    C0_lab = matrix([[I_plus, R_lab[0] * I_plus + R_lab[1] * I_minus, I_minus, R_lab[2] * I_plus + R_lab[3] * I_minus]])
    return C0_lab

# Note: 
# C0_sam = get_Uinv_sam_lab2(AGUIDE) * C0_lab
    

def get_U_sam_lab(AGUIDE):
    C = complex(cos(AGUIDE/2.0*pi/180.))
    IS = 1j * sin(AGUIDE/2.0*pi/180.)
    U = matrix([ [C , IS, 0 , 0 ],
                 [IS, C , 0 , 0 ],
                 [0 , 0 , C , IS],
                 [0 , 0 , IS, C ] ])
    return U

def get_Uinv_sam_lab(AGUIDE):
    C = complex(cos(AGUIDE/2.0*pi/180.))
    NS = 1j * -sin(AGUIDE/2.0*pi/180.)
    Uinv = matrix([ [C , NS, 0 , 0 ],
                    [NS, C , 0 , 0 ],
                    [0 , 0 , C , NS],
                    [0 , 0 , NS, C ] ])
    return Uinv
    
def get_U_sam_lab2(AGUIDE):
    C = complex(cos(AGUIDE/2.0*pi/180.))
    IS = 1j * sin(AGUIDE/2.0*pi/180.)
    U = matrix([ [C , IS],
                 [IS, C ] ])
    return U

def get_Uinv_sam_lab2(AGUIDE):
    C = complex(cos(AGUIDE/2.0*pi/180.))
    NS = 1j * -sin(AGUIDE/2.0*pi/180.)
    Uinv = matrix([ [C , NS],
                    [NS, C ] ])
    return Uinv

def unitary_LAB_SAM_LAB(A, AGUIDE):
    """ perform rotation of coordinate system on one side of matrix
    and perform inverse on the other side """
    U = get_U_sam_lab(AGUIDE)
    Uinv = get_Uinv_sam_lab(AGUIDE)
    CST =  (U * A) * Uinv
    return CST

def unitary_LAB_SAM_LAB2(A, AGUIDE):
    """ perform rotation of coordinate system on one side of matrix
    and perform inverse on the other side """
    U = get_U_sam_lab2(AGUIDE)
    Uinv = get_Uinv_sam_lab2(AGUIDE)
    CST =  (U * A) * Uinv
    return CST

def _test():
    rhoN_mult = 2e-4
    rhoM_mult = 2e-4
    dz_mult = 5
    rhoN = array([0.0, 2.0, 1.0,  2.0, 1.0, 2.0, 1.0, 0.0]) * rhoN_mult
    rhoM = array([0.0, 1.0, 0.0,  1.0, 0.0, 1.0, 0.0, 0.0]) * rhoM_mult
    my =   array([0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 0.0])
    mx =   array([0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0])
    mz =   array([0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0])
    dz =   array([1.0, 1.0, 1.0,  1.0, 1.0, 1.0, 1.0, 1.0]) * dz_mult
    AGUIDE = 0.0 # 90.0 would be along y
    kz_array = linspace(0.01, 1.5, 201)
    rpp = []
    rpm = []
    rmp = []
    rmm = []
    B = []
    for kzi in kz_array:
        result = list(calculateRB(kzi, dz, rhoN, rhoM, mx, my, mz, AGUIDE))
        rpp.append(result[0])
        rpm.append(result[1])
        rmp.append(result[2])
        rmm.append(result[3])
        B.append(result[4])
        
    rpp = array(rpp)
    rpm = array(rpm)
    rmp = array(rmp)
    rmm = array(rmm)
    B = array(B)
    
    from pylab import *
    plot(2*kz_array, abs(rpp)**2, label="r++")
    plot(2*kz_array, abs(rmp)**2, label="r-+")
    plot(2*kz_array, abs(rpm)**2, label="r+-")
    plot(2*kz_array, abs(rmm)**2, label="r--")
    legend()
    show()
    
    return rpp, rpm, rmp, rmm, B

def multiply4x4(A, B):
    C = matrix( [[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0]], dtype='complex' )
    for i in range(4):
        for j in range(4):
            for k in range(4):
                C[i,j] += (A[i,k] * B[k,j])
    return C

if __name__ == '__main__': 
    rpp, rpm, rmp, rmm, B = _test()
    
    
