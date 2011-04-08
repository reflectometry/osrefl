# -*- coding: utf-8 -*-
# a module for calculating the wavefunction inside layers of a film
# likely to be used in DWBA as the input function

# first let's try calculating the reflectivity...

from numpy import zeros, sqrt, pi, eye, array, dot, cos, sin, ndarray, ones, empty, exp, sum, indices

def calc_r(kz_in, array_of_sld):
  # array_of_sld is array of SLD, thickness, mu for each layer
  # (starting with the top layer, ending with the substrate)
  
  # generate matrix elements
  layer_num_total = array_of_sld.shape[0]
  M_l = zeros((layer_num_total,2,2), dtype='complex')
  M = eye(2)
  k0z = kz_in.astype('complex')
  for layer_num in range(layer_num_total):
    SLD,thickness,mu = array_of_sld[layer_num]
    nz = sqrt( 1 - 4 * pi * SLD / k0z**2 )
    kz = nz * k0z
    M_l[layer_num] = array([[cos(kz * thickness), 1/kz * sin(kz * thickness)],[-kz * sin(kz * thickness), cos(kz * thickness)]])
    M = dot(M_l[layer_num], M) # cumulative product moving right to left along M_j-1*M_j-2*...*M_1
   
  # assume substrate is vacuum - k_f = k_in
  kfz = k0z
  r = (-1j * kfz * M[0,0] + k0z * kfz * M[0,1] + M[1,0] + 1j * k0z * M[1,1]) / (-M[1,0] + 1j * k0z * M[1,1] + 1j * kfz * M[0,0] + k0z * kfz * M[0,1])
  return r

def calc_r_array(kz_in_array, array_of_sld):
  r_array = []
  for kz_in in kz_in_array:
    r_array.append(calc_r(kz_in, array_of_sld))

  return array(r_array)

def calc_r_born(kz_in, array_of_sld):
  z,sld_z = sld_discretize(array_of_sld)
  z.shape = (z.shape[0], 1)
  sld_z.shape = (sld_z.shape[0], 1)
  kz = kz_in.copy()
  kz.shape = (1, kz.shape[0])
  ft = exp(1j * 2.0 * kz * z) * sld_z
  intens = sum(ft, axis = 0)
  
  return intens, z[:,0], sld_z[:,0]
  
def calc_r_born_prediscretized(kz_in, z, sld_z):
  z.shape = (z.shape[0], 1)
  sld_z.shape = (sld_z.shape[0], 1)
  kz = kz_in.copy()
  kz.shape = (1, kz.shape[0])
  ft = exp(1j * 2.0 * kz * z) * sld_z
  intens = sum(ft, axis = 0)
  return intens
  
def calc_r_born_2d_discrete(Qz, Qy, dz, dy, sld_z_y):
  sld_z_y.shape = (sld_z_y.shape[0], sld_z_y.shape[1])
  z,y = indices(sld_z_y.shape, dtype=float)
  z *= dz
  y *= dy
  ft = sld_z_y * exp(1j * Qz * z) * exp(1j * Qy * y)
  ft_sum = sum(sum(ft))
  qy_mult = ( 1. - exp( -1j*Qy*dy ) )/(1j * Qy) if Qy!=0 else dy
  qz_mult = ( 1. - exp( -1j*Qz*dz ) )/(1j * Qz) if Qz!=0 else dz
  ft_sum *= qy_mult * qz_mult
  return ft_sum
  

class neutron_wavefunction:
  from numpy import zeros, sqrt, pi, eye, array, dot, cos, sin, exp, complex, complex128
  """object contains wavefunction psi (complex) for a given scattering-length-density,
  and incoming k0z 
  inputs: k0z, kfz, array_of_sld
  returns nothing
  member functions include:
  calling self (neutron_wavefunction_object(z)) which returns depth-dependent wavefunction
  r, which returns total reflectivity
  t, which returns transmission
  c[i] and d[i], which return coefficients of wavefunction in each layer
  """
  
  def __init__(self, kz_in, array_of_sld):
    from numpy import exp
    """initialize the class with necessary input variables.
    kz_in is wavevector for incoming neutrons (in vacuum)
    array of sld is array of [SLD, thickness, mu] for each layer
    (starting with incident medium, ending with the substrate)
    in the incident medium and substrate, only the SLD is used;
    mu and thickness are discarded"""
      
    self.kz_in = kz_in
    self.array_of_sld = array_of_sld
    
    layer_num_total = array_of_sld.shape[0]
    self.layer_num_total = layer_num_total
    self.total_thickness = sum(array_of_sld[1:-1,1])
    
    M_l = zeros((layer_num_total,2,2), dtype=complex)
    M = eye(2, dtype=complex)
    nz = zeros((layer_num_total), dtype=complex)
        
    SLD_incident = array_of_sld[0,0]
   
   #this next correction is only valid for side-entry of fronting layer:  
   # it assumes that kz_in is equal to kz_M (in laboratory frame)
   # and we create a fake k0z that multiplies away correctly
   # (kz[0] = k0z * nz[0] = kz_in)
    k0z = sqrt(complex(kz_in**2 + 4 * pi * SLD_incident))
    
    nz = sqrt( complex(1) - 4 * pi * array_of_sld[:,0] / k0z**2 )
    # calculate the M matrix
    for layer_num in range(1, layer_num_total-1):
      #leaving off the incident medium and substrate from sum
      SLD,thickness,mu = array_of_sld[layer_num]

      nz[layer_num] = sqrt(complex( 1 - 4 * pi * SLD/ k0z**2 ))
      kz = nz[layer_num] * kz_in
      n = nz[layer_num]
      M_l[layer_num] = array([[cos(kz * thickness), 1/n * sin(kz * thickness)],[-n * sin(kz * thickness), cos(kz * thickness)]])
      M = dot(M_l[layer_num], M) # cumulative product moving right to left along M_j-1*M_j-2*...*M_1
      
    self.nz = nz
    self.M = M
    self.M_l = M_l
    
    # calculate r
    SLD_substrate = array_of_sld[-1,0] # take the sld from the last element of array_of_sld
    
    # old r: don't know where I got it
    # r = (-1j * kfz * M[0,0] + k0z * kfz * M[0,1] + M[1,0] + 1j * k0z * M[1,1]) / (-M[1,0] + 1j * k0z * M[1,1] + 1j * kfz * M[0,0] + k0z * kfz * M[0,1])
    # r that I calculated myself for M above:
    r = (M[0,0] + 1j * nz[0] * M[0,1] + 1/(1j * nz[-1])*( -M[1,0] - 1j * nz[0] * M[1,1])) / (-M[0,0] + 1j * nz[0] * M[0,1] + 1/(1j * nz[-1])*( M[1,0] - 1j * nz[0] * M[1,1])) 
    self.r = r # make visible to the outside world
    if nz[-1].real == 0:
      self.t = complex(0)
    else:
      self.t = 1.0 + self.r
    
    self.kz_transmitted = nz[-1] * k0z
    
    # calculate c, d for each layer
    c = zeros((layer_num_total), dtype=complex)
    d = zeros((layer_num_total), dtype=complex)
    psi_l = zeros((layer_num_total), dtype=complex)
    psi_prime_l = zeros((layer_num_total), dtype=complex)
    c[0] = 1 # incident beam has intensity 1
    d[0] = r # reflected beam has intensity |r|**2
    
    psi_l[0] = 1 + r
    psi_prime_l[0] = 1j * nz[0] * (1 - r)
    z_interface = 0.
    p = complex(1 + r) #psi 
    pp = complex(1j * nz[0] * (1 - r)) #psi prime
    for l in range(1,layer_num_total):
      ## this algorithm works all the way into the substrate
      SLD,thickness,mu = array_of_sld[l]
      #print l, z_interface
      c[l] = 0.5 * ( p + ( pp / (1j * nz[l]) ) ) * exp(-1j * self.kz_in * nz[l] * z_interface)
      d[l] = 0.5 * ( p - pp/(1j * nz[l]) ) * exp(1j * self.kz_in * nz[l] * z_interface)
      z_interface += thickness

      p,pp = dot(M_l[l], array([[p],[pp]]))
      p = p[0]
      pp = pp[0]
      
    # fill final c,d
    self.c = c
    self.d = d
    self.d[-1] = 0.0
       
    return None

  def partial_layer_r(self, start_layer, end_layer):
    M = eye(2, dtype=complex)
    if start_layer > end_layer:
      M_l = flipud(self.M_l[end_layer:start_layer+1])
      k0z = self.nz[end_layer] * self.kz_in
      kfz = self.nz[start_layer+1] * self.kz_in
    else:
      M_l = self.M_l[start_layer:end_layer+1]
      k0z = self.nz[start_layer] * self.kz_in
      kfz = self.nz[end_layer+1] * self.kz_in
    for m in M_l:
      M = dot(m, M)
    
    r = (-1j * kfz * M[0,0] + k0z * kfz * M[0,1] + M[1,0] + 1j * k0z * M[1,1]) / (-M[1,0] + 1j * k0z * M[1,1] + 1j * kfz * M[0,0] + k0z * kfz * M[0,1])
    return r

  def __call__(self, z):
    from numpy import exp
    if z < 0: 
      # we're in the incident medium
      l = 0
      return self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) + self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z)
    elif z >= self.total_thickness:
      # we're in the substrate
      l = self.layer_num_total - 1
      print 'l = ', l
      print 1j * self.kz_in * self.nz[1] * z
      return self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) + self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z)
    elif (z >= 0 and z < self.total_thickness):
      z_next = 0.
      for layer_num in range(1, self.layer_num_total - 1):
        SLD, thickness, mu = self.array_of_sld[layer_num]
        z_next += thickness
        if z < z_next:
          l = layer_num
          return self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) + self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z)

  def prime(self, z):
    from numpy import exp
    if z < 0: 
      # we're in the incident medium
      l = 0
      return 1j * self.nz[l] * (self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) - self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z))
    elif z >= self.total_thickness:
      # we're in the substrate
      l = self.layer_num_total - 1
      return 1j * self.nz[l] * (self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) - self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z))
    elif (z >= 0 and z < self.total_thickness):
      z_next = 0.
      for layer_num in range(1, self.layer_num_total - 1):
        SLD, thickness, mu = self.array_of_sld[layer_num]
        z_next += thickness
        if z < z_next:
          l = layer_num
          return 1j * self.nz[l] * (self.c[l] * exp(1j * self.kz_in * self.nz[l] * z) - self.d[l] * exp(-1j * self.kz_in * self.nz[l] * z))

class escaping_neutron:
  from numpy import flipud
  """if a neutron wavefunction happens to pop up in the middle of a stack, what to do?
  how much of this neutron will escape out the top, and how much out the bottom?"""
  def __init__(self, kz_escapee, array_of_sld, layer_escaping_from):
    self.kz_in = kz_escapee
    self.array_of_sld = array_of_sld
    self.l = layer_escaping_from
    thick = self.array_of_sld[self.l][1]
    
    layer_num_total = array_of_sld.shape[0]
    self.layer_num_total = layer_num_total
    self.total_thickness = sum(array_of_sld[1:-1,1])

    # break SLD up into two blocks: one above and one below
    self.sld_below = array_of_sld[self.l:]
    self.sld_above = array_of_sld[:self.l + 1]
    
    if self.kz_in > 0:
      c = 1.0  / (1.0 - r_a * r_b * exp(2.0 * 1j * self.kz_in * thick))
      d = r_a * exp(2.0 * 1j * self.kz_in * thick) / (1.0 - r_a * r_b * exp(2.0 * 1j * self.kz_in * thick))
    else: # self.kz_in <= 0
      c = r_b / (1.0 - r_a * r_b * exp(2.0 * 1j * self.kz_in * thick))
      d = 1.0 / (1.0 - r_a * r_b * exp(2.0 * 1j * self.kz_in * thick))
    
    M_l = zeros((layer_num_total,2,2), dtype=complex)
    
    Ma = eye(2, dtype=complex)
    Mb = eye(2, dtype=complex)
    
    nz = zeros((layer_num_total), dtype=complex)
        
    SLD_incident = array_of_sld[0,0]
    nz[0] = 1.0
    n0z = nz[0]
    #k0z = kz_in * n0z # calculate kz in the incident medium (kz_in is in vacuum)
    
    k0z = kz_in
    


def sld_discretize(sld, min_feature_pixels = 10):
  """taking array of SLD that is in form [SLD, thickness, mu]
  and converting to discrete SLD as function of depth
  smallest layer is broken up into 10 pieces"""
  total_thickness = sum(sld[:,1])
  from numpy import cumsum, concatenate, linspace, repeat, arange
  from pylab import plot
  z_in = repeat(cumsum(sld[:,1]),2)
  z_in = concatenate(([0], z_in[:-1])) #
  sld_in = repeat(sld[:,0], 2)
  #plot(z_in,sld_in)
  
  min_feature_size = sld[:-1,1].min() # not including substrate
  
  z_out = arange(0,total_thickness, min_feature_size / 10.0)
  from scipy import interpolate
  f = interpolate.interp1d(z_in, sld_in)
  sld_out = f(z_out)
  return z_out, sld_out
  

        
def bottom_scatterer():
  from numpy import arange, array, zeros, sin, arcsin, arctan
  from numpy.linalg import inv
  qx = arange(-0.003, 0.003, 0.00005)
  qx.shape = (qx.shape[0],1)
  qz = arange(0,0.14,0.001)
  qz.shape = (qz.shape[0],1)
  sld = array([[0,0,0],[250,4.5e-6,0],[0,1.027e-6,0]])
  k0 = 2 * pi / 5.0
  p = zeros((qx.shape[0],qz.shape[0]), dtype='complex128')
  pp = zeros((qx.shape[0],qz.shape[0]), dtype='complex128')
  q = sqrt(qx*qx + (qz * qz).T)
  tilt = arctan( qx * (1.0/qz).T)
  A4 = 2. * arcsin(q / (2 * k0))
  th_in = A4/2.0 - tilt
  th_out = A4/2.0 + tilt
  ki = k0 * sin(th_in)
  kf = k0 * sin(th_out)
  for i in range(qx.shape[0]):
    for j in range(qz.shape[0]):
      psi_in = neutron_wavefunction(ki[i,j],sld)
      psi_out = neutron_wavefunction(kf[i,j],sld)
      Mi = psi_in.M
      Mf = psi_out.M
      output = dot( inv(Mf), array([psi_in(250),psi_in.prime(250)]) )
      #output = array([psi_out(0), psi_out.prime(0)])
      #print output
      p[i,j] = output[0]
      pp[i,j] = output[1]

  return p,pp,ki,kf

# from rebin import rebinned_data

# class DWBA_check(rebinned_data):
  