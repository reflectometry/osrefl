# Holder Class for Sample Objects
# Wayne Treible - June 4, 2012

from .sample_prep import *

class Sample(object):
    
    def empty(self):       
        return
    
    def getScene(self):
        
        self.scenelist = []
        if(self.MaterialC != None):
            self.scenelist.append(self.MaterialC)
        self.scenelist.append(self.MaterialB)
        for i in range(len(self.MaterialA)):
            self.scenelist.append(self.MaterialA[i])
        
        return Scene(self.scenelist)   
        
class CylinderSample(Sample):
    
    def __init__(self, shell_dim, core_dim, 
                 x_increment, y_increment,
                 base_depth = 0.0, 
                 offset = [0,0,0], 
                 offset2 = [0,0,0]):

        self.base_depth = base_depth
        self.shell_dim = shell_dim
        self.core_dim = core_dim

        self.offsetx = offset[0]
        self.offsety = offset[1]
        self.offsetz = offset[2]

        self.x_increment = x_increment
        self.y_increment = y_increment
        
        self.offset2x = offset2[0]
        self.offset2y = offset2[1]
        self.offset2z = offset2[2]
        
        self.offset_iterator = [self.offsetx, self.offsety, self.offsetz]   
       
    def Create(self, core_SLD, core_Ms, shell_SLD, shell_Ms, base_SLD = 0.0, base_Ms = 0.0):
    
        self.MaterialA = []
        self.MaterialB = Parallelapiped(SLD = shell_SLD, dim = self.shell_dim, Ms = shell_Ms)
        self.MaterialC = Parallelapiped(SLD = 2.5e-6,  dim = [self.shell_dim[0], self.shell_dim[1], 50.0], Ms = 0.0e-6 )
        self.MaterialC.on_top_of(self.MaterialB) 
 
        ## Creates 3 similar columns of cylinders
        for x in range(3):
            self.offset_iterator[1] = self.offsety
            for y in range(6):
                self.MaterialA.append(Ellipse(SLD = core_SLD, dim = self.core_dim, Ms = core_Ms))  
                self.MaterialA[len(self.MaterialA)-1].is_core_of(self.MaterialB, self.offset_iterator)
                self.offset_iterator[1] += self.y_increment 
            self.offset_iterator[0] += self.x_increment
        
        
        ## Sets new starting anchor position in the X-Y plane
        self.offsetx += self.offset2x
        self.offsety += self.offset2y
        
        self.offset_iterator = [self.offsetx, self.offsety, self.offsetz]        
         
        ## Creates the other 3 similar columns of cylinders    
        for x in range(3):  
            self.offset_iterator[1] = self.offsety
            for y in range(5):
                self.MaterialA.append(Ellipse(SLD = core_SLD, dim = self.core_dim, Ms = 2.162e-6))
                self.MaterialA[len(self.MaterialA)-1].is_core_of(self.MaterialB, self.offset_iterator)
                self.offset_iterator[1] += self.y_increment 
            self.offset_iterator[0] += self.x_increment 
            

class AlternatingSample(Sample):
    
    def __init__(self, shell_dim, core_dim, 
                 x_increment, y_increment,
                 base_depth = 0.0, 
                 offset = [0,0,0], 
                 offset2 = [0,0,0]):

        self.base_depth = base_depth
        self.shell_dim = shell_dim
        self.core_dim = core_dim

        self.offsetx = offset[0]
        self.offsety = offset[1]
        self.offsetz = offset[2]

        self.x_increment = x_increment
        self.y_increment = y_increment
        
        self.offset2x = offset2[0]
        self.offset2y = offset2[1]
        self.offset2z = offset2[2]
        
        self.offset_iterator = [self.offsetx, self.offsety, self.offsetz]   
       
    def Create(self, core_SLD, core_Ms, shell_SLD, shell_Ms, front_SLD = 0.0, front_Ms = 0.0, substrate_SLD=5.0e-6, substrate_Ms=0.0 ):
    
        self.StripeMaterial = []
        self.Matrix = Parallelapiped(SLD = shell_SLD, dim = self.shell_dim, Ms = shell_Ms)
        self.Substrate = Parallelapiped(SLD = substrate_SLD,  dim = [self.shell_dim[0], self.shell_dim[1], 100.0], Ms = substrate_Ms )
        self.Fronting = Parallelapiped(SLD = front_SLD,  dim = [self.shell_dim[0], self.shell_dim[1], 100.0], Ms = front_Ms )
        
        self.Matrix.on_top_of(self.Fronting)
        self.Substrate.on_top_of(self.Matrix)
        
          
        for i in range(6):  
            self.StripeMaterial.append(Parallelapiped(SLD = core_SLD, dim = self.core_dim, Ms = core_Ms))
            self.StripeMaterial[i].is_core_of(self.Matrix, self.offset_iterator)
            self.offset_iterator[1] += self.y_increment
            self.offset_iterator[0] += self.x_increment

    def getScene(self):
        
        self.scenelist = []
        if(self.Fronting != None):
            self.scenelist.append(self.Fronting)
        if(self.Substrate != None):
            self.scenelist.append(self.Substrate)
        self.scenelist.append(self.Matrix)
        for i in range(len(self.StripeMaterial)):
            self.scenelist.append(self.StripeMaterial[i])

        
        return Scene(self.scenelist)  
            

class TriPrismSample(Sample):
    
    def __init__(self, shell_dim, core_dim, 
                 x_increment, y_increment,
                 base_depth = 0.0, 
                 offset = [0,0,0], 
                 offset2 = [0,0,0]):

        self.base_depth = base_depth
        self.shell_dim = shell_dim
        self.core_dim = core_dim

        self.offsetx = offset[0]
        self.offsety = offset[1]
        self.offsetz = offset[2]

        self.x_increment = x_increment
        self.y_increment = y_increment
        
        self.offset2x = offset2[0]
        self.offset2y = offset2[1]
        self.offset2z = offset2[2]
        
        self.offset_iterator = [self.offsetx, self.offsety, self.offsetz]   
       
    def Create(self, core_SLD, core_Ms, shell_SLD, shell_Ms, base_SLD = 0.0, base_Ms = 0.0):
     
        self.MaterialA = []
        self.MaterialB = Parallelapiped(SLD = shell_SLD, dim = self.shell_dim, Ms = shell_Ms) 
        self.MaterialC = None
        
        for i in range(6):  
            self.MaterialA.append(TriangularPrism(SLD = core_SLD, dim = self.core_dim, Ms = core_Ms))
            self.MaterialA[i].is_core_of(self.MaterialB, self.offset_iterator)
            self.offset_iterator[1] += self.y_increment
            
            
            
