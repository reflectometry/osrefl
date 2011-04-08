"""
Some machines allow the sample to be rotated by 180
degrees, reversing the usual sense of the angles.

We provide two corrections, one which reverses the usual sense of the
angles, and the other which assumes that all reflection is from the
top of the film and there is no back reflectivity.

ReverseCorrection
FrontCorrection
"""

from osrefl.loaders.correction import Correction

class ReverseCorrection(Correction):
    """
    Reverse the sense of the reflection angles, making positive angles
    negative and vice versa
    """
    def apply(self,data):
        data.sample.angle_x = -data.sample.angle_x
        data.detector.angle_x = -data.detector.angle_x
    def __str__(self): return "ReverseCorrection()"

class FrontCorrection(Correction):
    """
    Assume all reflection is off the top surface, reversing the sense
    of negative angles.
    """
    def apply(self,data):
        data.sample.angle_x = abs(data.sample.angle_x)
        data.detector.angle_x = abs(data.detector.angle_x)
    def __str__(self): return "FrontCorrection()"
