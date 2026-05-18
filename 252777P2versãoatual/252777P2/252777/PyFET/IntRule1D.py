import numpy as np
import math
from .OOInterface import *

class IntRule1D(OOInterface):
  """Class used to generate integration rule for 1D elements
  """
#   ****************** 
#      INITIALIZER
#   ******************    
  def __init__(self, porder: int):
    """Initializer for IntRule1D class

    Args:
        porder (int): polynomial order to integrate
    """
    super().__init__()
    super().deactivateChecks()
    self.porder = porder

      

    # Calculate the number of points and the Gauss-Legendre quadrature points and weights
    self.num_points = math.ceil((porder + 1) / 2)

    # Generate integration rule with points organized as self.x = [[],[],[],...] and weights organizes as self.w = [,,,...]
    from numpy.polynomial.legendre import leggauss
    x, w = leggauss(self.num_points)
    self.x = x # Fill this attribute!
    self.w = w # Fill this attribute!
    
    super().activateChecks()
    
#   ****************** 
#        METHODS
#   ******************      
  def rule(self):
    # Calculate the Gauss-Legendre quadrature points and weights
    return self.x, self.w

  def numPoints(self):
    return self.num_points