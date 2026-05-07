from .TrussDSM import *

class TrussDSMBC(TrussDSM):
  """Class that implements boundary conditions for TrussDSM elements
  """
#   ****************** 
#      INITIALIZER
#   ******************
  def __init__(self, nodevec: list[Node] = [], type: str = "None", val: list[float] = [], E: float = 0., A: float = 0.):    
    """Initializer for TrussDSMBC class

    Args:
        nodevec (list[Node], optional): Vector with node objects. Defaults to [].
        type (str, optional): Type of BC (Displacement, Load, etc.). Defaults to "None".
        val (list[float], optional): Value to be imposed. Defaults to [].
        E (float, optional): Young modulus. Defaults to 0..
        A (float, optional): Section area. Defaults to 0..
    """
    super().__init__(E=E,A=A,nodevec=nodevec)
    super().deactivateChecks()
    self.type = type
    if len(val) != 2:
      DebugStop("Val must have 2 positions (x,y)")
    self.val = val
    self.dim = 0
    super().activateChecks()  

#   ****************** 
#        METHODS
#   ****************** 
  def calcstiff(self)->tuple[np.ndarray,np.ndarray]:  

    # Create kel and fel with the right sizes
    kel = np.zeros((self.ndofs(),self.ndofs()))
    fel = np.zeros((self.ndofs(),1))

    if self.type == "Displacement":
      kel[0,0] += self._bignumber
      kel[1,1] += self._bignumber
      fel[0,0] += self._bignumber * self.val[0]
      fel[1,0] += self._bignumber * self.val[1]

    elif self.type == "Load":
      fel[0,0] += self.val[0]
      fel[1,0] += self.val[1]

    elif self.type == "DisplacementX":        
      kel[0,0] += self._bignumber
      if np.isscalar(self.val):
        fel[0,0] += self._bignumber * self.val
      else:
        fel[0,0] += self._bignumber * self.val[0]

    elif self.type == "DisplacementY":        
      kel[1,1] += self._bignumber
      if np.isscalar(self.val):
        fel[1,0] += self._bignumber * self.val
      else:
        fel[1,0] += self._bignumber * self.val[1]

    else:
      DebugStop("Not defined or implemented boundary condition")

    return kel, fel

# ------------------------------------------------
# ------------------------------------------------

def vtktype(self) -> int:
  return 1 # 1 is the tag for points in vtk
  
# ------------------------------------------------
# ------------------------------------------------

def intrule(self):    
  DebugStop("Should not be called")

# ------------------------------------------------
# ------------------------------------------------    