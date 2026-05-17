from .Element import *
from .IntRule1D import *
import numpy as np
from numpy import linalg as LA

class Bar(Element):
  """Bar element class for 1D simulations.
  """
#   ****************** 
#      INITIALIZER
#   ******************
  def __init__(self, E: float = -1., A: float = -1., b: float = 0., nodevec: list[Node] = []):
    """Initializer for Bar element.

    Args:
        E (float, optional): Young modulus. Defaults to -1..
        A (float, optional): Section area. Defaults to -1..
        b (float, optional): Distribute load. Defaults to 0..
        nodevec (list[int], optional): Element nodes. Defaults to [].
    """
    if E < 0 or A < 0:
      DebugStop("Inconsistent material properties")
    if len(nodevec) > 0:
      if (False in [len(nod.coord) == 2 or len(nod.coord) == 3 for nod in nodevec]):
        DebugStop("Nodes should always have two coordinates (x,y) for Bar simulations")
        
    super().__init__(nodevec=nodevec,dim=1)
    super().deactivateChecks()
    self.E = E
    self.A = A
    self.Dmat = E * A
    self.b = b
    super().activateChecks()
 
#   ****************** 
#        METHODS
#   ******************   
  def __str__(self):    
    return f"Material Properties:\n E = {self.E}, A = {self.A}, b = {self.b}" + super().__str__()
    
# ------------------------------------------------
# ------------------------------------------------

  def physicalDerivatives(self, dNdqsi: list[float], invjac: float)->np.array:
    dNdx = [d*invjac for d in dNdqsi]
    return dNdx   

# ------------------------------------------------
# ------------------------------------------------

  def createB(self, dNdx: list[float])->np.array:
    B=np.array(dNdx).reshape(1,-1)
    return B
  
# ------------------------------------------------
# ------------------------------------------------

  def jacobian(self, dNdqsi: list[float])->tuple[float,list[float],list[float]]:
    # Assume equally spaced nodes for both Bar2Node and Bar3Node
    x1=self.nodevec[0].coord[0]
    x2=self.nodevec[1].coord[0]
    lel=x2-x1
    jac=lel/2
    invjac=1/jac
    return (jac,invjac,self.physicalDerivatives(dNdqsi, invjac)) 
  
# ------------------------------------------------
# ------------------------------------------------

  def postprocvar(self, var, qsivec: list[float], ue: list[float])->list[float]:
    # Compute shape functions and derivatives
    N, dNdqsi = self.shape(qsivec)
    jac, invjac, dNdx = self.jacobian(dNdqsi)
    # Do an if/else for the different tipes of post processing variables (var)
    if var == "Displacement":
      u=np.dot(N,ue)
      return [u]
    elif var == "Strain":
      B=self.createB(dNdx)
      strain=np.dot(B,ue)[0]
      return [strain]
    elif var == "Stress":
      B=self.createB(dNdx)
      strain=np.dot(B,ue)[0]
      stress=self.E*strain
      return [stress]
    # For instance, to compute strain you wull need to createB and do B@ue
    else:
      DebugStop("Invalid post processing variable")  


# ------------------------------------------------
# ------------------------------------------------
  
  def nstatevar(self):
    return 1 # only ux

# ------------------------------------------------
# ------------------------------------------------

  def intrule(self)->IntRule1D:    
    return IntRule1D(self.pord)

# ------------------------------------------------
# ------------------------------------------------
