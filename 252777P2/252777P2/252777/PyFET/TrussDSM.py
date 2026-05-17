from .Element import *
import numpy as np
from numpy import linalg as LA

class TrussDSM(Element):
  """
  TrussDSM element class
  Implements the truss element through a structural analysis perspective
  """
#   ****************** 
#      INITIALIZER
#   ******************
  def __init__(self, E: float = -1., A: float = -1., b: float = 0., nodevec: list[Node] = []):
    """Initializer for TrussDSM class

    Args:
        E (float, optional): Young modulus. Defaults to -1..
        A (float, optional): Section area. Defaults to -1..
        b (float, optional): Distributed load. Defaults to 0..
        nodevec (list[int], optional): Vector with node objects. Defaults to [].  
    """
    if E < 0 or A < 0:
      DebugStop("Inconsistent material properties")
    if len(nodevec) > 2 and len(nodevec) < 1:
      DebugStop("Inconsistent nodevec vector")
    if len(nodevec) > 0:
      if (False in [len(nod.coord) == 2 or len(nod.coord) == 3 for nod in nodevec]):
        DebugStop("Nodes should always have two coordinates (x,y) for TrussDSM simulations")
        
    super().__init__(nodevec=nodevec,dim=2)
    super().deactivateChecks()
    self.E = E
    self.A = A
    self.b = b
    super().activateChecks()
 
#   ****************** 
#        METHODS
#   ******************   
  def __str__(self):    
    return super().__str__() + f"Material Properties:\n E = {self.E}, A = {self.A}, b = {self.b}"

# ------------------------------------------------
# ------------------------------------------------

  def qsinod(self, nodind: int)->list[float]:
    if nodind < 0 or nodind > 1:
      DebugStop("Invalid node index")
    return [nodind]

# ------------------------------------------------
# ------------------------------------------------

  def calcstiff(self)->tuple[np.ndarray,np.ndarray]:
    x1=self.nodevec[0].coord[0]
    y1=self.nodevec[0].coord[1]
    x2=self.nodevec[1].coord[0]
    y2=self.nodevec[1].coord[1]
    dx=x2-x1
    dy=y2-y1
    l=np.sqrt(dx*dx+dy*dy)
    c=dx/l
    s=dy/l
    kel=((self.E*self.A)/(l))*np.array([[c*c, c*s, -c*c, -c*s],
                                  [c*s, s*s, -c*s, -s*s],
                                  [-c*c, -c*s, c*c, c*s],
                                  [-c*s, -s*s, c*s, s*s]])
    fel=np.array([[self.b*l/2*c],[self.b*l/2*s],[self.b*l/2*c],[self.b*l/2*s]])
    return kel, fel

    # Compute element length l
    # Compute direction cosines c and s
    # Compute element stiffness matrix kel
    # Compute element load vector fel
    # return kel, fel

# ------------------------------------------------
# ------------------------------------------------

  def nstatevar(self):
    return 2 # ux and uy (2d truss)

# ------------------------------------------------
# ------------------------------------------------

  def postprocvar(self, var, qsivec: list[float], ue: list[float]) -> list[float]:
    ue = [u.item() if hasattr(u, 'item') else float(u) for u in ue]

    x1, y1 = self.nodevec[0].coord
    x2, y2 = self.nodevec[1].coord
    l = ((x2 - x1)**2 + (y2 - y1)**2)**0.5
    c = (x2 - x1) / l
    s = (y2 - y1) / l

    u1_local = c * ue[0] + s * ue[1]
    u2_local = c * ue[2] + s * ue[3]
    eps = (u2_local - u1_local) / l
    N = float(self.E) * float(self.A) * eps

    if var == "displacement" or var == "solution":
        if qsivec[0] == 0:
            return [ue[0], ue[1]]
        elif qsivec[0] == 1:
            return [ue[2], ue[3]]
        else:
            raise ValueError("qsivec inválido")

    if var == "axialdisplacement":
        return [float(u1_local)] if qsivec[0] == 0 else [float(u2_local)]

    if var == "strain":
        return [float(eps)]

    if var == "axialforce" or var == "normal":
        print("Axial Force:", N)
        return [float(N)] if qsivec[0] == 0 else [float(N)]
# ------------------------------------------------
# ------------------------------------------------

  def shape(self, qsivec: list[float])->tuple[list[float],list[float]]:
    DebugStop("TrussDSM element does not have shape functions")
    return [None]

# ------------------------------------------------
# ------------------------------------------------

  def nodeorder(self) -> list[int]:
    return [0,1]
  
# ------------------------------------------------
# ------------------------------------------------

  def vtktype(self) -> int:
    return 3

# ------------------------------------------------
# ------------------------------------------------

  def jacobian(self, dNdqsi: list[float])->tuple[float,list[float],list[float]]:
    DebugStop("Should never be called in TrussDSM")
  
# ------------------------------------------------
# ------------------------------------------------
    
  def intrule(self):    
    DebugStop("Should never be called in TrussDSM")

# ------------------------------------------------
# ------------------------------------------------        