from .OOInterface import *
from .Element import *
from .Mesh import *
import numpy as np
from .TrussDSMBC import *
class Analysis(OOInterface):
  """
  Class responsible for assembling and solving the global system of equations
  """
#   ****************** 
#      INITIALIZER
#   ******************
  def __init__(self,mesh: Mesh):
    """Initializer for Analysis

    Args:
        mesh (Mesh): FEM mesh to be used in the analysis
    """
    super().__init__()
    self.mesh = mesh
    self.K = None
    self.F = None
    self.u = None
    super().activateChecks()

#   ****************** 
#        METHODS
#   ******************  
  def assemble(self):
    # CHALLENGE1: Implement the assembly of the global stiffness matrix and load vector using parallel computing
    # CHALLENGE2: Create the global stiffness matrix using sparse matrix paradigms
    print("--------- Starting assemble ---------")
    print(f"Number of elements (counting bcs) = {self.mesh.nel()}")
    print(f"Number of nodes = {self.mesh.nnodes()}")
    

    # Rough sketch of the assembly algorithm
    # Initialize global matrix and load vector with zeros
    ndofsglobal=self.mesh.nnodes()*self.mesh.dim
    self.K = np.zeros((self.mesh.nnodes()*self.mesh.dim,self.mesh.nnodes()*self.mesh.dim))
    self.F = np.zeros((self.mesh.nnodes()*self.mesh.dim,1))
    # Loop over all elements
    for el in self.mesh.elvec:

      el.buildEFT()
      if len(el.nodevec) == 2:
        kel, fel = TrussDSM.calcstiff(el)

      elif len(el.nodevec) == 1:
        kel, fel = TrussDSMBC.calcstiff(el)

      else:
        DebugStop(f"Tipo não reconhecido: {type(el)}")
    # Compute local stiffness matrix and load vector
    # Compute local stiffness matrix and load vector
    # Add to global matrix and global load vector
      for i in range(el.ndofs()):
        Ii=el.eft[i]
        self.F[Ii,0] += fel[i,0]
        for j in range(el.ndofs()):
          Jj=el.eft[j]
          self.K[Ii,Jj] += kel[i,j]
    print("--------- Finished assemble ---------")            

# ------------------------------------------------
# ------------------------------------------------

  def solve(self):
    print("--------- Starting solve ---------")

    # NOTE: Use np.linalg.solve to solve the system of equations
    self.u = np.linalg.solve(self.K, self.F)
    print("u = ",self.u)
    print("--------- Finished solve ---------")

# ------------------------------------------------
# ------------------------------------------------

  def run(self):
    self.assemble()
    self.solve()
    # Call the assemble and solve methods

# ------------------------------------------------
# ------------------------------------------------

  def generateSolutionVTK(self,filename: str, vecnames: list[str], scalnames: list[str]):

    print("--------- Starting generateSolutionVTK ---------")

    # First, generate the VTK file with just the mesh
    self.mesh.generatePostProcVTK(filename)

    # Second, write the solution at the nodes for vector variables
    ntotalnodes = np.sum([el.nnodes() for el in self.mesh.elvec if el.dim == self.mesh.dim])
    if ntotalnodes == 0:
      DebugStop("No nodes to write in the VTK file. Check dimension attribute of mesh and elements")
    f = open(filename,"a")
    f.write(f"\nPoint_data {ntotalnodes}\n")      
    for var in vecnames:
      f.write(f"\nVectors {var} float\n")  
      for el in self.mesh.elvec:
        if el.dim != self.mesh.dim: continue
        sol = self.u[el.eft]
        for inod in range(el.nnodes()):
          qsivec = el.qsinod(inod)
          solnod = el.postprocvar(var,qsivec,sol)
          for idim in range(len(solnod)):
            f.write(f"{solnod[idim]} ")         
          for bogus in range(el.dim+1,4):
            f.write("\t0.0") # Adding z coordinate
          f.write("\n")
            

    # Third, write the solution at the nodes for scalar variables
    for var in scalnames:
      f.write(f"\nScalars {var} float\n")
      f.write("LOOKUP_TABLE default\n")
      for el in self.mesh.elvec:
        if el.dim != self.mesh.dim: continue
        sol = self.u[el.eft]
        for inod in range(el.nnodes()):
          qsivec = el.qsinod(inod)
          solnod = el.postprocvar(var,qsivec,sol)
          f.write(f"{solnod[0]}\n")  # <- ÚNICO valor por nó (float limpo)
      f.write("\n")

    print("--------- Finished generateSolutionVTK ---------")

# ------------------------------------------------
# ------------------------------------------------

  def computeEnergyNormError(self, Uex: float)->float: 
    DebugStop("YOUR CODE GOES HERE")
    # Use the form 1/2 * a^T * K_glob * a.
    # In this class variables, it is self.u.T @ self.K @ self.u
    # Note that the stiffness matrix cannot have any boundary condition applied to it
    # so you will have to recompute it without the boundary conditions applied through big numbers
    # # The error is the sqrt of abs(1- Uh/Uex). 

# ------------------------------------------------
# ------------------------------------------------