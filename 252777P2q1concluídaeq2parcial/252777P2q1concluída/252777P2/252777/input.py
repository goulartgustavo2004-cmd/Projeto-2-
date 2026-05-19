import numpy as np
import matplotlib.pyplot as plt
from PyFET.Node import *
from PyFET.Bar2Node import *
from PyFET.Mesh import *
from PyFET.TrussDSM import *
from PyFET.TrussDSMBC import *
from PyFET.Analysis import *
def main():
    E=2.5e10
    A=0.01
    b=0.0
    #nós inferiores
    n1=Node(np.array([0.,0.]))
    n2=Node(np.array([2.,0.]))
    n3=Node(np.array([4.,0.]))
    n4=Node(np.array([7.,0.]))
    #nós superiores
    n5=Node(np.array([2.,2.]))
    n6=Node(np.array([4.,2.]))
    n7=Node(np.array([7.,2.]))
    n8=Node(np.array([9.,2.]))
    #barras horizontais inferiores
    b1=TrussDSM(E=E,A=A,b=b,nodevec=[n1,n2])
    b3=TrussDSM(E=E,A=A,b=b,nodevec=[n2,n3])
    b4=TrussDSM(E=E,A=A,b=b,nodevec=[n3,n4])
    #barras horizontais superiores
    b8=TrussDSM(E=E,A=A,b=b,nodevec=[n5,n6])
    b10=TrussDSM(E=E,A=A,b=b,nodevec=[n6,n7])
    b12=TrussDSM(E=E,A=A,b=b,nodevec=[n7,n8])
    #barras verticais
    b2=TrussDSM(E=E,A=A,b=b,nodevec=[n2,n5])
    b9=TrussDSM(E=E,A=A,b=b,nodevec=[n3,n6])
    b11=TrussDSM(E=E,A=A,b=b,nodevec=[n4,n7])
    #barras diagonais
    b6=TrussDSM(E=E,A=A,b=b,nodevec=[n1,n5])
    b7=TrussDSM(E=E,A=A,b=b,nodevec=[n3,n5])
    b5=TrussDSM(E=E,A=A,b=b,nodevec=[n4,n6])
    b13=TrussDSM(E=E,A=A,b=b,nodevec=[n4,n8])
    #bcs
    bc1=TrussDSMBC(nodevec=[n1], type="Displacement", val=[0.,0.], E=E, A=A)
    bc3y=TrussDSMBC(nodevec=[n3], type="DisplacementY", val=[None,0.], E=E, A=A)
    bc8=TrussDSMBC(nodevec=[n8], type="Displacement", val=[0.,0.], E=E, A=A)
    mesh=Mesh(dim=2)
    l1=TrussDSMBC(nodevec=[n5], type="Load", val=[1000.,-2000.], E=E, A=A)
    l2=TrussDSMBC(nodevec=[n7], type="Load", val=[0.,-1000.], E=E, A=A)
    l3=TrussDSMBC(nodevec=[n4], type="Load", val=[1000.,0.], E=E, A=A)
    mesh.addNode(n1)
    mesh.addNode(n2)
    mesh.addNode(n3)
    mesh.addNode(n4)
    mesh.addNode(n5)
    mesh.addNode(n6)
    mesh.addNode(n7)
    mesh.addNode(n8)
    mesh.addEl(b1)
    mesh.addEl(b2)
    mesh.addEl(b3)
    mesh.addEl(b4)
    mesh.addEl(b5)
    mesh.addEl(b6)
    mesh.addEl(b7)
    mesh.addEl(b8)
    mesh.addEl(b9)
    mesh.addEl(b10)
    mesh.addEl(b11)
    mesh.addEl(b12)
    mesh.addEl(b13)
    mesh.addEl(bc1)
    mesh.addEl(bc3y)
    mesh.addEl(bc8)
    mesh.addEl(l1)
    mesh.addEl(l2)
    mesh.addEl(l3)
    mesh.generateGeometryVTK("geometry.vtk")
    analysis= Analysis(mesh=mesh)
    analysis.run()
    analysis.generateSolutionVTK(filename="solution.vtk", vecnames=["displacement"], scalnames=["normal"])

if __name__ == "__main__":
    main()