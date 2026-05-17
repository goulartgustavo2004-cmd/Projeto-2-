from PyFET import Node,Bar2Node,Mesh,Analysis
from PyFET.Bar2Node import Bar2Node
from PyFET.Node import Node
from PyFET.Mesh import Mesh
from PyFET.BarBC import BarBC
from PyFET.Analysis import Analysis
import numpy as np
import matplotlib.pyplot as plt

L=1
A0=1
E=1
F=1
n=int(input("Enter the number of elements: "))
le=L/n
center=[]
aux=0
Area=[]
mesh=Mesh(dim=1)
while aux<=L:
    node=Node(coord=np.array([aux, 0.0]))
    mesh.addNode(node)
    c=aux+le/2
    A=A0*np.exp(-c/L)
    if c<=L:
        Area.append(A)
        center.append(c)
    aux=aux+le
nodes_list=mesh.nodevec
for i in range(n):
    node1=nodes_list[i]
    node2=nodes_list[i+1]
    nodevec=[node1,node2]
    cA=Area[i]
    element=Bar2Node(A=cA,E=E,nodevec=nodevec)
    mesh.addEl(element)
    print(element)
bc1=BarBC(nodevec=[nodes_list[0]],type="Displacement",val=[0.0],E=E,A=Area[0])
bc2=BarBC(nodevec=[nodes_list[-1]],type="Load",val=[F],E=E,A=Area[-1])
mesh.addEl(bc1)
mesh.addEl(bc2)
analysis=Analysis(mesh)
analysis.run()
print(f"\nDisplacements at the nodes:\n{analysis.u}")