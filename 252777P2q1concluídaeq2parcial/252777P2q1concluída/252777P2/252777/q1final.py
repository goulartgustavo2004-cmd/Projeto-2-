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
F=-1
enorm=[]
n_list=[1,2,4,8,16,32,64]

for n in n_list:
    print(f"\n{'='*50}")
    print(f"Processando: n = {n} elementos")
    print(f"{'='*50}")
    
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
    print(Area)
    for i in range(n):
        node1=nodes_list[i]
        node2=nodes_list[i+1]
        nodevec=[node1,node2]
        cA=Area[i]
        element=Bar2Node(A=cA,E=E,nodevec=nodevec)
        mesh.addEl(element)
    bc1=BarBC(nodevec=[nodes_list[0]],type="Displacement",val=[0.0],E=E,A=Area[0])
    bc2=BarBC(nodevec=[nodes_list[-1]],type="Load",val=[F],E=E,A=Area[-1])
    mesh.addEl(bc1)
    mesh.addEl(bc2)
    analysis=Analysis(mesh)
    analysis.run()
    print("\n--- VARIÁVEIS PÓS-PROCESSADAS ---")
    for i, el in enumerate(analysis.mesh.elvec):
        if type(el).__name__ == "Bar2Node":
            node1_idx = el.nodevec[0].index
            node2_idx = el.nodevec[1].index
            ue = np.array([analysis.u[node1_idx, 0], analysis.u[node2_idx, 0]])
            qsi = [0.0]
            displacement = el.postprocvar("Displacement", qsi, ue)[0]
            strain = el.postprocvar("Strain", qsi, ue)[0]
            stress = el.postprocvar("Stress", qsi, ue)[0]
            print(f"Elemento {i}:")
            print(f" Deslocamento (qsi={qsi}): {displacement:.6e}")
            print(f" Strain (qsi={qsi}): {strain:.6e}")
            print(f" Stress/Normal Force (qsi={qsi}): {stress:.6e}")
    Uex = (1.71828*F*1)/(E*A0*2)
    error = analysis.computeEnergyNormError(Uex)
    print(f"\n--- ERRO EM NORMA DE ENERGIA ---")
    print(f"Energia exata (Uex): {Uex:.10f}")
    print(f"Erro em norma de energia: {error:.10e}")
    enorm.append(error)

print(f"\n{'='*50}")
print(f"Erros calculados para todas as malhas:")
print(f"enorm = {enorm}")
print(f"{'='*50}")
print(f"\n{'='*50}")
print(f"Erros calculados para todas as malhas:")
print(f"enorm = {enorm}")
print(f"{'='*50}")

# Calcular comprimento de elemento h = L/n
h_list = [L/n for n in n_list]
# Plotar convergência em escala log-log (h vs erro)
fig, ax = plt.subplots(figsize=(10, 6))
ax.set_xscale('log')
ax.set_yscale('log')
ax.loglog(h_list, enorm, 'bo-', linewidth=2, markersize=8, label='Erro relativo na norma de energia')

ax.set_xlabel('Comprimento do elemento (h)', fontsize=12, fontweight='bold')
ax.set_ylabel('Erro relativo na norma de energia', fontsize=12, fontweight='bold')
ax.set_title('Convergência: Erro relativo vs Comprimento do Elemento (log-log)', fontsize=13, fontweight='bold')
ax.grid(True, which='both', alpha=0.3, linestyle='--')
ax.legend(fontsize=11)

plt.tight_layout()
plt.savefig('convergencia_loglog.png', dpi=300, bbox_inches='tight')
print("✅ Gráfico salvo como 'convergencia_loglog.png'")
plt.show()