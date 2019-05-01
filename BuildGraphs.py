import numpy as np
import pandas as pd
import os, glob, re, sys
import networkx as nx
from scipy.spatial import distance_matrix
import subprocess as sp
from sklearn.preprocessing import MinMaxScaler

"""       BUILDS A BINARY OR WEIGHTED GRAPH FROM XYZ FILE


        Possible Inputs: 3x3x3 Supercell, NeighShell (15 molecules)

        Graph Edges seen as contacts distance cutoff: sum(vdw_radii)+2 A
        Edges = Distances ≤ sum of the VdW radii + 2 Angstrom will be considered
        (CCDC_Compack inspired cutoff)

"""

# ‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹ USER INPUT ››››››››››››››››››››››››››››››››››››››››››››››››››
inp_folder = sys.argv[1]   # dir with xyz to convert into graphs


# read files without the 2 lines at the top
## variables which can be easily changed, hard coded for comparison of graphs
skip = 0#2
add_to_threshold = 2  # extra 2 A

atom_vdw_radii = {
              'Al': 2, 'Sb': 2, 'Ar': 1.88, 'As': 1.85, 'Ba': 2,
              'Be': 2, 'Bi': 2, 'B': 2, 'Br': 1.85, 'Cd': 1.58,
              'Cs': 2, 'Ca': 2, 'C': 1.7, 'Ce': 2, 'Cl': 1.75,
              'Cr': 2, 'Co': 2, 'Cu': 1.4, 'Dy': 2, 'Er': 2,
              'Eu': 2, 'F':  1.47, 'Gd': 2, 'Ga': 1.87, 'Ge': 2,
              'Au': 1.66, 'Hf': 2, 'He': 1.4, 'Ho': 2, 'H': 1.09,
              'In': 1.93, 'I': 1.98, 'Ir': 2, 'Fe': 2, 'Kr': 2.02,
              'La': 2, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2, 'Mg': 1.73,
              'Mn': 2, 'Hg': 1.55, 'Mo': 2, 'Nd': 2, 'Ne': 1.54,
              'Ni': 1.63, 'Nb': 2, 'N':  1.55, 'Os': 2, 'O':  1.52,
              'Pd': 1.63, 'P': 1.8, 'Pt': 1.72, 'K': 2.75, 'Pr': 2,
              'Pa': 2, 'Re': 2, 'Rh': 2, 'Rb': 2, 'Ru': 2, 'Sm': 2,
              'Sc': 2, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27,
              'Sr': 2, 'S': 1.8, 'Ta': 2, 'Te': 2.06, 'Tb': 2,
              'Tl': 1.96, 'Th': 2, 'Tm': 2, 'Sn': 2.17, 'Ti': 2,
              'W': 2, 'U':  1.86, 'V':  2, 'Xe': 2.16, 'Yb': 2,
              'Y': 2, 'Zn': 1.29, 'Zr': 2, 'X':  1.0, 'D':  1.0
                 }


class vdw_distance_graph:

    def __init__(self,inp_folder, name, labels, coords):
        self.inp_folder = inp_folder
        self.name = name
        self.labels = labels
        self.coords = coords
        self.dist_mat = distance_matrix(coords, coords)


    def build_graph(self, weighted=False, sparsify=True):
        # turns distance matrix into a simple 0 1 binary graph according to
        # sum of VdW distance + X A. see header
        # get vdw radii for each atom to compute sum of vdw radii
        labels_vdw= []
        for i in range(0,len(self.labels)):
            vdw = atom_vdw_radii.get(str(self.labels[i]))
            labels_vdw.append(vdw)
        # turning distance matrix into an unweighted graph - simplest version
        # decide what kind of distance matrix you want:
        counter=0

        if weighted==False:
            """ BINARY GRAPH """
            for i in range(0,int(np.shape(self.labels)[0])):
                for j in range(0,int(np.shape(self.labels)[0])):
                    sum_vdw_radii = labels_vdw[i]+labels_vdw[j]
                    threshold = float(sum_vdw_radii) +2

                    if float(self.dist_mat[i, j]) <= threshold:
                            self.dist_mat[i,j]= 1
                    else:
                        self.dist_mat[i,j]= 0

                    counter =counter + 1


        else:
            """ WEIGHTED GRAPH """
            for i in range(0,int(np.shape(self.labels)[0])):
                for j in range(0,int(np.shape(self.labels)[0])):
                    sum_vdw_radii = labels_vdw[i]+labels_vdw[j]
                    threshold = float(sum_vdw_radii) + add_to_threshold

                    if float(self.dist_mat[i, j]) <= threshold:
                            self.dist_mat[i,j]= threshold- self.dist_mat[i,j]
                    else:
                        self.dist_mat[i,j]= 0

                    counter =counter + 1

            scaler = MinMaxScaler(feature_range=(0, 1))
            self.dist_mat = scaler.fit_transform(self.dist_mat)

        # sparsify matrix nx
        if sparsify==True:
            G = nx.from_numpy_matrix(self.dist_mat)
            final_mat = nx.to_scipy_sparse_matrix(G)
        else:
            final_mat = G


        df = pd.DataFrame(final_mat.toarray())
        df.to_csv(f"./csv_{self.inp_folder}/{self.name}.csv", index=None, header=None)
        return self.dist_mat



#‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹‹››››››››››››››››››››››››››››››››››››››

def main():

    # create output folder if it doesn't already exist
    if os.path.isdir(f"csv_{inp_folder}") == True:
        pass
    else:
        sp.call(f"mkdir csv_{inp_folder}", shell=True)


    # loop through input folder
    for name in glob.glob(f"./{inp_folder}/*.xyz"):
        basename = re.findall("\d+", name)[-1] #last item get right number if res3_supercell_0.xyz => extracts 0
        print(basename)

        if os.path.exists(f"csv_{inp_folder}/{basename}.csv") == True:
            print("File already exists.")

        else:
            """ READ INPUT FILES """
            labels = np.genfromtxt(name,usecols=0,dtype=str, skip_header=skip)
            coords = np.genfromtxt(name, skip_header=skip)[:,1:]

            """ BUILD BINARY DISTANCE MATRIX """
            graph= vdw_distance_graph(inp_folder, basename, labels, coords)
            weighted_dist_mat = graph.build_graph(weighted=True)


if __name__=="__main__":
    main()
