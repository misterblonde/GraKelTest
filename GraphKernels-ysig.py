import numpy as np
import pandas as pd
import os, re, sys
import grakel
import time

##############################################################################

# INPUT DATA

# Folder To Compare Using Kernels:
inp_folder = str(sys.argv[1])  # eg. csv_inp

# Kernel Selection:
req_kernel = str(sys.argv[2])  #  eg. "nh"


# kernels sp kernel for now
if req_kernel == "sp": #or "shortest_path":
    kernel = grakel.GraphKernel(kernel={"name": "shortest_path"}, normalize=True)
elif req_kernel == "pm": #or "pyramid_match":
    kernel = grakel.GraphKernel(kernel={"name": "pyramid_match"}, normalize=True) # with_labels=True, L=2
elif req_kernel == "rw": #or "random_walk":
    kernel = grakel.GraphKernel(kernel={"name": "random_walk"}, normalize=True)
elif req_kernel == "nh": #or "neighborhood_hash":
    kernel = grakel.GraphKernel(kernel={"name": "neighborhood_hash"}, normalize=True)
else:
    print("Error: The kernel you're requesting doesn't exist.")
    print("Possible kernels are: shortest_path (sp), pyramid_match (pm), random_walk (rw), neighborhood_hash (nh)")
    exit()

#print(f"The selected kernel parameters are {kernel.get_params()}")


##############################################################################


def label_maker(n_atoms,shell_size, c_end):
    """
    Generate Hydrocarbon Node labels for GraKel library
    """
    node_labels = {}
    atom_no =0
    for n_mols in range(0, shell_size):

        for i in range(0,n_atoms):

            if i < c_end:
                node_labels[atom_no]= 'C'

            else:
                node_labels[atom_no] = 'H'

            atom_no += 1


    return node_labels



# def KernelSimilarity(new_graphs, node_labels):
#     """
#
#     Combines Node Labels and adds them onto graphs which are read in as .csv files
#     Builds a similarity matrix which is then saved as a .csv file
#
#     """
#
# # Initialise Empty DataFrame for Output
#     df_ =pd.DataFrame(index=new_graphs, columns=new_graphs)
#
#     ##############################################################################
#
#
#     for ref_index, ref_file in enumerate(new_graphs):
#
#         # Compute Reference Kernel all other Kernels will be compared across the row:
#         # Generate graph that all others are compared to in this row:
#         df = pd.read_csv(f"{current_path}/{inp_folder}/{ref_file}.csv", header=None)
#         G = df.values.tolist()
#         # join networkX graph data and fake labels created previously:
#         graph = [[G, node_labels]]
#         start = time.time()
#         # reference kernel in row:
#         ref_kernel = kernel.fit_transform(graph)
#         df_.loc[ref_file, ref_file]= ref_kernel.flat[0]
#         print(f"{ref_file}\t{ref_file}\t{ref_kernel.flat[0]}")
#         end = time. time()
#         print("Time needed for Reference Kernel computation (s): ", end - start)
#
#
#
#         for index, name in enumerate(new_graphs):
#         # compare similarity between reference kernel against all others in this row:
#             print("Current file is: ", name)
#             start2 = time. time()
#             graph_file = pd.read_csv(f"{current_path}/{inp_folder}/{name}.csv", header=None)#dtype=int,
#             graph_list= graph_file.values.tolist()
#             grakel_graph = [[graph_list, node_labels]]
#
#
#             # compare to reference kernel:
#             similarity = kernel.transform(grakel_graph)
#             print(f"{ref_file}\t{name}\t{similarity.flat[0]}")
#             end2 = time. time()
#             print("Time needed for Kernel comparison (s): ", end2 - start2)
#
#
#             # append similarity to Output DataFrame
#             df_.loc[ref_file,name] = similarity.flat[0]
#
#
#     # Write Matrix to file
#     df_.to_csv(f"{req_kernel}_{inp_folder}_SimilarityKernel.csv", header=True, index=True)

def KernelSimilarity(new_graphs, node_labels):
    """

    Combines Node Labels and adds them onto graphs which are read in as .csv files
    Builds a similarity matrix which is then saved as a .csv file

    """

    graphs = []
    for ref_index, ref_file in enumerate(new_graphs):

        # Compute Reference Kernel all other Kernels will be compared across the row:
        # Generate graph that all others are compared to in this row:
        df = pd.read_csv(f"{current_path}/{inp_folder}/{ref_file}.csv", header=None)
        G = df.values.tolist()
        # join networkX graph data and fake labels created previously:
        graphs.append([G, node_labels])


    start = time.time()
    out = kernel.fit_transform(graphs)
    end = time.time()
    print("Time needed for Reference Kernel computation (s): ", end - start)

    # Write Matrix to file
    pd.DataFrame(out, index=new_graphs, columns=new_graphs).to_csv(f"{req_kernel}_{inp_folder}_SimilarityKernel.csv", header=True, index=True)


#####################################################################

if __name__ == '__main__':
    # LABELS GENERATOR - create fake labels because GraKel only works with Labels:
    """
    Repeating atom labels (C and H )
    """
    node_labels = label_maker(n_atoms=26, shell_size=8, c_end=10)
    #n_atoms                # number of atoms in single molecule
    #shell_size             # number of neighbouring moleucles in shell
    #c_end                  # position of last carbon atom in XYZ file // c per mol


    current_path = os.getcwd()
    # List all files in folder (with extension):
    comp_crystals = os.listdir(f"{current_path}/{inp_folder}")
    # Gather File Numbers of files in folder into a List (lose extension):
    new_graphs = [x.strip(".csv") for x in comp_crystals]

    # Run Similarity Analysis
    KernelSimilarity(new_graphs, node_labels)
