# GraKelTest


## HOW TO RUN:

0. the _inp_ folder contains the dummy test xyz files to be compared. In this case, for simplicity all 3 files are identical.  

1.  ```
    python BuildGraphs.py inp
    ```
    
    - **INPUT:** _inp_      folder with XYZ of structures
    - **OUTPUT:** _csv_inp_  folder with graphs of structures
    
    - Takes a XYZ file and creates a weighted or unweighted graph for it. The new graph is saved in a separate folder in .csv format. 

2.  ```
    python GraphKernels.py csv_inp nh
    ```

    - **INPUT:** _csv_inp_
    - **OUTPUT:** _KernelSimilarity.csv_    A similarity matrix

    
 
    - creates dummy node labels (corresponding to possible atom labels)
    - performs a 1 vs 1 kernel comparison using 2 for loops (not ideally how GraKel should be used)
