# GraKelTest


## HOW TO RUN:

0. the _inp_ folder contains the dummy test xyz files to be compared. In this case they are identical.  

1.  ```
    python BuildGraphs.py inp*
    ```
    
    INPUT: _inp_.      folder with XYZ of structures
    OUTPUT: _csv_inp_  folder with graphs of structures

2.  ```
    python GraphKernels-bak.py csv_inp nh
    ```

    INPUT: _csv_inp_
    OUTPUT: _KernelSimilarity.csv_    A similarity matrix
