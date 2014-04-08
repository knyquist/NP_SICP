# Kalafut and Visscher SIC step fitting algorithm implementation  

I have the code organized into three files...these are in the 'Kalafut' folder  
     -- kal_main.c       the main file   
     -- SICalgorithm.c   file that does SIC algorithm for no step hypoth. and for adding steps  
     -- InputTrace.c     file that reads in a dataset from a file  
(there are corresponding header files where needed)  
To run the code:  
     1) open kal_main.c and type in the path/filename of the dataset and your desired output file(line 17/18)  
                -- sidenote: I have pre-generated datasets in the 'Data' folder  
     2) type "make -f KalafutMake" to compile  
     3) type "./kal" to run  