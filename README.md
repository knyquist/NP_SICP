# Nyquist and Presse SICP step fitting algorithm
(based on SIC implentation by Kalafut and Visscher but includes a prior for the expected noise of the system)

I have the code organized into three files...these are in the 'Kalafut' folder  
*   kal_main.c;       the main file   
*   SICalgorithm.c;   file that does SIC algorithm for no step hypoth. and for adding steps  
*   InputTrace.c;     file that reads in a dataset from a file  
(there are corresponding header files where needed)  

This code is very similar to my c-lang implementation of KV's SIC, but now includes two parameters which form the prior for the expected noise.
To run the code:  
1) type make to compile, kal becomes the executable  
       
2) to run: type
```
./kal dataset bias_params > output
```  
    
Notes: the dataset should have three columns:  
```
time  force  position  
```
the bias_params has a very particular required format:
```
nu num1
So num2
``` 
where num is a double and nu, So are just strings 'naming the numbers' (so to speak). The code uses the strings to confirm that the numbers are in the right order. See the provided example bias_params file to make sure you understand the format. The code will either break or not spit out the right result if you don't provide the proper inputs.
       
the output file has two columns: 
```
time  position
```

You can test the code on the provided fake dataset. Luckily, the contents in params.txt contain the proper nu and So to handle the data!

**PROTIP** if you want to generate a params.txt file on the fly you can one line it with something like:
```
touch params.txt; printf "nu num1\nSo num2" > params.txt
```
