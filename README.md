# MSciCode
IPython code for Classification

MSciCode\FinalCode is where the codes to be used can be found. The rest are the test codes used throughout this project to rest various sections of the codes - containing both faulty and working code.
Advice for using:

 #### Project_Functions.py contains all the functions needed
_For SDSS data only, Download the following_
* Storing.py  : This matches the Superset objects to the plate data'. 
* NeuralNetwork_Train.py : This is where a neural network model will be trained, and then saved in the working directory
* NeuralNetwork_Test.py : This is where the saved model is used to then predict spectra. 
###### optNN.pkl and scaler.save are the optimised neural network (and its feature scaling) used in this project

note:
Please enter the pixel rejection threshold and Bin size (for spectral binning)
Platedir - is the location of the plate folders; this code assumes the superset file is in the working directory
Bin_platedir - is where you will like to store the new FITS files

_For SDSS data and DESI, Download the following_:

* Storing.py  : This matches the Superset objects to the plate data'.
* DESI_NeuralNetwork.py: Trains a neural network and prints out the classification. This can easily be changed to a text file by saving 
                         the printed variables in a text file. 
