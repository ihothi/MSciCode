# MSciCode
IPython code for Classification

Advice for using:
 #### Project_Functions.py contains all the functions needed
Ignore any .ipynb files, as they were used for testing sections of the script on a smaller dataset.

To run the neural network on the (binned and stored) data set, use:  NeuralNetwork_Bin.py.   <---- In Neural_Networks
To run this in splinter, the shell script is: NeuralNetworkBin.sh

To test the performance on various amounts of testing and training plates use: PerformanceTest.py <---- In Neural_Networks
To run this in splinter, the shell script is: Performance.sh

Change TrainP = [1,3,10,30,100,300] into however many plates you want. Note: The plate numbers used are the same for testing, to vary testing plates differently change TestP.

For 
