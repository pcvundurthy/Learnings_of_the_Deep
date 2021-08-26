#The following are the instructions to 1. Generate Virtual Training data for different layouts and 2. Building a Neural Network Model to solve the inverse problem.

The folder contains five files, namely:
1. Generate_training_data_for_rectangular_layout.ipynb
2. Generate_training_data_for_hexagonal_layout.ipynb
3. Neural_Network_Model.ipynb
4. X_hex_data_5_lat.csv
5. y_hex_data_5_lat.csv

Step 1: 

To generate virtual training data, two different layouts can be choosen, namely, rectangular layouts and hexagonal layouts.

A) For Rectangular layouts:

open file (Generate_training_data_for_rectangular_layout.py)

a regular (like 4x4 or 3x3 etc..) arrangements can be chosen or an irregular arrangement (like 8x4) can be chosen.

Corresponding crack locations based on layout size and distance between each sensor, the dimensions need to be altered in the program where they are specifically commented.

The generated data file needs to be saved separately as input and output file.

B) For Hexagonal Layout:

open file (Generate_training_data_for_hexagonal_layout.py)

The size of the hexagon can be altered where commented, if necessary.

The generated data file needs to be saved separately as input and outut file.


Step 2:

To build a Neural Network, open file (Neural_Network_Model.py).

Load the input and output files generated earlier (using pandas), as commented in the program.

Change the number of input neurons based on number of sensors in the layout:
if Rectangular (4X4) = 16*3=48 neurons
if Rectangular (8X4) = 32*3=96 neurons
if Hexagonal = 13*3= 39 neurons

to change learning rate scheduler to a set of epochs, change epochs in the function to 5, 20 etc.. as commented.

The number of epochs and batch size could be altered, if necessary.

After running the program:
1.'model' variable contains the NN model
2. a graph is plotted for Model performance(Training and Validation Loss vs. Epochs)
3. a relation between prediction error in sin(beta) & cos(beta) and Validation loss is also obtained as a constant.

As a reference, the model is run for Hexagon layout with 5p for optimized hyperparameters(epochs=30, batch_size=62). 

The data files, X_hex_data_5_lat.csv and y_hex_data_5_lat.csv are used for running this model.

##########################################################THE END#####################################################