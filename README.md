# NS_MoI
NS_MoI is a python code to solve the stellar structure equations for a slowly rotating Neutron Star. It gives the mass, radius, moment of inertia and the dimensionless tidal deformability parameter ($\Lambda$).

The main code is contained in the 'code' directory and is divided into two parts, driver has the methods to calculate the required quantities and the main module implements them with multiprocessor support. 

There is no need for installation but the code requires a working python 3 installation and few other packages which can be easily installed by the pip command:

1. Numpy
2. Scipy
3. Pandas
4. Matplotlib (only for the plots)

After cloning the repository, go to the params.py file and change the paramters according to your use case and specify the equation of state. Sample code for 'ALF2' is provided.

Run the code using the following command:

```
python MI_Main.py
```
