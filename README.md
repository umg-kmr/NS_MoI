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

The code will output a 'EOS_Name.dat' file in the output directory and the execution time on the screen. The output file structure is as follows:

> Mass ($M_\odot$), Radius(Km), $\Lambda$, Moment of Inertia ($kg \ m^2$)

You can now either write your own code to visualise the data or use the provided jupyter notebooks in the 'notebooks' directory.

## Expample with Anaconda

If you a working anaconda install, you can create a new environment and use the code inside that.

Create a new python 3 environment with:

```
conda create -n NS_MoI python=3.11 
```
Install the required packages:

```
conda install -c conda-forge numpy scipy pandas matplotlib
```

Optionally install jupyter notebook:
```
conda install -c anaconda jupyter
``` 

## Acknowledgments
We would like to acknowledge the Ashoka High Performance Computing Cluster (The Chanakya@Ashoka) for providing computational resources and support. The program was tested and all calculations were computed on the Chanakya@Ashoka.

## Contact
For any inquiries or questions, please contact umang.kumar_phd21@ashoka.edu.in
