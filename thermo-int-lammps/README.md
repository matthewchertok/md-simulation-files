# Thermodynamic Integration
This folder contains an example of running thermodynamic integration simulation in LAMMPS, in which a water molecule is inserted into bulk water (N=500) to calculate the free energy difference between these two states (N=500 and N=501).

There are many exisitng useful resources about the theory details of thermodynamic integration(for example, [Wiki](https://en.wikipedia.org/wiki/Thermodynamic_integration) and [Alchemistry](https://www.alchemistry.org/wiki/Thermodynamic_Integration)). Here I will focus on how to perform the molecule insertion in LAMMPS.
## Files

### ```500``` : The simulations files
Thermodynamic integration requires a numerical integration along a coupling parameter $\lambda$. 

$\Delta F \approx \sum_{k=1}^N w_k\left\langle\frac{dU}{d\lambda}\right\rangle_{\lambda_k} $

Here $N$ is the total number of points between 0 and 1 to perfrom the numerical integration, $\lambda_k$ and $w_k$ are the $\lambda$ position and coresponding weights. In this example, $N=12$ is used and I set the 12 $\lambda$ to a 12-point [Gaussian quadrature](https://en.wikipedia.org/wiki/Gaussian_quadrature). 

To calculate the $\left\langle\frac{dU}{d\lambda}\right\rangle_{\lambda_k}$, we will do another numerical calculation

$\left\langle\frac{dU}{d\lambda}\right\rangle_{\lambda_k} \approx \left\langle\frac{U(\lambda+\delta)-U(\lambda-\delta)}{2\delta}\right\rangle_{\lambda_k}$

Here $\delta$ is a small number that allows the numerical differentiation.

In this example, we want to insert a water molecule. To achieve that, we use $\lambda$ to describe a [soft core potential](https://docs.lammps.org/pair_fep_soft.html) between the inserted water molecule and other water molecules. When $\lambda=0$, there is no interaction bewteen them and the system is effectively a 500 water molecules system. When $\lambda=1$, there is full, normal interaction bewteen them and the system is a 501 water molecules system.

```lmp.data``` and ```lmp.settings``` are the LAMMPS system data and force field. The inserted water has atomtype **3** and **4**, other waters have atomtype **1** and **2**.

```b.in``` is the LAMMPS script to perform a simulation with a particular $\lambda$. It first equilibrates the system, then dumps a trajectory file for the following calculation of $\left\langle\cdot\right\rangle_{\lambda_k}$.

```c.in``` is the LAMMPS script to rerun the trajectory from ```b.in``` to calculate $U(\lambda)$, $U(\lambda+\delta)$, and $U(\lambda-\delta)$.

### ```slurm``` : The slurm jobs submit files

First perform all 12 ```b.in``` scripts to collect 12 trajectories at different $\lambda$, then perform all 12 ```c.in``` scripts to calculate the potential energies.

### ```TI.py``` : The results analysis script

Read in results from ```c.in``` and calculate the free energy from thermodynamic integration. Three numerical integrations are performed: $\sum_{k=1}^N w_k\left\langle\frac{U(\lambda+\delta)-U(\lambda-\delta)}{2\delta}\right\rangle_{\lambda_k}$, $\sum_{k=1}^N w_k\left\langle\frac{U(\lambda+\delta)-U(\lambda)}{\delta}\right\rangle_{\lambda_k}$, and $\sum_{k=1}^N w_k\left\langle\frac{U(\lambda)-U(\lambda-\delta)}{\delta}\right\rangle_{\lambda_k}$. These three integrations show give similar results. Large differences indicates the $\delta$ is too large. However, with smaller $\delta$, longer simulations is likely required to converage the results.