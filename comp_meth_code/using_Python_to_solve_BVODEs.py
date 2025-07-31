"""Calculations for Using Python to Solve BVODEs """

# import libraries
import numpy as np
from scipy.integrate import solve_bvp
import matplotlib.pyplot as plt

# set the dpi for figures
plt.rc("savefig", dpi = 300)

# global constants
Dax = 4.8E-4 # dm^2 /min
k = 0.72 # /min
K = 1.0
L = 1.25 # dm
D = 0.07 # dm
Vdot = 0.0023 # dm^3 /min
nA_feed = 0.0023 # mol /min
nZ_feed = 0.0 # mol /min

# derivatives function
def derivatives(ind, dep):
    # extract the mesh size and dependent variables
    nMesh = ind.size
    nA = dep[0,:]
    nZ = dep[1,:]
    wA = dep[2,:]
    wZ = dep[3,:]

    # allocate storage for dydz
    ddz = np.zeros((dep.shape))

    # evaluate the derivatives at each mesh point
    for i in range(0,nMesh):
        # evaluate the derivatives
        dnAdz = wA[i]
        dnZdz = wZ[i]
        dwAdz = 1/Dax*4*Vdot/np.pi/D**2*wA[i] + k/Dax*nA[i] - k/Dax/K*nZ[i]
        dwZdz = 1/Dax*4*Vdot/np.pi/D**2*wZ[i] - k/Dax*nA[i] + k/Dax/K*nZ[i]
        ddz[:,i] = [dnAdz, dnZdz, dwAdz, dwZdz]

    # return the derivatives
    return ddz

# boundary conditions residuals function
def residuals(dep_lb, dep_ub):
    # extract the boundary values needed to evaluate the residuals
    nAlb = dep_lb[0]
    nZlb = dep_lb[1]
    wAlb = dep_lb[2]
    wZlb = dep_lb[3]
    wAub = dep_ub[2]
    wZub = dep_ub[3]

    # evaluate the residuals
    epsilon_1 = nAlb - nA_feed - Dax*np.pi*D**2/4/Vdot*wAlb
    epsilon_2 = nZlb - nZ_feed - Dax*np.pi*D**2/4/Vdot*wZlb
    epsilon_3 = wAub
    epsilon_4 = wZub

    # return the residuals
    epsilon = np.array([epsilon_1, epsilon_2, epsilon_3, epsilon_4])
    return epsilon

# axial dispersion reactor function
def reactor_variables():
	# set the initial mesh with 20 mesh points
    ind = np.linspace(0, L, 20)

    # set the guess
    depGuess = np.zeros((4,20))

    # solve the bvodes
    soln = solve_bvp(derivatives, residuals, ind, depGuess)

    # check if a solution was found
    if not soln.success:
        print("")
        print(f"The solver was not successful: {soln.message}")

    # extract and return the profiles
    z = soln.x
    nA = soln.y[0,:]
    nZ = soln.y[1,:]
    wA = soln.y[2,:]
    wZ = soln.y[3,:]
    return z, nA, nZ, wA, wZ

# deliverables function
def deliverables():
	# get the reactor profiles
    z, nA, nZ, wA, wZ = reactor_variables()

    # report the initial and final mesh sizes
    print("")
    print(f"Size of initial mesh: 20, size of final mesh: {len(z)}")
    print("")
    
    # plot the results
    plt.figure(1) 
    plt.plot(z, nA, color = 'k', label = 'A') 
    plt.plot(z, nZ, color = 'b', label = 'Z')
    plt.xlabel("z (dm)")
    plt.ylabel("Molar Flow Rate (mol/min)")
    plt.legend()

    # save and show the figure
    plt.savefig("results.png")
    plt.show()
    return

if __name__=="__main__":
    deliverables()
