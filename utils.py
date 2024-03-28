import numpy as np
import cvxpy as cp
from scipy.linalg import block_diag
import control
import scipy as sp

#----------------------------------------
# Dynamics
#----------------------------------------

class PendulumSimulator:
    def __init__(self, params_dic, Delta):
        # Initialize parameters
        self.m = params_dic['m']
        self.g = params_dic['g']
        self.l = params_dic['l']
        self.b = params_dic['b']
        self.Delta = Delta

    def sim(self, x, u):
        u = u[0]  # Unpack the control input
        # Integrate the dynamics with scipy solve_ivp
        x_next = sp.integrate.solve_ivp(self.dynamics, [0, self.Delta], x, args=(u,))
        return x_next.y[:, -1] # Return the final state

    def dynamics(self, t, x, u):
        # Unpack states
        th = x[0]
        th_dot = x[1]

        # Compute the acceleration
        th_ddot = (-self.b*th_dot - self.m*self.g*self.l*np.sin(th) + u)/(self.m*self.l**2)
        return (th_dot, th_ddot)

def linearize_dyn(x_lin, u_lin, params_dic, Delta):
    # Extract parameters.
    m = params_dic['m']
    g = params_dic['g']
    l = params_dic['l']
    b = params_dic['b']

    theta_lin = x_lin[0]

    A = np.array([[0, 1],
                  [-g/l*np.cos(theta_lin), -b/(m*l**2)]])
    
    B = np.array([[0],
                  [1/(m*l**2)]])
    return A, B

#----------------------------------------
# Pi-MPC, Periodic Disturbance Augmentation Functions
#----------------------------------------

def create_Sd(nd, Nperiod):
    # Create a block diagonal matrix with Nperiod copies of the identity matrix of size nd
    # The matrix should be of size (Nperiod*nd, Nperiod*nd)
    assert isinstance(nd, int)
    assert isinstance(Nperiod, int)

    Sd = block_diag(*[np.eye(nd) for _ in range(Nperiod)])

    # Roll the matrix by nd to the right
    Sd = np.roll(Sd, nd, axis=1)

    return Sd

def create_Ssel(nd, Nperiod):
    # Create a selection matrix that selects the first nd elements of the state vector
    Ssel = np.block([np.eye(nd), np.zeros((nd, (Nperiod-1)*nd))])

    return Ssel

def periodic_augmented(A, B, C, Bd, Cd, Nperiod):
    # Get dimensions
    nx = A.shape[0]
    nu = B.shape[1]
    nd = Bd.shape[1]

    Sd = create_Sd(nd, Nperiod)
    Ssel = np.block([np.eye(nd), np.zeros((nd, (Nperiod-1)*nd))])

    # Augmented state space model
    A_aug = np.block([[A,                          Bd @ Ssel],
                   [np.zeros((Nperiod*nd, nx)), Sd]])
    B_aug = np.block([[B],
                   [np.zeros((Nperiod*nd, nu))]])
    C_aug = np.block([[C, Cd @ Ssel]])

    return A_aug, B_aug, C_aug

#----------------------------------------
# MPC Functions
#----------------------------------------

def cp_block_diag(A, n):
    Afull = []
    dim1 = A.shape[0]
    dim2 = A.shape[1]
    for i in range(n):
        cur = [np.zeros((dim1, dim2))] * n
        cur[i] = A
        Afull.append(cur)
    return cp.bmat(Afull)

def MPC(simulator, dims_dic, MPC_dic, dyn_dic):
    # Baseline MPC function
    
    # Unpack dictionaries
    # dims_dic = {'nx':nx,'nu':nu,'ny':ny,'nr':nr,'nd':nd}
    nx = dims_dic['nx']
    nu = dims_dic['nu']
    ny = dims_dic['ny']
    nr = dims_dic['nr']
    nd = dims_dic['nd']

    # MPC_dic = {'Nsim':Nsim,'Nhor':Nhor,'Q':Q,'R':R,'x0':x0,'u0':u0,'ref_z':ref_z}
    Nsim = MPC_dic['Nsim']
    Nhor = MPC_dic['Nhor']
    Q = MPC_dic['Q']
    R = MPC_dic['R']
    x0 = MPC_dic['x0']
    u0 = MPC_dic['u0']
    ref_z = MPC_dic['ref_z']
    u_max = MPC_dic['u_max']

    # dyn_dic = {'A':A,'B':B,'C':C,'H':H,'Delta':Delta}
    A = dyn_dic['A']
    B = dyn_dic['B']
    C = dyn_dic['C']
    H = dyn_dic['H']
    Delta = dyn_dic['Delta']

    # Design observer
    Q_kalman = block_diag(np.eye(nx))
    R_kalman = np.eye(ny)
    K, S, E = control.dlqr(A.T, C.T, Q_kalman, R_kalman)
    Lx = -K.T

    # Construct the cvxpy problem
    cvx_usp = cp.Parameter((Nhor * nu))
    cvx_usp.value = np.zeros((Nhor * nu))
    cvx_zsp = cp.Parameter(((Nhor+1) * nr))
    cvx_zsp.value = np.zeros(((Nhor+1) * nr))
    cvx_x0 = cp.Parameter(nx)
    cvx_x0.value = x0
    cvx_u = cp.Variable((Nhor * nu))
    cvx_x = cp.Variable(((Nhor+1) * nx))

    # Define cost function and constraints
    Qfull = block_diag(*[Q]*(Nhor+1))
    Rfull = block_diag(*[R]*(Nhor))
    HCfull = cp_block_diag(H@C, Nhor+1)
    obj = cp.Minimize(cp.quad_form(HCfull @ cvx_x- cvx_zsp, Qfull) + cp.quad_form(cvx_u - cvx_usp, Rfull))
    # Constraints
    Afull = cp_block_diag(A, Nhor)
    Bfull = cp_block_diag(B, Nhor)
    constr = [cvx_x[0:nx] == cvx_x0] # Initial condition
    constr += [cvx_x[nx:] == Afull @ cvx_x[:-nx] + Bfull @ cvx_u] # Dynamics
    constr += [cvx_u >= -u_max, cvx_u <= u_max] # Input constraints
    # Build problem
    prob = cp.Problem(obj, constr)

    # Functions to generate reference trajectory for the next Nhor steps
    def set_zsp(k):
        for i in range(Nhor+1):
            t = (k + i)*Delta
            cvx_zsp.value[i*nr:(i+1)*nr] = ref_z(t)
    def set_usp(k):
        for i in range(Nhor):
            cvx_usp.value[i*nu:(i+1)*nu] = 0

    # Initialize
    x = np.zeros((Nsim, nx))
    x[0,:] = x0
    u = np.zeros((Nsim, nu))
    u[0,:] = u0
    y = np.zeros((Nsim, ny))
    y[0,:] = C @ x0
    z = np.zeros((Nsim, nr))
    z[0,:] = H @ y[0,:]

    ref = np.zeros((Nsim, nr))

    # Observer variables
    xhat = np.zeros((Nsim, nx))
    xhat[0,:] = x0
    innov = np.zeros((Nsim, nx))

    x_planned = np.zeros((Nsim, nx))

    # Solve the problem
    for k in range(Nsim-1):
        # Set parameters
        cvx_x0.value = xhat[k,:]
        set_zsp(k)
        ref[k,:] = cvx_zsp.value[:nr]
        set_usp(k)
        # Solve the problem
        try:
            prob.solve(solver=cp.OSQP, warm_start=True)
        except:
            print("***Error during optimization!")
            break

        # Store the solution
        if k < Nsim-Nhor and k % Nhor == 0:
            for i in range(Nhor):
                x_planned[k+i, :] = cvx_x.value[(i)*nx:(i+1)*nx]
        u[k,:] = cvx_u.value[:nu]
        if u[k,:] > u_max+1:
            print("***Warning: Input saturation! k = ", k)
            break
        # Simulate with nonilnear model.
        try:
            x[k+1,:] = simulator.sim(x[k,:], u[k,:])
        except:
            print("***Error during simulation!")
            break
        y[k+1,:] = C @ x[k+1,:]
        z[k+1,:] = H @ y[k+1,:]

        # Simulate observer
        innov[k,:] = Lx @ ( C @ xhat[k,:] - y[k,:] )
        xhat[k+1,:] = A @ xhat[k,:] + B @ u[k,:] + innov[k,:nx]
    
    print("MPC: Done")
    return x, u, y, z, xhat, ref, x_planned, innov