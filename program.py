import matplotlib.pyplot as plt
import numpy as np

def Euler(f,x0,T,h,*args,**kwargs) :
    '''
    To perform integration stepwise using the Euler method
    
    Parameters
    ----------
    f : function 
        to calculate the RHS of the differential equation.
        
    x0 : numpy array
        initial values of the variables
        
    T : float
        the time interval
        
    h : float
        the time step
        
    *args, **kwargs : arguments and keyword arguments are passed down the function f
    
    Returns
    -------
    x : numpy array
        calculated values of the variables at each time step

    '''
    x = np.zeros((2,int((T/h)+1)))
    x[:,0]=x0
    for i in range(int(T/h)):
        x[:,i+1] = x[:,i] + h*f(x[:,i],*args, **kwargs)
    return x

def heun(f,g,x0,T,h,*args,**kwargs) :
    '''
    To perform integration stepwise using the Heun Method

    Parameters
    ----------
    f : function
        to calculate the RHS of the ode
    x0 : numpy array
        initial condition of the system
    T : float
        the time interval
    h : float
        the time step
    *args : additional arguments passed down to functions g and f
    **kwargs : additional arguments passed down to functions g and f
    

    Returns
    -------
    x : numpy array
        calculated values of the system at each time step

    '''
    x=np.zeros((2,int((T/h)+1)))
    x[:,0]=x0
    y=g(f,x[:,0],T,h,*args,**kwargs)
    for i in range(int(T/h)):
        x[:,i+1] = x[:,i] + (h/2)*(f(x[:,i],*args,**kwargs)+f(y[:,i+1],*args,**kwargs))
    return x

def rk4(f,x0,T,h,*args,**kwargs) :
    '''
    To perform integration stepwise using the Runge-Kutta stage 4 Method

    Parameters
    ----------
    f : function
        to calculate the RHS of the ode
    x0 : numpy array
        initial condiiton of the system
    T : float
        the time interval
    h : float
        the time step
    *args : additional arguments to pass down to the function f
    **kwargs : additional keyword arguments to pass down to the function f

    Returns
    -------
    x : numpy array
        calculated values of the system at each step

    '''
    x=np.zeros((2,int((T/h)+1)))
    x[:,0]=x0
    for i in range(int(T/h)):
        k1=f(x[:,i],*args,**kwargs)
        k2=f(x[:,i]+(h/2)*k1,*args,**kwargs)
        k3=f(x[:,i]+(h/2)*k2,*args,**kwargs)
        k4=f(x[:,i]+h*k3,*args,**kwargs)
        x[:,i+1] = x[:,i] + (h/6)*(k1+2*k2+2*k3+k4)
    return x
        

def ode(x,d,u,al,gm,b) :
    '''
    

    Parameters
    ----------
    x : numpy array
        DESCRIPTION.
    d : float
        DESCRIPTION.
    u : float
        DESCRIPTION.
    al : float
        DESCRIPTION.
    gm : float
        DESCRIPTION.
    b : float
        DESCRIPTION.

    Returns
    -------
    numpy array
        DESCRIPTION.

    '''
    return np.array([d*x[0]*-1.0+u*np.tanh(al*x[0]+gm*x[1]) +b ,-d*x[1]+u*np.tanh(al*x[1]+gm*x[0])+b])
       
      


T = 7.

#Part A : Euler Method 1
t1 = np.arange(0.0,T+0.001,0.001)
xe1 = Euler(ode, np.array([-1,-1]), T, 0.001, d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0 )
#Part B : Euler Method 2
t2= np.arange(0.0,T+0.5,0.5)
xe2 = Euler(ode,np.array([-1,-1]),T,0.5,d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0)
#Part C : Heun Method
xh = heun(ode,Euler,np.array([-1,-1]),T,0.5,d=1.0,u=0.31,al=1.2,gm=-1.3,b=0)
#Part D : Runge-Kutta 4 Method
xrk = rk4(ode,np.array([-1,-1]),T,0.5,d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0)

#Comparison for x_1

plt.plot(t1,xe1[0,:],label="$Euler 1 : x_1$")
plt.plot(t2,xe2[0,:],label="$Euler 2 : x_1$")
plt.plot(t2,xh[0,:], label="$Heun Method : x_1$")
plt.plot(t2,xrk[0,:], label="$Runge-Kutta 4 : x_1$")
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.title("$Comparison 1$")
plt.legend()
plt.show()

#Comparison for x_2

plt.plot(t1,xe1[1,:],label="$Euler 1 : x_2$")
plt.plot(t2,xe2[1,:],label="$Euler 2 : x_2$")
plt.plot(t2,xh[1,:], label="$Heun Method : x_2$")
plt.plot(t2,xrk[1,:], label="$Runge-Kutta 4 : x_2$")
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.title("$Comparison 2$")
plt.legend()
plt.show()

#Phase Plots

#Euler 1
plt.plot(xe1[0,:],xe1[1,:],label="$Euler 1$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Euler Method 1 phase plot$")
plt.show()

#Euler 2
plt.plot(xe2[0,:],xe2[1,:],label="$Euler 2$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Euler Method 2 phase plot$")
plt.show()

#Heun
plt.plot(xh[0,:],xh[1,:],label="$Heun$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Heun Method phase plot$")
plt.show()

#Runge-Kutta 4 
plt.plot(xrk[0,:],xrk[1,:],label="$Runge Kutta 4$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Runge-Kutta 4 Method phase plot$")
plt.show()

#AGENTS START WITH OPPOSING OPINIONS 

#Part A : Euler Method 1
t1 = np.arange(0.0,T+0.001,0.001)
xe1 = Euler(ode, np.array([1,-1]), T, 0.001, d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0 )
#Part B : Euler Method 2
t2= np.arange(0.0,T+0.5,0.5)
xe2 = Euler(ode,np.array([1,-1]),T,0.5,d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0)
#Part C : Heun Method
xh = heun(ode,Euler,np.array([1,-1]),T,0.5,d=1.0,u=0.31,al=1.2,gm=-1.3,b=0)
#Part D : Runge-Kutta 4 Method
xrk = rk4(ode,np.array([1,-1]),T,0.5,d=1.0, u=-0.31,al=1.2,gm=-1.2, b=0)

#Comparison for x_1

plt.plot(t1,xe1[0,:],label="$Euler 1 : x_1$")
plt.plot(t2,xe2[0,:],label="$Euler 2 : x_1$")
plt.plot(t2,xh[0,:], label="$Heun Method : x_1$")
plt.plot(t2,xrk[0,:], label="$Runge-Kutta 4 : x_1$")
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.title("$Comparison_1 Opposing $")
plt.legend()
plt.show()

#Comparison for x_2

plt.plot(t1,xe1[1,:],label="$Euler 1 : x_2$")
plt.plot(t2,xe2[1,:],label="$Euler 2 : x_2$")
plt.plot(t2,xh[1,:], label="$Heun Method : x_2$")
plt.plot(t2,xrk[1,:], label="$Runge-Kutta 4 : x_2$")
plt.xlabel("$t$")
plt.ylabel("$x(t)$")
plt.title("$Comparison_2 Opposing$")
plt.legend()
plt.show()

#Phase Plots

#Euler 1
plt.plot(xe1[0,:],xe1[1,:],label="$Euler 1$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Euler Method 1 phase plot opp$")
plt.show()

#Euler 2
plt.plot(xe2[0,:],xe2[1,:],label="$Euler 2$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Euler Method 2 phase plot opp$")
plt.show()

#Heun
plt.plot(xh[0,:],xh[1,:],label="$Heun$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Heun Method phase plot opp$")
plt.show()

#Runge-Kutta 4 
plt.plot(xrk[0,:],xrk[1,:],label="$Runge Kutta 4$")
plt.xlabel("$x_1 (t)$")
plt.ylabel("$x_2 (t)$")
plt.title("$Runge-Kutta 4 Method phase plot opp$")
plt.show()