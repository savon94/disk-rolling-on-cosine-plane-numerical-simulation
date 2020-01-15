
import numpy as np
import matplotlib.pyplot as plt
from math import sin 
from math import cos
from matplotlib.animation import FuncAnimation


#parameter
m1 = 100      #Kg
g = 9.8     #Gravity
r = 0.2  #meter

a=0
tperiod =3
n = 2000
t= a; dt = (tperiod/n)

h = 0.8 #slopeheight in meter
omg = 1# frequency or amount of mountain or hil...


#initialization
dq = np.zeros ((2,3))   # (xv,yv,jv)    velocity
                        # (xa,ya,ja)    acceleration

q= np.zeros ((2,3))     # (x,y,j)       distance
                        # (xv,yv,jv)    velocity

lam=np.zeros ((2))      # constraint force      
result =np.zeros ((5,1))# second order derrivative result 



def eom(t,q,dq,lam): #loop for equation of motion
    

    slope=-omg*h*sin(omg*q[0][0])
    teta=abs(np.arctan (slope)) #slopeangle
   
    

    

# mass matrix
    M=np.zeros( (5,5) )
    M[0][0]=m1
    M[1][1]=m1
    M[2][2]=(m1*(r**2))/2
    M[3][0]=1
    M[3][2]=-r
    M[4][1]=1
    M[0][3]=1
    M[2][3]=-r
    M[1][4]=1
#invers
    inv_M = np.linalg.inv(M)


#Force
    if slope<0:
        F = [[m1*g*sin(teta)],
             [-m1*g*cos(teta)],
             [r*m1*g*sin(teta)],
             [0],
             [0]]

    elif slope>0: 
        F = [[-m1*g*sin(teta)],
              [-m1*g*cos(teta)],
              [-r*m1*g*sin(teta)], 
              [0],
              [0]]

    else:
        F = [[0],
             [-m1*g],
             [0],
             [0],
             [0]]
    
# second derivative
    result = [[0],
              [0],
              [0],
              [0],
              [0]]

    lam =   [[0],
             [0]]

# iterate through rows of M
    for i in range(len(inv_M)):
   # iterate through columns of F
        for j in range(len(F[0])):
       # iterate through rows of F
            for k in range(len(F)):
                result[i][j] += inv_M[i][k] * F[k][j]
                
                asd= np.transpose(result)
                dq[1][0:3]=asd[0][0:3].T
                lam[0]=result[3][0]
                lam[1]=result[4][0]
                
                dq[0][0:3]=q[1][0:3]
               
    return dq
eom(t,q,dq,lam)

#Runge Kutta 4th order
def rk4(t,dt,q,dq,lam,n):
    
    k1 = np.zeros ((2,3),float)
    k2 = np.zeros ((2,3),float)
    k3 = np.zeros ((2,3),float)
    k4 = np.zeros ((2,3),float)
    b  = np.zeros ((2,3),float)
    fr = np.zeros ((2,3),float)
    fr = eom(t,q,dq,lam)
      
    for i in range (0,n):
        fr[i]
        k1[i] = dt*dq[i]
    for i in range (0,n):
        b[i] = q[i]+0.5*k1[i]
   
    eom(t+dt/2,q,dq,lam)
    k2 = dt*dq
    for i in range (0,n):
        b[i]=q[i]+0.5*k2[i]
    
    eom(t+dt/2,q,dq,lam)
    k3 = dt*dq
    for i in range (0,n):
        b[i] = q[i]+0.5*k3[i]
    
    eom(t+dt,q,dq,lam)
    k4 = dt*dq
    for i in range (0,2):
        q[i]= q[i] + (k1[i] + 2*(k2[i] + k3[i]) +k4[i])/6                 
    return q

#initial state
    

q[0][0]=0            #x position starting point
q[0][1]=r            #y position starting point
q[0][2]=0            #teta whell

q[1][0]=2            #x INITIAL velocity
q[1][1]=0            #y velocity must be ZERO
q[1][2]=q[1][0]/r    #angular velocity(omega=v/r)

#list
List = [] 

bb=[]
cc=[]
dd=[]
ee=[]
ff=[]
gg=[]
hh=[]
jj=[]
kk=[]
ll=[]




#time loop
while (t<tperiod):   
    if((t+dt)>tperiod): 
        dt= tperiod -t
    q = rk4(t,dt,q,dq,lam,2)

    
#graph 
    List.append(t)  #time
    bb.append(q[0][0])  #X displacement
    cc.append(q[0][2])
    dd.append(q[0][0]-r*sin(q[0][2])) #x direction center trajectory
    ee.append((r+h*cos(omg*q[0][0]))-r*cos(q[0][2]))#x direction point p trajectory 
    hh.append(r+h*cos(omg*q[0][0])) #y direction center trajectory
    jj.append(h*cos(omg*q[0][0])) #surface
    ff.append(-r*sin(q[0][2]))
    gg.append(-r*cos(q[0][2]))
    kk.append(dq[0][0])
    ll.append(dq[1][0])

    t=t+dt

#Output Result

plt.figure(200) #Velocity
plt.subplot(111)
plt.plot(List,kk,label='Disk velocity')
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('velocity (m/s)', fontsize=12)
plt.legend()
plt.xlim([0, t])

plt.show()

plt.figure(300)#Acceleration
plt.subplot(111)
plt.plot(List,ll,label='Disk acceleration')    # Disk accelertion
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('acceleration  $(m/s^{2} )$', fontsize=12)
plt.legend()
plt.xlim([0, t])
plt.show()

plt.figure(400)
plt.subplot(111)
plt.plot(List,bb,label='Travelling Distance')
plt.xlabel('Time (s)', fontsize=12)    # Traveling Distance
plt.ylabel('Traveling distance (m)', fontsize=12)
plt.legend()
plt.xlim([0, t])
plt.ylim([0, q[0][0]])
plt.show()

plt.figure(100) #Trajectory
plt.subplot(111)
plt.plot(bb,hh,label='Disk center trajectory')  # Disk center trajectory
plt.plot(dd,ee,label='Disk point P trajectory') # Disk point P trajectory (an arbitary point on a circufmference circle)
plt.plot(bb,jj,label='Surface')          # Surface
plt.xlabel('X Displacement', fontsize=12)
plt.ylabel('y Displacement', fontsize=12)
plt.legend()
plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.00) ,ncol=1)
plt.xlim([0, q[0][0]])
plt.ylim([-q[0][0]/2,q[0][0]/2])
plt.gca().set_aspect('equal', adjustable='box')
plt.show()



# set up figure and animation
fig = plt.figure(500)
ax = plt.subplot(111)
skip = 30



def init_func():
    ax.clear()
    plt.xlim((bb[0], bb[-1]))
    plt.gca().set_aspect('equal', adjustable='box')


    
def animate(i):

    ax.plot(bb[i:i+skip], hh[i:i+skip], color='blue')
    ax.plot(dd[i:i+skip], ee[i:i+skip], color='red')
    
    ax.plot(bb,jj, color='black')

    

anim = FuncAnimation(fig,
                     animate,
                     frames=np.arange(0, len(List), skip),
                     init_func=init_func,
                     interval=2)
plt.show()



