#fully Fuzzy linear systems in Python


#including modules for matrix operations
import numpy as np
from numpy.linalg import inv
np.set_printoptions(suppress=True)
#imput Fuzzy Matrix
a=[[6,5,3],[12,14,8],[24,32,20]]
m=[[1,2,2],[8,12,8],[10,30,19]]
n=[[4,2,1],[20,15,10],[34,30,24]]
b=[[58],[142],[316]]
h=[[30],[139],[297]]
g=[[60],[257],[514]]
#matrix for ST decomposition [A=ST]
s=[[0,0,0],[0,0,0],[0,0,0]]
t=[[0,0,0],[0,0,0],[0,0,0]]
#Algorithm for ST Decomposition
s[0][0]=a[0][0]
s[0][1]=a[1][0]
s[0][2]=a[2][0]
s[1][0]=a[1][0]
s[1][1]=((a[1][0]*a[1][0])+(a[0][0]
*a[1][1])-(a[1][0]*a[0][1]))/a[0][0]
s[1][2]=((a[2][0]*a[1][0])+(a[0][0]
*a[2][1])-(a[2][0]*a[0][1]))/a[0][0]
s[2][0]=a[2][0]
s[2][1]=s[1][2]
s[2][2]=((a[2][2]*(a[1][0]*a[1][0]))-
(a[2][2]*a[0][0]*s[1][1])-
(a[2][0]*a[1][0]*a[1][2])+(a[2][0]*a[0][2]*s[1][1])+(2*a[2][
0]*s[1][2]*a[1][0])-((a[2][0]*a[2][0])*s[1][1])-
(s[1][2]*a[0][2]*a[1][0])+(a[1][2]
*a[0][0]*s[1][2])-(a[0][0]*(s[1][2]*s[1][2])))/((a[1][0]
*a[1][0])-(a[0][0]*s[1][1]))
t[0][0]=1
t[0][1]=(a[0][1]-a[1][0])/a[0][0]
t[0][2]=((a[1][0]*a[1][2])-(a[0][2]*s[1][1])-
(a[1][0]*s[1][2])+(a[2][0]*s[1][1]))/((a[1][0]*a[1][0])-
(a[0][0]*s[1][1]))
t[1][0]=0
t[1][1]=1
t[1][2]=((a[0][2]*a[1][0])-(a[0][0]*a[1][2])-
(a[2][0]*a[1][0])+(a[0][0]*s[1][2]))/((a[1][0]*a[1][0])-
(a[0][0]*s[1][1]))
t[2][0]=0
t[2][1]=0
t[2][2]=1
#Displaying S and T matrix
print(s) # Symmetrical matrix
print(t) #Upper Triangular matrix
#Finding Inverse matrix of S and T
sinv=np.matrix(inv(s))
tinv=np.matrix(inv(t))
#Displaying Inverse matrix of S and T
print(sinv)
print(tinv)
#Multiplying inverse(T) and inverse(S)
tinv_sinv=tinv*sinv
#mtinv_sinv=np.matrix(tinv_sinv) Storing as matrix
print(tinv_sinv)
#Computation of x,y,z values using ST decomposition Algorithm
x=tinv_sinv*b
y=tinv_sinv*(h-m*x)
z=tinv_sinv*(g-n*x)
#Displaying Solutions of ST Decompostion
print("\n")
print("x=")
print(x)
print("y=")
print(y)
print("z=")
print(z)