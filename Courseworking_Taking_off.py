import numpy as np
import math as m
from scipy.integrate import ode
import matplotlib.pyplot as plt

h=0
Re=6371000
G=6.67e-11
Me=5.97e24
coeff_Cx=0.1 #В описанни же стоит аэродинамический коэффициент или нет?
totalTime=0

S3_Diametr=6.6
TotalMass=5500 #нужно учесть всю массу!
V0=0
angle=0

def fout(t, y):# обработчик шага 
        ts.append(t)
        ys.append(list(y.copy()))
        y1, y2, y3, y4 = y
        h = m.sqrt(y1*y1 +y3*y3) - Re
        if (h>=185000):
          return -1
  
    
def rho(x,y): #нужно ли это учитывать при взлёте?
    '''с помощью линейной апроксимации определяет давление на необходимой нам высоте'''
    Space = [0, 0, 1.85 * 0.00001, 1.5*0.0001, 3 * 0.0001, 1.03 * 0.001, 4 * 0.001, 7.26 * 0.001, 0.0136, 0.0251, 0.0469, 0.0889, 0.1216, 0.1665, 0.2279, 0.3119, 0.3648, 0.4135, 0.4671, 0.5258, 0.59, 0.6601, 0.7365, 0.8194, 0.9093, 1,1.1]
            # плотность для разных высот
    Space_lst = [185000, 100000, 80000, 70000, 60000, 50000, 40000, 36000, 32000, 28000, 24000, 20000, 18000, 16000, 14000, 12000, 11000, 10000, 9000, 8000, 7000, 6000, 5000, 4000, 3000,2000,1000]    
   
    h = m.sqrt(x*x +y*y) - Re
    
    i = 25
    while h > Space_lst[i]:
        i-=1
        if i < 0:
            return 0
           
    delta = h - Space_lst[i]
# разница между высотой и ближайшим значением
    delta_h = Space_lst[i-1] - Space_lst[i]
# разница между ближайшими соседями
    otn = delta/delta_h
# относительное отклонение
    p = Space[i] + ((Space[i-1] - Space[i]) * otn)
    return p


        
# функция правых частей системы ОДУ
def f(t, y):
             
         y1, y2, y3, y4 = y
               
         vv=m.sqrt(y2*y2+y4*y4)
         coeff_Cy=0.34*coeff_Cx
         Sm=4/(m.pi*S3_Diametr*S3_Diametr) #площадь поверхности
         angle=40
         resistant_a=coeff_Cx*rho(y1,y3)*vv*vv*Sm/2
         lift_a=coeff_Cy*rho(y1,y3)*vv*vv*Sm/2*m.cos(angle*m.pi/180)
         ax=y1*G*Me/((y1*y1+y3*y3)**1.5)-resistant_a*(y2/vv)-lift_a*(y4/vv)
         ay=y3*G*Me/((y1*y1+y3*y3)**1.5)-resistant_a*(y4/vv)+lift_a*(y2/vv)
         #a=m.sqrt(ax*ax+ay*ay)
         #if a>120:
             #print('oops',a)
        
         return [y2,ax, y4,ay] 

tmax=900  
       

x_start=0
Vx_start=0
y_start=Re
Vy_start=0.1
#Vz_start=0
#z_start=0

xc,yc=[],[]
for i in range(0, 630):
    xc.append(Re*m.cos(i/100))
    yc.append(Re*m.sin(i/100))


y0,t0=[x_start,  Vx_start, y_start, Vy_start ], 0 
ODE=ode(f)
ODE.set_integrator('dopri5')
ODE.set_solout(fout)
ts, ys = [ ],[ ]
ODE.set_initial_value(y0, t0) 
ODE.integrate(tmax)      
Y=np.array(ys)

H=m.sqrt(Y[-1:,0]*Y[-1:,0] +Y[-1:,2]*Y[-1:,2])-Re
Vx=Y[-1:,1]
Vy=Y[-1:,3]
V0=m.sqrt(G*Me/(Re+H)) # Орбитальная скорость

plt.plot(Y[:,0],Y[:,2],linewidth=3)
plt.axis('equal')
plt.plot(xc,yc,linewidth=2)

plt.grid(True)
plt.show()

print('Высота =',H, 'Vx =', Vx, 'Vy =', Vy)
print('Время =',ts[-1]);