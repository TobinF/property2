# Standard library imports
from math import sqrt
import numpy as np
import json
import cmath

# 从json中读取系数
def load_coefficient():
    with open('Temperature-Density\\Data\\D0_coefficient.json','r',encoding='utf-8') as f:
        data = json.load(f)
    popt_D = np.array(data['T-D'])
    u0 = data['u0']
    T0 = data['T0']
    # popt_Cp = np.array(data['p-Cp'])
    # popt_eta = np.array(data['p-eta'])
    # popt_tcx = np.array(data['p-tcx'])
    return popt_D,u0,T0

# 从json中读取u0
def load_fitting_nodes():
    with open('Temperature-Density\\Data\\fitting_nodes.json','r',encoding='utf-8') as f:
       data = json.load(f)
    T = np.array(data['温度(K)'])
    D = np.array(data['密度(kg/m3)'])
    # Cp = np.array(data['比热容[J/(kg·K)]'])
    # eta = np.array(data['粘度[microPa.s (10^-6 Pa.s)]'])
    # tcx = np.array(data['导热系数[W/(m·K)]'])
    return T,D

# 载入转换初值
def load_transfer_data():
    with open('Temperature-Density\\Data\\transfer_data.json','r',encoding='utf-8') as f:
       data = json.load(f)

    D0 = np.array(data['密度(kg/m3)'])
    # Cp0 = np.array(data['比热容[J/(kg·K)]'])
    # eta0 = np.array(data['粘度[microPa.s (10^-6 Pa.s)]'])
    # tcx0 = np.array(data['导热系数[W/(m·K)]'])
    return D0

def transfer_D(ut0,y1,y2,y3,filepath):
    # y1,y2,y3 = solve_z()
    # 需要修改成对应转换关系
    # data = {
    # 'y1':((y1.real + u0)/0.15 +100).tolist(),
    # 'y2':((y2.real + u0)/0.15 +100).tolist(),
    # 'y3':((y3.real + u0)/0.15 +100).tolist(),
    # # 'y':y.real.tolist(),
    # }
    data = {
    'y1':((y1.real + ut0)).tolist(),
    'y2':((y2.real + ut0)).tolist(),
    'y3':((y3.real + ut0)).tolist(),
    # 'y':y.real.tolist(),
    }
    json_data = json.dumps(data,ensure_ascii=False)
    with open(filepath,'w',encoding='utf-8') as f:       
            f.write(json_data)

def solve_z(T,popt,u0,T0,ut0,filepath):
    '''
        显式求解方法
    ''' 
    print (popt)
    y1 = np.zeros(len(T),dtype=complex)
    y2 = np.zeros(len(T),dtype=complex)
    y3 = np.zeros(len(T),dtype=complex)  

    for i in range(len(T)):
        # 显式求解系数
        # a = popt[2]
        # b = popt[1]*p[i] + popt[5]
        # c = popt[0]*((p[i])**2) + popt[4]*p[i] + popt[7]
        # d = (p[i])**3 + popt[3]*((p[i])**2) + popt[6]*p[i]
        # 显式求解系数
        # a = -(1/u0 + popt[0]*(T0/u0) + popt[1]*(T0**2/u0) + popt[2]*(T0**3/u0) + popt[3]*(1/u0**2)+ popt[4]*(T0/u0**2) + popt[5]*(T0**2/u0**2) + popt[6]*(T0**3/u0**2) + popt[7]*(1/u0**3) + popt[8]*(T0/u0**3) + popt[9]*(T0**2/u0**3) +popt[10]*(T0**3/u0**3))
        # b = 1.0 + popt[0]*T[i] + popt[1]*T[i]**2 + popt[2]*T[i]**3
        # c = popt[3] + popt[4]*T[i] + popt[5]*T[i]**2 + popt[6]*T[i]**3
        # d = popt[7] + popt[8]*T[i] + popt[9]*T[i]**2 + popt[10]*T[i]**3
        # 显式求解系数
        a = popt[10] + popt[6]*T[i] + popt[2]*T[i]**2       
        b = popt[9] + popt[5]*T[i] + popt[1]*T[i]**2
        c = popt[8] + popt[4]*T[i] + popt[0]*T[i]**2
        d = popt[7] + popt[3]*T[i] + T[i]**2 - (1/u0 + popt[0]*(T0/u0) + popt[1]*(T0**2/u0) + popt[2]*(T0**3/u0) + popt[3]*(1/u0**2)+ popt[4]*(T0/u0**2) + popt[5]*(T0**2/u0**2) + popt[6]*(T0**3/u0**2) + popt[7]*(1/u0**3) + popt[8]*(T0/u0**3) + popt[9]*(T0**2/u0**3) +popt[10]*(T0**3/u0**3))*T[i]**3


        # 将方程转化为x3+rx=q的形式
        r = float((c/a - (b**2)/(3*a**2)))
        q = float((d/a + (2*b**3)/(27*a**3) - (b*c)/(3*a**2)))
        # 判别式Δ>0有一个实根；Δ<0有三个实根
        # 尝试显式计算 D
       
        k = 2j/2
        w1 = (1/2)*(-1 + sqrt(3)*k)
        w2 = (1/2)*(-1 - sqrt(3)*k)
        
        # 显式三次方程解法
        y1[i] = ((-1/2)*q + cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) + ((-1/2)*q - cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) - b/(3*a)
        
        y2[i] = w1*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w2*(((-1/2)*q - cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)
        
        y3[i] = w2*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w1*(((-1/2)*q-cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)
        
    # 判断适用于哪一种转换函数,并存储拟合后的数据
                         
    transfer_D(ut0,y1,y2,y3,filepath)
    # elif mark ==2:
    #     transfer_Cp(u0,y1,y2,y3,filepath)
    # elif mark == 3:
    #     transfer_eta(u0,y1,y2,y3,filepath)
    # elif mark == 4:
    #     transfer_tcx(u0,y1,y2,y3,filepath)
    
    # return y1,y2,y3

def main():
    popt_D,D0,T0= load_coefficient()
    Dt0 = load_transfer_data()
    T,D = load_fitting_nodes()  
    solve_z(T,popt_D,D0,T0,Dt0,'Temperature-Density\\Data\\fitting_D.json')
    # solve_z(p,popt_Cp,Cp0,'Pressure-Density\\Data\\u_fitting_Cp.json',2)
    # solve_z(p,popt_eta,eta0,'Pressure-Density\\Data\\u_fitting_eta.json',3)
    # solve_z(p,popt_tcx,tcx0,'Pressure-Density\\Data\\u_fitting_tcx.json',4)
        

if __name__ == '__main__':
    main()
    # print(__name__)