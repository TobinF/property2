from math import sqrt
from matplotlib.pyplot import plot
import numpy as np
import matplotlib.pylab as plt
import cmath
import json

def get_data():
    p = input('请输入压力（MPa）:')
    T_start = float(input('请输入起始温度（K）:'))
    T_end = float(input('请输入终止温度（K）:'))
    num = int(input('请输入要获取的数据点数:'))
    return np.linspace(T_start,T_end,num),np.zeros(num),float(p)
  
def select(mark):

    if mark == 1:
        with open('Pressure-Density\\Data\\D0_supercritical_coefficient.json','r',encoding='utf-8') as f:
            data = json.load(f)
        popt_D = np.array(data['p-D'])
        with open('Pressure-Density\\Data\\transfer_supercritical_data.json','r',encoding='utf-8') as f:
            data = json.load(f)
        D0 = data['密度(kg/m3)']
               
    elif mark == 2:
        with open('Pressure-Density\\Data\\D0_superheat_coefficient.json','r',encoding='utf-8') as f:
            data = json.load(f)
        popt_D = np.array(data['p-D'])
        with open('Pressure-Density\\Data\\transfer_superheat_data.json','r',encoding='utf-8') as f:
            data = json.load(f)
        D0 = data['密度(kg/m3)']
                
    return popt_D,D0

def solve_P(p,popt,u0):
    '''
        显式求解方法
    ''' 
    # popt_supercritical_D,popt_superheat_D,supercritical_D0,superheat_D0 = load_data()
    # print (popt)
    # y1 = np.zeros(len(p),dtype=complex)
    # y2 = np.zeros(len(p),dtype=complex)
    # y3 = np.zeros(len(p),dtype=complex)  

    # for i in range(len(p)):
        # 显式求解系数
    a = popt[2]
    b = popt[1]*p + popt[5]
    c = popt[0]*((p)**2) + popt[4]*p + popt[7]
    d = (p)**3 + popt[3]*((p)**2) + popt[6]*p

    # 将方程转化为x3+rx=q的形式
    r = float((c/a - (b**2)/(3*a**2)))
    q = float((d/a + (2*b**3)/(27*a**3) - (b*c)/(3*a**2)))
    # 判别式Δ>0有一个实根；Δ<0有三个实根
    # 尝试显式计算 D
    
    k = 2j/2
    w1 = (1/2)*(-1 + sqrt(3)*k)
    w2 = (1/2)*(-1 - sqrt(3)*k)
    
    # 显式三次方程解法
    # y1[i] = ((-1/2)*q + cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) + ((-1/2)*q - cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) - b/(3*a)
    
    # y2[i] = w1*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w2*(((-1/2)*q - cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)
    
    y3 = w2*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w1*(((-1/2)*q-cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)

    D = y3.real + u0
    
    return D

def solve_T(T,popt,u0,T0,ut0):
    # T[0]=0.00001
    # 显式求解系数
    a = popt[10] + popt[6]*T + popt[2]*T**2       
    b = popt[9] + popt[5]*T + popt[1]*T**2
    c = popt[8] + popt[4]*T + popt[0]*T**2
    d = popt[7] + popt[3]*T + T**2 - (1/u0 + popt[0]*(T0/u0) + popt[1]*(T0**2/u0) + popt[2]*(T0**3/u0) + popt[3]*(1/u0**2)+ popt[4]*(T0/u0**2) + popt[5]*(T0**2/u0**2) + popt[6]*(T0**3/u0**2) + popt[7]*(1/u0**3) + popt[8]*(T0/u0**3) + popt[9]*(T0**2/u0**3) +popt[10]*(T0**3/u0**3))*T**3


    # 将方程转化为x3+rx=q的形式
    r = float((c/a - (b**2)/(3*a**2)))
    q = float((d/a + (2*b**3)/(27*a**3) - (b*c)/(3*a**2)))
    # 判别式Δ>0有一个实根；Δ<0有三个实根
    # 尝试显式计算 D
    
    k = 2j/2
    w1 = (1/2)*(-1 + sqrt(3)*k)
    w2 = (1/2)*(-1 - sqrt(3)*k)
    
    # 显式三次方程解法
    y1 = ((-1/2)*q + cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) + ((-1/2)*q - cmath.sqrt((q/2)**2 + (r/3)**3))**(1/3) - b/(3*a)
    
    y2 = w1*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w2*(((-1/2)*q - cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)
    
    y3 = w2*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w1*(((-1/2)*q-cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)

    D = y3.real + ut0

    return D

def main():
    T,D,p = get_data()
    with open('Pressure-Density\\Data\\original_p_D_data.json','r',encoding='utf-8') as f:
        data = json.load(f)
    # p =data["参考压力(MPa)"]
    Pc = data["临界压力(MPa)"]
    with open('Temperature-Density\\Data\\transfer_data.json','r',encoding='utf-8') as f:
            data = json.load(f)
    T_Dt0 = data['密度(kg/m3)']   
    with open('Temperature-Density\\Data\\D0_coefficient.json','r',encoding='utf-8') as f:
            data = json.load(f)
    popt_T_D = np.array(data['T-D'])
    T0 = data['T0']
    # with open('Temperature-Density\\Data\\original_T_D_data.json','r',encoding='utf-8') as f:
    #         data = json.load(f)
    # p = data["参考压力(MPa)"]
    # u0 = data['u0']  
    if p > Pc:
        mark = 1
        p = p - Pc
    else: 
        mark = 2
        p = p - 2
    popt_p_D,p_Dt0 = select(mark)
    # T,popt_D,D0,T0,Dt0
    
    D0 = solve_P(p,popt_p_D,p_Dt0)
    for i in range(len(T)):
        D[i] = solve_T(T[i]-T0,popt_T_D,D0,T0,T_Dt0)

    # print('方程系数为：',popt_T_D)
    print('该压力下对应密度为',D)
    plt.plot(T,D)
    plt.show()
    
        

if __name__ == '__main__':
    main()
    # print(__name__)