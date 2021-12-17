from math import sqrt
from matplotlib.pyplot import plot
import numpy as np
import matplotlib.pylab as plt
import cmath
import json

def get_data():
    p = input('请输入压力（MPa）:')
    return float(p)
  
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

def solve_z(p,popt,u0,filepath):
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
    # transfer_D(u0,y3,filepath)
        # y1,y2,y3 = solve_z()
    # 需要修改成对应转换关系
    # data = {
    # # 'y1':((y1.real + u0)).tolist(),
    # # 'y2':((y2.real + u0)).tolist(),
    # 'D':((y3.real + u0)).tolist(),
    # }
    # json_data = json.dumps(data,ensure_ascii=False)
    # with open(filepath,'w',encoding='utf-8') as f:       
    #         f.write(json_data)
    return D

def transfer_D(u0,y,filepath):
    # y1,y2,y3 = solve_z()
    # 需要修改成对应转换关系
    data = {
    # 'y1':((y1.real + u0)).tolist(),
    # 'y2':((y2.real + u0)).tolist(),
    'D':((y.real + u0)).tolist(),
    }
    json_data = json.dumps(data,ensure_ascii=False)
    with open(filepath,'w',encoding='utf-8') as f:       
            f.write(json_data)
    return y.real + u0

def main():
    p = get_data() - 2
    with open('Pressure-Density\\Data\\original_p_D_data.json','r',encoding='utf-8') as f:
        data = json.load(f)
    Pc = data["临界压力(MPa)"]

    if p > Pc:mark = 1
    else: mark = 2

    popt_D,D0 = select(mark)

    D = solve_z(p,popt_D,D0,'Pressure-Density\\Data\\fitting_D.json')

    print('方程系数为：',popt_D)
    print('该压力下对应密度为',D)
    # select()
    
        

if __name__ == '__main__':
    main()
    # print(__name__)