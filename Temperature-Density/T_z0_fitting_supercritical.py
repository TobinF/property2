# Standard library imports
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import numpy as np
import json
import cmath
from math import sqrt

def get_data():
    # p = input('请输入压力（MPa）:')
    p = 2
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
    
    y3 = w2*(((-1/2)*q+cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) + w1*(((-1/2)*q-cmath.sqrt((q/2)**2+(r/3)**3))**(1/3)) - b/(3*a)

    D = y3.real + u0
    
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

def get_D0(p):
    # p = get_data()
    with open('Pressure-Density\\Data\\original_p_D_data.json','r',encoding='utf-8') as f:
        data = json.load(f)
    Pc = data["临界压力(MPa)"]

    if p > Pc:
        mark = 1
        p = p - Pc
    else: 
        mark = 2
        p = p - 2

    popt_D,D0 = select(mark)

    D = solve_z(p,popt_D,D0,'Pressure-Density\\Data\\fitting_D.json')

    print('方程系数为：',popt_D)
    print('该压力下对应密度为',D)
    # select()
    
        

# 隐函数拟合形式
def implicit_func(X,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10):
    u,T,mark = X
    with open('data\\T_data\\original_T_z0_data.json','r',encoding='utf-8') as f:
       data = json.load(f)
    u0 = get_D0(p)
    T0 = 273.15-80

    a = (1/u0 - 1/u) + a1*(T0/u0 - T/u) + a2*(T0**2/u0 - T**2/u) + a3*(T0**3/u0 - T**3/u) + a4*(1/u0**2 - 1/u**2) + a5*(T0/u0**2 - T/u**2) + a6*(T0**2/u0**2 - T**2/u**2) + a7*(T0**3/u0**2 - T**3/u**2) + a8*(1/u0**3 - 1/u**3) + a9*(T0/u0**3 - T/u**3) + a10*(T0**2/u0**3 - T**2/u**3) + a10*(T0**3/u0**3 - T**2/u**3)
    return a
    
# 读入原始数据
def load_z0():
    with open('Pressure-Density\\Data\\supercritical.json','r',encoding='utf-8') as f:
       data = json.load(f)
    p = np.array(data['压力(MPa)'])
    D = np.array(data['密度(kg/m3)'])
    # Cp = np.array(data['比热容[J/(kg·K)]'])
    # eta = np.array(data['粘度[microPa.s (10^-6 Pa.s)]'])
    # tcx = np.array(data['导热系数[W/(m·K)]'])
    Tc = data['临界温度(K)']
    Pc = data['临界压力(MPa)']
    return p,D,Tc,Pc

def trans_func():
    """
    转换函数，执行z*=f(z)\n
    做差：u = z* - z*0
    """ 
    p,D,Tc,Pc = load_z0()
    # 自变量p,无需转换，但是需要减去初值 v = p - p0，记作p_u
    p_u = p - p[0]
    # 对自变量z执行一个简单的转换
    # D_t = (D - 100) * 0.15
    D_t = D
    # Cp_t = (Cp + 200) * 10
    # eta_t = eta - 10
    # tcx_t = tcx*100
    # 存储转换后的物性参数初值，用以后边的还原计算
    udata_storage(p[0],D_t[0],'Pressure-Density\\Data\\transfer_supercritical_data.json')
    # 对自变量减去初值,作为隐式拟合模型中的u，记作z_u
    D_u = D_t - D_t[0]
    # Cp_u = Cp_t[0] - Cp_t
    # eta_u = eta_t - eta_t[0]
    # tcx_u = tcx_t - tcx_t[0]
    # 用于作为拟合节点
    # plt.plot(p,Cp_u)
    # plt.show()
    udata_storage(p_u,D_u,'Pressure-Density\\Data\\supercritical_fitting_nodes.json')

    return p_u,D_u

# 用最小二乘法进行拟合，获取获取参数 a1 - a8
def fitting_func(p,D):
    z = 0 # 隐式方程的解，及z=F(x,y)=0
    popt_D, pcov1 = curve_fit(implicit_func, (p, D), z)
    # popt_Cp, pcov2 = curve_fit(implicit_func, (p, Cp), z)
    # popt_eta, pcov3 = curve_fit(implicit_func, (p, eta), z)
    # popt_tcx, pcov4 = curve_fit(implicit_func, (p, tcx), z)
   
    # print('最大误差%.3f\n'%(max(abs(implicit_func((p, D),*popt_D)))))
    # print('最大误差%.3f\n'%(max(abs(implicit_func((p, Cp),*popt_Cp)))))
    

    return popt_D

def coefficient_storage(popt_D):
    data = {
        'p-D': popt_D.tolist(),
        # 'p-Cp':popt_Cp.tolist(),
        # 'p-eta':popt_eta.tolist(),
        # 'p-tcx':popt_tcx.tolist(),
        }
    json_data = json.dumps(data,ensure_ascii=False)
    with open('Pressure-Density\\Data\\D0_supercritical_coefficient.json','w',encoding='utf-8') as f:       
            f.write(json_data)
    pass

def udata_storage(p,D,filepath):
    data = {
        '压力(MPa)': p.tolist(),
        '密度(kg/m3)':D.tolist(),
        # '比热容[J/(kg·K)]':Cp.tolist(),
        # '粘度[microPa.s (10^-6 Pa.s)]':eta.tolist(),
        # '导热系数[W/(m·K)]':tcx.tolist(),
        }
    json_data = json.dumps(data,ensure_ascii=False)
    with open(filepath,'w',encoding='utf-8') as f:       
            f.write(json_data)
    pass

def main():
    p_u,D_u = trans_func()
    popt_D= fitting_func(p_u,D_u,)
    
    coefficient_storage(popt_D)
    
    return 0

if __name__ == '__main__':
    main()
    # print(__name__)