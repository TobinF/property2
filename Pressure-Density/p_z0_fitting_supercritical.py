# Standard library imports
from scipy.optimize import curve_fit
import matplotlib.pylab as plt
import numpy as np
import json



# 隐函数拟合形式
# def implicit_func(X,a1,a2,a3,a4,a5,a6,a7,a8):
#     u,v = X
#     return (u**3) + a1*(u**2)*v + a2*u*(v**2) + a3*(v**3) + a4*(u**2) + a5*u*v + a6*(v**2) + a7*u + a8*v 
def implicit_func(X,a1,a2,a3,a4,a5,a6,a7,a8):
    u,v = X
    f = (u**3) + a1*(u**2)*v + a2*u*(v**2) + a3*(v**3) + a4*(u**2) + a5*u*v + a6*(v**2) + a7*u + a8*v
    return f
    
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