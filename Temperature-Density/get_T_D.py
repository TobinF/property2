# Standard library imports
import os
import numpy as np
import ctREFPROP.ctREFPROP as ct
import json
import matplotlib.pylab as plt



# 加载64位的refprop  dll文件
# r = ct.REFPROPFunctionLibrary(os.environ['RPPREFIX'],'dll')  #自动加载，需要配置环境变量RPPREFIX
r = ct.REFPROPFunctionLibrary('D:\\Others\\Tools\\Refprop\\REFPRP64.DLL', 'dll')  #需要有64位dll文件
r.SETPATHdll(os.environ['RPPREFIX'])

# 定义初试数据
def init_data():
    p0 = 2*1000
    
    # 组分比例
    mix_ratio = [0.2,0.8] 
    # 节点数
    num = 200
    T_start = 213.15 # 2MPa
    T_end = 313.15 # 40MPa
    T_range = np.linspace(T_start,T_end,num)

    # 存放各物性参数
    wm = np.zeros(num) # 摩尔质量比
    D = np.zeros(num) # 密度
    return p0,T_range,mix_ratio,wm,D

# 定义get_z0函数获取T0下不同压力下的物性

def get_z0(p0,T,mix_ratio,wm,D):
    ''''
    参数：
        T0:一大于临界温度的温度值
        p:一数组，包含所有要获取物性的压力节点
        mix_ratio:混合气体、液体的摩尔比
        wm,D,Cp,eta,tcx：分别为摩尔质量比、密度、比热容、粘度、导热系数
    '''
    # 使用循环读取T0下对应的热物性
    # 指定fulid,配置摩尔组分
    
    r.SETUPdll(len(mix_ratio), 'HYDROGEN.FLD|METHANE.FLD', 'HMX.BNC', 'DEF')
    r.SETREFdll("DEF",1,mix_ratio,0,0,0,0)
    # 获取临界参数
    Tc,Pc,Dc,ierr,trim = r.CRITPdll(mix_ratio)
    # 存放数据
    # stor_z0(wm,D,Cp,eta,tcx,len(p))

    for i in range(len(T)):
        # 通过REFPROP获取物性
        D[i], Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = r.TPFLSHdll(T[i],p0,mix_ratio)
        wm[i] = r.WMOLdll(mix_ratio)
        eta, tcx, ierr, herr = r.TRNPRPdll(T[i],D[i],x)
    data_storage(T,p0,wm,D,Tc,Pc)
    # return wm,D,Cp,eta,tcx       

def data_storage(T,p0,wm,D,Tc,Pc):
    data = {
        '温度(K)': (T).tolist(),
        '密度(kg/m3)':(D*wm).tolist(),
        # '比热容[J/(kg·K)]':(Cp/wm).tolist(),
        # '粘度[microPa.s (10^-6 Pa.s)]':(eta).tolist(),
        # '导热系数[W/(m·K)]':tcx.tolist(),
        '参考压力(MPa)':p0/1000,
        '临界温度(K)':Tc,
        '临界压力(MPa)':Pc/1000
        }
    # plt.plot(T,D*wm)
    # plt.show()    
    json_data = json.dumps(data,ensure_ascii=False)
    with open('Temperature-Density\\Data\\original_T_D_data.json','w',encoding='utf-8') as f:       
            f.write(json_data)
    

# 用以获取在T0下的物性
def main():
    p0,T_range,mix_ratio,wm,D=init_data()
    # wm,D,Cp,eta,tcx,Tc,Pc 
    get_z0(p0,T_range,mix_ratio,wm,D) 
    # 存储原始数据
    # data_storage(p_range,wm,D,Cp,eta,tcx,Tc,Pc) 

    return 0


if __name__ == '__main__':
    main()
    # print(__name__)

