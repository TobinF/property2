from math import sqrt
from matplotlib.pyplot import plot
import numpy as np
import matplotlib.pylab as plt
import json

#解决中文显示问题
plt.rcParams['font.sans-serif']=['SimHei']
plt.rcParams['axes.unicode_minus'] = False

def error(z0,z):
       e = z0-z
       e0= e/z0
       # print(abs(e).max())
       print(abs(e0).max())
       return abs(e0)

def plot_fitting(p,U,u1,u2,u3,mark):
       # 初始数据
       e1 = error(U,u1)
       e2 = error(U,u2)
       e3 = error(U,u3)
       
       
       plt.subplot(2,2,1)
       plt.plot(p,e1,'r')
       plt.plot(p,e2,'g')
       plt.plot(p,e3,'b')
       plt.xlabel('p',position=(1,0))
       plt.ylabel('e',rotation=0,position=(0,1))
       plt.title('误差')

       plt.subplot(2,2,2)
       plt.plot(p,U,'r')
       plt.plot(p,u1,'g')
       plt.xlabel('p',position=(1,0))
       plt.ylabel('y1',rotation=0,position=(0,1))
       plt.title(mark+'2')

       # plt.show()
       plt.subplot(2,2,3)
       plt.plot(p,U,'r')
       plt.plot(p,u2,'g')
       plt.xlabel('p',position=(1,0))
       plt.ylabel('y2',rotation=0,position=(0,1))
       plt.title(mark+'3')

       plt.subplot(2,2,4)
       plt.plot(p,U,'r')
       plt.plot(p,u3,'g')
       plt.xlabel('p',position=(1,0))
       plt.ylabel('y3',rotation=0,position=(0,1))
       plt.title(mark+'4')
       
       plt.tight_layout()
       plt.show()

def load_fitting_data(filepath):
       with open(filepath,'r',encoding='utf-8') as f:
              data = json.load(f)
       U1 = np.array(data['y1'])
       U2 = np.array(data['y2'])
       U3 = np.array(data['y3'])
       return U1,U2,U3

# with open('data\\p_data\\supercritical.json','r',encoding='utf-8') as f:
with open('Pressure-Density\\Data\\superheat.json','r',encoding='utf-8') as f:
              data = json.load(f)
p = np.array(data['压力(MPa)'])
D = np.array(data['密度(kg/m3)'])
# Cp = np.array(data['比热容[J/(kg·K)]'])
# eta = np.array(data['粘度[microPa.s (10^-6 Pa.s)]'])
# tcx = np.array(data['导热系数[W/(m·K)]'])


# d = load_fitting_data('data\\p_data\\u_fitting_D.json')
# D1,D2,D3= load_fitting_data('data\\p_data\\u_fitting_supercritical_D.json')
# Cp1,Cp2,Cp3 = load_fitting_data('data\\p_data\\u_fitting_supercritical_Cp.json')
# tcx1,tcx2,tcx3= load_fitting_data('data\\p_data\\u_fitting_supercritical_tcx.json')
# eta1,eta2,eta3= load_fitting_data('data\\p_data\\u_fitting_supercritical_eta.json')
D1,D2,D3= load_fitting_data('Pressure-Density\\Data\\u_fitting_D.json')
# Cp1,Cp2,Cp3 = load_fitting_data('data\\p_data\\u_fitting_Cp.json')
# tcx1,tcx2,tcx3= load_fitting_data('data\\p_data\\u_fitting_tcx.json')
# eta1,eta2,eta3= load_fitting_data('data\\p_data\\u_fitting_eta.json')
# etaf = load_fitting_data('data\\p_data\\u_fitting_D.json')
# tcxf = load_fitting_data('data\\p_data\\u_fitting_D.json')

# plot_error(p,D,d)
plot_fitting(p,D,D1,D2,D3,'密度')
# plot_fitting(p,Cp,Cp1,Cp2,Cp3,'比热容')
# plot_fitting(p,eta,eta1,eta2,eta3,'粘度')
# plot_fitting(p,tcx,tcx1,tcx2,tcx3,'导热系数')

# plot_error(p,tcx,tcxf)
