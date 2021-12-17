import json
import numpy as np

def load_data():
    with open('Pressure-Density\\Data\\original_p_D_data.json','r',encoding='utf-8') as f:
       data = json.load(f)
    p = np.array(data['压力(MPa)'])
    D = np.array(data['密度(kg/m3)'])
    # Cp = np.array(data['比热容[J/(kg·K)]'])
    # eta = np.array(data['粘度[microPa.s (10^-6 Pa.s)]'])
    # tcx = np.array(data['导热系数[W/(m·K)]'])
    Tc = data['临界温度(K)']
    Pc = data['临界压力(MPa)']
    
    return p,D,Tc,Pc

def stronge(p,D,Tc,Pc,filepath):
    data = {
    '压力(MPa)': p.tolist(),
    '密度(kg/m3)':D.tolist(),
    # '比热容[J/(kg·K)]':Cp.tolist(),
    # '粘度[microPa.s (10^-6 Pa.s)]':eta.tolist(),
    # '导热系数[W/(m·K)]':tcx.tolist(),
    '临界温度(K)':Tc,
    '临界压力(MPa)':Pc
    }
    json_data = json.dumps(data,ensure_ascii=False)
    with open(filepath,'w',encoding='utf-8') as f:       
            f.write(json_data)
    
    pass


# data\\p_data\\z0_data.json
def main():
    p,D,Tc,Pc = load_data()
    if p.min() > Pc:
        stronge(p,D,Tc,Pc,'Pressure-Density\\Data\\supercritical.json')
        pass
    elif p.min() <= Pc:
        for i in range(len(p)):
            if p[i] > Pc:
                n = i    
                break

        p1 = p[0:n]
        p2 = p[n:]
        D1 = D[0:n]
        D2 = D[n:]
        # Cp1 = Cp[0:n]
        # Cp2 = Cp[n:]
        # eta1 = eta[0:n]
        # eta2 = eta[n:]
        # tcx1 = tcx[0:n]
        # tcx2 = tcx[n:]  

        # 下标1 表示过热区域；下标2 表示超临界区域；
        stronge(p1,D1,Tc,Pc,'Pressure-Density\\Data\\superheat.json')
        stronge(p2,D2,Tc,Pc,'Pressure-Density\\Data\\supercritical.json')
        
    return 0


if __name__ == '__main__':
    main()
    # print(__name__)
