import coolprop as rp

# import refprop
#定义组分
rp.setup('def','CO2','nitrogen')
x=[0.5,0.5] #定义组分，摩尔分数
#prop=rp.xmole(x) #质量分数转摩尔分数
prop=rp.xmass(x) #摩尔分数转质量分数
print(prop)
x=prop['x'] #质量分数
wmix=prop['wmix'] #摩尔质量
#基本物性，t,p,D,cp,cv,wm
prop=rp.flsh('tp',250,220*wmix,x)
print(prop)
#传输特性，ets，tcx
prop=rp.trnprp(t=250,D=prop['D'],x=x)
print(prop)