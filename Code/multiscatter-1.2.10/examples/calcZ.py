from math import *

lamb=0.00318979*1e3
kext1=0.0003766
S=15.9156
back1=kext1/S
Z1=log10(back1*lamb**4/pi**5*10.**6/0.176*pi)
salb1=0.9835
g1=0.2259
print Z1

kext2=0.00130145 
S=23.3279
back2=kext2/S
g2=0.0559
salb2=0.41738
Z2=log10(back2*lamb**4/pi**5*10.**6/0.176*pi)
print back1,back2, log10(back2/back1)*10.
print Z1,Z2

print 0.00113203*715*2
print 0.0
theta=0.00113203*180/pi
for i in range(0,80):
    z=(i+0.5)*.250
    if z<4:
        print kext2, salb2, g2, Z2*10.
    else:
        if z<8:
            print kext1, salb1, g1, Z1*10.
        else:
            print 0.1*kext2, 0., 0., -99.
print theta*2
