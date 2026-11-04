# PC-SAFT (Gross & Sadowski) yogunluk hesabi: su + 1-propanol, 2B asosiyasyon
# Kullanim: python3 pcsaft_density.py   (T, P, x asagida tanimli)
import numpy as np
def brentq(f,a,b):
    fa=f(a)
    for _ in range(60):
        c=(a+b)/2; fc=f(c)
        if fa*fc<=0: b=c
        else: a,fa=c,fc
    return (a+b)/2
kB=1.380649e-23; NA=6.02214076e23
a0=[0.9105631445,0.6361281449,2.6861347891,-26.547362491,97.759208784,-159.59154087,91.297774084]
a1=[-0.3084016918,0.1860531159,-2.5030047259,21.419793629,-65.255885330,83.318680481,-33.746922930]
a2=[-0.0906148351,0.4527842806,0.5962700728,-1.7241829131,-4.1302112531,13.776631870,-8.6728470368]
b0=[0.7240946941,2.2382791861,-4.0025849485,-21.003576815,26.855641363,206.55133841,-355.60235612]
b1=[-0.5755498075,0.6995095521,3.8925673390,-17.215471648,192.67226447,-161.82646165,-165.20769346]
b2=[0.0976883116,-0.2557574982,-9.1558561530,20.642075974,-38.804430052,93.626774077,-29.666905585]
# water(2B), 1-propanol(2B)  Gross & Sadowski 2002
m=np.array([1.0656,2.9997]); s=np.array([3.0007,3.2522]); e=np.array([366.51,233.40])
eAB=np.array([2500.7,2276.8]); kAB=np.array([0.034868,0.015268]); M=np.array([18.015,60.096])
x=np.array([0.5,0.5]); T=360+273.15; P=100*101325.0; kij=0.0
def ares(rho):  # rho in 1/A^3
    d=s*(1-0.12*np.exp(-3*e/T))
    z=[np.pi/6*rho*np.sum(x*m*d**n) for n in range(4)]
    z0,z1,z2,z3=z; mb=np.sum(x*m)
    ahs=(3*z1*z2/(1-z3)+z2**3/(z3*(1-z3)**2)+(z2**3/z3**2-z0)*np.log(1-z3))/z0
    dij=np.outer(d,d)/np.add.outer(d,d)
    g=1/(1-z3)+dij*3*z2/(1-z3)**2+dij**2*2*z2**2/(1-z3)**3
    ahc=mb*ahs-np.sum(x*(m-1)*np.log(np.diag(g)))
    sij=(s[:,None]+s[None,:])/2; eij=np.sqrt(np.outer(e,e))*(1-kij*(1-np.eye(2)))
    xx=np.outer(x,x)*np.outer(m,m)
    m2e1=np.sum(xx*(eij/T)*sij**3); m2e2=np.sum(xx*(eij/T)**2*sij**3)
    eta=z3; ai=[a0[i]+(mb-1)/mb*a1[i]+(mb-1)*(mb-2)/mb**2*a2[i] for i in range(7)]
    bi=[b0[i]+(mb-1)/mb*b1[i]+(mb-1)*(mb-2)/mb**2*b2[i] for i in range(7)]
    I1=sum(ai[i]*eta**i for i in range(7)); I2=sum(bi[i]*eta**i for i in range(7))
    C1=1/(1+mb*(8*eta-2*eta**2)/(1-eta)**4+(1-mb)*(20*eta-27*eta**2+12*eta**3-2*eta**4)/((1-eta)*(2-eta))**2)
    adisp=-2*np.pi*rho*I1*m2e1-np.pi*rho*mb*C1*I2*m2e2
    eABij=(eAB[:,None]+eAB[None,:])/2
    kABij=np.sqrt(np.outer(kAB,kAB))*(np.sqrt(np.outer(s,s))/sij)**3
    D=sij**3*g*kABij*(np.exp(eABij/T)-1)
    X=np.ones(2)
    for _ in range(500):
        Xn=1/(1+rho*(D@(x*X))); 
        if np.max(abs(Xn-X))<1e-13: X=Xn;break
        X=0.5*X+0.5*Xn
    aas=np.sum(x*2*(np.log(X)-X/2+0.5))
    return ahc+adisp+aas, X
def Z(rho):
    h=rho*1e-6
    return 1+rho*(ares(rho+h)[0]-ares(rho-h)[0])/(2*h)
def Pcalc(rho): return Z(rho)*rho*1e30*kB*T
if __name__ == "__main__":
    d=s*(1-0.12*np.exp(-3*e/T)); rmax=0.74/(np.pi/6*np.sum(x*m*d**3))
    grid=np.linspace(1e-7,rmax*0.999,4000); f=[Pcalc(r)-P for r in grid]
    roots=[brentq(lambda r:Pcalc(r)-P,grid[i],grid[i+1]) for i in range(len(grid)-1) if f[i]*f[i+1]<0]
    Mm=np.sum(x*M)
    for r in roots:
        a,X=ares(r); z=Z(r); G=a+z-1-np.log(z)
        print(f"rho={r*1e30/NA/1e3:.4f} kmol/m3  {r*1e30/NA*Mm/1e6*1e3:.2f} kg/m3  Z={z:.4f}  eta={np.pi/6*r*np.sum(x*m*d**3):.4f}  X={X}  g_res/RT={G:.4f}")
