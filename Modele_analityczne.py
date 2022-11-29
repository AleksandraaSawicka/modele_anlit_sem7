import math
import numpy as np
import pyttsx3


class ConstParam:
    def __init__(self, glebokosc, natezenie, opornosc1, opornosc2, r_kuli):
        self.glebokosc = glebokosc
        self.natezenie = natezenie
        self.opornosc1 = opornosc1
        self.opornosc2 = opornosc2
        self.r_kuli = r_kuli

    def __str__(self):
        return 'Glebokosc:' + str(self.glebokosc) + '\tNatezenie:' + str(self.natezenie) +\
               '\tOpornosc1:' + str(self.opornosc1) + '\tOpornosc2:' + str(self.opornosc2) +\
               '\tPromien kuli:' + str(self.r_kuli)


def geometry(param, x, a):
    h = param.glebokosc
    b = param.r_kuli

    da = math.sqrt(x**2+(h+b)**2)
    dm = math.sqrt((x+a)**2+(h+b)**2)
    dn = math.sqrt((x+2*a)**2+(h+b)**2)
    db = math.sqrt((x+3*a)**2+(h+b)**2)

    theta_am = (da**2+dm**2-a**2)/(2*da*dm)
    theta_an = (da**2+dn**2-(2*a)**2)/(2*da*dn)
    theta_bm = (db**2+dm**2-(2*a)**2)/(2*db*dm)
    theta_bn = (db**2+dn**2-a**2)/(2*db*dn)

    return x, a, da, db, dm, dn, theta_am, theta_an, theta_bm, theta_bn


def legendre(number, f):
    pn = np.zeros(number)
    pn[0] = 1
    pn[1] = f
    for n in range(2, number):
        pn[n] = 1/n * (2*(n-1)+1)*pn[n-1] - (n-1)*pn[n-2]
    return pn


def coef_bn(param, number):
    h = param.glebokosc
    b = param.r_kuli
    I = param.natezenie
    ro1 = param.opornosc1
    ro2 = param.opornosc2

    bn = np.zeros(number)
    for n in range(number):
        bn[n] = I*n/(2*math.pi) * b**(2*n+1)/d**(n+1) * ro1*(ro2-ro1)/(n*ro1+(n+1)*ro2)
    return bn


def potential(param, geom, number, f):
    h = param.glebokosc
    b = param.r_kuli
    I = param.natezenie
    ro1 = param.opornosc1
    ro2 = param.opornosc2
    # h, I, ro1, ro2, b = param.values()

    x, a, da, db, dm, dn, theta_am, theta_an, theta_bm, theta_bn = geom
    pn = legendre(number, f)
    bn = coef_bn(param, number)
    v1 = 0
    # Potencjał od elektrody A(+), dm!=da
    for n in range(pn.size):
        if dm < da:
            v1 += (I*ro1/(2*math.pi) * dm**n/da**(n+1) + bn*dm**(-n-1))*pn
        elif dm > da:
            v1 += (I*ro1/(2*math.pi) * da**n/dm**(n+1) + bn*dm**(-n-1))*pn
    # Potencjał od elektrody B(-)
    for n in range(pn.size):
        if dm < db:
            v1 += (-I*ro1/(2*math.pi) * dm**n/db**(n+1) + bn*dm**(-n-1))*pn
        elif dm > db:
            v1 += (-I*ro1/(2*math.pi) * db**n/dm**(n+1) + bn*dm**(-n-1))*pn
    return v1


def change(param, x, dx, pn, bn):
    m = x/dx
    v = np.zeros(m)
    for i in range(m):
        v[i] = pot_zew(param, x+i*dx, dx, pn, bn)
    return v


# engine = pyttsx3.init()
# engine.say("Proszę podać parametry")
# engine.runAndWait()

# parametry = stale_param()
# print(parametry)
# print(parametry.values())

param = ConstParam(10, 2, 20, 50, 1)
print(param)

x, a, dA, dB, dM, dN, cos_am, cos_an, cos_bm, cos_bn = geom = geometry(param, -30, 1)
print(geom)

# def const_param():
#     glebokosc = input("Podaj glebokosc: ")
#     natezenie = input("Podaj nateznie: ")
#     opornosc1 = input("Podaj opornosc 1: ")
#     opornosc2 = input("Podaj opornosc 2: ")
#     r_kuli = input("Podaj promien kuli: ")
#     param = {'Glebokosc': glebokosc, 'Natezenie': natezenie,
#              'Opornosc1': opornosc1, 'Opornosc2': opornosc2, 'Promien kuli': r_kuli}
#     return param
