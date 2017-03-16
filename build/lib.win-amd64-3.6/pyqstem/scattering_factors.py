from ase.data import chemical_symbols
from scipy.special import kn
import numpy as np
import csv

def get_prefactor(units):
    if units=='SI':
        e=1.60217662*10**(-19)
        a0=0.52917721067
        epsilon0=8.854187817*10**(-12)
        prefactor=a0*e/(2*epsilon0)
    elif units=='ASE':
        e=1
        a0=0.52917721067
        epsilon0=0.07957747154594767
        prefactor=a0*e/(2*epsilon0)
    elif units=='QSTEM':
        prefactor=1
    else:
        raise RuntimeError('Units not recognized')

    return prefactor

def read_lobato():
    with open('lobato.csv', 'rb') as  csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        parameters={}
        for i,row in enumerate(reader):
            cs=chemical_symbols[i/2+1]
            if i%2==0:
                parameters[cs]={'a':[float(r) for r in row]}
            else:
                parameters[cs].update({'b':[float(r) for r in row]})
    return parameters

def fe_lobato(g,element):
    parameters=read_lobato()

    aa=parameters[element]['a']
    bb=parameters[element]['b']

    fe=np.sum([a*(2+b*g**2)/(1+b*g**2)**2 for a,b in zip(aa,bb)],axis=0)

    return fe

def V_lobato(r,element,units='QSTEM'):
    parameters=read_lobato()

    aa=parameters[element]['a']
    bb=parameters[element]['b']

    prefactor=get_prefactor(units)

    V=prefactor*np.sum([np.pi**2*a/b**(3/2.)*(b**(1/2.)/(np.pi*r)+1)*
                      np.exp(-2*np.pi*r/b**(1/2.)) for a,b in zip(aa,bb)],axis=0)

    return V

def V_projected_lobato(r,element,units='QSTEM'):
    parameters=read_lobato()

    aa=parameters[element]['a']
    bb=parameters[element]['b']

    prefactor=get_prefactor(units)

    V=prefactor*np.sum([2*np.pi**2*a/b**(3/2.)*(kn(0,2*np.pi*r/b**(1/2.))+
                                              r*kn(1,2*np.pi*r/b**(1/2.))) for a,b in zip(aa,bb)],axis=0)
    return V

def read_kirkland():
    with open('kirkland.csv', 'rb') as  csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        parameters={}
        for i in range(103):
            cs=chemical_symbols[i+1]
            parameters[cs]={}
            element_data=[]
            for j in range(3):
                element_data+=[float(r) for r in reader.next()]
            parameters[cs]['a']=element_data[0:6:2]
            parameters[cs]['b']=element_data[1:6:2]
            parameters[cs]['c']=element_data[6::2]
            parameters[cs]['d']=element_data[7::2]
    return parameters

def fe_kirkland(g,element):
    parameters=read_kirkland()

    aa=parameters[element]['a']
    bb=parameters[element]['b']
    cc=parameters[element]['c']
    dd=parameters[element]['d']

    fe=np.sum([a/(g**2+b) for a,b in zip(aa,bb)],axis=0)
    fe+=np.sum([c*np.exp(-d*g**2) for c,d in zip(cc,dd)],axis=0)

    return fe

def V_kirkland(r,element,units='QSTEM'):
    parameters=read_kirkland()

    aa=parameters[element]['a']
    bb=parameters[element]['b']
    cc=parameters[element]['c']
    dd=parameters[element]['d']

    prefactor=get_prefactor(units)

    V=prefactor*np.sum([np.pi*a/r*np.exp(-2*np.pi*r*np.sqrt(b)) for a,b in zip(aa,bb)],axis=0)
    V+=prefactor*np.sum([np.pi**(3/2.)*c/d**(3/2.)*np.exp(-np.pi**2*r**2/d) for c,d in zip(cc,dd)],axis=0)

    return V

def V_projected_kirkland(r,element,units='QSTEM'):
    parameters=read_kirkland()

    aa=parameters[element]['a']
    bb=parameters[element]['b']
    cc=parameters[element]['c']
    dd=parameters[element]['d']

    prefactor=get_prefactor(units)

    V=prefactor*np.sum([2*np.pi*a*kn(0,2*np.pi*r*np.sqrt(b)) for a,b in zip(aa,bb)],axis=0)
    V+=prefactor*np.sum([np.pi*c/d*np.exp(-np.pi**2*r**2/d) for c,d in zip(cc,dd)],axis=0)

    return V

def hydrogen_analytic_projected(r):
    e=1.60217662*10**(-19)
    a0=0.52917721067
    epsilon0=8.854187817*10**(-12)

    return e/(4*np.pi*epsilon0)*(2*kn(0,2.*r/a0)+2*r/a0*kn(1,2.*r/a0))
