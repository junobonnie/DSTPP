import numpy as np
import matplotlib.pyplot as plt
import units
from atom import *

def binary_search(list_, data):
    low = 0
    high = len(list_) - 2
    while low <= high:
        mid = (low + high) // 2
        if list_[mid] < data < list_[mid+1]:
            return mid
        elif list_[mid] > data:
            high = mid - 1
        else:
            low = mid + 1
    return low

def pseudo_isothermal_profile(rho_0, R, r_):
    return rho_0/(1+(r_/R)**2)

def NFW_profile(rho_0, R, r_):
    return rho_0/((r_/R)*(1+r_/R)**2)
    

def center_of_mass(atoms):
    result = Vector(0, 0)
    for atom in atoms:
        result += atom.pos
    return result / len(atoms)

def draw_density_profile(atoms):
    cm = center_of_mass(atoms)
    bins = 100
    maximum_exponent = 2.6
    r = np.logspace(-1, maximum_exponent, bins)
    r_ = []
    
    for i in range(bins-1):
        r_.append((r[i]+r[i+1])/2)
    r_ = np.array(r_)
    rho = np.zeros((bins-1))
    
    for atom in atoms:
        d = abs(atom.pos - cm)
        if d < 10**maximum_exponent:
            index = binary_search(r, d)
            rho[index] += atom.element.mass
            
    yerr = 1/np.sqrt(rho/5e8)
    
    for i in range(bins-1):
        rho[i] /= np.pi*(r[i+1]**2-r[i]**2)
        
    yerr *= rho
    
    plt.figure(dpi=200)
    plt.plot(r_,pseudo_isothermal_profile(2e9, 6, r_), '--', label='Pseudo Isothermal')
    plt.plot(r_, NFW_profile(5e7, 110, r_), '--', label='NFW Profile')
    plt.errorbar(r_, rho, yerr=yerr , fmt='r+', label='Simulation')
    plt.legend()
    plt.xlabel('r[kpc](distance form center of mass)')
    plt.ylabel(r'$\rho[M_⊙/kpc^2]$')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
        
if __name__ == "__main__":

    width = 1000
    height = 800
    
    screen = pg.display.set_mode((width, height))
    render = Render(screen, width, height)
    clock = pg.time.Clock()
    
    black = pg.Color('black')
    white = pg.Color('white')
    red = pg.Color('red')
    green = pg.Color('green')
    blue = pg.Color('blue')
        
    walls = []
    atoms = []
    
    gravity = Vector(0, 0)
    world = World(0, atoms, walls, gravity)
    
    simulator = Simulator(0.01, world, render)
    
    for i in range(0, 501, 100):
        simulator.load_snapshot('snapshots/v-dep-BH_demo_1/snapshot_%08d.txt'%(i))
        atoms = simulator.world.atoms
        draw_density_profile(atoms)