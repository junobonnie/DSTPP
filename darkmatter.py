# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 17:58:38 2023

@author: replica
"""

from vectortools import *
from atom import *
import random
import math
import sys
import units
limit_number = 15000
sys.setrecursionlimit(limit_number)

random.seed(42)
softening_length = 5#5.5
theta_ = 0.5
class Atom(Atom):
    def __init__(self, element, pos, vel = Vector(0, 0)):
        self.element = element
        self.pos = pos
        self.vel = vel
        
    def kinetic_energy(self):
        return (1/2)*self.element.mass*self.vel.dot(self.vel)
        
    def potential_energy(self, other):
        r = self.pos - other.pos
        if not self == other:# and r.dot(r) > self.element.radius+other.element.radius:
            return -self.element.mass*other.element.mass/m.sqrt(r.dot(r)+softening_length**2)
        else:
            return 0
        
    def fusion(self, other_Atom):
        new_Atom = None
        if not self == other_Atom:
            d = self.pos - other_Atom.pos
            if (d.dot(d) < (self.element.radius + other_Atom.element.radius)**2):
                new_element = Element(name = 'New atom', mass = self.element.mass + other_Atom.element.mass, 
                                      radius = m.sqrt(self.element.radius**2 + other_Atom.element.radius**2),
                                      color = self.element.color + other_Atom.element.color)
                new_Atom = Atom(element = new_element, 
                                pos = (self.element.mass*self.pos + other_Atom.element.mass*other_Atom.pos)/(self.element.mass + other_Atom.element.mass),
                                vel = (self.element.mass*self.vel + other_Atom.element.mass*other_Atom.vel)/(self.element.mass + other_Atom.element.mass))
        return new_Atom
    
class World(World):
    def __init__(self, t, atoms, walls, G, gravity = Vector(0, 0)):
        self.t = t
        self.atoms = atoms
        self.walls = walls
        self.G = G
        self.gravity = gravity
            
class Simulator(Simulator):
    def __init__(self, dt, world, render, grid_size = 100):
        self.dt = dt
        self.world = world
        self.render = render
        self.count_screen = 0
        self.count_snapshot = 0
        self.grid_size = grid_size
        self.grid = None
        
    # Function to build the Barnes-Hut tree
    def build_tree(self, atoms, x, y, width, height):
        # draw box
        draw_box = False
        if draw_box:
            positions = [Vector(x, y),
                         Vector(x+width, y),
                         Vector(x+width, y+height),
                         Vector(x, y+height)]
            self.render.polygon(positions, 'red')
        
        if len(atoms) == 0:
            return None
    
        if len(atoms) == 1:
            return atoms[0]
    
        # Calculate center of mass and total mass for the combined atoms
        total_mass = 0
        center = Vector(0, 0)
        for atom in atoms:
            total_mass += atom.element.mass
            center += atom.element.mass*atom.pos
        center /= total_mass
    
        # Divide the atoms into quadrants
        nw_atoms, ne_atoms, sw_atoms, se_atoms = [], [], [], []
        for atom in atoms:
            if atom.pos.x < x+width/2:
                if atom.pos.y < y+height/2:
                    nw_atoms.append(atom)
                else:
                    sw_atoms.append(atom)
            else:
                if atom.pos.y < y+height/2:
                    ne_atoms.append(atom)
                else:
                    se_atoms.append(atom)
    
        # Recursively build the tree
        nw = self.build_tree(nw_atoms, x, y, width / 2, height / 2)
        ne = self.build_tree(ne_atoms, x + width / 2, y, width / 2, height / 2)
        sw = self.build_tree(sw_atoms, x, y + height / 2, width / 2, height / 2)
        se = self.build_tree(se_atoms, x + width / 2, y + height / 2, width / 2, height / 2)
    
        return [nw, ne, sw, se, center, total_mass, (width**2+height**2)**(1/2)]
 
    def calculate_force(self, atom1, atom2):
        r = atom1.pos - atom2.pos
        f = -self.world.G*atom2.element.mass*r/((r.dot(r)+softening_length**2)**(3/2))
        return f
  
    def calculate_net_force(self, atom, tree):
        f = Vector(0, 0)
        if tree is None:
            return f

        if isinstance(tree, Atom):
            if tree != atom:
                force = self.calculate_force(atom, tree)
                f += force
                
        else:
            d = abs(tree[4] - atom.pos)
            if  tree[6] < theta_*d:
                force = self.calculate_force(atom, Atom(pos = tree[4], element=Element('tree', tree[5], 1, 'red')))
                f += force
            else:
                f_nw = self.calculate_net_force(atom, tree[0])
                f_ne = self.calculate_net_force(atom, tree[1])
                f_sw = self.calculate_net_force(atom, tree[2])
                f_se = self.calculate_net_force(atom, tree[3])

                f = f_nw + f_ne + f_sw + f_se

        return f
    
    def delete_outer_atom(self):
        for atom in self.world.atoms:
            if not ((-self.render.width/2 < atom.pos.x < self.render.width/2) and 
                    (-self.render.height/2 < atom.pos.y < self.render.height/2)):
                self.world.atoms.remove(atom)

    def v_dep_collision(self):
        sigma_0 = 5e8*units.M_sun*2.17*(units.cm**2)/units.g
        w = 180*units.km/units.s
        self.make_grid()
        for atom in self.world.atoms:
            result = 0
            atoms = self.get_near_atoms(atom)
            for other_atom in atoms:
                #self.render.polygon([atom.pos, other_atom.pos], 'red')
                r = atom.pos - other_atom.pos
                if not atom == other_atom:# and r.dot(r) > self.element.radius+other.element.radius:
                    v = atom.vel - other_atom.vel
                    v_ = (atom.vel + other_atom.vel)/2
                    if abs(r) < 5 and random.random() < sigma_0/(math.pi*25)*((1+v.dot(v)/(w**2))**(-2)):
                        theta = 2*math.pi*random.random()
                        atom.vel = v_ + SO2(theta).dot(v/2)
                        other_atom.vel = v_ + SO2(theta).dot(-v/2)
                        #print('collision')
                
    def main(self):
        self.delete_outer_atom()
        self.v_dep_collision()
        tree = self.build_tree(self.world.atoms, -self.render.width/2, -self.render.height/2, self.render.width, self.render.height)
        x_ = []
        v_ = []
        for atom in self.world.atoms:
            force = self.calculate_net_force(atom, tree)
            new_v = atom.vel + force*self.dt + self.world.gravity*self.dt
            v_.append(new_v)
            x_.append(atom.pos + new_v*self.dt)
        
        count = 0
        for atom in self.world.atoms:
            atom.pos = x_[count]
            atom.vel = v_[count]
            count = count + 1
        
def recursive_safety(atoms):
    for atom in atoms:
        for other_atom in atoms:
            if not atom == other_atom:
                d = atom.pos - other_atom.pos
                if d.dot(d) < 0.01:
                    atoms.remove(other_atom)
                    
if __name__ == '__main__':
    DEBUG = False
        
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

    wall1 = Wall(1000, 50, 0, Vector(-500, -400), red)
    wall2 = Wall(50, 800, 0, Vector(-500, -400), blue)
    wall3 = Wall(50, 800, 0, Vector(450,-400), blue)
    wall4 = Wall(1000, 50, 0, Vector(-500, 350), blue)
    wall5 = Wall(100, 50, m.pi/4, Vector(-300, 0), blue)

    e1 = Element(name = 'Helium', mass = 5e8*units.M_sun, radius = 2, color = blue)
    e2 = Element(name = 'Uranium', mass = 10, radius = 2, color = red)
   
    # atom1 = Atom(e1, Vector(-200, 0), Vector(50, 0))
    # atom2 = Atom(e1, Vector(0, 0))
    # atom3 = Atom(e1, Vector(25, -10))
    # atom4 = Atom(e1, Vector(25, 10))
    # atom5 = Atom(e1, Vector(50, -20))
    # atom6 = Atom(e1, Vector(50, 0))
    # atom7 = Atom(e1, Vector(50, 20))

    walls = [] # [wall1, wall2, wall3, wall4]#, wall5]
    atoms = [] # [atom1, atom2, atom3, atom4, atom5, atom6, atom7]
    
    # atom1 = Atom(e1, Vector(-400, 0), Vector(500, 0))
    # atom2 = Atom(e1, Vector(400, 0), Vector(-500, 0))
    # atoms = [atom1, atom2]
    
    for i in range(5000):
        theta = random.random()*2*math.pi
        r_theta = random.random()*2*math.pi
        rV = SO2(theta).dot(Vector(random.randrange(0, 200, 2*e1.radius) ,0)) - 0*Vector(250, 0)
        atoms.append(Atom(e1, rV, SO2(theta).dot(Vector(5, 0)))) #abs((rV + 0*Vector(250, 0))/200)*10*3*SO2(theta).dot(Vector(0, 20))))
    
    # for i in range(10000):
    #     theta = r.random()*2*m.pi
    #     rV = SO2(theta).dot(Vector(r.randrange(0, 200, 2*e2.radius) ,0)) + Vector(250, 0)
    #     atoms.append(Atom(e2, rV, abs((rV - Vector(250, 0))/200)*10*3*SO2(theta).dot(Vector(0, 30))))
        
    recursive_safety(atoms)
    
    G = 6.67384e-11*(units.m**3)*(units.kg**-1)*(units.s**-2)
    gravity = Vector(0, 0)

    world = World(0, atoms, walls, G, gravity)

    simulator = Simulator(0.01, world, render, 5)
    if DEBUG:  
        t_list = []
        K_list = []
        P_list = []
        TOT_E_list = []
    
    while True:
        t = simulator.clock()
        simulator.draw_background(white)
        simulator.draw_grid(100)
        simulator.draw_wall()
        simulator.main()
        # simulator.atom_wall_collision()
        # simulator.atom_atom_collision()
        # simulator.atom_atom_fusion()
        simulator.draw_atom()

        # render.text('pos = (%.2f, %.2f)'%(atoms[0].pos.x, atoms[0].pos.y) , None, 30, Vector(atoms[0].pos.x -100, atoms[0].pos.y - 30), black)
        # render.text('vel = (%.2f, %.2f)'%(atoms[0].vel.x, atoms[0].vel.y) , None, 30, Vector(atoms[0].pos.x -100, atoms[0].pos.y - 50), black)

        # render.text('pos = (%.2f, %.2f)'%(atoms[50].pos.x, atoms[50].pos.y) , None, 30, Vector(atoms[50].pos.x -100, atoms[50].pos.y - 30), blue)
        # render.text('vel = (%.2f, %.2f)'%(atoms[50].vel.x, atoms[50].vel.y) , None, 30, Vector(atoms[50].pos.x -100, atoms[50].pos.y - 50), blue)
        
        if DEBUG:   
            K = 0
            P = 0
         
            for atom in atoms:
                K = K + atom.kinetic_energy()
                for other_atom in atoms:
                    P = P + world.G*atom.potential_energy(other_atom)
            P = P/2
                
            t_list.append(t)
            K_list.append(K)
            P_list.append(P)
            TOT_E_list.append(K+P)
            
            render.text('t = %.2f'%(t) , None, 30, Vector(-480, -270), red)
            render.text('K_E = %.2f'%(K) , None, 30, Vector(-480, -300), red)
            render.text('P_E = %.2f'%(P) , None, 30, Vector(-480, -330), red)
            render.text('TOT_E = %.2f'%(K+P) , None, 30, Vector(-480, -360), red)
    
        for event in pg.event.get():
            if event.type == pg.QUIT:
                sys.exit()
        clock.tick(100)
        pg.display.update()
        
        # you need 'images/Barnes_Hut_demo_1' directory path
        #simulator.save_screen('images/v-indep-BH_demo_1')
        #simulator.save_snapshot('snapshots/v-indep-BH_demo_1')
        
        if t > 5:
            break
        
    if DEBUG:      
        import matplotlib.pyplot as plt
        
        plt.figure(figsize = (10,10))
        plt.plot(t_list, K_list, color='blue', label = 'Kinetic energy')
        plt.plot(t_list, P_list, color='orange', label = 'Potential energy')
        plt.plot(t_list, TOT_E_list, label = 'Total energy')
        plt.xlabel('time')
        plt.ylabel('energy')
        plt.axhline(sum(K_list)/len(t_list), 0, max(t_list), color='blue', linestyle='--', linewidth='1', label = 'Kinetic energy avg')
        plt.axhline(sum(P_list)/len(t_list), 0, max(t_list), color='orange', linestyle='--', linewidth='1', label = 'Potential energy avg')
        plt.legend(loc = 'best')
        plt.savefig('energy.png')
        plt.show()
        plt.close()
        
        mass_list = []
        kinetic_energy_list = []
        speed_list = []
        e1_speed_list = []
        e2_speed_list = []
        r_list = []
        e1_r_list = []
        e2_r_list = []
        
        for atom in atoms:
            mass_list.append(atom.element.mass)
            kinetic_energy_list.append(atom.kinetic_energy())
            speed_list.append(abs(atom.vel))
            r_list.append(abs(atom.pos))
            
            if atom.element.name == e1.name:
                e1_speed_list.append(abs(atom.vel))
                e1_r_list.append(abs(atom.pos))
                
            elif atom.element.name == e2.name:
                e2_speed_list.append(abs(atom.vel))
                e2_r_list.append(abs(atom.pos))
            
        plt.hist(mass_list, bins = 50)
        plt.xlabel('mass')
        plt.savefig('mass.png')
        plt.show()
        plt.close()
        
        plt.hist(kinetic_energy_list, bins = 50)
        plt.xlabel('kinetic energy')
        plt.savefig('kinetic.png')
        plt.show()
        plt.close()
        
        plt.hist(speed_list, bins = 50, label = 'Total', alpha = 0.5)
        plt.hist(e1_speed_list, bins = 50, label = 'e1', alpha = 0.5)
        plt.hist(e2_speed_list, bins = 50, label = 'e2', alpha = 0.5)
        plt.xlabel('speed')
        plt.legend(loc = 'best')
        plt.savefig('speed.png')
        plt.show()
        plt.close()
        
        plt.hist(r_list, bins = 50, label = 'Total', alpha = 0.5)
        plt.hist(e1_r_list, bins = 50, label = 'e1', alpha = 0.5)
        plt.hist(e2_r_list, bins = 50, label = 'e2', alpha = 0.5)
        plt.xlabel('distance')
        plt.legend(loc = 'best')
        plt.savefig('distance.png')
        plt.show()
        plt.close()
