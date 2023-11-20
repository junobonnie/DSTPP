import numpy as np
import matplotlib.pyplot as plt
from atom import *

class CIC:
    def __init__(self, atoms, grid_size, width, height):
        self.atoms = atoms
        self.grid_size = grid_size
        self.width = width
        self.height = height

    def density_map(self):
        x = np.linspace(-self.width/2, self.width/2, self.width//self.grid_size+1)
        y = np.linspace(-self.height/2, self.height/2, self.height//self.grid_size+1)
        density = np.zeros((len(y), len(x))) + 1
        for atom in self.atoms:
            i = np.searchsorted(x, atom.pos.x, side='right') - 1
            j = np.searchsorted(y, atom.pos.y, side='right') - 1
            if (0 <= i < len(x)-1) and (0 <= j < len(y)-1):
                x0, x1 = x[i], x[i+1]
                y0, y1 = y[j], y[j+1]
                dx = (atom.pos.x - x0) / self.grid_size
                dy = (atom.pos.y - y0) / self.grid_size
                density[j][i] += 5 * (1-dx) * (1-dy)
                density[j][i+1] += 5 * dx * (1-dy)
                density[j+1][i] += 5 * (1-dx) * dy
                density[j+1][i+1] += 5 * dx * dy
        return density

    def draw_density_map(self,path, i):
        density = np.log(self.density_map())
        plt.figure(figsize=(self.width/self.grid_size/4, self.height/self.grid_size/4))
        plt.imshow(density, cmap = 'hot')
        plt.savefig(path + '/%08d.png'%(i), bbox_inches='tight')
        plt.close()
        
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
    
    
    
    for i in range(0, 501):
        simulator.load_snapshot('snapshots/v-indep-BH_demo_1/snapshot_%08d.txt'%(i))
        atoms = simulator.world.atoms
        cic = CIC(atoms, 5, width, height)
        cic.draw_density_map('images/v-indep-BH_density_map', i)
        