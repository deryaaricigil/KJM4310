# Author: Audun Skau Hansen (a.s.hansen@kjemi.uio.no), 2022

import numpy as np

import matplotlib.pyplot as plt

square_lattice_basis = np.array([[ 1, 0],
                                 [-1, 0],
                                 [ 0, 1],
                                 [ 0,-1]], dtype = int)

# representations

class polymer():
    def __init__(self, config, conformation = None, basis = square_lattice_basis):
        self.config = config
        self.conformation = conformation
        self._edict = {"H":1, "P":0}
        self.basis = basis
        self.conformations = []
        
    def determine_conformations(self):
        """
        Explore and determine possible conformations
        
        **Note** The generated conformations contain a lot of 
                 redundancy due to symmetries.
                 However, we simply want a set which is compact 
                 enough for a somewhat efficient exploration of
                 small polymers.
        """
        self.advance_conformation([np.array([1,0])], [np.array([0,0]), np.array([1,0])])
        
        # generate transformed conformations
        transformed_conformations = []
        for i in self.conformations:
            transformed_conformations.append(i)
            transformed_conformations.append([np.array([j[0], -j[1]]) for j in i])
            transformed_conformations.append([np.array([-j[0], j[1]]) for j in i])
            transformed_conformations.append([np.array([-j[0], -j[1]]) for j in i])
            
            transformed_conformations.append([np.array([j[1],  j[0]]) for j in i])
            transformed_conformations.append([np.array([-j[1],  j[0]]) for j in i])
            transformed_conformations.append([np.array([j[1], -j[0]]) for j in i])
            transformed_conformations.append([np.array([-j[1], -j[0]]) for j in i])
            
        self.conformations = transformed_conformations

            
        
        
    def advance_conformation(self, conformation, points):
        """
        Recurrent function for discovering 
        possible conformations.
        """
        if len(conformation) >= len(self.config)-1:
            self.conformations.append(conformation)
        else:
            for i in self.basis:
                #if not np.all(i == -1*conformation[-1]):
                if np.sum((i + conformation[-1])**2)>0:
                    new_point = points[-1] + i
                    if not np.any(np.sum(np.abs(new_point - points), axis = 1)==0): # and np.sum(new_point**2)>0:
                        new_points = points + [new_point]
                        new_conformation = conformation + [i]
                        self.advance_conformation(new_conformation, new_points)

        
        
    def set_conformation(self, conformation):
        self.conformation = conformation
        
    def energy(self):
        # compute self energy of polymer
        e = np.array([self._edict[i] for i in self.config])
        return np.sum(np.abs(e[1:]-e[:-1]))
       
        
        


class lattice():
    def __init__(self, Nx = 10, Ny = 10):
        
        self.lattice = -1*np.ones((Nx,Ny), dtype = int) # initially populated with water
        
        self._edict = {"W":-1, "H":0, "P" :1}
        
        self.conformations = []
    
    def validate_placement(self, i,j,conf):
        """
        Only sites with water are replaced
        """
        if self.lattice[i,j]!=-1:
            return False
        di,dj = i,j
        for c in conf:
            di, dj = di+c[0], dj+c[1]
            if self.lattice[di,dj]!=-1 or di<0 or dj<0:
                return False
        
        return True
            
            
        
        
        
    def find_all_possible_placements(self, polym, uniques = True):
        """
        Returns all (unique) possible conformations of polym
        that fit in the lattice
        """
        config = [self._edict[k] for k in polym.config]
        
        valid_conformations = [] 
        
        unique_lattices = []
        for i in range(self.lattice.shape[0]):
            for j in range(self.lattice.shape[1]):
                if self.lattice[i,j] == -1:
                    for ci in range(len(polym.conformations)):
                        c = polym.conformations[ci]
                        if self.validate_placement(i,j,c):
                            conf = [np.array([i,j])] + c
                            conf = np.cumsum(np.array(conf), axis = 0)
                            
                            
                            # determine if conformation is unique                            
                            conf_unique = True
                            for k in valid_conformations:
                                if np.sum((k-conf)**2)==0:
                                    conf_unique = False
                            
                            if conf_unique:
                                valid_conformations.append(conf) #, self.place_polymer_at(i,j,polym, c)])
                                #unique_lattices.append(lat)
        return valid_conformations
                    

    def place_polymer_at(self, i,j,polym, conformation = None):
        """
        Returns a lattice where polymer is placed
        at coordinates i,j in lattice
        If no conformation is specified the function
        will chose the first available conformation
        which fits.
        
        """
        assert(i<self.lattice.shape[0]), "invalid placement"
        assert(j<self.lattice.shape[1]), "invalid placement"
        
        config = [self._edict[k] for k in polym.config]
        
        lattice = self.lattice*1
        
        if conformation is None:
            placed = False
            for c in range(len(polym.conformations)):
                pc = polym.conformations[c]
                if self.validate_placement(i,j,pc):
                    placed = True
                    print(pc)
                    # place polymer
                    di,dj = i,j
                    lattice[di,dj] = self._edict[polym.config[0]]
                    for k in range(len(polym.conformations[c])):
                        di += polym.conformations[c][k][0]
                        dj += polym.conformations[c][k][1]
                        lattice[di,dj] = self._edict[polym.config[k+1]]
                        
            if not placed:
                print("Unable to fit polymer in lattice")
        else:
            if self.validate_placement(i,j,conformation):
                di,dj = i,j
                lattice[di,dj] = self._edict[polym.config[0]]
                for k in range(len(conformation)):
                    di += conformation[k][0]
                    dj += conformation[k][1]
                    lattice[di,dj] = self._edict[polym.config[k+1]]
                    
            else:
                print("Conformation does not fit in lattice")
        return lattice
        
    
    def energy(self, lattice = None):
        """
        Computes the energy of the lattice 
        (if no lattice provided, the energy of self.lattice is computed)
        """
        if lattice is None:
            lattice = self.lattice
        energy = 0
        #for i in range(lattice.shape[0]-1):
        #    for j in range(lattice.shape[1]-1):
        #        if np.abs(lattice[i,j] - lattice[i,j+1]) == 1:
        #            energy += 1
        #        if np.abs(lattice[i,j] - lattice[i+1,j]) == 1:
        #            energy += 1
        return np.sum(np.abs(lattice[1:,:]  - lattice[:-1,:])==1) + np.sum(np.abs(lattice[:, 1:]  - lattice[:, :-1])==1)


    
    



# Visualization


def show_lattice_placement(l, p, c):

    """
    Show how the lattice l looks when polymer p is placed 
    as defined by positions c
    """



    if type(c) is not list:
        c = [c]
        
    n = int(np.sqrt(len(c)) + 1)
    
    counter = 0
    
    plt.figure(figsize=(2*l.lattice.shape[0], 2*l.lattice.shape[1]))
    
    dx, dy = l.lattice.shape[1]+3, l.lattice.shape[0]+3
    
    unique_lattices = [] #list to hold unique lattice configurations
    
    
    
    
    
    
    for i in range(n):
        for j in range(n):
            if counter<len(c):

                #compute energy of polymer
                energy_polymer = p.energy()

                #compute energy of empty cavity
                energy_cavity = l.energy() 


                pts = np.array(np.meshgrid(np.arange(l.lattice.shape[1])+dx*i, np.arange(l.lattice.shape[0])+dy*j)).reshape(2,-1).T
                lat = l.lattice*1

                config = [l._edict[k] for k in p.config]

                lat[c[counter][:,0], c[counter][:, 1]] = config

                # compute energy of filled cavity
                energy_filled = l.energy(lat)

                lat = lat.ravel()

                lw = lat==-1
                lh = lat== 0
                lp = lat== 1

                plt.plot(pts[lw,1], pts[lw,0], "o", markersize = 5, color = (.3,.3,.8))
                plt.plot(pts[lh,1], pts[lh,0], "o", markersize = 5, color = (.9,.9,.3))
                plt.plot(pts[lp,1], pts[lp,0], "o", markersize = 5, color = (.8,.4,.3))
                plt.plot(pts[lw,1], pts[lw,0], "o", markersize = 7, color = (0,0,0), zorder = -1)
                plt.plot(pts[lh,1], pts[lh,0], "o", markersize = 7, color = (0,0,0), zorder = -1)
                plt.plot(pts[lp,1], pts[lp,0], "o", markersize = 7, color = (0,0,0), zorder = -1)
                plt.plot(c[counter][:,0]+dy*j, c[counter][:,1]+dx*i, "-", color = (0,0,0), zorder = -1, linewidth = 2)


                #print("Polymer self-energy    :", energy_polymer)
                #print("Energy of empty cavity :", energy_cavity)
                #print("Energy of filled cavity:", energy_filled)
                #print("Total energy.          :", energy_filled - energy_cavity - energy_polymer)
                plt.text(dy*(j+.25), dx*i-1, "$\epsilon$ = %i" % (energy_filled - energy_cavity - energy_polymer), ha = "center", va = "center", fontsize = 8)
              
                plt.text(dy*(j+.25), dx*i-2, "Filled E = %i" % (energy_filled), ha = "center", va = "center", fontsize = 8)
                
                counter += 1

    plt.axis("off")
    
    #myedit 
    plt.savefig('exercise4.png')
    
    plt.show()



def show_conformations(polym):
    """
    Show all possible conformations of polym(er)
    """

    n = int(np.sqrt(len(polym.conformations))+1)
    sep = 6.5
    plt.figure(figsize = (10,10))
    c = 0
    lp = len(polym.conformations)
    hi = np.array([i for i in polym.config])=="H"
    pi = np.array([i for i in polym.config])=="P"
    
    for i in range(n):
        for j in range(n):
            if c<lp:
                conf = [np.array([0,0])] + polym.conformations[c]
                conf = np.cumsum(np.array(conf), axis = 0)

                conf = conf - np.mean(conf, axis = 0)[None, :]

                plt.plot(conf[:,0] + i*sep, conf[:,1]-j*sep, "-", color = (0,0,0))

                plt.plot(conf[hi,0] + i*sep, conf[hi,1]-j*sep, "o", color = (.8,.3,0), markersize = 2)
                plt.plot(conf[pi,0] + i*sep, conf[pi,1]-j*sep, "o", color = (0 ,.3,.8), markersize = 2)
                c +=1

    plt.xlim(-sep, sep*n+1)
    plt.ylim(-sep*n, sep)
    plt.axis("off")
    plt.show()
    
    
def remove_rotational_redundance(conf, remove_reversed = False):
    """
    Remove rotational redundancies from the set of conformations conf
    """
    nonredundant_set = []  #nonredundant conformations to be added here
    
    for c in conf:
        nonredundant = True # if nonredundant, add to nonredundant_set
        cc = np.array(c)
        
        # check for redundancies
        for m in nonredundant_set:
            mc = np.array(m)
            for j in range(4):
                if np.sum((mc-cc)**2) < 1e-10:
                    nonredundant = False
                cc = cc.dot(np.array([[0,1],[-1,0]])) #rotate polymer 90 degrees
            
            if remove_reversed:
                cc = cc[::-1] # reverse polymer
                for j in range(4):
                    if np.sum((mc-cc)**2) < 1e-10:
                        nonredundant = False
                    cc = cc.dot(np.array([[0,1],[-1,0]]))

                
        if nonredundant:
            nonredundant_set.append(c)
    return nonredundant_set


#my edit
def show_conformations_redundant(polym, conformations):
    """
    Show all possible conformations of polym(er)
    """

    n = int(np.sqrt(len(conformations))+1)
    sep = 6.5
    plt.figure(figsize = (10,10))
    c = 0
    lp = len(conformations)
    hi = np.array([i for i in polym.config])=="H"
    pi = np.array([i for i in polym.config])=="P"
    
    for i in range(n):
        for j in range(n):
            if c<lp:
                
                conf = [np.array([0,0])] + conformations[c]
                conf = np.cumsum(np.array(conf), axis = 0)

                conf = conf - np.mean(conf, axis = 0)[None, :]

                plt.plot(conf[:,0] + i*sep, conf[:,1]-j*sep, "-", color = (0,0,0))

                plt.plot(conf[hi,0] + i*sep, conf[hi,1]-j*sep, "o", color = (.8,.3,0), markersize = 2)
                plt.plot(conf[pi,0] + i*sep, conf[pi,1]-j*sep, "o", color = (0 ,.3,.8), markersize = 2)
                
                
                
                c +=1

    plt.xlim(-sep, sep*n+1)
    plt.ylim(-sep*n, sep)
    plt.axis("off")
    plt.savefig("Plot.png")
    plt.show()
   

def show_conformations_redundant_in_oil(polym, conformations):
    """
    Show all possible conformations of polym(er)
    """

    n = int(np.sqrt(len(conformations))+1)
    sep = 6.5
    plt.figure(figsize = (10,10))
    c = 0
    lp = len(conformations)
    hi = np.array([i for i in polym.config])=="H"
    pi = np.array([i for i in polym.config])=="P"
    
    for i in range(n):
        for j in range(n):
            if c<lp:
                
                conf = [np.array([0,0])] + conformations[c]
                conf = np.cumsum(np.array(conf), axis = 0)

                conf = conf - np.mean(conf, axis = 0)[None, :]

                plt.plot(conf[:,0] + i*sep, conf[:,1]-j*sep, "-", color = (0,0,0))

                plt.plot(conf[hi,0] + i*sep, conf[hi,1]-j*sep, "o", color = (0 ,.3,.8), markersize = 2)
                plt.plot(conf[pi,0] + i*sep, conf[pi,1]-j*sep, "o", color = (.8,.3,0), markersize = 2)
                
                
                
                c +=1

    plt.xlim(-sep, sep*n+1)
    plt.ylim(-sep*n, sep)
    plt.axis("off")
    plt.show()
    
    
    