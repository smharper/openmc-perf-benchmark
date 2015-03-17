from __future__ import print_function

import xml.etree.ElementTree as ET

import numpy as np


class Particle():
    def __init__(self, xyz, radius):
        self.xyz = np.array(xyz, dtype=np.float64)
        assert self.xyz.shape == (3,)
        self.radius = radius
        assert self.radius >= 0.0
        self.lat_inds = None


class Lattice():
    def __init__(self, lower_left, pitch, n_cells):
        self.lower_left = np.array(lower_left, dtype=np.float64)
        assert self.lower_left.shape == (3,)
        self.pitch = np.array(pitch, dtype=np.float64)
        assert self.pitch.shape == (3,)
        self.n_cells = np.array(n_cells, dtype=int)
        assert self.n_cells.shape == (3,)
        self.dic = dict()

        for i_z in range(1, self.n_cells[2]+1):
            for i_y in range(1, self.n_cells[1]+1):
                for i_x in range(1, self.n_cells[0]+1):
                    self.dic[(i_x, i_y, i_z)] = []

    def add_particle(self, particle):
        xyz = particle.xyz - self.lower_left
        xyz_low = xyz - particle.radius
        xyz_hi = xyz + particle.radius

        i_xyz_low = np.floor(xyz_low / self.pitch).astype(int)
        i_xyz_hi = np.floor(xyz_hi / self.pitch).astype(int)

        i_xyz_low = np.maximum(i_xyz_low, (0, 0, 0))
        i_xyz_hi = np.minimum(i_xyz_hi, self.n_cells-1)

        for i_z in range(i_xyz_low[2], i_xyz_hi[2]+1):
            for i_y in range(i_xyz_low[1], i_xyz_hi[1]+1):
                for i_x in range(i_xyz_low[0], i_xyz_hi[0]+1):
                    origin = (np.array((i_x, i_y, i_z))
                         + np.array((0.5, 0.5, 0.5))) * self.pitch
                    translated = Particle(xyz - origin, particle.radius)
                    self.dic[(i_x+1, i_y+1, i_z+1)].append(translated)


def read_locs(fname):
    particles = []

    with open(fname) as fin:
        for line in fin.readlines():
            line = [x for x in line.split(' ') if x != '']
            xyz = (line[0], line[1], line[2])
            radius = 0.0565
            particles.append(Particle(xyz, radius))

    return particles


if __name__ == '__main__':
    #### Set-up environment. ####
    # Declare lattice parameters.
    lower_left = (-0.9, -0.9, 0.0)
    pitch = (0.12, 0.12, 0.135)
    n_unique_cells = (15, 15, 50)
    n_axial_copies = 20
    lower_left_translation = (0.0, 0.0, -67.5)

    #pitch = (0.06, 0.06, 0.0675)
    #n_unique_cells = (30, 30, 100)
    #n_axial_copies = 20

    #pitch = (0.3, 0.3, 0.3375)
    #n_unique_cells = (6, 6, 20)
    #n_axial_copies = 20

    #pitch = (0.3, 0.3, 0.3)
    #n_unique_cells = (6, 6, 150)
    #n_axial_copies = 3

    # Declare XML element id parameters.
    lat_id = 51
    particle_universe = 11
    surrounding_universe = 3
    surrounding_material = 604
    starting_particle_id = 5001

    # Read the base geometry file.
    tree = ET.parse('base_geo.xml')
    root = tree.getroot()

    # Read the particle locations.
    particles = read_locs('locs_R0p8_b0p0_t6p75_pf0p35_r0p0565.txt')
    #particles = read_locs('locs_R0p8_b0_t45_pf0p35_r0p0565.txt')

    # Generate a lattice mapping of the particles.
    lat = Lattice(lower_left, pitch, n_unique_cells)
    #for p in particles: lat.add_particle(p)
    i = 1
    for p in particles:
        lat.add_particle(p)

    # Adjust the lattice position
    lat.lower_left = lat.lower_left + np.array(lower_left_translation)

    #### Add the particle definitions to the XML tree. ####
    universe_id = starting_particle_id
    cell_id = starting_particle_id
    for i_z in range(1, lat.n_cells[2]+1):  # OpenMC is 1-indexed
        for i_y in range(1, lat.n_cells[1]+1):
            for i_x in range(1, lat.n_cells[0]+1):
                surfs = []  # Keep a list of surface id's in this universe
                for p in lat.dic[(i_x, i_y, i_z)]:
                    # Add surface definition (outer edge of particles).
                    coeffs = [str(x) for x in p.xyz]
                    coeffs.append(str(p.radius))
                    surf = ET.Element('surface',
                         {'id':str(cell_id),
                          'type':'sphere',
                          'coeffs':' '.join(coeffs)})
                    root.append(surf)
                    surfs.append(cell_id)

                    # Add cell definition (inside of each particle).
                    cell = ET.Element('cell',
                         {'id':str(cell_id),
                          'universe':str(universe_id),
                          'fill':str(particle_universe),
                          'surfaces':str(-cell_id),
                          'translation':' '.join(str(s) for s in p.xyz)})
                    root.append(cell)
                    cell_id += 1

                # Add cell definition (outside of each particle).
                cell = ET.Element('cell',
                     {'id':str(cell_id),
                      'universe':str(universe_id),
                      'material':str(surrounding_material),
                      'surfaces':' '.join(str(x) for x in surfs)})
                root.append(cell)
                cell_id += 1
                universe_id += 1

    #### Add the lattice definition to the tree. ####
    # Make list of universe id's with y-axis inverted.
    universes = []
    universe_id = starting_particle_id
    for i_z in range(0, lat.n_cells[2]):
        plane = []
        for i_y in range(0, lat.n_cells[1])[::-1]:
            offset = (i_z*lat.n_cells[0]*lat.n_cells[1] + i_y*lat.n_cells[0]
                 +starting_particle_id)
            row = [str(i_x + offset) for i_x in range(0, lat.n_cells[0])]
            row = ' '.join(row)
            plane.append(row)
        plane = '  '.join(plane)
        universes.append(plane)
    universes = '\n'.join(universes)
    universes = '\n\n'.join([universes]*n_axial_copies)

    # Make an XML definition of the universes.
    lat_unis = ET.Element('universes')
    lat_unis.text = universes

    # Make an XML definition of the lattice.
    real_n_cells = lat.n_cells * np.array([1, 1, n_axial_copies])
    lat = ET.Element('lattice',
         {'id':str(lat_id),
          'type':'rect',
          'outer':str(surrounding_universe),
          'dimension':' '.join([str(s) for s in real_n_cells]),
          'width':' '.join([str(s) for s in lat.pitch]),
          'lower_left':' '.join([str(s) for s in lat.lower_left])})
    lat.append(lat_unis)
    root.append(lat)

    # Write the XML file.
    tree.write('geometry.xml')
