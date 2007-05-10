from MMTK import *
from MMTK.Geometry import Box
from Scientific.Functions.Interpolation import InterpolatingFunction
from Scientific.indexing import index_expression
import Numeric as N
from LinearAlgebra import eigenvectors
import os, string, struct, tempfile

class DensityMap:

    def __init__(self, filename):
        if filename is None:
            return
        filetype = os.path.splitext(filename)[1]
        if filetype.lower() == '.ezd':
            self.readEZD(filename)
        elif filetype.lower() == '.ccp4':
            self.readCCP4(filename)
        else:
            raise ValueError("Unknown file type %s" % filetype)
        self.map = None
        self.box = Box(Vector(self.x_axis[0],self.y_axis[0],self.z_axis[0]),
                       Vector(self.x_axis[-1],self.y_axis[-1],self.z_axis[-1]))
        self.normalize()

    def readEZD(filename):
        file = open(filename)
        while 1:
            line = file.readline()
            if not line:
                raise IOError, "unexpected end of file"
            words = string.split(line)
            if words[0] == 'MAP':
                break
            if words[0] == 'CELL':
                cell_x, cell_y, cell_z, alpha, beta, gamma = \
                        map(float, words[1:])
                if alpha != 90. or beta != 90. or gamma != 90.:
                    raise ValueError, "cell must be rectangular"
            if words[0] == 'EXTENT':
                self.nx, self.ny, self.nz = map(float, words[1:])
            if words[0] == 'GRID':
                gnx, gny, gnz = map(float, words[1:])

        data = []
        while 1:
            line = file.readline()
            if not line:
                raise IOError, "unexpected end of file"
            line = string.join(string.split(line, '.-'), '. -')
            words = string.split(line)
            if words[0] == 'END':
                break
            for value in map(float, words):
                data.append(value)

        data = N.array(data)
        data.shape = (self.nz, self.ny, self.nx)
        self.data = N.transpose(data)
        self.x_axis = N.arange(self.nx)*(cell_x/gnx)*Units.Ang
        self.y_axis = N.arange(self.ny)*(cell_y/gny)*Units.Ang
        self.z_axis = N.arange(self.nz)*(cell_z/gnz)*Units.Ang

    def readCCP4(self, filename):
        mapfile = file(filename)
        header_data = mapfile.read(1024)
        NC, NR, NS, MODE, NCSTART, NRSTART, NSSTART, NX, NY, NZ, X, Y, Z, \
            ALPHA, BETA, GAMMA, MAPC, MAPR, MAPS, AMIN, AMAX, AMEAN, \
            ISPG, NSYMBT, LSKFLG = struct.unpack('=10l6f3l3f3l',
                                                 header_data[:4*25])
        if MODE == 2:
            byte_order = '='
        elif MODE == 33554432:
            NC, NR, NS, MODE, NCSTART, NRSTART, NSSTART, NX, NY, NZ, X, Y, Z, \
                ALPHA, BETA, GAMMA, MAPC, MAPR, MAPS, AMIN, AMAX, AMEAN, \
                ISPG, NSYMBT, LSKFLG = struct.unpack('>10l6f3l3f3l',
                                                     header_data[:4*25])
            byte_order = '>'
            if MODE == 33554432:
                NC, NR, NS, MODE, NCSTART, NRSTART, NSSTART, NX, NY, NZ, \
                    X, Y, Z, ALPHA, BETA, GAMMA, MAPC, MAPR, MAPS, \
                    AMIN, AMAX, AMEAN, ISPG, NSYMBT, LSKFLG \
                    = struct.unpack('<10l6f3l3f3l', header_data[:4*25])
                byte_order = '<'

        else:
            raise IOError("Not a mode 2 CCP4 map file")

        symmetry_data = mapfile.read(NSYMBT)
        map_data = mapfile.read(4*NS*NR*NC)

        if byte_order == '=':
            array = N.fromstring(map_data, N.Float32, NC*NR*NS)
        else:
            array = N.zeros((NS*NR*NC,), N.Float32)
            index = 0
            while len(map_data) >= 4*10000:
                values = struct.unpack(byte_order + '10000f',
                                       map_data[:4*10000])
                array[index:index+10000] = N.array(values, N.Float32)
                index += 10000
                map_data = map_data[4*10000:]
            values = struct.unpack(byte_order + '%df' % (len(map_data)/4),
                                   map_data)
            array[index:] = N.array(values, N.Float32)

        del map_data

        array.shape = (NS, NR, NC)
        self.data = N.transpose(array)

        resolution_x = X*Units.Ang/NX
        resolution_y = Y*Units.Ang/NY
        resolution_z = Z*Units.Ang/NZ

        self.x_axis = (NCSTART+N.arange(NC))*resolution_x
        self.y_axis = (NRSTART+N.arange(NR))*resolution_y
        self.z_axis = (NSSTART+N.arange(NS))*resolution_z

    def __getitem__(self, item):
        if not isinstance(item, tuple) or len(item) != 3:
            raise ValueError("indexation requires three slices")
        sx, sy, sz = item
        if not (isinstance(sx, slice) and isinstance(sy, slice) \
                and isinstance(sz, slice)):
            raise ValueError("indexation requires three slices")
        new_map = DensityMap(None)
        new_map.data = self.data[sx, sy, sz]
        new_map.x_axis = self.x_axis[sx]
        new_map.y_axis = self.y_axis[sy]
        new_map.z_axis = self.z_axis[sz]
        new_map.map = None
        new_map.box = Box(Vector(new_map.x_axis[0], new_map.y_axis[0],
                                 new_map.z_axis[0]),
                          Vector(new_map.x_axis[-1], new_map.y_axis[-1],
                                 new_map.z_axis[-1]))
        return new_map

    def normalize(self):
        self.data /= N.sum(N.ravel(self.data))
        
    def makePositive(self):
        min = N.minimum.reduce(N.ravel(self.data))
        if min < 0:
            nonzero_mask = self.data != 0
            self.data = (self.data - min)*nonzero_mask

    def _makeMapObjects(self):
        if self.map is None:
            self.map = InterpolatingFunction((self.x_axis, self.y_axis,
                                              self.z_axis), self.data, 0.)
            self.map_gx = self.map.derivative(0)
            self.map_gy = self.map.derivative(1)
            self.map_gz = self.map.derivative(2)

    def center(self):
        self.map = None
        x_center = N.sum(N.ravel(self.x_axis[:, N.NewAxis, N.NewAxis]
                                 * self.data))
        y_center = N.sum(N.ravel(self.y_axis[N.NewAxis, :, N.NewAxis]
                                 * self.data))
        z_center = N.sum(N.ravel(self.z_axis[N.NewAxis, N.NewAxis, :]
                                 * self.data))
        self.x_axis = self.x_axis - x_center
        self.y_axis = self.y_axis - y_center
        self.z_axis = self.z_axis - z_center
        self.box = Box(Vector(self.x_axis[0],self.y_axis[0],self.z_axis[0]),
                       Vector(self.x_axis[-1],self.y_axis[-1],self.z_axis[-1]))

    def principalAxes(self):
        r = self._rGrid()
        cm = N.sum(N.sum(N.sum(self.data[..., N.NewAxis]*r)))
        r = r - cm[N.NewAxis, N.NewAxis, N.NewAxis, :]
        nx, ny, nz = self.data.shape
        t = 0.
        for i in range(nx): # make loops explicit to conserve memory
            for j in range(ny):
                for k in range(nz):
                    t = t + self.data[i, j, k] * r[i, j, k, N.NewAxis, :] * \
                        r[i, j, k, :, N.NewAxis]
        ev, axes = eigenvectors(t)
        return map(lambda a, b: (Vector(a), b), axes, ev)

    def principalPoints(self):
        (ex, vx), (ey, vy), (ez, vz) = self.principalAxes()
        axes = N.array([ex.array, ey.array, ez.array])
        r = self._rGrid()
        cm = N.sum(N.sum(N.sum(self.data[..., N.NewAxis]*r)))
        r = r - cm[N.NewAxis, N.NewAxis, N.NewAxis, :]
        regions = N.greater(N.dot(r, N.transpose(axes)), 0)
        regions = N.sum(N.array([[[[1, 2, 4]]]])*regions, -1)
        points = []
        for i in range(8):
            mask = N.equal(regions, i)
            weight = N.sum(N.sum(N.sum(mask*self.data)))
            cmr = N.sum(N.sum(N.sum((mask*self.data)[..., N.NewAxis]*r))) \
                  / weight + cm
            points.append(Vector(cmr))
        return points

    def overlap(self, object):
        self._makeMapObjects()
        sum = 0.
        for a in object.atomList():
            x, y, z = a.position()
            sum = sum + self.map(x, y, z)
        return sum

    def gradient(self, object):
        self._makeMapObjects()
        g = ParticleVector(object.universe())
        for a in object.atomList():
            x, y, z = a.position()
            g[a] = Vector(self.map_gx(x, y, z),
                          self.map_gy(x, y, z),
                          self.map_gz(x, y, z))
        return g

    def atomMap(self, object, r0=0.3):
        r = self._rGrid()
        atom_map = N.zeros(self.data.shape, N.Float)
        cutoff = 4.*r0
        for a in object.atomList():
            # An over-eager optimization: it should use
            # an enlarged box
            #if not self.box.enclosesPoint(a.position()):
            #    continue
            ra = a.position().array
            xi1 = N.sum(self.x_axis < ra[0]-cutoff)
            xi2 = N.sum(self.x_axis < ra[0]+cutoff)
            yi1 = N.sum(self.y_axis < ra[1]-cutoff)
            yi2 = N.sum(self.y_axis < ra[1]+cutoff)
            zi1 = N.sum(self.z_axis < ra[2]-cutoff)
            zi2 = N.sum(self.z_axis < ra[2]+cutoff)
            if xi2 > xi1 and yi2 > yi1 and zi2 > zi1:
                dr = r[xi1:xi2, yi1:yi2, zi1:zi2] - \
                     ra[N.NewAxis, N.NewAxis, N.NewAxis, :]
                w = N.exp(-0.5*N.sum(dr**2, axis=-1)/r0**2)
                lmap = atom_map[xi1:xi2, yi1:yi2, zi1:zi2]
                N.add(lmap, w, lmap)
        N.divide(atom_map, N.sum(N.sum(N.sum(atom_map))), atom_map)
        return atom_map

    def fit(self, object, r0=0.4, w_neg = 1.):
        diff = self.atomMap(object, r0)-self.data
        N.multiply(diff, 1.+(w_neg-1.)*N.less(diff, 0.), diff)
        return N.sum(N.sum(N.sum(diff**2)))

    def fitWithGradient(self, object, r0=0.4, w_neg = 1.):
        r = self._rGrid()
        atom_map = N.zeros(self.data.shape, N.Float)
        cutoff = 4.*r0
        for a in object.atomList():
            ra = a.position().array
            xi1 = N.sum(self.x_axis < ra[0]-cutoff)
            xi2 = N.sum(self.x_axis < ra[0]+cutoff)
            yi1 = N.sum(self.y_axis < ra[1]-cutoff)
            yi2 = N.sum(self.y_axis < ra[1]+cutoff)
            zi1 = N.sum(self.z_axis < ra[2]-cutoff)
            zi2 = N.sum(self.z_axis < ra[2]+cutoff)
            if xi2 > xi1 and yi2 > yi1 and zi2 > zi1:
                dr = r[xi1:xi2, yi1:yi2, zi1:zi2] - \
                     ra[N.NewAxis, N.NewAxis, N.NewAxis, :]
                w = N.exp(-0.5*N.sum(dr**2, axis=-1)/r0**2)
                lmap = atom_map[xi1:xi2, yi1:yi2, zi1:zi2]
                N.add(lmap, w, lmap)
        norm_factor = 1./N.sum(N.sum(N.sum(atom_map)))
        N.multiply(atom_map, norm_factor, atom_map)
        N.subtract(atom_map, self.data, atom_map)
        weight = 1.+(w_neg-1.)*N.less(atom_map, 0.)
        N.multiply(atom_map, weight, atom_map)
        g = ParticleVector(object.universe())
        for a in object.atomList():
            ra = a.position().array
            xi1 = N.sum(self.x_axis < ra[0]-cutoff)
            xi2 = N.sum(self.x_axis < ra[0]+cutoff)
            yi1 = N.sum(self.y_axis < ra[1]-cutoff)
            yi2 = N.sum(self.y_axis < ra[1]+cutoff)
            zi1 = N.sum(self.z_axis < ra[2]-cutoff)
            zi2 = N.sum(self.z_axis < ra[2]+cutoff)
            if xi2 > xi1 and yi2 > yi1 and zi2 > zi1:
                dr = r[xi1:xi2, yi1:yi2, zi1:zi2] - \
                     ra[N.NewAxis, N.NewAxis, N.NewAxis, :]
                lmap = atom_map[xi1:xi2, yi1:yi2, zi1:zi2]
                lw = weight[xi1:xi2, yi1:yi2, zi1:zi2]
                w = N.exp(-0.5*N.sum(dr**2, axis=-1)/r0**2)
                v = N.sum(N.sum(N.sum(lmap[..., N.NewAxis]
                                      * lw[..., N.NewAxis] * dr
                                      * w[..., N.NewAxis])))
                g[a] = Vector(v)
        return N.sum(N.sum(N.sum(atom_map**2))), -2*norm_factor*g/r0**2

    def _rGrid(self):
        x = N.add.outer(N.add.outer(self.x_axis, 0*self.y_axis),
                        0*self.z_axis)[..., N.NewAxis]
        y = N.add.outer(N.add.outer(0*self.x_axis, self.y_axis),
                        0*self.z_axis)[..., N.NewAxis]
        z = N.add.outer(N.add.outer(0*self.x_axis, 0*self.y_axis),
                        self.z_axis)[..., N.NewAxis]
        return N.concatenate((x, y, z), axis=-1)
        

    def clipBelow(self, cutoff):
        mask = self.clipMask(cutoff)
        self.data *= mask
        self.data /= N.add.reduce(N.ravel(self.data))
        self.map = None

    def clipAbove(self, cutoff):
        mask = self.clipMask(cutoff)
        self.data *= (1-mask)
        self.data /= N.add.reduce(N.ravel(self.data))
        self.map = None

    def clipMask(self, cutoff):
        max_value = N.maximum.reduce(N.ravel(self.data))
        mask = N.greater_equal(self.data, cutoff*max_value)
        gradient_masks = []
        for i in range(len(mask.shape)):
            upper_index = i*index_expression[::] + index_expression[1::]
            lower_index = i*index_expression[::] + index_expression[:-1:]
            gmask1 = N.greater(self.data[upper_index], self.data[lower_index])
            gmask2 = N.greater(self.data[lower_index], self.data[upper_index])
            gradient_masks.append((gmask1, gmask2))
        while 1:
            new_mask = 0*mask
            for i in range(len(mask.shape)):
                upper_index = i*index_expression[::] + index_expression[1::]
                lower_index = i*index_expression[::] + index_expression[:-1:]
                upper_mask, lower_mask = gradient_masks[i]
                N.logical_or(new_mask[upper_index],
                             N.logical_and(mask[lower_index], lower_mask),
                             new_mask[upper_index])
                N.logical_or(new_mask[lower_index],
                             N.logical_and(mask[upper_index], upper_mask),
                             new_mask[lower_index])
            N.logical_and(new_mask, N.logical_not(mask), new_mask)
            N.logical_or(mask, new_mask, mask)
            if N.sum(N.ravel(new_mask)) == 0:
                break
        return mask

    def dataDump(self, filename):
        outfile = file(filename, 'w')
        data = N.ravel(N.transpose(self.data))
        data = 10000.*data/N.maximum.reduce(data)
        data.shape = (len(data)/7, 7)
        for line in data:
            for number in line:
                outfile.write('%.1f ' % number)
            outfile.write('\n')
        outfile.close()

    def writeVTK(self, filename):
        import pyvtk
        origin = N.array([self.x_axis[0], self.y_axis[0], self.z_axis[0]])
        spacing = N.array([self.x_axis[1], self.y_axis[1], self.z_axis[1]]) \
                  - origin
        values = pyvtk.Scalars(N.ravel(N.transpose(self.data)),
                               'electron density')
        data = pyvtk.VtkData(pyvtk.StructuredPoints(self.data.shape,
                                                    origin, spacing),
                             'Density map',
                             pyvtk.PointData(values))
        data.tofile(filename, format='binary')

    def view(self, gmodule, lower, upper, file=None):
        scene = gmodule.Scene()
        scale = gmodule.ColorScale(1.)
        size = 0.5*min(self.x_axis[1], self.y_axis[1], self.z_axis[1])
        for i in range(self.nx):
            for j in range(self.ny):
                for k in range(self.nz):
                    value = self.data[i, j, k]
                    if value > lower and value < upper:
                        center = Vector(self.x_axis[i], self.y_axis[j],
                                        self.z_axis[k])
                        m = gmodule.Material(diffuse_color=scale(value))
                        object = gmodule.Sphere(center, size,
                                                material=m)
                        scene.addObject(object)
        if file is None:
            scene.view()
        else:
            scene.writeToFile(file)

    def viewVMD(self, filename=None, pdb_filename = None):
        run_vmd = 0
        if filename is None:
            filename = tempfile.mktemp()
            run_vmd = 1
        file = open(filename, 'w')
        if pdb_filename is not None:
            file.write('mol load pdb %s\n' % pdb_filename)
        file.write('mol volume top "Electron Density" \\\n')
        x_origin = self.x_axis[0]/Units.Ang
        y_origin = self.y_axis[0]/Units.Ang
        z_origin = self.z_axis[0]/Units.Ang
        x_length = self.x_axis[-1]/Units.Ang-x_origin
        y_length = self.y_axis[-1]/Units.Ang-y_origin
        z_length = self.z_axis[-1]/Units.Ang-z_origin
        file.write('  {%f %f %f} \\\n' % (x_origin, y_origin, z_origin))
        file.write('  {%f 0. 0.} \\\n' % x_length)
        file.write('  {0. %f 0.} \\\n' % y_length)
        file.write('  {0. 0. %f} \\\n' % z_length)
        file.write('  %d %d %d \\\n' % self.data.shape)
        file.write('  {')
        factor = 1./N.maximum.reduce(N.ravel(self.data))
        for iz in range(self.data.shape[2]):
            for iy in range(self.data.shape[1]):
                for ix in range(self.data.shape[0]):
                    file.write(str(factor*self.data[ix, iy, iz]) + ' ')
        file.write('}\n')
        file.write('mol addrep top\nmol modstyle 1 top isosurface\n')
        if run_vmd:
            file.write('file delete %s\n' % filename)
        file.close()
        if run_vmd:
            os.system('vmd -e ' + filename + ' 1> /dev/null 2>&1')

    def writeXPlor(self, filename):
        from Scientific.IO.FortranFormat import FortranFormat, FortranLine
        file = open(filename, 'w')
        file.write('\n       1 !NTITLE\n')
        file.write('REMARKS Electronic density map\n')
        data = [self.data.shape[0], 1, self.data.shape[0],
                self.data.shape[1], 1, self.data.shape[1],
                self.data.shape[2], 1, self.data.shape[2]]
        file.write(str(FortranLine(data, '9I8')) + '\n')
        x = (self.x_axis[-1]-self.x_axis[0])/Units.Ang
        y = (self.y_axis[-1]-self.y_axis[0])/Units.Ang
        z = (self.z_axis[-1]-self.z_axis[0])/Units.Ang
        data = [x, y, z] + 3*[90.]
        file.write(str(FortranLine(data, '6E12.5')) + '\n')
        file.write('ZYX\n')
        map_data = N.ravel(self.data)
        map_data = map_data / N.maximum.reduce(map_data)
        average = N.sum(map_data)/len(map_data)
        sd = N.sqrt(N.sum((map_data-average)**2)/len(map_data))
        map_data.shape = self.data.shape
        for i in range(self.data.shape[2]):
            file.write(str(FortranLine([i], 'I8')) + '\n')
            data = list(N.ravel(N.transpose(map_data[:, :, i])))
            while data:
                file.write(str(FortranLine(data[:6],
                                           '%dE12.5' % min(6, len(data)))) \
                           + '\n')
                data = data[6:]
        file.write(str(FortranLine([-9999], 'I8')) + '\n')
        file.write(str(FortranLine([average, sd], '2(E12.4,1X)')) + '\n')
        file.close()


if __name__ == '__main__':

    from Scientific.Visualization import VMD; Graphics = VMD
    #map = DensityMap('d/dtg8.ezd')
    map = DensityMap('d/e1-p.ezd')
    map.center()
    #map.view(Graphics, 0.6, 1.)
    map.data = map.data*10000.
    map.viewVMD('map_e1-p.vmd', 'initial_e1-p_1iwob_1.pdb')
    map.writeVTK('map_e1-p.vtk')

