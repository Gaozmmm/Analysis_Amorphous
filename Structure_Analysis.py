import numpy as np
from ase.io import read
from ase.formula import Formula
from scipy.spatial import distance_matrix, Voronoi
import time
import matplotlib.pyplot as plt
from re import compile, findall


class RadialAnalysis(object):
    def __init__(self, structure: "POSCAR", r_cutoff=9, bin_size=0.1,  # support only for POSCAR
                 step_length=1, overall_density=0, smooth=2, supercell=np.eye(3, 3), gaussian_smearing_sigma=0):
        self.supercell = supercell
        if len(self.supercell) != 3:
            raise ValueError("You have to put in a 3x3 numpy array")
        self._frame = make_supercell(read(structure), self.supercell)
        self.coor_f = self._frame.positions
        self.cell_f = self._frame.cell.array
        self.r_cut = r_cutoff
        self.bin_size = bin_size
        self.step_l = step_length
        self.atom_num = len(self.coor_f)
        self.NDF = {}
        self.PDF = {}
        self.dis_matrix = np.zeros(self.atom_num)
        self.density = overall_density
        self.original_dis_matrix = distance_matrix(self.coor_f, self.coor_f) // self.bin_size * self.bin_size
        print("Distance Matrix Finished!")
        self.n_decimal = (str(self.bin_size))[::-1].find('.')
        self.dis_matrix = np.around(self.original_dis_matrix, self.n_decimal)
        self._bins = self.n_bins()
        self.smooth = smooth
        self.sigma = gaussian_smearing_sigma
        print(self.atom_num)

        if self.density == 0:
            self.density = self.atom_num / self._frame.cell.volume


    def n_bins(self):
        _bins = int(self.r_cut // self.bin_size)
        if _bins < 2:
            raise ValueError("Increase r_cutoff or Decrease bin_size")
        return _bins

    def set_system_NDF(self):  # use this to initialize the distance dictionary for the whole system
        for i in range(1, self._bins):
            dis = round(i * self.bin_size, self.n_decimal)
            self.NDF[dis] = np.count_nonzero(self.dis_matrix == dis)
        print("Counting distance Finished!")
        if self.smooth > 0:
            self.get_smooth()
        print("Smoothing Finished!")

    def get_system_NDF(self):
        self.set_system_NDF()
        return self.NDF

    def set_pair_NDF(self, element_1, element_2):
        print(element_2, element_1)
        # if element_1 or element_2 == 'What?':
        #     print('OK')
        #     raise ValueError('If you DONT tell me which atom, Why you use this method!!!')
        symbol = str(self._frame.get_chemical_formula())
        pattern = compile(r'([A-Za-z]+)([0-9]+)')  # find the specific rows of element
        start, end = 0, 0
        index_frame = []
        for (letters, numbers) in findall(pattern, symbol):
            end += int(numbers)
            index_frame.append([letters, start, end])
            start = end

        s_1, s_2 = 0, 0
        e_1, e_2 = 0, 0
        for i in range(len(index_frame)):
            if index_frame[i][0] == element_1:
                s_1, e_1 = index_frame[i][1], index_frame[i][2]
            if index_frame[i][0] == element_2:
                s_2, e_2 = index_frame[i][1], index_frame[i][2]
        pair_dis_matrix = self.dis_matrix[s_1:e_1, s_2:e_2]

        for i in range(1, self._bins):
            dis = round(i * self.bin_size, self.n_decimal)
            self.NDF[dis] = np.count_nonzero(pair_dis_matrix == dis)
        print("Counting distance Finished!")
        if self.smooth > 0:
            self.get_smooth()
        print("Smoothing Finished!")

        return (e_1 - s_1), (e_2 - s_2)

    def get_pair_NDF(self, element_1='What?', element_2='What?'):
        self.set_pair_NDF(element_1, element_2)
        return self.NDF

    def set_system_PDF(self):
        self.set_system_NDF()
        print("start to convert RDF to PDF")
        print(self.density, self.atom_num, self.bin_size)
        for radii in self.NDF.keys():
            self.PDF[radii] = self.NDF[radii] / (4 * 3.14159 * (radii ** 2) * self.density * self.atom_num
                                                 * self.bin_size)

    def get_system_PDF(self):
        self.set_system_PDF()
        return self.PDF

    def set_pair_PDF(self, element1, element2, element2_density=0):
        n_1, n_2 = self.set_pair_NDF(element1, element2)
        print(n_1, n_2)
        if element2_density == 0:
            element2_density = n_2 / self._frame.cell.volume
        print("start to convert RDF to PDF")
        for radii in self.NDF.keys():
            # print(n_2, self.bin_size, element2_density, radii)
            self.PDF[radii] = self.NDF[radii] / (4 * 3.14159 * (radii ** 2) * element2_density
                                                 * n_2 * self.bin_size)

    def get_pair_PDF(self, element_a, element_b):
        self.set_pair_PDF(element_a, element_b)
        return self.PDF

    def plot_pair_PDF(self, element_1, element_2, title=''):
        self.set_pair_PDF(element_1, element_2)
        x, y = [], []
        for j in self.PDF.keys():
            x.append(j)
            y.append(self.PDF[j])
        plt.plot(x, y)
        plt.xlabel("$r$ / Angstrom)")
        plt.ylabel(element_1+'-'+element_2+" PDF($r$)")
        if title == '':
            title = element_1+'-'+element_2+' PDF'
        plt.title(title)
        return plt


    def plot_system_PDF(self, title='PDF'):
        self.set_system_PDF()
        x, y = [], []
        for j in self.PDF.keys():
            x.append(j)
            y.append(self.PDF[j])

        y = smear(y, self.sigma)
        plt.plot(x, y)
        plt.xlabel(r"$r (\AA)$", fontsize=16)
        plt.ylabel(r"$g(r)$", fontsize=16)
        plt.title(title)
        return plt

    def plot_system_NDF(self, title='NDF'):
        self.set_system_NDF()
        x, y = [], []
        for j in self.NDF.keys():
            x.append(j)
            y.append(self.NDF[j])
        plt.plot(x, y)
        plt.xlabel("$r$ / Angstrom)")
        plt.ylabel("NDF($r$)")
        plt.title(title)
        return plt

    def plot_pair_NDF(self, element1='What', element2='What', title='NDF'):
        self.set_pair_NDF(element1, element2)
        x, y = [], []
        for j in self.NDF.keys():
            x.append(j)
            y.append(self.NDF[j])
        plt.plot(x, y)
        plt.xlabel("$r$ / Angstrom)")
        plt.ylabel(element1+'-'+element2+" NDF($r$)")
        plt.title(title)
        return plt

    def get_smooth(self):
        """
        Helper function to recursively smooth RDFs using a 5-parameter Savitzky-Golay filter.
        Args:
            RDFs: A dictionary of partial radial distribution functions
            with pairs as keys and RDFs as values.
            passes: number of times the filter is applied during smoothing.
        Returns
            RDFs dictionary with with each RDF smoothed.
        inspired from Github --- mpmorph
        """
        while self.smooth != 0:
            print("smooth still has: ", self.smooth, " attempts to go")
            j_i = [0, 0, 0, 0, 0]
            for j in range(3, self._bins - 2):
                for k in -2, -1, 0, 1, 2:
                    j_i[k+2] = round(j * self.bin_size - k * self.bin_size, self.n_decimal)
                self.NDF[j_i[2]] = (-3 * self.NDF[j_i[0]] + 12 * self.NDF[j_i[1]] + 17 * self.NDF[j_i[2]] +
                                    12 * self.NDF[j_i[3]] - 3 * self.NDF[j_i[4]]) / 35.0
            self.smooth -= 1


class VoronoiAnalysis(object):
    def __init__(self, structure: "POSCAR", r_cutoff = 5):
        self._frame = read(structure)
        self.coor_f = self._frame.positions

    def set_voronoi(self):
        points = self.coor_f
        print(points)
        print(type(points))
        vor = Voronoi(points)
        print(vor)
        print(vor.points)
        print(vor.ridge_dict)
        from scipy.spatial import voronoi_plot_2d


class BondAngleAnalysis(object):
    def __init__(self):
        pass

def smear(data, sigma):
    from scipy.ndimage.filters import gaussian_filter1d
    """
    Apply Gaussian smearing to spectrum y value.
    Args:
        sigma: Std dev for Gaussian smear function
    """
    if sigma == 0:
        return data
    else:
        return gaussian_filter1d(data, sigma)



# a = VoronoiAnalysis('C:\Gao\lmp_MC_relax_8.vasp')
# a.set_voronoi()




from ase.build import make_supercell

a0 = read('C:\Gao\lNaCl_mp_22862_primitive.vasp')

s = np.array([[5,0,0],[0,5,0],[0,0,5]])

S1 = RadialAnalysis('C:\Gao\Ar.vasp', bin_size=0.08, supercell=s, smooth=0, r_cutoff=15, gaussian_smearing_sigma=2)
F = S1.plot_system_PDF()
F.show()
# S1 = RadialAnalysis('C:\Gao\lmp_MC_relax_8.vasp', smooth=1)
# # S4 = StructureAnalysis('C:\Gao\lmp_MC_relax_8.vasp', smooth=4)
# F = S1.get_pair_PDF(element_a='Nb', element_b='Nb')
# F = S1.plot_pair_PDF(element_1='Nb', element_2='Nb')
# F.show()
# # F = S4.plot_RDF()
# print(time.process_time())