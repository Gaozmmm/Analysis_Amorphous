import numpy as np
from ase.io import read
from scipy.spatial import distance_matrix, Voronoi, ConvexHull
import time
import matplotlib.pyplot as plt
from re import compile, findall
from ase.build import make_supercell
from math import pi


class RadialAnalysis(object):
    def __init__(self, structure: "POSCAR", r_cutoff=9, bin_size=0.1,  # support only for POSCAR
                 step_length=1, overall_density=0, supercell=np.eye(3, 3), gaussian_smearing_sigma=0):
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
        self.RDF = {}
        self.n_decimal = (str(self.bin_size))[::-1].find('.')
        self.dis_matrix = np.zeros(self.atom_num)
        self.center_coor_f = self._frame.positions
        self.original_dis_matrix = np.zeros(self.atom_num)
        self.density = overall_density
        # self.dis_matrix_cal()
        self.original_dis_matrix = distance_matrix(self.coor_f, self.coor_f) // self.bin_size * self.bin_size
        print("Distance Matrix Finished!, it's shape is:" + str(np.shape(self.dis_matrix)))
        self.dis_matrix = np.around(self.original_dis_matrix, self.n_decimal)
        self._bins = self.n_bins()
        self.sigma = gaussian_smearing_sigma
        print(self.atom_num)
        self.Sq = {}

        if self.density == 0:
            self.density = self.atom_num / self._frame.cell.volume

    def dis_matrix_cal(self):
        count = 0
        for i in range(len(self.coor_f)):
            if self.coor_f[i][0] <= self.r_cut or self.coor_f[i][0] >= self.cell_f[0][0] - self.r_cut\
                or self.coor_f[i][1] <= self.r_cut or self.coor_f[i][1] >= self.cell_f[1][1] - self.r_cut\
                    or self.coor_f[i][1] <= self.r_cut or self.coor_f[i][2] >= self.cell_f[2][2] - self.r_cut:
                i -= count
                self.center_coor_f = np.delete(self.center_coor_f, i, 0)
                count += 1
        self.original_dis_matrix = distance_matrix(self.center_coor_f, self.coor_f)
        self.dis_matrix = np.around(self.original_dis_matrix, self.n_decimal)

    def n_bins(self):
        _bins = int(self.r_cut // self.bin_size)
        if _bins < 2:
            raise ValueError("Increase r_cutoff or Decrease bin_size")
        return _bins

    def set_system_NDF(self):  # use this to initialize the distance dictionary for the whole system
        for i in range(1, self._bins):
            dis = round(i * self.bin_size, self.n_decimal)
            self.NDF[dis] = np.count_nonzero(self.dis_matrix == dis)
        print(self.NDF)
        print("Counting distance Finished!")

    def get_system_NDF(self):
        self.set_system_NDF()
        return self.NDF

    def set_pair_NDF(self, element_1, element_2):
        print(element_1, element_2)
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
        print(index_frame)
        for i in range(len(index_frame)):
            if index_frame[i][0] == element_1:
                s_1, e_1 = index_frame[i][1], index_frame[i][2]
            if index_frame[i][0] == element_2:
                s_2, e_2 = index_frame[i][1], index_frame[i][2]
        print(s_1, s_2, e_1, e_2)
        pair_dis_matrix = self.dis_matrix[s_1:e_1, s_2:e_2]
        print(np.shape(pair_dis_matrix))

        for i in range(1, self._bins):
            dis = round(i * self.bin_size, self.n_decimal)
            self.NDF[dis] = np.count_nonzero(pair_dis_matrix == dis)
        print("Counting distance Finished!")

        return (e_1 - s_1), (e_2 - s_2)

    def get_pair_NDF(self, element_1='What?', element_2='What?'):
        self.set_pair_NDF(element_1, element_2)
        return self.NDF

    def set_system_PDF(self):
        self.set_system_NDF()
        print("start to convert NDF to PDF")
        print("System density is : " + str(self.density))
        print("System number of atoms is : " + str(self.atom_num))
        print("System bin size is : " + str(self.bin_size))
        print(np.shape(self.dis_matrix))
        for radii in self.NDF.keys():
            self.PDF[radii] = self.NDF[radii] / (4 * 3.14159 * (radii ** 2) * self.density * len(self.dis_matrix)
                                                 * self.bin_size)

    def get_system_PDF(self):
        self.set_system_PDF()
        return self.PDF

    # def set_system_RDF(self):
    #     self.set_system_NDF()
    #     print("start to convert NDF to RDF")
    #     sum_of_atom = 0
    #     for radii in self.NDF.keys():
    #         sum_of_atom += self.NDF[radii]
    #         self.RDF[radii] = sum_of_atom / (4 * 3.14159 * (radii ** 2) * self.density)
    #     print(self.RDF)

    # def get_system_RDF(self):
    #     self.set_system_RDF()
    #     return self.RDF

    def set_pair_PDF(self, element1, element2, element2_density=0):
        n_1, n_2 = self.set_pair_NDF(element1, element2)
        print(n_1, n_2)
        if element2_density == 0:
            element2_density = n_2 / self._frame.cell.volume
        print("start to convert NDF to PDF")
        for radii in self.NDF.keys():
            # print(n_2, self.bin_size, element2_density, radii)
            print(np.shape(self.dis_matrix))
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
        y = smear(y, self.sigma)
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

        y = smear(y, self.sigma)
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
        y = smear(y, self.sigma)
        plt.plot(x, y)
        plt.xlabel("$r$ / Angstrom)")
        plt.ylabel(element1+'-'+element2+" NDF($r$)")
        plt.title(title)
        return plt

    # def plot_system_RDF(self, title='RDF'):
    #     self.set_system_RDF()
    #     x, y = [], []
    #     for j in self.RDF.keys():
    #         x.append(j)
    #         y.append(self.RDF[j])
    #     y = smear(y, self.sigma)
    #     plt.plot(x, y)
    #     plt.xlabel("$r$ / Angstrom)")
    #     plt.ylabel("RDF($r$)")
    #     plt.title(title)
    #     return plt

    def set_structure_factor(self, k_limit=1.4, is_2d=False):

        kxint = 2 * 3.1415926535 / self.cell_f[0][0]
        kyint = 2 * 3.1415926535 / self.cell_f[1][1]
        if self.cell_f[2][2] != 0 or is_2d is False:
            kzint = 2 * 3.1415926535 / self.cell_f[2][2]
        else:
            kzint = 19980807
        print(kxint, kyint, kzint)
        k_num = int(k_limit // (2 * 3.1415926535 / max(self.cell_f[0][0], self.cell_f[1][1], self.cell_f[2][2])))
        if k_num <= 1:
            print("Please increase k_limit")
        k_matrix = np.zeros(shape=(k_num**3, 3), dtype=float)
        count = 0
        for kx_num in range(0, k_num):
            for ky_num in range(0, k_num):
                for kz_num in range(0, k_num):
                    k_matrix[count][0] = kxint * kx_num
                    k_matrix[count][1] = kyint * ky_num
                    k_matrix[count][2] = kzint * kz_num
                    count += 1
        k_matrix = np.delete(k_matrix, 0, axis=0)
        print("start plot1")
        Sq_count = {}
        coor_f_t = np.transpose(self.coor_f)
        for k in k_matrix:
            q = np.linalg.norm(k)
            q = np.round(q, 2)
            if q in self.Sq.keys():
                # b = (np.abs(np.sum(np.cos((np.dot(k, coor_f_t)))))) ** 2
                # c = (np.abs(np.sum(np.sin((np.dot(k, coor_f_t)))))) ** 2
                b = (np.abs(np.sum(np.exp(1j * (np.dot(k, coor_f_t)))))) ** 2
                # self.Sq[q] += (b + c)
                self.Sq[q] += b
                Sq_count[q] += 1
            else:
                # b = (np.abs(np.sum(np.cos((np.dot(k, coor_f_t)))))) ** 2
                # c = (np.abs(np.sum(np.sin((np.dot(k, coor_f_t)))))) ** 2
                b = (np.abs(np.sum(np.exp(1j * (np.dot(k, coor_f_t)))))) ** 2
                # self.Sq[q] = (b + c)
                self.Sq[q] = b
                Sq_count[q] = 1
        print("start plot")
        for key in self.Sq.keys():
            self.Sq[key] = self.Sq[key] / (self.atom_num * Sq_count[key])
        self.Sq = {k: self.Sq[k] for k in sorted(self.Sq.keys())}  # sort the dict

        # self.set_system_PDF()
        # for i in range(1, 105):
        #     k = i * 2 * pi / self.cell_f[0][0]
        #     self.Sq[k] = 1
        #     for r in self.PDF:
        #         inte = (4 * 3.1415926 * self.density) * (r * (self.PDF[r]-1) * np.sin(k * r)) \
        #                * self.bin_size / k
        #         if r == 0 or r == self.r_cut:
        #             self.Sq[k] += 0.5 * inte
        #         else:
        #             self.Sq[k] += inte


    def get_structure_factor_1d(self):
        self.set_structure_factor()
        return self.Sq

    def plot_system_Sq_1d(self, title='Structure Factor'):
        print("start plot")
        self.set_structure_factor()
        x, y = [], []
        for k in self.Sq.keys():
            # kxint = 2 * 3.1415926535 / self.cell_f[0][0]
            x.append(k)
            y.append(self.Sq[k])

        y = smear(y, self.sigma)
        plt.plot(x, y)
        plt.xlabel("$k$ / $2{\pi}/a^{-1}$")
        plt.ylabel("S($k$)")
        plt.title(title)
        return plt


    def set_3d_structure_factor(self):
        from math import cos, sin
        k = [1, 0, 0]
        Sq_3d = {}
        temp = 0
        co, si = 0, 0
        dis_m_f = distance_matrix(self.coor_f, self.coor_f)
        countk = 0
        for k in range():
            for i in range(self.atom_num):
                co += cos(np.dot(k, self.coor_f[i]))
                si += sin(np.dot(k, self.coor_f[i]))
            temp += ((abs(co)) ** 2 + (abs(si)) ** 2)
            countk += 1
            Sq_3d[k] = temp
            pass



    def num_var(self):
        va = {}
        f1 = 0
        for i in range(30, 200):
            i *= 0.1
            count1 = 0
            for j in range(len(self.dis_matrix)):
                num1 = np.count_nonzero(self.dis_matrix[j] == i)
                count1 += num1
                f1 += (num1 ** 2)
            if count1 == 0:
                continue
            f1 = f1 / count1
            f2 = (np.count_nonzero(self.dis_matrix == i) / count1) ** 2
            va[i] = (f1 - f2) / (i ** 2)
        xx, yy = [], []
        for key in va.keys():
            xx.append(key)
            yy.append(va[key])
        yy = smear(yy, self.sigma)
        plt.plot(xx, yy)
        plt.xlabel('R / $\AA$')
        plt.ylabel('$\sigma^2_N(R)/R^2$')
        plt.show()



class VoronoiAnalysis(object):
        def __init__(self):
        self.vor_ensemble = None

    @staticmethod
    def voronoi_analysis(filename: "POSCAR,cif,vasp", n=0, cutoff=5.0, qhull_options="Qbb Qc Qz"):
        structure = mg.Structure.from_file(filename)
        center = structure[n]
        neighbors = structure.get_sites_in_sphere(center.coords, cutoff)
        neighbors = [i[0] for i in sorted(neighbors, key=lambda s: s[1])]
        qvoronoi_input = np.array([s.coords for s in neighbors])
        voro = Voronoi(qvoronoi_input, qhull_options=qhull_options)
        vor_index = np.array([0, 0, 0, 0])

        for key in voro.ridge_dict:
            if 0 in key:
                "This means if the center atom is in key"
                if -1 in key:
                    "This means if an infinity point is in key"
                    print("Cutoff too short. Exiting.")
                    return None
                else:
                    try:
                        vor_index[len(voro.ridge_dict[key]) - 3] += 1
                    except IndexError:
                        # If a facet has more than 7 edges, it's skipped here
                        pass
        return vor_index

    def from_structures(self, structures, cutoff=4.0, step_freq=10, qhull_options="Qbb Qc Qz"):
        print("This might take a while...")
        voro_dict = {}
        step = 0
        # for structure in structures:
        #     step += 1
        #     if step % step_freq != 0:
        #         continue

        v = []
        for n in range(len(structures)):
            v.append(str(self.voronoi_analysis(structures, n=n, cutoff=cutoff,
                                               qhull_options=qhull_options).view()))
        for voro in v:
            if voro in voro_dict:
                voro_dict[voro] += 1
            else:
                voro_dict[voro] = 1
        self.vor_ensemble = sorted(voro_dict.items(), key=lambda x: (x[1], x[0]), reverse=True)[:15]
        return self.vor_ensemble

    @property
    def plot_vor_analysis(self):
        t = list(zip(*self.vor_ensemble))
        labels = t[0]
        val = list(t[1])
        tot = np.sum(val)
        val = [float(j) / tot for j in val]
        pos = np.arange(len(val)) + .5  # the bar centers on the y axis
        plt.figure(figsize=(4, 4))
        plt.barh(pos, val, align='center', color='navy', alpha=0.5)
        plt.yticks(pos, labels)
        plt.xlabel('Frequency')
        plt.title('Voronoi Spectra')
        plt.grid(axis='x')
        return plt

    def voronoi_volumes(points):
        v = Voronoi(points)
        vol = np.zeros(v.npoints)
        for i, reg_num in enumerate(v.point_region):
            indices = v.regions[reg_num]
            if -1 in indices:  # some regions can be opened
                vol[i] = np.inf
            else:
                vol[i] = ConvexHull(v.vertices[indices]).volume
        return vol


class BondAngleAnalysis(object):
    def __init__(self):
        pass

def smear(data, sigma):    # In order to smooth the curve
    from scipy.ndimage.filters import gaussian_filter1d
    if sigma == 0:
        return data
    else:
        return gaussian_filter1d(data, sigma)



# a = VoronoiAnalysis('C:\Gao\lmp_MC_relax_8.vasp')
# a.set_voronoi()





# a0 = read('C:\Gao\lNaCl_mp_22862_primitive.vasp')


S1 = RadialAnalysis('C:\Gao\lmp_MC_relax_8.vasp', bin_size=0.1, r_cutoff=8)
a = S1.plot_system_Sq_1d()
a.show()
# print(a)
# print(time.process_time())
# f = VoronoiAnalysis('C:\Gao\lNaCl_mp_22862_primitive.vasp', supercell=s)
# f.set_voronoi()
