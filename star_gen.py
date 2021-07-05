import time as t
import math
import random as r
import os

global sss
sss = {'d':0, 'f':0, 'c':0}


# Class STAR which get point in space and mass to generate object STAR
# Coordinates are given in pc, but inside the class they transformed to METERS, mass are given im M0(mass of sun), and
# also transformed into KG inside the class
class Star:
    def __init__(self, crd, m):  #
        self.G = 6.674 * 10 ** (-11)
        self.pc = 3.085677581 * 10 ** 16  # meters
        self.M0 = 1.989 * (10 ** 30)  # kg
        self.c = 299792458  # m\s
        self.dm = 0

        self.mass = self.m = m  # kg
        self.x = crd[0]  # meters
        self.y = crd[1]  # meters
        self.h = crd[2]  # meters


class Generator:
    def __init__(self):
        self.path = os.path.dirname(os.path.abspath(__file__)) + "\\"
        self.stars = []
        self.pc = 3.085677581 * 10 ** 16
        self.G = 6.674 * 10 ** (-11)
        self.M0 = 1.989 * (10 ** 30)
        self.c = 299792458
        self.r = 100
        self.random_list = []
        self.get_good_random_mass_prep()

    def einasto_distibution(self, height, population, center=False):
        if population == 'all':
            population = ['flat', 'disk', 'corona']
        r0 = 8.5
        a = (r0 ** 2 + height ** 2) ** 0.5

        # [e, a0, M, N, h, k], for corona + a00=xka0, lgx = 1.4
        consts = {'disk+': [0.1, 4.177, 6.872, 1.221, 6.2771, 0.3195],
                  'disk-': [0.47, 0.893, -0.314, 1.221, 6.2771, 0.3195],
                  'flat+': [0.02, 5.372, 0.637, 0.546, 1.6990, 1.0626],
                  'flat-': [0.04, 3.209, -0.227, 0.546, 1.6990, 1.0626],
                  'corona': [1, 60, 200, 0.5, 8.3206, 0.25817, 400],
                  'halo': [0.64, 1.096, 0.457, 7.072, 4.155 * 10 ** 6, 1.925 * 10 ** -9],
                  'buldge': [0.44, 0.211, 0.340, 0.417, 1.3727, 1.2442]}
        p = 0
        if 'corona' in population:
            e, a0, M, N, h, k, a00 = consts['corona']
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            if a > a00:
                p += 0
            else:
                _ = (1 + (a / (k * a0)) ** 2) ** (-N) - (1 + (a00 / (k * a0)) ** 2) ** (-N)
                p += p0 * (_ ** (1 / N))
                # print('corona', p0, _, p0 * _)
        if 'disk' in population:
            e, a0, M, N, h, k = consts['disk+']
            c = a0 * (1 - e ** 2) ** 0.5
            a_disk = 0.5 * (((r0 - c) ** 2 + height ** 2) ** 0.5 + ((r0 + c) ** 2 + height ** 2) ** 0.5)
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a_disk / (k * a0)) ** (1 / N))
            # print('disk+', p0, _, p0 * _)
            p += p0 * _
            e, a0, M, N, h, k = consts['disk-']
            c = a0 * (1 - e ** 2) ** 0.5
            a_disk = 0.5 * (((r0 - c) ** 2 + height ** 2) ** 0.5 + ((r0 + c) ** 2 + height ** 2) ** 0.5)
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a_disk / (k * a0)) ** (1 / N))
            # print('disk-', p0, _, p0 * _)
            p += p0 * _
        if 'flat' in population:
            e, a0, M, N, h, k = consts['flat+']
            c = a0 * (1 - e ** 2) ** 0.5
            a_flat = 0.5 * (((r0 - c) ** 2 + height ** 2) ** 0.5 + ((r0 + c) ** 2 + height ** 2) ** 0.5)
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a_flat / (k * a0)) ** (1 / N))
            # print('flat+', p0, _, p0 * _)
            p += p0 * _
            e, a0, M, N, h, k = consts['flat-']
            c = a0 * (1 - e ** 2) ** 0.5
            a_flat = 0.5 * (((r0 - c) ** 2 + height ** 2) ** 0.5 + ((r0 + c) ** 2 + height ** 2) ** 0.5)
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a_flat / (k * a0)) ** (1 / N))
            # print('flat-', p0, _, p0 * _)
            p += p0 * _
        if center:
            e, a0, M, N, h, k = consts['halo']
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a / (k * a0)) ** (1 / N))
            # print('halo', p0, _, p0 * _)
            p += p0 * _
            e, a0, M, N, h, k = consts['buldge']
            p0 = (h * M) / (4 * math.pi * e * (a0 ** 3))
            _ = math.exp(-(a / (k * a0)) ** (1 / N))
            # print('buldge', p0, _, p0 * _)
            p += p0 * _
        return p

    def create_realisitc_set_of_stars(self, filename):
        def fill_thin_cylinder(level, file_stars):
            #conc = self.einasto_distibution(level / 1000, 'all')
            disk_conc = self.einasto_distibution(level / 1000, 'disk')
            flat_conc = self.einasto_distibution(level / 1000, 'flat')
            corona_conc = self.einasto_distibution(level / 1000, 'corona')

            sum_mass = 0
            volume = math.pi * (1) * self.r * self.r
            cnt_local = 0
            ct = 0
            star_type = 'd'
            lst = [(disk_conc, 'd'), (flat_conc, 'f'), (corona_conc, 'c')]
            lst.sort(reverse=True)


            conc = 0
            for i in range(3):
                conc += lst[i][0]
                star_type = lst[i][1]
                while sum_mass < conc * volume:
                    ct += 1
                    x = self.r - r.random() * 2 * self.r
                    y = self.r - r.random() * 2 * self.r
                    if self.is_point_in((x, y)):
                        mass = self.get_good_random_mass()
                        h = (level + r.random())
                        # x = self.short_num(x, 4)
                        # y = self.short_num(y, 4)
                        # h = self.short_num(h, 4)
                        # mass = self.short_num(mass, 4)

                        #if level % 1000 == 0 and star_type == 'c': print('yes')

                        print(x, y, h, mass, star_type, file=file_stars)

                        sum_mass += mass
                        cnt_local += 1
                        sss[star_type] += 1
            return cnt_local, sum_mass

        data_name = self.path + filename

        data = open(data_name, 'w')
        cnt = 0
        # print('Creating of stars started')
        start_time = t.time()
        height = 0
        mass_of_part = 0
        n = 390000

        while height < n:
            cntd, mass_of_partd = fill_thin_cylinder(height, data)
            mass_of_part += mass_of_partd
            cnt += cntd
            if height % 25 == 0:
                pass
                # print('h:', height, ', stars now:', cnt, ', time:', int((t.time() - start_time) * 10) / 10, 's')
            height += 1
            if height % 1000 == 0:
                data.close()
                data = open(data_name, 'a')
            if height % (n // 100) == 0:
                mass_of_part = 0
        #print(sss['d'], sss['f'], sss['c'])

        print(cnt, 'stars there', file=data)
        data.close()
        f = open(self.path + "_stars_config.txt", "w")
        print(data_name, file=f)
        print(cnt, file=f)
        f.close()
        # print('Creating of stars ended in', int(t.time() - start_time), 'seconds,', cnt, 'stars was generated')

    def is_point_in(self, point):
        return point[0] ** 2 + point[1] ** 2 <= self.r ** 2

    def get_good_random_mass_prep(self):
        size = 10 ** 10
        dm = 0.01
        mass_now = 0.01
        lst = []  # хранит пару (начиная с какого числа мы будем брать эту массу, масса, к которой будем прибавять dm*r)
        value = 0
        sum_mass = 0
        while mass_now < 63:
            lst.append([value, mass_now])
            value += self.get_probability(mass_now) * size
            sum_mass += self.get_probability(mass_now) * mass_now
            mass_now += dm
        self.random_list = [lst, value, dm, sum_mass]

    def get_good_random_mass(self):
        lst, value, dm, sm = self.random_list
        rand_ind = int(r.random() * (value - 1)) + 1
        mass = self.get_elem(lst, rand_ind) + r.random() * dm
        return mass

    @staticmethod
    def get_elem(lst, x):  # binary search in lst to get at which mass will be this value(aka x), returning mass
        right, left = len(lst) - 1, 0
        while right - left > 1:
            mid = (right + left) // 2
            if lst[mid][0] > x:
                right = mid
            else:
                left = mid
        return lst[left][1]

    @staticmethod
    def get_probability(mass):
        if 0.01 < mass < 0.08:
            return 0.158
        elif 0.08 <= mass < 1:
            return 0.158 * math.exp((math.log(mass / 0.079, 10) ** 2) / (-0.9577))
        elif 1 <= mass < 3.47:
            return 0.044 * (mass ** (-4.37))
        elif 3.47 <= mass < 18.2:
            return 0.015 * (mass ** (-3.53))
        elif 18.2 <= mass < 63:
            return 2.5 * 10 ** -4 * (mass ** (-2.11))
        elif mass <= 0.01 or mass >= 63:
            return 10 ** (-8)

    @staticmethod
    def short_num(num, digit=2):
        mod = 10 ** digit
        num = int(num * mod) / mod
        return num


gen = Generator()
gen.create_realisitc_set_of_stars("stars.txt")
