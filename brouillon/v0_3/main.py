#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from time import time
from utils import hydrogrammeLavabre, read_hecras_data, reverse_data
from irregularSection import IrregularSection
from rectangularSection import RectangularSection
from trapezoidalSection import TrapezoidalSection
from profile import Profile
from granulometry import Granulometry
from sedimentTransport.sedimentTransportLaw import SedimentTransportLaw 
from sedimentTransport.rickenmann1990 import Rickenmann1990 
from sedimentTransport.rickenmann1991 import Rickenmann1991


def main():
    ############################################################################################## FOR HECRAS COMP

    ############################################################################################# steep contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [20, 6, 20]
    # slopes = [0.05, 0.05, 0.05]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile(0, 500, 5, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################# mild contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [20, 18, 20]
    # slopes = [0.001, 0.001, 0.001]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1200, 5, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################# very very very mild
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8]
    # slopes = [0.05, 0.0001]
    # granulometries = [granulometry, granulometry]
    # profile = build_profile(0, 800, 5, widths, slopes, granulometries, z_min_diff=2)
    ############################################################################################# Adverse slope
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 10, 10]
    # slopes = [0.03, -0.01, 0.03]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile(0, 600, 5, widths, slopes, granulometries, manning=0.03, z_min_diff=0)
    #############################################################################################

    ############################################################################################## 

    ############################################################################################## Irregular section (profile simple)
    # x = np.arange(0, 1000, 10)
    # points = [(0, 100), (5, 0), (15, 0), (20, 100)] # \_/
    # section_list = []
    # s0 = 0.05
    # z0 = 1000
    # z = z0-s0*x
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # # granulometry = Granulometry(dm=0.041, d30=0.0229, d50=0.034, d90=0.071, d84tb=0.063, d84bs=0.063, Gr=2)
    # # granulometry = Granulometry(dm=0.065, d30=0.01, d50=0.03, d90=0.09, d84tb=0.04, d84bs=2.2, Gr=2.5)
    # for i in range(len(x)):
    #     section_new = IrregularSection(points, x[i], z[i], z_min=z[i], granulometry=granulometry, s0=s0, manning=0.013)
    #     section_list.append(section_new)
    # profile = Profile(section_list)
    ############################################################################################## Rectangular section (profile simple)
    # x = np.arange(0, 500, 10)
    # section_list = []
    # s0 = 0.1
    # z0 = 1000
    # z = z0-s0*x
    # b0 = 8
    # b = b0*np.ones(len(x))
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # # granulometry = Granulometry(dm=0.041, d30=0.0229, d50=0.034, d90=0.071, d84tb=0.063, d84bs=0.063, Gr=2)
    # # granulometry = Granulometry(dm=0.065, d30=0.01, d50=0.03, d90=0.09, d84tb=0.04, d84bs=2.2, Gr=2.5)
    # for i in range(len(x)):
    #     section_list.append(RectangularSection(x[i], z[i], b[i], z_min=z[i], granulometry=granulometry, s0=s0))
    # profile = Profile(section_list)
    ############################################################################################## Rectangular section, profile complex, S - S - M -S
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.02, 0.05, 0.001, 0.05]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1000, 10, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################## M - M - S - M
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.001, 0.002, 0.05, 0.001]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1000, 10, widths, slopes, granulometries, z_min_diff=1)
    ############################################################################################# Contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [20, 4, 20]
    # slopes = [0.02, 0.02, 0.02]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile(0, 900, 20, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################# with drops
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.05, 0.001, 0.001, 0.05]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # height_drops = [10, 0, 0]
    # profile = build_profile(0, 800, 10, widths, slopes, granulometries, height_drops=height_drops, z_min_diff=1)
    ############################################################################################# double jump
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.05, 0.001, 0.05, 0.001]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1000, 10, widths, slopes, granulometries, z_min_diff=2)
    ############################################################################################# very very very mild
    granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    widths = [8, 8]
    slopes = [0.05, 0.0001]
    granulometries = [granulometry, granulometry]
    profile = build_profile(0, 300, 20, widths, slopes, granulometries, z_min_diff=10)
    ############################################################################################# Null or Negative slope test
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.03, 0.01, -0.005, 0.05]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1000, 20, widths, slopes, granulometries, z_min_diff=2)
    ############################################################################################## Adverse slope 1
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.02, 0.05, 0.001, 0.05] # 2% 5% 0.1% 5%
    # height_drops = [0, 0, -3]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 1000, 10, widths, slopes, granulometries, z_min_diff=1, height_drops=height_drops)
    # ############################################################################################## big one
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 10, 10, 10, 10, 10, 10]
    # slopes = [0.03, 0.03, 0.001, -0.01, 0.001, 0.02, 0.07]
    # height_drops = [0, 10, -1, 0, 0, -6]
    # granulometries = [granulometry, granulometry, granulometry, granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile(0, 7000, 20, widths, slopes, granulometries, manning=0.03, z_min_diff=2, height_drops=height_drops)
    ############################################################################################# Bin
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 10]
    # slopes = [0.05, 0.001]
    # height_drops = [0, 0]
    # granulometries = [granulometry, granulometry]
    # profile = build_profile(0, 600, 10, widths, slopes, granulometries, manning=0.03, z_min_diff=0, height_drops=height_drops, section="Trapezoidal")
    ############################################################################################## 

    # x_compare, y_compare = read_hecras_data("../../data_hecras/clean/adverse_slope.txt")
    # x_compare, y_compare = reverse_data(x_compare, y_compare)
    # print([(s.get_x(), s.get_z(), s.get_b()) for s in profile.get_section_list()])

    profile.complete(10)
    # profile.plot3D()
    # Q = 30
    # y_list = profile.compute_depth(Q, plot=True, hydraulic_jump_analysis=True, compare=None)

    dt = 10
    t = np.arange(0, 6*60*60, dt)
    Qmax = 50
    tm = 1*60*60
    alpha = 0
    Qbase = 30
    hydrogram = hydrogrammeLavabre(Qmax,tm,alpha,Qbase,t)
    law = Rickenmann1991()
    ani = profile.compute_event(hydrogram, dt, law, backup=True, debug=False, smooth=False)

    plt.show()
    return


def build_profile(x_begin, x_end, dx, widths, slopes, granulometries, height_drops=None, z_min_diff=0, manning=0.013, section="Rectangular"):
    """
    Made by and for devs.
    """
    x = np.arange(x_begin, x_end, dx)
    n = len(x)
    part_nb = len(widths)
    part_len = int(np.ceil(n/part_nb))
    height_drops = [0 for _ in range(part_nb-1)] if height_drops==None else height_drops
    z = [200]
    b = []
    granulometry = []
    old_part_index = 0
    for i in range(n):
        part_index = i // part_len
        b.append(widths[part_index])
        granulometry.append(granulometries[part_index])
        if i == 0:
            old_part_index = part_index
            continue
        else:
            z_new = z[-1]-dx*slopes[part_index]
            if part_index > old_part_index:
                z_new -= height_drops[old_part_index]
            z.append(z_new)
            old_part_index = part_index

    section_list = []
    for i in range(n):
        if section == "Rectangular":
            s = RectangularSection(x[i], z[i], b[i], y_max=50, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning)
        elif section == "Trapezoidal":
            s = TrapezoidalSection(x[i], z[i], b[i], 1, y_max=50, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning)
        elif section == "Irregular":
            points = [(0, 50), (0, 0), (b[i], 0), (b[i], 50)] # rectangular, but it will be treated as an irregular section
            s = IrregularSection(points, x[i], z[i], z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning)
        section_list.append(s)
    profile = Profile(section_list)
    return profile























if __name__=="__main__":
    main()
