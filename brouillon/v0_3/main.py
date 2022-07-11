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
from perf import Performance
from sedimentTransport.sedimentTransportLaw import SedimentTransportLaw 
from sedimentTransport.rickenmann1990 import Rickenmann1990 
from sedimentTransport.rickenmann1991 import Rickenmann1991
from sedimentTransport.lefort2015 import Lefort2015
from sedimentTransport.lefort1991 import Lefort1991


def main():
    ############################################################################################## FOR HECRAS COMP

    ############################################################################################# steep contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 5, 10]
    # slopes = [0.02, 0.02, 0.02]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile([0, 101, 201], [100, 200, 300], 5, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################# mild contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [20, 18, 20]
    # slopes = [0.001, 0.001, 0.001]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile([0, 405, 805], [400, 800, 1200], 5, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
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
    # profile = build_profile([0, 205, 405], [200, 400, 600], 5, widths, slopes, granulometries, manning=0.03, z_min_diff=0)
    ############################################################################################# new hecras comp
    # profile = Profile.import_profile("./results/hecras_comp/hecras_mix.pkl")
    # profile = Profile.import_profile("./results/hecras_comp/hecras_mild_contraction.pkl")
    # profile = Profile.import_profile("./results/hecras_comp/hecras_steep_contraction.pkl")
    # profile = Profile.import_profile("./results/hecras_comp/hecras_drops.pkl")
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
    # profile = build_profile([0, 251, 501, 751], [250, 500, 750, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################## M - M - S - M
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.001, 0.002, 0.05, 0.001]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 251, 501, 751], [250, 500, 750, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=2)
    # profile_copy = build_profile([0, 251, 501, 751], [250, 500, 750, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=2)
    ############################################################################################# Contraction
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 3, 10]
    # slopes = [0.03, 0.03, 0.001]
    # granulometries = [granulometry, granulometry, granulometry]
    # height_drops = [0, 0]
    # profile = build_profile([0, 301, 601], [300, 600, 900], 10, widths, slopes, granulometries, height_drops=height_drops , manning=0.013, z_min_diff=5)
    # profile_copy = build_profile([0, 301, 601], [300, 600, 900], 10, widths, slopes, granulometries, height_drops=height_drops , manning=0.013, z_min_diff=0) 
    ############################################################################################# Modele reduit
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [0.075, 0.0375, 0.075]
    # slopes = [0.0261, 0.0261, 0.0261]
    # granulometries = [granulometry, granulometry, granulometry]
    # profile = build_profile([0, 0.72, 0.82], [0.71, 0.81, 1.22], 0.01, widths, slopes, granulometries, manning=0.013, z_min_diff=0) # 71cm, 9cm, 40cm
    ############################################################################################# with drops
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.05, 0.001, 0.001, 0.05]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # height_drops = [10, 0, 0]
    # profile = build_profile([0, 251, 501, 751], [250, 500, 750, 1000], 10, widths, slopes, granulometries, height_drops=height_drops, manning=0.013, z_min_diff=0)
    ############################################################################################# triple jump
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8, 8, 8]
    # slopes = [0.05, 0.001, 0.05, 0.001, 0.05, 0.001]
    # granulometries = [granulometry, granulometry, granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 303, 603, 903, 1203, 1503], [300, 600, 900, 1200, 1500, 1800], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=0, z0=1000)
    ############################################################################################# very very very mild
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8]
    # slopes = [0.15, 0.06] # 5% 0.01% vs 15%, 6% pour laves torrentielles
    # granulometries = [granulometry, granulometry]
    # profile = build_profile([0,  501], [500, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=10)
    # profile_copy = build_profile([0,  501], [500, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=10)
    ############################################################################################# Null or Negative slope test
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.03, 0.01, -0.005, 0.05]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 251, 501, 751], [250, 500, 750, 1000], 10, widths, slopes, granulometries, manning=0.013, z_min_diff=0)
    ############################################################################################## Adverse slope 1
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [8, 8, 8, 8]
    # slopes = [0.02, 0.05, 0.001, 0.05] # 2% 5% 0.1% 5%
    # height_drops = [0, 0, -3]
    # granulometries = [granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 260, 510, 760], [250, 500, 750, 1000], 10, widths, slopes, granulometries, height_drops=height_drops, manning=0.013, z_min_diff=0)
    ############################################################################################## drops streak
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 10, 5, 5, 5, 10, 10]
    # slopes = [0.2, 0.2, 0.2, 0.1, 0.1, 0.1, 0.1]
    # height_drops = [5, 6, 7, 4, 8, 5]
    # granulometries = [granulometry, granulometry, granulometry, granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 201, 401, 601, 801, 1001, 1201], [200, 400, 600, 800, 1000, 1200, 1400], 10, widths, slopes, granulometries,height_drops=height_drops, manning=0.013, z_min_diff=0)
    ############################################################################################## drops streak reduction
    granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    widths = [10, 10, 6, 6, 10]
    slopes = [0.2, 0.002, 0.002, 0.001, 0.001]
    height_drops = [0, -6, 6, 6]
    granulometries = [granulometry, granulometry, granulometry, granulometry, granulometry]
    profile = build_profile([0, 101, 201, 401, 601], [100, 200, 400, 600, 800], 10, widths, slopes, granulometries,height_drops=height_drops, manning=0.013, z_min_diff=2)
    ############################################################################################## Mix
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10, 10, 5, 5, 5, 5, 10]
    # slopes = [0.2, 0.001, 0.001, 0.05, 0.001, 0.01, 0.01]
    # height_drops = [0, 0, 0, 2, -8, 0]
    # granulometries = [granulometry, granulometry, granulometry, granulometry, granulometry, granulometry, granulometry]
    # profile = build_profile([0, 201, 401, 601, 801, 1001, 1201], [200, 400, 600, 800, 1000, 1200, 1400], 20, widths, slopes, granulometries, height_drops=height_drops, manning=0.013, z_min_diff=0)
    ############################################################################################# debug examples
    # profile = Profile.import_profile("./results/profile_test_1.pkl")
    # profile = Profile.import_profile("./results/profile_test_2.pkl")
    # profile = Profile.import_profile("./results/test_boucle_infinie_1.pkl")
    # profile = Profile.import_profile("./results/seuil_bizarre.pkl")
    # profile = Profile.import_profile("./results/seuil_bizarre_2.pkl")
    # profile = Profile.import_profile("./results/seuil_bizarre_3.pkl")
    # profile = Profile.import_profile("/home/samuel/Documents/3a/pfe/ONF/stage_evofond/evofond/brouillon/v0_3/results/wtf_dig.pkl")
    # profile = Profile.import_profile("./results/seuil_bizarre_reduit_1.pkl")
    # profile = Profile.import_profile("./results/pic_first_method.pkl")
    # profile = Profile.import_profile("./results/increased_head.pkl") # Q = 40 !
    # profile = Profile.import_profile("./results/instability_amplification.pkl")
    ##############################################################################################

    # x_compare, y_compare = read_hecras_data("../../data_hecras/new/steep_contraction.txt")
    # x_compare, y_compare = reverse_data(x_compare, y_compare)
    # print([(s.get_x(), s.get_z(), s.get_b()) for s in profile.get_section_list()])

    # profile.complete(10)
    # Q_list = [30 for _ in range(profile.get_nb_section())]
    # Q = 20 # 10/3600
    # y_list = profile.compute_depth_bis(Q, plot=True, hydraulic_jump_analysis=False, method="RungeKutta", friction_law="Manning-Strickler", compare=(x_compare, y_compare))
    # y_list = profile.compute_depth(Q, plot=True, hydraulic_jump_analysis=True, compare=(x_compare, y_compare), method="ImprovedEuler", friction_law="Manning-Strickler") 
    # profile.plot(y=y_list)
    # profile.export("./hecras_steep_contraction.pkl")

    # for _ in range(40):
    #     if (_+1)%10!=0:
    #         profile.update_bottom(Q, y_list, QsIn0=0, dt=1, law=Rickenmann1991(), plot=False)
    #     else:
    #         profile.update_bottom(Q, y_list, QsIn0=0, dt=1, law=Rickenmann1991(), plot=True)
    #     y_list = profile.compute_depth_bis(Q, plot=False, hydraulic_jump_analysis=False, method="ImprovedEuler")
    #     # y_list = profile.compute_depth(Q, plot=False, hydraulic_jump_analysis=True, compare=None, method="ImprovedEuler")
    # for _ in range(0):
    #     profile.update_bottom(Q, y_list, QsIn0=0, dt=1, law=Rickenmann1991(), plot=True)
    #     y_list = profile.compute_depth_bis(Q, plot=False, hydraulic_jump_analysis=False, method="ImprovedEuler")
    #     # y_list = profile.compute_depth(Q, plot=False, hydraulic_jump_analysis=True, compare=None, method="ImprovedEuler")

    dt = 10
    t = np.arange(0, 6*60*60+dt, dt)
    Qmax = 60
    tm = 2*60*60
    alpha = 2
    Qbase = 30
    hydrogram = hydrogrammeLavabre(Qmax,tm,alpha,Qbase,t)
    law = Rickenmann1991()
    ani = profile.compute_event_bis(hydrogram, t, law, backup=True, debug=False, method="ImprovedEuler", friction_law="Ferguson", cfl=20)
    return ani

    ############################################################################################## lave torrentielle
    # # fixed parameters
    # granulometry = Granulometry(dm=0.065, d30=0.03, d50=0.03, d90=0.065, d84tb=0.04, d84bs=0.08, Gr=2.6)
    # widths = [10]
    # height_drops = []
    # granulometries = [granulometry]
    # x_begin_list = [0]
    # x_end_list = [100]

    # # variables
    # slopes = np.arange(0.07, 0.15, 0.001)
    # Q_list = [100, 100, 100, 180, 180, 180, 300, 300, 300]
    # tauc_over_rho_list = [0.5, 1, 1.5, 0.5, 1, 1.5, 0.5, 1, 1.5]

    # nb_test = len(Q_list)
    # plt.figure()
    # for i in range(nb_test):
    #     y_list = []
    #     for s in slopes:
    #         profile = build_profile(x_begin_list, x_end_list, 5, widths, [s], granulometries, tauc_over_rho=tauc_over_rho_list[i])
    #         yn_list = profile.get_yn_list(Q_list[i], friction_law="Coussot")
    #         y_list.append(yn_list[0])
    #     plt.plot(slopes, y_list, label=fr"Q = {Q_list[i]}$m^3/s$ - $\tau_c/\rho = ${tauc_over_rho_list[i]}$m^2/s^2$")
    # plt.legend()

    return



def build_profile(x_begin_list, x_end_list, dx, widths, slopes, granulometries, z0=200, height_drops=None, z_min_diff=0, manning=0.013, tauc_over_rho=1, K_over_tauc=0.3, section="Rectangular"):
    """
    Made by and for devs.
    """
    x = []
    limit_index = []
    for i in range(len(x_begin_list)):
        x += list(np.arange(x_begin_list[i], x_end_list[i]-0.5*dx, dx)) + [x_end_list[i]]
        limit_index.append(len(x))
    n = len(x)
    height_drops = [0 for _ in range(len(x_begin_list)-1)] if height_drops==None else height_drops
    z = [z0]
    b = []
    granulometry = []
    old_part_index = 0
    for i in range(n):
        part_index = old_part_index if i < limit_index[old_part_index] else old_part_index+1
        b.append(widths[part_index])
        granulometry.append(granulometries[part_index])
        if i == 0:
            old_part_index = part_index
            continue
        else:
            z_new = z[-1]-(x[i]-x[i-1])*slopes[old_part_index]
            if part_index > old_part_index:
                z_new -= height_drops[old_part_index]
            z.append(z_new)
            old_part_index = part_index

    section_list = []
    for i in range(n):
        if section == "Rectangular":
            s = RectangularSection(x[i], z[i], b[i], y_max=50, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Trapezoidal":
            s = TrapezoidalSection(x[i], z[i], b[i], 1, y_max=10, z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        elif section == "Irregular":
            points = [(0, 50), (0, 0), (b[i], 0), (b[i], 50)] # rectangular, but it will be treated as an irregular section
            s = IrregularSection(points, x[i], z[i], z_min=z[i]-z_min_diff, granulometry=granulometry[i], manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        section_list.append(s)
    profile = Profile(section_list)
    return profile























if __name__=="__main__":
    Performance.start()
    ani = main()
    Performance.stop()
    Performance.print_perf()
    Performance.save_perf("./results/perf/perf.txt")
    plt.show()