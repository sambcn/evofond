import numpy as np
import copy
import scipy.optimize as op
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from time import time
from irregularSection import IrregularSection
from utils import Y_MIN, G, get_matrix_max, real_roots_cubic_function, read_hecras_data, reverse_data

class Profile():

    def __init__(self, section_list):
        if len(section_list) < 2:
            raise(ValueError("You need at least 2 sections to initialize a profile"))
        self.__upstream = section_list[0]
        self.__downstream = section_list[-1]
        self.__section_list = section_list
        self.setup_section_list()

    def complete(self, dx, starting_index=0, ending_index=-1):
        """
        Add new sections to the profile between the starting_index and the ending_index, such that there is at least one section every dx in this interval.
        This methods do not change upstream et downstream attributes.
        It needs the section_list to be in a coherent state (see more details in setup_section_list method).
        """
        if ending_index==-1:
            ending_index = len(self.__section_list)-1
        section_sublist = self.__section_list[starting_index:ending_index+1]
        new_list = [section_sublist.pop(0)]
        x_current = self.get_section(starting_index).get_x()
        x_down = self.get_section(ending_index).get_x()
        while x_current < x_down:
            if x_current + dx >= section_sublist[0].get_x():
                new_section = section_sublist.pop(0)
                x_current = new_section.get_x()
                new_section.set_up_section(new_list[-1])
                new_list.append(new_section)
            else:
                x_current += dx
                up_section = new_list[-1]
                down_section = section_sublist[0]
                new_section = up_section.interp_as_up_section(down_section, x_current)
                new_list.append(new_section)
                up_section.set_down_section(new_section)
        self.__section_list[starting_index:ending_index+1] = new_list

    def setup_section_list(self):
        """
        Sort sections by x ascending, setup the upSections / downSections / is_upstream / is_downstream attributes.
        """
        self.__section_list.sort(key=lambda section: section.get_x()) 
        section_list = self.__section_list
        if len(section_list) > 1:
            for i in range(len(section_list)):
                section = section_list[i]
                if i==0:
                    section.set_is_upstream(True)
                    section.set_is_downstream(False)
                    # section.set_up_section(section) # This is automatically done when setting is_upstream on True
                    section.set_down_section(section_list[i+1])
                elif i == len(section_list) - 1:
                    section.set_is_upstream(False)
                    section.set_is_downstream(True)
                    section.set_up_section(section_list[i-1])
                    # section.set_down_section(section) # This is automatically done when setting is_upstream on True
                else:
                    section.set_is_upstream(False)
                    section.set_is_downstream(False)
                    section.set_up_section(section_list[i-1])
                    section.set_down_section(section_list[i+1])

    def smooth(self):
        """
        used to erase growing instabilities.
        When it detects a pattern such like "z increase, z decrease, z increase" then it keeps every following section which keep this alternation, and replace their z by the "mean" value.
        Look at documentation for more details.
        """
        i = 0
        while i < self.get_nb_section()-3:
            j = i
            while j < self.get_nb_section()-3:
                s0_1 = self.get_section(j).get_S0()
                s0_2 = self.get_section(j+1).get_S0()
                s0_3 = self.get_section(j+2).get_S0()
                if (s0_1-s0_2)*(s0_2-s0_3) < 0:
                    j += 1
                else:
                    break
            if j > i:
                x0 = self.get_section(i).get_x()
                z0 = self.get_section(i).get_z()
                xf = self.get_section(j+1).get_x()
                zf = self.get_section(j+1).get_z()
                for k in range(i+1, j+1):
                    s = self.get_section(k)
                    try:
                        s.set_z(np.interp(s.get_x(), [x0, xf], [z0, zf]))
                    except ValueError as e:
                        s.set_z(s.get_z_min())
            i = j+1
        return

    def copy(self):
        """return a safe copy of this section"""
        section_list = []
        for s in self.__section_list:
            section_list.append(s.copy())
        copied_profile = Profile(section_list)
        return copied_profile

    # hydraulic methods

    def compute_depth(self, Q, plot=False, hydraulic_jump_analysis=False, compare=None): # Trying to accelerate computations
        """
        Compute the water depth at every section of the profile, for an upstream water discharge Q. Use standard step method.
        """
        k = 0.1 # Borda coef 
        # yn_list = self.get_yn_list(Q)
        yc_list = self.get_yc_list(Q)
        n = len(self.__section_list)

        # supercritical computation
        y_list_supercritical = [self.get_upstream_boundary_condition(Q)]
        hs_list_supercritical = [self.__upstream.get_Hs(Q, y_list_supercritical[0])]
        i_current = 0
        while i_current < n-1:
            current_section = self.get_section(i_current)
            next_section = self.get_section(i_current+1)

            slope = current_section.get_S0() - current_section.get_Sf(Q, y_list_supercritical[-1])
            current_hs = hs_list_supercritical[-1]
            dx = next_section.get_x() - current_section.get_x()
            hs_next = current_hs + dx*slope
            if hs_next + next_section.get_z() > hs_list_supercritical[-1] + current_section.get_z():
                hs_next = hs_list_supercritical[-1] + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, yc_list[i_current+1]):
                hs_next = next_section.get_Hs(Q, yc_list[i_current+1])

            # roots = np.roots([1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2)]) # ONLY FOR RECTANGULAR
            # positive_real_roots = []
            # for r in roots:
            #     if r.imag < 1e-3 and r.real > 0:
            #         positive_real_roots.append(r.real)
            # positive_real_roots.sort()
            positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            y_supercritical = min(positive_real_roots)

            # improved euler

            next_slope = next_section.get_S0(up_direction=True) - next_section.get_Sf(Q, y_supercritical)
            hs_next = current_hs + 0.5*dx*(slope+next_slope)

            if hs_next + next_section.get_z() > hs_list_supercritical[-1] + current_section.get_z():
                hs_next = hs_list_supercritical[-1] + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, yc_list[i_current+1]):
                hs_next = next_section.get_Hs(Q, yc_list[i_current+1])

            # roots = np.roots([1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2)]) # ONLY FOR RECTANGULAR
            # positive_real_roots = []
            # for r in roots:
            #     if r.imag < 1e-5 and r.real > 0:
            #         positive_real_roots.append(r.real)
            # positive_real_roots.sort()
            positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            y_supercritical = min(positive_real_roots)

            # TEST BORDA
            # hs_next = hs_next - k*abs(next_section.get_V(Q, y_supercritical)**2-current_section.get_V(Q, y_list_supercritical[-1]))/(2*G)
            # if hs_next + next_section.get_z() > hs_list_supercritical[-1] + current_section.get_z():
            #     hs_next = hs_list_supercritical[-1] + current_section.get_z() - next_section.get_z()
            # if hs_next < next_section.get_Hs(Q, yc_list[i_current+1]):
            #     hs_next = next_section.get_Hs(Q, yc_list[i_current+1])
            # positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            # y_supercritical = min(positive_real_roots)



            hs_list_supercritical.append(hs_next)
            y_list_supercritical.append(y_supercritical)
            i_current += 1

        # subcritical computation
        y_list_subcritical = [self.get_downstream_boundary_condition(Q)]
        hs_list_subcritical = [self.__downstream.get_Hs(Q, y_list_subcritical[-1])]
        i_current = n-1
        while i_current > 0:
            current_section = self.get_section(i_current)
            next_section = self.get_section(i_current-1)

            slope = current_section.get_S0(up_direction=True) - current_section.get_Sf(Q, y_list_subcritical[0])
            current_hs = hs_list_subcritical[0]
            dx = next_section.get_x() - current_section.get_x()
            hs_next = current_hs + dx*slope
            if hs_next + next_section.get_z() < hs_list_subcritical[0] + current_section.get_z():
                hs_next = hs_list_subcritical[0] + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, yc_list[i_current-1]):
                hs_next = next_section.get_Hs(Q, yc_list[i_current-1])

            # roots = np.roots([1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2)]) # ONLY FOR RECTANGULAR
            # positive_real_roots = []
            # for r in roots:
            #     if r.imag < 1e-5 and r.real > 0:
            #         positive_real_roots.append(r.real)
            # positive_real_roots.sort()
            positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            y_subcritical = max(positive_real_roots)

            # improved euler

            next_slope = next_section.get_S0() - next_section.get_Sf(Q, y_subcritical)
            hs_next = current_hs + 0.5*dx*(slope+next_slope)

            if hs_next + next_section.get_z() < hs_list_subcritical[0] + current_section.get_z():
                hs_next = hs_list_subcritical[0] + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, yc_list[i_current-1]):
                hs_next = next_section.get_Hs(Q, yc_list[i_current-1])

            # roots = np.roots([1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2)]) # ONLY FOR RECTANGULAR
            # positive_real_roots = []
            # for r in roots:
            #     if r.imag < 1e-5 and r.real > 0:
            #         positive_real_roots.append(r.real)
            # positive_real_roots.sort()
            positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            y_subcritical = max(positive_real_roots)

            # TEST BORDA
            # hs_next = hs_next + k*abs(next_section.get_V(Q, y_supercritical)**2-current_section.get_V(Q, y_list_subcritical[0]))/(2*G)
            # if hs_next + next_section.get_z() < hs_list_subcritical[0] + current_section.get_z():
            #     hs_next = hs_list_subcritical[0] + current_section.get_z() - next_section.get_z()
            # if hs_next < next_section.get_Hs(Q, yc_list[i_current-1]):
            #     hs_next = next_section.get_Hs(Q, yc_list[i_current-1])
            # positive_real_roots = real_roots_cubic_function(1, -hs_next, 0, Q**2 / (2*G*next_section.get_b()**2))
            # y_supercritical = max(positive_real_roots)


            hs_list_subcritical.insert(0, hs_next)
            y_list_subcritical.insert(0, y_subcritical)
            i_current -= 1

        y_list = []
        hydraulic_index = [] # list of the index of the section where an hydraulic jump occurs

        # flow_supercritical = None
        # for i in range(self.get_nb_section()):
        #     if hs_list_subcritical[i] >=  hs_list_supercritical[i]:
        #         j = i-1
        #         y_list.append(y_list_subcritical[i])
        #         while j>0 and (self.get_section(j).get_Fs(Q, y_list_subcritical[j]) > self.get_section(j).get_Fs(Q, y_list_supercritical[j])):
        #             y_list[j] = (y_list_subcritical[j])
        #             j -= 1 
        #         if flow_supercritical != None and flow_supercritical:
        #             hydraulic_index.append(j)
        #         flow_supercritical = False
        #     else:
        #         y_list.append(y_list_supercritical[i])
        #         flow_supercritical = True

        # flow_supercritical = None   
        # for i, s in enumerate(self.__section_list):
        #     if flow_supercritical != None and not(flow_supercritical):
        #         y_list.append(y_list_subcritical[i])
        #         flow_supercritical_new = not(y_list_subcritical[i] - yc_list[i] > 1e-4)
        #     else:
        #         Fs_supercritical = s.get_Fs(Q, y_list_supercritical[i])
        #         Fs_subcritical = s.get_Fs(Q, y_list_subcritical[i])
        #         if Fs_subcritical <= Fs_supercritical:
        #             y_list.append(y_list_supercritical[i])
        #             flow_supercritical_new = True 
        #         else:
        #             y_list.append(y_list_subcritical[i])
        #             flow_supercritical_new = False        
        #     if flow_supercritical != None and flow_supercritical and not flow_supercritical_new: #ressaut hydraulique
        #         hydraulic_index.append(i-1)
        #     flow_supercritical = flow_supercritical_new

        flow_supercritical = None   
        for i, s in enumerate(self.__section_list):
            Fs_supercritical = s.get_Fs(Q, y_list_supercritical[i])
            Fs_subcritical = s.get_Fs(Q, y_list_subcritical[i])
            if flow_supercritical == None:
                flow_supercritical = Fs_supercritical >= Fs_subcritical
                y_list.append(y_list_supercritical[i] if flow_supercritical else y_list_subcritical[i])
                continue
            changing_flow = (Fs_supercritical >= Fs_subcritical)^(flow_supercritical) # XOR operator. True if flow changes
            y_new = y_list_supercritical if (changing_flow and not(flow_supercritical)) or (not(changing_flow) and flow_supercritical) else y_list_subcritical
            y_current = y_list_supercritical if flow_supercritical else y_list_subcritical 
            if changing_flow:
                if s.get_H(Q, y_new[i]) <= s.get_up_section().get_H(Q, y_list[-1]):
                    y_list.append(y_new[i])
                    if flow_supercritical:
                        hydraulic_index.append(i-1)
                    flow_supercritical = not(flow_supercritical)
                else:
                    y_list.append(y_current[i])
            else:
                y_list.append(y_new[i])

        if plot:
            fig, axs = plt.subplots(2)
            axs[0].plot(self.get_x_list(), 0 + np.array(hs_list_supercritical), label="hs supercritical")
            axs[0].plot(self.get_x_list(), 0 + np.array(hs_list_subcritical), label="hs subcritical")
            axs[0].plot(self.get_x_list(), 0 + np.array([s.get_Hs(Q, yc_list[i]) for i, s in enumerate(self.__section_list)]), 'r--', label="hs critical")            
            axs[0].set_title("Comparaison de l'énergie spécifique entre le calcul torrentiel et fluvial")
            axs[0].legend()
            axs[1].plot(self.get_x_list(), [s.get_Fs(Q, y_list_supercritical[i]) for i, s in enumerate(self.__section_list)], label="Fs supercritical")
            axs[1].plot(self.get_x_list(), [s.get_Fs(Q, y_list_subcritical[i]) for i, s in enumerate(self.__section_list)], label="Fs subcritical")
            axs[1].set_title("Comparaison des forces spécifiques entre le régime torrentiel et fluvial")
            axs[1].legend()
            fig = self.plot(y=y_list_supercritical, Q=Q)
            fig = self.plot(y=y_list_subcritical, Q=Q)
            fig = self.plot(y=y_list, Q=Q, compare=compare)
            # self.plot3D(y=y_list)
            if hydraulic_jump_analysis:
                ax1 = fig.get_axes()[0]
                string = ""
                for i, i_section in enumerate(hydraulic_index):
                    upSection = self.get_section(i_section)
                    downSection = self.get_section(i_section+1)
                    y_up = y_list[i_section]
                    y_down = y_list[i_section+1]
                    Fr_up = upSection.get_Fr(Q, y_up)
                    Fr_down = upSection.get_Fr(Q, y_down)
                    deltaH = upSection.get_H(Q, y_up) - downSection.get_H(Q, y_down)
                    string += f"\nHydraulic jump n°{i+1} analysis : \n"
                    string += f"occurs at x = {upSection.get_x()}\n"
                    string += f"result : y2 = {y_down}\n"
                    string += f"theory : y2 = {-y_up*0.5+y_up*0.5*np.sqrt(1+8*Fr_up**2)}\n"
                    string += f"result : deltaH = {deltaH} \n"
                    string += f"theory : deltaH = {(y_down-y_up)**3 / (4*y_up*y_down)} \n"
                ax1.annotate(string, (self.get_x_min(), min(self.get_z_min_list())))
            
        return y_list

    def update_bottom(self, Q, y, QsIn0, dt, law, plot=False):
        """
        Q and dt are respectively the water discharge and the time step. 
        Compute the sediment discharge at the downstream and change the attribute __z for every section
        """
        z0 = self.get_z_list()
        V0 = self.get_stored_volume()
        QsIn = QsIn0
        if plot:
            fig = self.plot(Q=Q, y=y)
        y_up = y[0]
        for i, section in enumerate(self.__section_list):
            if section.is_downstream():
                y_down = y[i]
            else:
                y_down = y[i+1]
            QsIn = section.update_bottom(Q, y_up, y[i], y_down, QsIn, dt, law)
            y_up = y[i]

        if plot:
            x = self.get_x_list()
            ax1 = fig.get_axes()[0]
            ax1.plot(x, self.get_z_list(), "s-", label="new z")
            ax1.annotate(f"Sediment : \nIn = {QsIn0*dt}\nOut = {QsIn*dt}\nStored = {self.get_stored_volume()-V0}\nSum = {QsIn0*dt-QsIn*dt-(self.get_stored_volume()-V0)}", (self.get_x_min(), min(self.get_z_min_list())))
            plt.legend()
            ax2 = ax1.twinx()
            ax2.set_ylabel("difference of height between the start and the end")
            ax2.plot(x, np.array(self.get_z_list())-np.array(z0), color="pink", linestyle="dashdot")

        return QsIn

    def compute_event(self, hydrogram, dt, law, backup=False, debug=False, smooth=False):
        start_computation = time()
        y_matrix = [] # list of the water depth during the event
        z_matrix = [self.get_z_list()] # list of the bottom height during the event
        h_matrix = [] # list of head during the event
        V_in = 0 # solid volume gone into the profile during the event
        V_out = 0 # solid volume gone out of the profile
        stored_volume_start = self.get_stored_volume() # stored volume of sediment at the start of the event 
        initial_profile = self.copy()
        log_string = f"[backup={backup}, debug={debug}, smooth={smooth}]"

        if debug:
            profile_list = [initial_profile]
            current_stored_volume = stored_volume_start
            total_volume_difference = [] 
            one_step_volume_difference = []

        for i, Q in enumerate(hydrogram):
            print(f"{i}/{len(hydrogram)} "+log_string)

            # hydraulic computations
            try:
                y = self.compute_depth(Q)
            except ValueError as e:
                print("ERROR IN COMPUTING THIS EVENT : COULD NOT FINISH\n plotting the last state ...\n")
                self.plot(Q=Q, y=y)
                plt.show()
                break

            y_matrix.append(y)
            h_matrix.append([s.get_H(Q, y[i]) for i, s in enumerate(self.__section_list)])

            # solid transport
            QsIn0 = law.compute_Qs(initial_profile.get_upstream_section(), Q, y[0], y[1]) # Gonna change, it is a given parameter, chosen by users
            V_in += QsIn0*dt
            QsOut = self.update_bottom(Q, y, QsIn0, dt, law)
            V_out += QsOut*dt
            if smooth:
                self.smooth() # to smooth instabilities
            z_matrix.append(self.get_z_list())

            # debug
            if debug:
                total_volume_difference.append(V_in - V_out - (self.get_stored_volume() - stored_volume_start))
                one_step_volume_difference.append(QsIn0*dt - QsOut*dt - (self.get_stored_volume() - current_stored_volume))
                current_stored_volume = self.get_stored_volume()
                if abs(one_step_volume_difference[-1]) > 0.1:
                    print(f"WARNING : HUGE SEDIMENT CREATION/DISAPPEARANCE ON THE STEP OF TIME i={i}")
                profile_list.append(self.copy())

        try:       
            y_matrix.append(self.compute_depth(hydrogram[-1]))
            h_matrix.append([s.get_H(Q, y_matrix[-1][i]) for i, s in enumerate(self.__section_list)])
        except ValueError:
            pass
        stored_volume_end = self.get_stored_volume()
        end_computation = time()
        print(f"computation time = {end_computation-start_computation}")
        x = self.get_x_list()

        title = 'Sediment transport :\n' + \
            f'Volume gone in :  {V_in}\n' + \
            f'Volume gone out : {V_out}\n' + \
            f'Stored volume : {stored_volume_end - stored_volume_start}\n' + \
            f'Sum : {V_in - V_out - (stored_volume_end - stored_volume_start)}'
        fig0, axs = plt.subplots(2)
        fig0.suptitle(title)
        axs[0].plot(x, z_matrix[0], "r", label="z start")
        axs[0].plot(x, z_matrix[-1], "orange", label="z end")
        axs[0].plot(x, self.get_z_min_list(), "g--", label="zmin")
        axs[0].set(xlabel="x", ylabel="height (m)")
        # axs[0].plot(x, [s.get_H(Q, y_matrix[-1][i]) for i, s in enumerate(self.__section_list)], label="Energy grade line")
        axs[0].annotate(f"event of {len(hydrogram)*dt}s\ndt = {dt}s\nQmax = {max(hydrogram):.3f}m3/s\nQmean = {np.mean(hydrogram):.3f}m3/s\nsediment transport law = {str(law)}\ntime of computation = {end_computation-start_computation:.4f}s", (self.get_x_min(), min(self.get_z_min_list())))
        axs[1].plot(x, np.array(z_matrix[-1])-np.array(z_matrix[0]), "b", label="newz-z")
        axs[1].set(xlabel="x", ylabel="difference of height (m)")
        fig0.legend()
        fig0.set_size_inches(10.5, 9.5)
        if backup:
            print("saving result plot...")
            fig0.savefig("./results/result.png", dpi=200, format="png")
            print("plot saved")

        if debug:
            plt.figure()
            plt.plot([i*dt for i in range(len(total_volume_difference))], total_volume_difference, label="acumulated solid creation or disappearance")
            plt.plot([i*dt for i in range(len(one_step_volume_difference))], one_step_volume_difference, label="solid creation or disappearance during this step of time")
            plt.legend()
            plt.title('sediment creation or diseppearance due to numerical errors')

        fig, ax1 = plt.subplots() 
        line, = ax1.plot(x, np.array(y_matrix[0])+np.array(z_matrix[0]), label="water line")
        ax1.set_ylabel("height (m)", color="blue")
        annotation = ax1.annotate(f"Q={hydrogram[0]}", (self.get_x_min(), min(self.get_z_min_list())))
        ax1.plot(x, z_matrix[0], "r", label="z at the begining")
        # ax1.vlines(200, 180, 210)
        line2, = ax1.plot(x, z_matrix[0], "orange", label="z")
        ax1.plot(x, self.get_z_min_list(), "g--", label="zmin")
        line3, = ax1.plot(x, h_matrix[0], color="pink", label="energy grade line")
        ax1.set_ylim(min(self.get_z_min_list()), get_matrix_max(z_matrix)+get_matrix_max(y_matrix)) 

        plt.xlim(self.get_x_min(), self.get_x_max())
        plt.xlabel("x")
        ax1.set_title("Water depth and bottom evolution")
        fig.set_size_inches(9.5, 5.5)
        plt.legend()
        def animate(i, flag=True): 
            if debug and flag:
                i = int(input("chose the frame of the animation : "))
                return animate(i, flag=False)
            y = y_matrix[i%(len(y_matrix))]
            z = z_matrix[i%(len(y_matrix))] # y_matrix because in case of error, there is one more element in z_matrix and we need y and z to be synchronized
            h = h_matrix[i%(len(y_matrix))] # len(y_matrix) = len(h_matrix) so whatever 
            annotation.set_text(f"Q={hydrogram[i%(len(hydrogram))]:.2f}\n {(i%(len(hydrogram)))*100/(len(hydrogram)):.1f}%\n i={i%len(y_matrix)}\n")
            line.set_data(x, np.array(y)+np.array(z))
            line2.set_data(x, z)
            line3.set_data(x, h)
            return line, line2, line3, annotation, 
        time_ani = 15 # seconds
        # print(f"interval = {time_ani*1000/len(hydrogram)}, frames = {len(hydrogram)}")
        ani = animation.FuncAnimation(fig, animate, frames=len(hydrogram), interval=time_ani*1000/len(hydrogram), repeat=True, repeat_delay=3000)
        if backup:
            print("saving animation...")
            ani.save("./results/animation.mp4", fps=len(hydrogram)/time_ani, dpi=100)
            print("animation saved.")
        if debug:
            print("[DEBUG] STARTING DEBUG")
            plt.show()
            while input("[DEBUG] write \"stop\" to leave the debug loop : ") != "stop":
                index = int(input("[DEBUG]\t CHOOSE THE INDEX : "))%len(profile_list)
                test_profile = profile_list[index]
                test_Q = hydrogram[index]
                test_y = test_profile.compute_depth(test_Q, plot=True)
                plt.show()



        return ani

    # getters and setters

    def get_x_list(self):
        """
        Return the list of curvilinear abscissa.
        """
        x_list = []
        for section in self.__section_list:
            x_list.append(section.get_x())
        return x_list

    def get_x_max(self):
        return self.__section_list[-1].get_x()

    def get_x_min(self):
        return self.__section_list[0].get_x()

    def get_z_list(self):
        """
        Return the list of height above the datum.
        """
        return [section.get_z() for section in self.__section_list]

    def get_yc_list(self, Q):
        return [section.get_yc(Q) for section in self.__section_list]        

    def get_yn_list(self, Q):
        return [section.get_yn(Q) for section in self.__section_list]

    def get_z_min_list(self):
        """
        Return the list of height above the datum.
        """
        z_list = []
        for section in self.__section_list:
            z_list.append(section.get_z_min())
        return z_list

    def set_z_list(self, z_list):
        for i, section in enumerate(self.__section_list):
            section.set_z(z_list[i])

    def get_nb_section(self):
        return len(self.__section_list)

    def get_section(self, i):
        try:
            return self.__section_list[i]
        except IndexError as e:
            raise IndexError(f"index out of range for get_section. {e}")    

    def get_section_list(self):
        return self.__section_list

    def get_upstream_section(self):
        return self.__upstream

    def get_downstream_section(self):
        return self.__downstream

    def get_stored_volume(self):
        """
        Compute the potential amount of solid stored in the profile.
        """
        s = 0
        section_list = self.__section_list
        for section in section_list:
            s += section.get_stored_volume()
        return s

    def get_upstream_boundary_condition(self, Q):
        s = self.get_upstream_section()
        s0 = s.get_S0(up_direction=False)
        yc = s.get_yc(Q)
        if s0 <= 0:
            return yc
        else:
            return min(yc, s.get_yn(Q))

    def get_downstream_boundary_condition(self, Q):
        s = self.get_downstream_section()
        s0 = s.get_S0(up_direction=True)
        yc = s.get_yc(Q)
        if s0 <= 0:
            return yc
        else:
            return max(yc, s.get_yc(Q))

    # plot methods

    def plot(self, y=None, Q=None, compare=None):
        fig, ax1 = plt.subplots()
        x = self.get_x_list()
        z = self.get_z_list()
        z_min = self.get_z_min_list()
        ax1.plot(x, z, 'r', label="z")
        ax1.plot(x, z_min, '--', color="black", label="z_min")
        ax1.set_xlabel('x abscissa')
        ax1.set_ylabel('height from datum (m)')
        if y != None:
            ax1.plot(x, np.array(z)+np.array(y), "b-", label="water depth")
            ax1.fill_between(x, np.array(z), np.array(z)+np.array(y), color="cyan")
        if Q != None:
            yn_list = self.get_yn_list(Q)
            yc_list = self.get_yc_list(Q)
            ax1.plot(x, np.array(z)+np.array(yn_list), color="violet" , linestyle="dashed", label="normal depth")
            ax1.plot(x, np.array(z)+np.array(yc_list), color="yellowgreen", linestyle="dashdot", label="critical depth")
        if y!=None and Q!=None:
            ax1.plot(x, [self.__section_list[i].get_H(Q, y[i]) for i in range(len(self.__section_list))], linestyle=(0, (5, 1)), color="orange", label="energy grade line")
            #ax1.plot(x, [self.__section_list[i].get_H(Q, yc_list[i]) for i in range(len(self.__section_list))], linestyle=(0, (5, 1)), color="purple", label="critic energy grade line")
        if y != None and compare != None:
            x_compare, y_compare = compare
            ax1.plot(x_compare, y_compare, "-x", color="green", label="compared data")
            y_min, y_max = ax1.get_ylim()
            ax1.set_ylim(y_min, y_max+0.2*(y_max-y_min))
            plt.legend(loc=3)
            ax2 = ax1.twinx()
            ax2.set_ylabel("relative error between water elevation", color="tomato")
            y_diff = (np.array(z)+np.array(y) - np.interp(x, x_compare, y_compare)) / (np.array(z)+np.array(y))
            max_diff = max(abs(y_diff))            
            ax2.set_ylim(-6*max_diff, 1.05*max_diff)
            ax2.plot(x, y_diff, color="tomato", label="relative error")
            ax2.plot(x, [0 for _ in range(len(x))], color="darkcyan", label="e=0")
            ax2.tick_params(axis='y', colors='tomato')
            

        plt.legend(loc=1)
        plt.title("Profile")
        return fig

    def plot3D(self, y=None):
        try:
            y0 = float(y)
            y = [y0 for _ in range(len(self.__section_list))]
        except ValueError:
            pass
        except TypeError:
            pass

        fig = plt.figure()
        ax = fig.add_subplot(projection="3d")
        data1 = [[], [], []]
        data2 = [[], [], []]
        alternate = -1
        for i, section in enumerate(self.__section_list):
            points = section.get_points()
            xi = section.get_x()
            zi = section.get_z()
            if y != None:
                wet_points = section.get_wet_section(y[i])
                data2[0].append(xi)
                data2[0].append(xi)
                data2[1].append(wet_points[0][0])
                data2[1].append(wet_points[-1][0])
                data2[2].append(zi + wet_points[0][1])
                data2[2].append(zi + wet_points[-1][1])
            if alternate==1: # trick for make a beautiful plot :)
                points = points[-1::-1]
            for point_x, point_y in points:
                data1[0].append(xi)
                data1[1].append(point_x)
                data1[2].append(zi+point_y)    
            alternate = - alternate
        ax.plot(data1[0], data1[1], data1[2], color="brown", label="profile")
        ax.plot(data2[0], data2[1], data2[2], color="blue", label="water depth")
        ax.set_xlabel('x profile abscissa')
        ax.set_ylabel('x section absicssa')
        ax.set_zlabel('z height (m)')
        plt.legend()
        return fig