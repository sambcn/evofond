import numpy as np
import copy
import scipy.optimize as op
import pickle as pkl
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from time import time
from irregularSection import IrregularSection
from perf import Performance
from utils import Y_MIN, G, get_matrix_max, read_hecras_data, reverse_data, time_to_string

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

    # def smooth(self):
    #     """
    #     used to erase growing instabilities.
    #     When it detects a pattern such like "z increase, z decrease, z increase" then it keeps every following section which keep this alternation, and replace their z by the "mean" value.
    #     Look at documentation for more details.
    #     """
    #     i = 0
    #     while i < self.get_nb_section()-3:
    #         j = i
    #         while j < self.get_nb_section()-3:
    #             s0_1 = self.get_section(j).get_S0()
    #             s0_2 = self.get_section(j+1).get_S0()
    #             s0_3 = self.get_section(j+2).get_S0()
    #             if (s0_1-s0_2)*(s0_2-s0_3) < 0:
    #                 j += 1
    #             else:
    #                 break
    #         if j > i:
    #             x0 = self.get_section(i).get_x()
    #             z0 = self.get_section(i).get_z()
    #             xf = self.get_section(j+1).get_x()
    #             zf = self.get_section(j+1).get_z()
    #             for k in range(i+1, j+1):
    #                 s = self.get_section(k)
    #                 try:
    #                     s.set_z(np.interp(s.get_x(), [x0, xf], [z0, zf]))
    #                 except ValueError as e:
    #                     s.set_z(s.get_z_min())
    #         i = j+1
    #     return

    def copy(self):
        """return a safe copy of this section"""
        section_list = []
        for s in self.__section_list:
            section_list.append(s.copy())
        copied_profile = Profile(section_list)
        return copied_profile

    def export(self, filename="exported_profile.pkl"):
        """export this profile using pickle"""
        pkl.dump(self, open("./results/"+filename, "wb"))
        print("profile exported successfully.")
        
    # resolution methods

    @Performance.measure_perf
    def compute_depth(self, Q, plot=False, hydraulic_jump_analysis=False, compare=None, method="ImprovedEuler", friction_law="Ferguson"): # Trying to accelerate computations
        """
        Compute the water depth at every section of the profile, for an upstream water discharge Q. Use standard step method.
        """
        k = 0.1 # Borda coef 
        method_set = {"Euler", "ImprovedEuler", "RungeKutta"}
        if not(method in method_set):
            print(f"WARNING : chosen method not in the available list : {method_set}, it has been set by default on ImprovedEuler")
            method = "ImprovedEuler"
        yc_list = self.get_yc_list(Q)
        n = len(self.__section_list)

        # supercritical computation
        y_list_supercritical = [self.get_upstream_boundary_condition(Q)]
        hs_list_supercritical = [self.get_upstream_section().get_Hs(Q, y_list_supercritical[0])]
        i_current = 0
        while i_current < n-1:
            current_section = self.get_section(i_current)
            next_section = self.get_section(i_current+1)

            y_supercritical = self.__compute_next_y(Q, current_section, next_section, y_list_supercritical[-1], hs_list_supercritical[-1], yc_list[i_current+1], supercritical=True, method=method, friction_law=friction_law)
            hs_next = next_section.get_Hs(Q, y_supercritical)

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
        y_list_subcritical = [self.get_downstream_boundary_condition(Q, friction_law=friction_law)]
        hs_list_subcritical = [self.get_downstream_section().get_Hs(Q, y_list_subcritical[-1])]
        i_current = n-1
        while i_current > 0:
            current_section = self.get_section(i_current)
            next_section = self.get_section(i_current-1)

            y_subcritical = self.__compute_next_y(Q, current_section, next_section, y_list_subcritical[0], hs_list_subcritical[0], yc_list[i_current-1], supercritical=False, method=method, friction_law=friction_law)
            hs_next = next_section.get_Hs(Q, y_subcritical)

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

        flow_supercritical = None
        for i in range(self.get_nb_section()):
            if hs_list_subcritical[i] >=  hs_list_supercritical[i]:
                j = i-1
                y_list.append(y_list_subcritical[i])
                while j>0 and (self.get_section(j).get_Fs(Q, y_list_subcritical[j]) >= self.get_section(j).get_Fs(Q, y_list_supercritical[j])):
                    y_list[j] = (y_list_subcritical[j])
                    j -= 1 
                if flow_supercritical != None and flow_supercritical:
                    hydraulic_index.append(j)
                flow_supercritical = False
            else:
                y_list.append(y_list_supercritical[i])
                flow_supercritical = True
        
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

        # flow_supercritical = None   
        # for i, s in enumerate(self.__section_list):
        #     Fs_supercritical = s.get_Fs(Q, y_list_supercritical[i])
        #     Fs_subcritical = s.get_Fs(Q, y_list_subcritical[i])
        #     if flow_supercritical == None:
        #         flow_supercritical = Fs_supercritical >= Fs_subcritical
        #         y_list.append(y_list_supercritical[i] if flow_supercritical else y_list_subcritical[i])
        #         continue
        #     changing_flow = (Fs_supercritical >= Fs_subcritical)^(flow_supercritical) # XOR operator. True if flow changes
        #     y_current = y_list_supercritical if flow_supercritical else y_list_subcritical 
        #     if changing_flow:
        #         y_new = y_list_subcritical if flow_supercritical else y_list_supercritical
        #         if s.get_H(Q, y_new[i]) <= s.get_up_section().get_H(Q, y_list[-1]):
        #             y_list.append(y_new[i])
        #             if flow_supercritical:
        #                 hydraulic_index.append(i-1)
        #             flow_supercritical = not(flow_supercritical)
        #         else:
        #             y_list.append(y_current[i])
        #     else:
        #         y_list.append(y_current[i])

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
            fig = self.plot(y=y_list_supercritical, Q=Q, friction_law=friction_law)
            fig = self.plot(y=y_list_subcritical, Q=Q, friction_law=friction_law)
            fig = self.plot(y=y_list, Q=Q, compare=compare, friction_law=friction_law)
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

    def find_best_dt(self, Q, y_list):
        dt_list = []
        for i, s in enumerate(self.get_section_list()[:-1]):
            dx = s.get_down_section().get_x() - s.get_x()
            v = 0.5*(s.get_V(Q, y_list[i]) + s.get_down_section().get_V(Q, y_list[i+1]))
            dt_list.append(dx/v)
        # print(dt_list)
        return min(dt_list)

    @Performance.measure_perf
    def update_bottom(self, Q, y, QsIn0, dt, law, plot=False, friction_law="Ferguson"):
        """
        Q and dt are respectively the water discharge and the time step. 
        Compute the sediment discharge at the downstream and change the attribute __z for every section
        """
        z0 = self.get_z_list()
        V0 = self.get_stored_volume()
        QsIn = QsIn0
        if plot:
            fig = self.plot(Q=Q, y=y, friction_law=friction_law)
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

    def compute_event(self, hydrogram, dt, law, backup=False, debug=False, smooth=False, method="ImprovedEuler", friction_law="Ferguson"):
        start_computation = time()
        y_matrix = [] # list of the water depth during the event
        z_matrix = [self.get_z_list()] # list of the bottom height during the event
        h_matrix = [] # list of head during the event
        V_in = 0 # solid volume gone into the profile during the event
        V_out = 0 # solid volume gone out of the profile
        stored_volume_start = self.get_stored_volume() # stored volume of sediment at the start of the event 
        initial_profile = self.copy()
        method_set = {"Euler", "ImprovedEuler", "RungeKutta"}
        if not(method in method_set):
            print(f"WARNING : chosen method not in the available list : {method_set}, it has been set by default on ImprovedEuler")
            method = "ImprovedEuler"
        log_string = f"[backup={backup}, debug={debug}, smooth={smooth}, method={method}]"

        if debug:
            profile_list = [initial_profile]
            current_stored_volume = stored_volume_start
            total_volume_difference = [] 
            one_step_volume_difference = []

        for i, Q in enumerate(hydrogram):
            if i%(1 if debug else 10)==0 or i==len(hydrogram)-1:
                print(f"{i}/{len(hydrogram)} "+log_string)

            # hydraulic computations
            try:
                y_list = self.compute_depth_bis(Q, method=method, friction_law=friction_law)
            except Exception as e:
                print("ERROR IN COMPUTING THIS EVENT : COULD NOT FINISH\n plotting the last state ...\n")
                try:
                    self.plot(Q=Q, y=y_list, friction_law=friction_law)
                except Exception:
                    pass
                plt.show()
                break

            y_matrix.append(y_list)
            h_matrix.append([s.get_H(Q, y_list[i]) for i, s in enumerate(self.__section_list)])

            # solid transport
            QsIn0 = law.compute_Qs(initial_profile.get_upstream_section(), Q, y_list[0], y_list[1]) # Gonna change, it is a given parameter, chosen by users
            V_in += QsIn0*dt
            QsOut = self.update_bottom(Q, y_list, QsIn0, dt, law, friction_law=friction_law)
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
            y_matrix.append(self.compute_depth_bis(hydrogram[-1]))
            h_matrix.append([s.get_H(Q, y_matrix[-1][i]) for i, s in enumerate(self.__section_list)])
        except ValueError:
            pass
        stored_volume_end = self.get_stored_volume()
        end_computation = time()
        print(f"computation time = {end_computation-start_computation}s")
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
        axs[0].annotate(f"event of {time_to_string(len(hydrogram)*dt)}\ndt = {dt}s\nQmax = {max(hydrogram):.3f}m3/s\nQmean = {np.mean(hydrogram):.3f}m3/s\nsediment transport law = {str(law)}\ntime of computation = {end_computation-start_computation:.4f}s", (self.get_x_min(), min(self.get_z_min_list())))
        axs[1].plot(x, np.array(z_matrix[-1])-np.array(z_matrix[0]), "b", label="newz-z")
        axs[1].set(xlabel="x", ylabel="difference of height (m)")
        self.__plot_width_background(axs[0], label=False)
        self.__plot_width_background(axs[1])
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
        ax1.set_ylabel("height (m)")
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
        self.__plot_width_background(ax1)
        plt.legend()
        def animate(i): 
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
            answer = str(input("[DEBUG] \t Do you want to see the animation ? [yes/no] :"))
            while answer != "yes" and answer != "no":
                answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
            if answer == "yes":
                plt.show()
            else:
                plt.close()
            while input("[DEBUG] write \"stop\" to leave the debug loop, else please press ENTER : ") != "stop":
                index = (int(input(f"[DEBUG]\t CHOOSE THE INDEX (<{len(profile_list)}) : ")))%(len(profile_list))
                test_profile = profile_list[index]
                if index == len(profile_list)-1:
                    print("[DEBUG] \t\t\t last profile chosen.")
                    test_profile.plot()
                    plt.show()
                else:
                    test_Q = hydrogram[index]
                    test_y = test_profile.compute_depth_bis(test_Q, plot=True)
                    plt.show()
                    answer = str(input("[DEBUG] \t Do you want to complete the profile ? [yes/no] :"))
                    while answer != "yes" and answer != "no":
                        answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
                    if answer=="yes":
                        dx = int(input("[DEBUG] \t\t chose a dx : "))
                        test_profile.complete(dx)
                        test_y = test_profile.compute_depth_bis(test_Q, plot=True)
                        plt.show()
                answer = str(input("[DEBUG] \t Do you want to save this profile ? [yes/no] :"))
                while answer != "yes" and answer != "no":
                    answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
                if answer=="yes":
                    filename = str(input("[DEBUG]\t please chose a filename : "))
                    test_profile.export(filename)
                    print("[DEBUG] profile saved.")


        return ani

    def compute_event_bis(self, hydrogram, t_hydrogram, law, backup=False, debug=False, method="ImprovedEuler", friction_law="Ferguson"):
        start_computation = time()
        y_matrix = [] # list of the water depth during the event
        z_matrix = [self.get_z_list()] # list of the bottom height during the event
        h_matrix = [] # list of head during the event
        V_in = 0 # solid volume gone into the profile during the event
        V_out = 0 # solid volume gone out of the profile
        Q_list = [] # list of water discharge at each step of the computation
        t_list = [] # list of the t time of each step of the commputation
        dt_list = [] # list of the time steps used at each iteration
        t = t_hydrogram[0]
        stored_volume_start = self.get_stored_volume() # stored volume of sediment at the start of the event 
        initial_profile = self.copy()
        method_set = {"Euler", "ImprovedEuler", "RungeKutta"}
        if not(method in method_set):
            print(f"WARNING : chosen method not in the available list : {method_set}, it has been set by default on ImprovedEuler")
            method = "ImprovedEuler"
        log_string = f"[backup={backup}, debug={debug}, method={method}, friction_law={friction_law}]"

        if debug:
            profile_list = [initial_profile]
            current_stored_volume = stored_volume_start
            total_volume_difference = [] 
            one_step_volume_difference = []

        while t <= t_hydrogram[-1]:
            Q = np.interp(t, t_hydrogram, hydrogram)
            t_list.append(t)
            Q_list.append(Q)

            # hydraulic computations
            try:
                y_list = self.compute_depth_bis(Q, method=method, friction_law=friction_law)
            except Exception as e:
                print("ERROR IN COMPUTING THIS EVENT : COULD NOT FINISH\n plotting the last state ...\n")
                try:
                    self.plot(Q=Q, y=y_list, friction_law=friction_law)
                except Exception:
                    pass
                plt.show()
                break

            y_matrix.append(y_list)
            h_matrix.append([s.get_H(Q, y_list[i]) for i, s in enumerate(self.__section_list)])
            dt_opti = self.find_best_dt(Q, y_list)
            t_aux = list(np.sort(abs(np.array(t_hydrogram) - t)))
            dt_hydrogram = t_aux[1] + t_aux[0]
            print(f"dt_opti={dt_opti:.3f}, dt_hydrogram={dt_hydrogram:.3f}")
            dt = min(dt_hydrogram, dt_opti)
            print(f"{t:.3f}/{t_hydrogram[-1]} (dt={dt:.3f}) "+log_string)
            t += dt
            dt_list.append(dt)

            # solid transport
            QsIn0 = law.compute_Qs(initial_profile.get_upstream_section(), Q, y_list[0], y_list[1]) # Gonna change, it is a given parameter, chosen by users
            V_in += QsIn0*dt
            QsOut = self.update_bottom(Q, y_list, QsIn0, dt, law, friction_law=friction_law)
            V_out += QsOut*dt
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
            Q = np.interp(t, t_hydrogram, hydrogram)
            t_list.append(t)
            Q_list.append(Q)
            y_matrix.append(self.compute_depth_bis(hydrogram[-1]))
            h_matrix.append([s.get_H(Q, y_matrix[-1][i]) for i, s in enumerate(self.get_section_list())])
        except Exception:
            pass
        stored_volume_end = self.get_stored_volume()
        end_computation = time()
        print(f"computation time = {end_computation-start_computation}s")
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
        axs[0].annotate(f"event of {time_to_string(t)}\ndt $\in$ [{min(dt_list):.3f}s, {max(dt_list):.3f}s]\nQmax = {max(hydrogram):.3f}m3/s\nQmean = {np.mean(hydrogram):.3f}m3/s\nfriction law = {friction_law}\nsediment transport law = {str(law)}\ntime of computation = {end_computation-start_computation:.4f}s", (self.get_x_min(), min(self.get_z_min_list())))
        axs[1].plot(x, np.array(z_matrix[-1])-np.array(z_matrix[0]), "b", label="newz-z")
        axs[1].set(xlabel="x", ylabel="difference of height (m)")
        self.__plot_width_background(axs[0], label=False)
        self.__plot_width_background(axs[1])
        fig0.legend()
        fig0.set_size_inches(10.5, 9.5)
        if backup:
            print("saving result plot...")
            fig0.savefig("./results/result.png", dpi=400, format="png")
            print("plot saved.")

        if debug:
            plt.figure()
            plt.plot([i*dt for i in range(len(total_volume_difference))], total_volume_difference, label="acumulated solid creation or disappearance")
            plt.plot([i*dt for i in range(len(one_step_volume_difference))], one_step_volume_difference, label="solid creation or disappearance during this step of time")
            plt.legend()
            plt.title('sediment creation or diseppearance due to numerical errors')

        fig, ax1 = plt.subplots() 
        line, = ax1.plot(x, np.array(y_matrix[0])+np.array(z_matrix[0]), label="water line")
        ax1.set_ylabel("height (m)")
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
        self.__plot_width_background(ax1)
        plt.legend()
        def animate(i): 
            y = y_matrix[i%(len(y_matrix))]
            z = z_matrix[i%(len(y_matrix))] # y_matrix because in case of error, there is one more element in z_matrix and we need y and z to be synchronized
            h = h_matrix[i%(len(y_matrix))] # len(y_matrix) = len(h_matrix) so whatever 
            annotation.set_text(f"Q={Q_list[i%(len(y_matrix))]:.2f}\n {(i%(len(y_matrix)))*100/(len(y_matrix)):.1f}%\n t={t_list[i%len(y_matrix)]:.3f}/{t_list[-1]:.3f}\n")
            line.set_data(x, np.array(y)+np.array(z))
            line2.set_data(x, z)
            line3.set_data(x, h)
            return line, line2, line3, annotation, 
        time_ani = 60 # seconds
        nb_frames = 1000 
        frames=(range(len(y_matrix)) if len(y_matrix)<nb_frames else range(0, len(y_matrix), len(y_matrix)//nb_frames))
        ani = animation.FuncAnimation(fig, animate, frames=frames, interval=time_ani*1000/len(y_matrix), repeat=True, repeat_delay=3000)
        if backup:
            print("saving animation...")
            print(f"number of frames : {len(frames)}")
            t0 = time()
            ani.save("./results/animation.mp4", fps=len(frames)/time_ani, dpi=150)
            print(f"animation saved ({time()-t0}).")
        if debug:
            print("[DEBUG] STARTING DEBUG")
            answer = str(input("[DEBUG] \t Do you want to see the animation ? [yes/no] :"))
            while answer != "yes" and answer != "no":
                answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
            if answer == "yes":
                plt.show()
            else:
                plt.close()
            while input("[DEBUG] write \"stop\" to leave the debug loop, else please press ENTER : ") != "stop":
                index = (int(input(f"[DEBUG]\t CHOOSE THE INDEX (<{len(profile_list)}) : ")))%(len(profile_list))
                test_profile = profile_list[index]
                if index == len(profile_list)-1:
                    print("[DEBUG] \t\t\t last profile chosen.")
                    test_profile.plot()
                    plt.show()
                else:
                    test_Q = hydrogram[index]
                    test_y = test_profile.compute_depth_bis(test_Q, plot=True)
                    plt.show()
                    answer = str(input("[DEBUG] \t Do you want to complete the profile ? [yes/no] :"))
                    while answer != "yes" and answer != "no":
                        answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
                    if answer=="yes":
                        dx = int(input("[DEBUG] \t\t chose a dx : "))
                        test_profile.complete(dx)
                        test_y = test_profile.compute_depth_bis(test_Q, plot=True)
                        plt.show()
                answer = str(input("[DEBUG] \t Do you want to save this profile ? [yes/no] :"))
                while answer != "yes" and answer != "no":
                    answer = str(input("[DEBUG]\t please write 'yes' or 'no' : "))
                if answer=="yes":
                    filename = str(input("[DEBUG]\t please chose a filename : "))
                    test_profile.export(filename)
                    print("[DEBUG] profile saved.")


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

    def get_yn_list(self, Q, friction_law="Ferguson"):
        return [section.get_yn(Q, friction_law=friction_law) for section in self.__section_list]

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

    def get_upstream_boundary_condition(self, Q, friction_law="Ferguson"):
        s = self.get_upstream_section()
        s0 = s.get_S0(up_direction=False)
        yc = s.get_yc(Q)
        if s0 <= 0:
            return yc
        else:
            return s.get_yn(Q, friction_law=friction_law)

    def get_downstream_boundary_condition(self, Q, friction_law="Ferguson"):
        s = self.get_downstream_section()
        s0 = s.get_S0(up_direction=True)
        yc = s.get_yc(Q)
        if s0 <= 0:
            return yc
        else:
            return max(yc, s.get_yc(Q))
    
    def has_only_rectangular_section(self):
        for s in self.get_section_list():
            if not(s.is_rectangular()):
                return False
        return True

    # computational stuff

    @Performance.measure_perf
    def __compute_next_y(self, Q, current_section, next_section, current_y, current_hs, next_yc, supercritical, method="ImprovedEuler", friction_law="Ferguson"):
        """
        private methods only used in compute_depth for computing at each step the new specific head thanks to Euler methods.
        method can be "Euler", "ImprovedEuler" or "RungeKutta". 
        """
        if method=="RungeKutta":
            up_direction = next_section.get_x() - current_section.get_x() < 0
            inter_section = next_section.interp_as_up_section(current_section) if up_direction else next_section.interp_as_down_section(current_section)
            inter_yc = inter_section.get_yc(Q)
            s1 = current_section.get_S0(up_direction=up_direction) - current_section.get_Sf(Q, current_y, friction_law=friction_law)
            dx = inter_section.get_x() - current_section.get_x()
            hs_next = current_hs + dx*s1
            if (hs_next + inter_section.get_z() - (current_hs + current_section.get_z()))*dx > 0:
                hs_next = current_hs + current_section.get_z() - inter_section.get_z()
            if hs_next < inter_section.get_Hs(Q, inter_yc):
                hs_next = inter_section.get_Hs(Q, inter_yc)
            y_inter = inter_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=inter_yc)
            s2 = inter_section.get_S0() - inter_section.get_Sf(Q, y_inter, friction_law=friction_law)
            hs_next = current_hs + dx*s2
            if (hs_next + inter_section.get_z() - (current_hs + current_section.get_z()))*dx > 0:
                hs_next = current_hs + current_section.get_z() - inter_section.get_z()
            if hs_next < inter_section.get_Hs(Q, inter_yc):
                hs_next = inter_section.get_Hs(Q, inter_yc)
            y_inter = inter_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=inter_yc)
            s3 = inter_section.get_S0() - inter_section.get_Sf(Q, y_inter, friction_law=friction_law)
            hs_next = current_hs + 2*dx*s3
            if (hs_next + next_section.get_z() - (current_hs + current_section.get_z()))*dx > 0:
                hs_next = current_hs + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, next_yc):
                hs_next = next_section.get_Hs(Q, next_yc)
            y_inter = next_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=inter_yc)
            s4 = next_section.get_S0(up_direction=not(up_direction)) - next_section.get_Sf(Q, y_inter, friction_law=friction_law)
            hs_next = current_hs + (s1+2*s2+2*s3+s4)*dx*(2/6)
            if (hs_next + next_section.get_z() - (current_hs + current_section.get_z()))*dx > 0:
                hs_next = current_hs + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, next_yc):
                hs_next = next_section.get_Hs(Q, next_yc)
            return next_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=next_yc)
        else:
            up_direction = next_section.get_x() - current_section.get_x() < 0
            s1 = current_section.get_S0(up_direction=up_direction) - current_section.get_Sf(Q, current_y, friction_law=friction_law)
            dx = next_section.get_x() - current_section.get_x()
            hs_next = current_hs + dx*s1
            if (hs_next + next_section.get_z() - (current_hs + current_section.get_z()))*dx > 0:
                hs_next = current_hs + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, next_yc):
                hs_next = next_section.get_Hs(Q, next_yc)
            next_y = next_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=next_yc)
            if method=="Euler":
                return next_y
            s2 = next_section.get_S0(up_direction=not(up_direction)) - next_section.get_Sf(Q, next_y, friction_law=friction_law)
            hs_next = current_hs+ 0.5*(s1+s2)*dx
            if (hs_next + next_section.get_z() - (current_hs + current_section.get_z()))*dx>0:
                hs_next = current_hs + current_section.get_z() - next_section.get_z()
            if hs_next < next_section.get_Hs(Q, next_yc):
                hs_next = next_section.get_Hs(Q, next_yc)
            return next_section.get_y_from_Hs(Q, hs_next, supercritical=supercritical, yc=next_yc)

    # plot methods

    def plot(self, y=None, Q=None, title="profile", compare=None, friction_law="Ferguson"):
        fig, ax1 = plt.subplots()
        x = self.get_x_list()
        z = self.get_z_list()
        z_min = self.get_z_min_list()
        ax1.plot(x, z, 'r', label="z")
        ax1.plot(x, z_min, '--', marker='+', color="black", label="z_min")
        ax1.set_xlabel('x abscissa')
        ax1.set_ylabel('height from datum (m)')
        if y != None:
            ax1.plot(x, np.array(z)+np.array(y), "b-", label="water depth")
            ax1.fill_between(x, np.array(z), np.array(z)+np.array(y), color="cyan")
        if Q != None:
            yn_list = self.get_yn_list(Q, friction_law=friction_law)
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
        self.__plot_width_background(ax1)
                        
        plt.legend(loc=1)
        plt.title(title)
        return fig

    def __plot_width_background(self, ax, label=True):
        """add a grey transparent background which depends on section width. It is relevant when ax has an x_axis which is the profile abscissa"""
        if self.has_only_rectangular_section():
            b_list = [s.get_b() for s in self.get_section_list()]
            b_max = max(b_list)
            b_min = min(b_list)
            b_diff = b_max-b_min if b_max != b_min else 1
            alpha_max = 0.8
            alpha_min = 0.2
            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            b_current = b_list[0]
            x_current = self.get_upstream_section().get_x()
            x_final = self.get_downstream_section().get_x()
            b_seen = []
            for i in range(1, self.get_nb_section()-1):
                if b_list[i] != b_current:
                    new_x = (self.get_section(i).get_x() + self.get_section(i-1).get_x())*0.5
                    if not(label) or b_current in b_seen:
                        ax.fill_betweenx([ymin, ymax], x_current, new_x, color='grey', alpha=alpha_min+(alpha_max-alpha_min)*((b_max-b_current)/(b_diff)))
                    else:
                        ax.fill_betweenx([ymin, ymax], x_current, new_x, color='grey', alpha=alpha_min+(alpha_max-alpha_min)*((b_max-b_current)/(b_diff)), label=f"width = {b_current}m")
                    x_current = new_x
                    b_seen.append(b_current)
                    b_current = b_list[i]
            if not(label) or b_current in b_seen:
                ax.fill_betweenx([ymin, ymax], x_current, x_final, color='grey', alpha=alpha_min+(alpha_max-alpha_min)*((b_max-b_current)/(b_diff)))
            else:    
                ax.fill_betweenx([ymin, ymax], x_current, x_final, color='grey', alpha=alpha_min+(alpha_max-alpha_min)*((b_max-b_current)/(b_diff)), label=f"width = {b_current}m")
        return

    def plot3D(self, y=None, title="profile"):
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
        ax.plot(data1[0], data1[1], data1[2], color="brown", label="profile", alpha=0.5)
        ax.plot(data2[0], data2[1], data2[2], color="blue", label="water depth")
        ax.set_xlabel('x profile abscissa')
        ax.set_ylabel('x section absicssa')
        ax.set_zlabel('z height (m)')
        plt.title(title)
        plt.legend()
        return fig

    # static methods

    def import_profile(filename):
        return pkl.load(open(filename, "rb"))

    # TEST

    @Performance.measure_perf
    def compute_depth_bis(self, Q, plot=False, hydraulic_jump_analysis=False, method="ImprovedEuler", friction_law="Ferguson", compare=None):
        
        method_set = {"Euler", "ImprovedEuler", "RungeKutta"}
        if not(method in method_set):
            print(f"WARNING : chosen method not in the available list : {method_set}, it has been set by default on ImprovedEuler")
            method = "ImprovedEuler"

        yc_list = self.get_yc_list(Q)
        y_list = yc_list[:]
        y_list[0] = self.get_upstream_boundary_condition(Q, friction_law=friction_law)
        y_list[-1] = self.get_downstream_boundary_condition(Q, friction_law=friction_law)
        hs_list = [s.get_Hs(Q, y_list[i]) for i, s in enumerate(self.get_section_list())]
        hsc_list = [s.get_Hs(Q, yc_list[i]) for i, s in enumerate(self.get_section_list())]
        Fs_list = [s.get_Fs(Q, y_list[i]) for i, s in enumerate(self.get_section_list())]
        hydraulic_index = []

        i_current = 0
        i_memory_1 = self.get_nb_section()-1
        i_memory_2 = self.get_nb_section()-1
        i_memory_3 = self.get_nb_section()-1        
        down_direction = True

        y_list_memory = y_list[:]
        while i_current > 0 or (down_direction and i_current == 0):
            # print(f"start at x = {self.get_section(i_current).get_x()} toward {'down' if down_direction else 'up'} direction")
            
            if down_direction:
                while i_current < self.get_nb_section()-1:#i_memory_3:
                    current_section = self.get_section(i_current)
                    next_section = self.get_section(i_current+1)
                    if y_list[i_current] > yc_list[i_current]:                    
                        y_next = self.__compute_next_y(Q, current_section, next_section, yc_list[i_current], hsc_list[i_current], yc_list[i_current+1], supercritical=True, method=method, friction_law=friction_law)
                    else:
                        y_next = self.__compute_next_y(Q, current_section, next_section, y_list[i_current], hs_list[i_current], yc_list[i_current+1], supercritical=True, method=method, friction_law=friction_law)
                    hs_next = next_section.get_Hs(Q, y_next)
                    Fs_next = next_section.get_Fs(Q, y_next)
                    if i_current != self.get_nb_section()-2 and (Fs_next > Fs_list[i_current+1] or hs_list[i_current+1] + next_section.get_z() > hs_list[i_current] + current_section.get_z()):
                        y_list[i_current+1] = y_next
                        hs_list[i_current+1] = hs_next
                        Fs_list[i_current+1] = Fs_next 
                    i_current += 1
                down_direction = False
                i_memory_3 = i_memory_2
            else:
                # i_current = i_memory_2
                update_flag = False
                while i_current > 0:
                    current_section = self.get_section(i_current)
                    next_section = self.get_section(i_current-1)
                    if y_list[i_current] < yc_list[i_current]:                    
                        y_next = self.__compute_next_y(Q, current_section, next_section, yc_list[i_current], hsc_list[i_current], yc_list[i_current-1], supercritical=False, method=method, friction_law=friction_law)
                    else:
                        y_next = self.__compute_next_y(Q, current_section, next_section, y_list[i_current], hs_list[i_current], yc_list[i_current-1], supercritical=False, method=method, friction_law=friction_law)
                    hs_next = next_section.get_Hs(Q, y_next)
                    Fs_next = next_section.get_Fs(Q, y_next)
                    if i_current != 1 and Fs_next > Fs_list[i_current-1]:
                        # if hs_next + next_section.get_z() >= hs_list[i_current] + current_section.get_z():
                        if not(update_flag):
                            i_memory_1 = i_current-1
                            update_flag = True
                        y_list[i_current-1] = y_next
                        hs_list[i_current-1] = hs_next
                        Fs_list[i_current-1] = Fs_next
                    elif update_flag: # end of an updated section : we need to refresh supercritical computation
                        i_memory_2 = i_current-1
                        down_direction = True
                        if not(i_current-1 in hydraulic_index):
                            hydraulic_index.insert(0, i_current-1)
                        i_current = i_memory_1
                        break
                    i_current -= 1
            
            # print(f"stoped at x={self.get_section(i_current).get_x()}")
            # self.plot(Q=Q, y=y_list, title="new profile")
            # self.plot(Q=Q, y=y_list_memory, title="previous profile")
            # y_list_memory = y_list[:]
            # plt.show()

        # y_list = yc_list

        if plot:
            fig = self.plot(y=y_list, Q=Q, friction_law=friction_law, compare=compare)
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