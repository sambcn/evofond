import logging as lg
import matplotlib.pyplot as plt
import numpy as np
import json
import os
import datetime as dt

from src.utils import parse_datafile, write_datafile, inter_xy, hydrogrammeLavabre
from src.rectangularSection import RectangularSection
from src.trapezoidalSection import TrapezoidalSection
from src.profile import Profile
from src.granulometry import Granulometry
from src.sedimentTransport.lefortsogreah1991 import LefortSogreah1991
from src.sedimentTransport.lefort2015 import Lefort2015
from src.sedimentTransport.rickenmann1990 import Rickenmann1990
from src.sedimentTransport.rickenmann1991 import Rickenmann1991
from src.sedimentTransport.meunier1989 import Meunier1989
from src.sedimentTransport.meyerpeter1948 import MeyerPeter1948
from src.sedimentTransport.pitonrecking2017 import PitonRecking2017
from src.sedimentTransport.piton2016 import Piton2016

def run_back(args, hydrau=None):
    """
    Interface between front and back : read args, initialize objects, then launch computations 
    """
    # granulometry
    granulometry_files = args["GRANULOMETRY_FILES"]
    granulometry_list = []
    for path in granulometry_files:
        granulo = json.load(open(path, 'r'))
        granulometry_list.append(Granulometry(dm=granulo["dm"], d30=granulo["d30"], d50=granulo["d50"], d90=granulo["d90"], d84tb=granulo["d84tb"], d84bs=granulo["d84bs"], Gr=granulo["Gr"]))
    
    # hydrogram
    if args["LAVABRE"]:
        if args["DURATION"] > 0 and args["DT"] > 0:
            t = np.arange(0, args["DURATION"]+args["DT"], args["DT"])
        else:
            print("ERROR : while building lavabre hydrogram time_list, DURATION and DT must be positive (>0)")
            return
        try:
            Q = hydrogrammeLavabre(args["QM"],args["TM"],args["ALPHA"],args["QB"],t)
        except ValueError:
            print("ERROR : simulation aborted while building lavabre hydrogram")
            return
    else:
        hydrogram_data = parse_datafile(args["HYDROGRAM_PATH"])
        try:
            t = hydrogram_data[0]
            Q = hydrogram_data[1]
        except IndexError:
            print("ERROR : hydrogram data must have exactly 2 columns : t Q")
            return
        except FileNotFoundError:
            print(f"ERROR : {args['HYDROGRAM_PATH']} does not exist")

    # profile
    section = args["SECTION"]
    try:
        data = parse_datafile(args["PROFILE_PATH"])
    except FileNotFoundError:
        print(f"ERROR : {args['PROFILE_PATH']} does not exist")
    try:
        first_line = open(args["PROFILE_PATH"], 'r').readlines()[0].split()
    except IndexError:
        print(f"ERROR : no data in {args['PROFILE_PATH']}")
        return
    try:
        x = data[first_line.index('x')]
        z = data[first_line.index('z')]
        if x[z.index(min(z))] < x[z.index(max(z))]:
            print("INFO : reversed abscissa (computations are done with upstream for lower x)")
            x_mini = min(x)
            x_maxi = max(x)
            x = [x_mini + x_maxi - xi for xi in x]
        b = data[first_line.index('b')]
    except ValueError:
        print("ERROR : the first line of your data file must be column titles with at least x, z, b")
        return
    try:
        z_min = data[first_line.index('zmin')]
    except ValueError:
        print("WARNING : zmin set equal to z because there is no column 'zmin' found")
        z_min = [None for _ in range(len(z))]
    try:
        y_max = data[first_line.index('ymax')]
    except ValueError:
        print("WARNING : ymax set on its default value because there is no column 'ymax' found")
        y_max = [None for _ in range(len(z))]
    try:
        granulo_index = data[first_line.index('granulometry')]
        granulo_index = [int(_)-1 for _ in granulo_index]
    except ValueError:
        print("WARNING : granulometry set on the first granulometry by default because there is no column 'granulometry' found")
        granulo_index = [0 for _ in range(len(z))]
    try:
        manning = data[first_line.index('manning')]
    except ValueError:
        print("WARNING : manning coefficient computed with granulometry because there is no column 'manning' found")
        manning = [None for _ in range(len(z))]        
    if section=="trapezoidal":
        try:
            s = data[first_line.index('s')]
        except ValueError:
            print("ERROR : the first line of your data file must be column titles. Trapezoidal section need one column called 's' for slope of the sides")
            return
    elif section=="rectangular":
        list_of_section = []
        for i in range(len(x)):
            list_of_section.append(RectangularSection(x[i], z[i], b[i], z_min=z_min[i], y_max=y_max[i], granulometry=granulometry_list[granulo_index[i]], manning=manning[i]))
        if len(list_of_section)<2:
            print(f"ERROR : you need at least 2 sections to build a profile (currently {len(list_of_section)})")
            return
        profile = Profile(list_of_section, name=args["NAME"])
        if args["INTERPOLATION"]:
            if args["DX"] == None:
                print("ERROR : interpolation set on true but no dx was given (dx=null)")
            else:
                profile.complete(args["DX"])
        # profile.export("./profile_export.pkl")
    else:
        print(f"unknown section : {section}. Execution aborted.")

    # sedimentogram
    I = float(args["UPSTREAM_SLOPE"])*0.01 # slope is given in %
    transport_law_value = args["TRANSPORT_LAW"]
    if transport_law_value == "Lefort2015":
        transport_law = Lefort2015()
    elif transport_law_value == "LefortSogreah1991":
        transport_law = LefortSogreah1991()
    elif transport_law_value == "Meunier1989":
        transport_law = Meunier1989()
    elif transport_law_value == "Rickenmann1991":
        transport_law = Rickenmann1991()
    elif transport_law_value == "Rickenmann1990":
        transport_law = Rickenmann1990()
    elif transport_law_value == "MeyerPeter1948":
        transport_law = MeyerPeter1948()
    elif transport_law_value == "PitonRecking2017":
        transport_law = PitonRecking2017()
    elif transport_law_value == "Piton2016":
        transport_law = Piton2016()
    else:
        print(f"ERROR : unknown sediment transport law (= {transport_law_value})")
        return
    i_upstream = x.index(min(x))
    QsIn = [transport_law.compute_Qs_formula(args["UPSTREAM_WIDTH"], granulometry_list[granulo_index[i_upstream]], Q[i], I) for i, ti in enumerate(t)]

    # result
    if hydrau == None:
        result = profile.compute_event(Q, t, transport_law, sedimentogram=QsIn, backup=False, debug=False, method="ImprovedEuler", friction_law=args["FRICTION_LAW"], cfl=args["CFL"], critical=args["CRITICAL"])
        try:
            os.chdir("./results")
        except FileNotFoundError:
            os.mkdir("./results")
            os.chdir("./results")
            print("./results folder created.")
        e = dt.datetime.now()
        try:
            os.mkdir(f"./{e.year}_{e.month}_{e.day}_{e.hour}_{e.minute}_{e.second}")
        except FileExistsError:
            pass
        os.chdir(f"./{e.year}_{e.month}_{e.day}_{e.hour}_{e.minute}_{e.second}")
        
        json.dump(args, open("./conf_reminder.json", 'w'), indent=6)
        step = args["BACKUP_TIME_STEP"]
        ani = result["animation"]
        x = result["abscissa"]
        y_matrix = result["water_depth"]
        z_matrix = result["bottom_height"]
        t_list = result["time"]
        h_matrix = result["energy"]
        try:
            os.mkdir("./np_files")
        except FileExistsError:
            pass
        os.chdir("./np_files")
        np.save("./x_list", np.array(x))
        np.save("./t_list", np.array(t_list))
        np.save("./y_matrix", np.array(y_matrix))
        np.save("./z_matrix", np.array(z_matrix))
        np.save("./h_matrix", np.array(h_matrix))
        os.chdir("..")
        try:
            os.mkdir("./txt_files")
        except FileExistsError:
            pass
        os.chdir("./txt_files")
        next_t = 0
        column_name_list = ["x", "y", "z", "h"]
        for i, ti in enumerate(t_list):
            if ti >= next_t:
                next_t += step
                filename = f"{ti:.0f}.txt"
                # save
                data = [x, y_matrix[i], z_matrix[i], h_matrix[i]]
                write_datafile(filename, column_name_list, data)
        
    else:
        y_list = profile.compute_depth(hydrau, plot=True, method="ImprovedEuler", friction_law=args["FRICTION_LAW"])
    plt.show()

    return