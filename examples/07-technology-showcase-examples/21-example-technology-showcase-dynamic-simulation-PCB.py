""".. _ref_dynamic_simulation_printed_circuit_board:

Modal analysis of a PCB  modelled as shell,solid and solid-shell
------------------------------------------------------

This examples shows how to use PyMAPDL to import an existing FE model and to
run a modal analysis. PyDPF modules are also used for post-processing.

The existing FE model (CDB file) contains imported traces too.

This example uses 3 modelling approaches :shell(SHELL181),solid(SOLID185) and solid-shell(SLOSH190) to model PCB board 

Additional Packages Used
~~~~~~~~~~~~~~~~~~~~~~~~
* `Matplotlib <https://matplotlib.org>`_ is used for plotting purposes.

"""
###############################################################################
# Starting MAPDL as a service and importing an external model
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The original FE model is given in the Ansys Mechanical APDL Technology
# Showcase Manual.  The .cdb contains a FE model of a single circuit board. The
# model is meshed with SOLID186, SHELL181 and BEAM188 elements. All components
# of the PCB model is assigned with linear elastic isotropic materials. Bonded
# and flexible surface-to-surface contact pairs are used to define the contact
# between the IC packages and the circuit board.

from ansys.mapdl.core import launch_mapdl
from ansys.mapdl.core.examples import download_tech_demo_data

# read model of single circuit board
# download the cdb file
pcb_shell_file = download_tech_demo_data("td-21", "pcb_shell_file.cdb")
pcb_solid_file = download_tech_demo_data("td-21", "pcb_solid_file.cdb")
pcb_solidshell_file = download_tech_demo_data("td-21", "pcb_solidshell_file.cdb")

#Creating paths
import os
myrunpath = r"C:\Users\vvenkata\OneDrive - ANSYS, Inc\Desktop\TRR\SD_Micro\apdlrun"
myjobname = "PCB"
myrstpath =os.path.join(myrunpath,myjobname+".rst")

#Launch mapdl
mapdl = launch_mapdl(jobname ="PCB",run_location =myrunpath,override =True)

#Verify if your mapdl is Listening: Verify the versions of MAPDL and PYMAPDL
print(mapdl)

#reads the cdb from my directory
mapdl.cdread("DB", pcb_shell_file)#Use shell_cdb_path or solid_cdb_path or solidshell_cdb_path based on what you want to run

# Prints the global status : Elements Nodes Loads ANTYPE from the imported cdb
print(mapdl.run("GSTAT"))

#Prints Element Types and Keyoptions 
print(mapdl.etlist("ALL"))

#Plot Elements  
mapdl.eplot(background='w',show_edges=True,notebook=True)

#Solve the model
mapdl.run("/SOLU")
print(mapdl.run("solve"))

#Exit Mapdl
mapdl.exit()

#Use DPF for Post processing 
from ansys.dpf import post
solution_path = 'file.rst'
solution = post.load_solution(solution_path)
print(solution)

#Plot First 6 Mode shapes
solution.displacement(time_scoping=1).vector.plot_contour()
solution.displacement(time_scoping=2).vector.plot_contour()
solution.displacement(time_scoping=3).vector.plot_contour()
solution.displacement(time_scoping=4).vector.plot_contour()
solution.displacement(time_scoping=5).vector.plot_contour()
solution.displacement(time_scoping=6).vector.plot_contour()

#Get the number of frequencies and frequency values and Plot Barchart using matplotlib
nset = solution.time_freq_support.n_sets
freq = solution.time_freq_support.time_frequencies.data
import matplotlib.pyplot as plt; plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
objects = range(1,nset+1,1)
y_pos = np.arange(len(objects))
performance = freq

plt.bar(y_pos, performance, align='center', alpha=0.5)
plt.xticks(y_pos, objects)
plt.ylabel('Frequency Hz')
plt.xlabel('Modes')
plt.title('Modes vs Frequency')
for index,data in enumerate(performance):
    plt.text(x=index , y =data+1 , s=f"{round(data,1)}")
plt.tight_layout()
plt.show()