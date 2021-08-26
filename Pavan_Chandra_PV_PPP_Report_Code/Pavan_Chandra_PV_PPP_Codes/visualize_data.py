'''
=========================************PERSONAL PROGRAMMING PROJECT****************===========================================!!!
CODE AUTHOR: PADMANABHA PAVAN CHANDRA VUNDURTHY(Mat.Nr.62750)(MSc. COMPUTATIONAL MATERIALS SCIENCE)
UNDER THE SUPERVISION OF DR. ARUNA PRAKASH
TITLE: BIG DATA VISUALIZATION OF LARGE AND ULTRA LARGE SCALE ATOMISTIC SIMULATION (UPTO 20 MILLION ATOMS) 
FILES NECESSARY AS INPUT: PARAMETER FILE, SORTED DATA FILES CONTAINING THE PROPERTIES NECESSARY FOR ANALYSIS THAT ARE MENTIONED IN PARAMETER FILE
OUTPUT FILES: CONTAINS THE AVERAGE PROPERTIES OF THE GRAINS WITH DIFFERENT STRUCTURE TYPES
THE FOLLOWING PROGRAM WORKS FOR 'N' NUMBER OF FILES AND ANY AMOUNT OF DATA PROVIDED FOR ANALYSIS
'''
#Importing Libraries
import numpy as np
import matplotlib.pyplot as plt
import argparse     #useful for parsing data from user for analysis
import sys

#Plot Set 1: Volume_Plots
def Volume_Plots(f,gr_vol_frac,vol_list,tot_cong_stress,nng,vol_mean):    
    count=f 
    plt.figure()
    d = list()
    d = gr_vol_frac[count-1]
    N=len(d)
    width=d*2*np.pi
    theta=np.zeros(N)
    #Assign angle between polar plots
    for j in range(1,N):
        theta[j]=theta[j-1]+width[j-1]/2+width[j]/2
    #
    radii=vol_list[count-1]  
    radii_stress = tot_cong_stress[count-1]
    ax = plt.subplot(111, projection='polar')
    bars = ax.bar(theta, radii, width=width, bottom=0.0,linewidth=1,edgecolor='black') 
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_theta_zero_location('N')  
    # Use custom colors and opacity
    c=nng[count-1]
    m=np.mean(radii_stress)
    #for j in range(0,N):
        #bars[j].set_color(colors[j])
        #plt.text(theta[j],radii[j]+0.3,c[j])
    for r, bar in zip(radii_stress, bars) :
        print (r)
        #print (plt.cm.bwr(r / 10.))
        bar.set_facecolor(plt.cm.bwr((r-m)/10.0))
        bar.set_alpha(0.5)
    for j in range(0,N):
        plt.text(theta[j],3*max(radii)/4,c[j],fontsize='20')
    for i in range(0,360):
        if i==359:
            ax.plot([i,i+np.pi/180],[vol_mean,vol_mean],linewidth=1.5, color='black',label='avg. volume of polycrystal')
        else:
            ax.plot([i,i+np.pi/180],[vol_mean,vol_mean],linewidth=1.5, color='black')            
    ax.legend(fontsize='20',loc='upper center',bbox_to_anchor=(0.55, 1.25))
    plt.savefig("plotvol_"+str(count)+".jpg",bbox_inches='tight', pad_inches=0.1)
    plt.show()
    return 
#Plot2 Stress Plots    
def Stress_Plots(q,gr_vol_frac,vol_list,tot_cong_stress,nng,Macro_VM,stress_in_cong,hydstress_cong):
    count=q
    plt.figure()
    d = list()
    d = gr_vol_frac[count-1]
    N=len(d)
    width=d*2*np.pi
    theta=np.zeros(N)
    for j in range(1,N):
        theta[j]=theta[j-1]+width[j-1]/2+width[j]/2
    radii=vol_list[count-1]  
    radii_stress = tot_cong_stress[count-1]
    ax = plt.subplot(111, projection='polar')
    bars = ax.bar(theta, radii_stress, width=width, bottom=0.0,linewidth=1,edgecolor='black')   
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_theta_zero_location('N') 
        # Use custom colors and opacity
    c=nng[count-1]
    m = np.mean(radii_stress)
    for r, bar in zip(radii_stress, bars) :
        print (r)
        bar.set_facecolor(plt.cm.bwr((r-m)/10))
        bar.set_alpha(0.5)
    #writing the grain numbers beyond plot
    for j in range(0,N):
        plt.text(theta[j],3*max(radii_stress)/4,c[j],fontsize='20')
    #Average stress in a Polycrystal plot
    for i in range(0,360):
        if i==359:
            ax.plot([i,i+np.pi/180],[Macro_VM,Macro_VM],linewidth=1.5, color='black',label='Average stress in a polycrystal')
        else:
            ax.plot([i,i+np.pi/180],[Macro_VM,Macro_VM],linewidth=1.5, color='black')            
    ax.legend(fontsize='20',loc='upper center',bbox_to_anchor=(0.55, 1.25))
    # FOR CONGLOMERATE USE VM_STRESS
    for i in range(0,360):
        if i==359:
            ax.plot([i,i+np.pi/180],[stress_in_cong[q-1],stress_in_cong[q-1]],linewidth=1.5, color='magenta',label='Average Stress in a conglomerate')
        else:
            ax.plot([i,i+np.pi/180],[stress_in_cong[q-1],stress_in_cong[q-1]],linewidth=1.5, color='magenta')            
    ax.legend(fontsize='20',loc='upper center',bbox_to_anchor=(0.55, 1.25))
    plt.savefig("plotstress_"+str(count)+".jpg",bbox_inches='tight', pad_inches=0.1)
    plt.show()
    return 
#Plot 3 Distribution Plots
def Distribution_Plots(r,gr_size,polarized_stress_hyd,hyd_stress,triaxiality,polarized_stress_vm,vm_stress):
    fig, host = plt.subplots()
    fig.subplots_adjust(right=0.75)
    par1 = host.twinx()
    if r ==1:

        par2 = host.twinx()
        par2.spines["right"].set_position(("axes", 1.2))
        make_patch_spines_invisible(par2)
        par2.spines["right"].set_visible(True)
        p1, = host.plot(gr_size,polarized_stress_hyd, "b*", label="Polarized Hrdrostatic stress")
        p2, = par1.plot(gr_size,hyd_stress, "r*", label="Hydrostatic stress")
        p3, = par2.plot(gr_size,triaxiality, "g*", label="Triaxiality")
        host.set_ylim(min(polarized_stress_hyd)-polarized_stress_hyd[1], max(polarized_stress_hyd)+polarized_stress_hyd[0])
        par2.set_ylim(min(triaxiality)-triaxiality[0],max(triaxiality)+triaxiality[2])
        host.set_xlabel("Grain Size")
        host.set_ylabel("Polarized Values")
        par1.set_ylabel("Hydrostatic Stress")
        par2.set_ylabel("Triaxiality")

        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())
        par2.yaxis.label.set_color(p3.get_color())

        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        par2.tick_params(axis='y', colors=p3.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)
        lines = [p1, p2, p3]
        host.legend(lines, [l.get_label() for l in lines])
        plt.title("Grain size vs polarized hydrostatic stress, hydrostatic_stress, triaxiality")

    if r==2:
        p1, = host.plot(gr_size,polarized_stress_vm, "b*", label="Polarized Von_Mises stress")
        p2, = par1.plot(gr_size,vm_stress, "r*", label="Von_Mises stress")
        host.set_ylim(min(polarized_stress_vm),max(polarized_stress_vm))
        par1.set_ylim(min(vm_stress)-vm_stress[1],max(vm_stress)+vm_stress[2])
        host.set_xlabel("Grain Size")
        host.set_ylabel("Polarized VM Values")
        par1.set_ylabel("Von_Mises Stress")
        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())
        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)
        lines = [p1, p2]
        host.legend(lines, [l.get_label() for l in lines])
        plt.title("Grain size vs polarized Von Mises,Von Mises stress")
   
    if r==3:
        p1, = host.plot(gr_size,hyd_stress, "b*", label="Hrdrostatic stress")
        p2, = par1.plot(gr_size,triaxiality, "r*", label="Triaxiality")
        host.set_ylim(min(hyd_stress),max(hyd_stress))
        host.set_xlabel("Grain Size")
        host.set_ylabel("Hydro_stress")
        par1.set_ylabel("Triaxiality")
        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())
        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)
        lines = [p1, p2]
        host.legend(lines, [l.get_label() for l in lines])
        plt.title("Grain Size vs hydrostatic_stress,triaxiality")
    plt.show()
    return 

##Below function is obtained from source: https://matplotlib.org/gallery/ticks_and_spines/multiple_yaxis_with_spines.html
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plt_circ(x_cen,y_cen,rad,plt_style):
    x_circ=np.zeros(360)
    y_circ=np.zeros(360)
    for i in range(0,360):
        x_circ[i]=x_cen+rad*np.cos(i*np.pi/180)
        y_circ[i]=y_cen+rad*np.sin(i*np.pi/180)       
    plt.plot(x_circ,y_circ,plt_style)
'''
Implementation of Argparse Module for taking User Input from Command Line Interface
'''
parser = argparse.ArgumentParser(description='Big Data Visulization of Large and Ultra Large Scale Atomistic Simulations')
#SET: ranges from 1 to 4 types of visualization Plots
parser.add_argument('set',type=int,help='Select the set of visualization:Set1:Average Volume Polar Plots;Set2:Average Stress Polar Plots;Set3:Distribution Plots;Set4:All the plots')
parser.add_argument('Max_GrainNumber', type=int, help='Specify the Maximum grains in the polycrystal')
parser.add_argument('Max_Struct_Typ', type=int, help='Specify the Maximum Structure type in the polycrystal')
parser.add_argument('start_grnnum', type=int, help='Specify the start of grain number the user is interested in..')
parser.add_argument('iteration', type=int, help='an integer specifying the iteration of grain numbers at interest')
parser.add_argument('end_grnnum', type=int, help='Specify the end of grain number the user is interested in..')
parser.add_argument('structure_type', help='Mention the interested structure type/types to be considered for visualization\nA:Other defects\nB:Fcc\nC:Hcp\nD:Fcc+hcp\nE:fcc+hcp+other defects\nF:all structure types')
parser.add_argument('Avg_prop_file_name', default = [], nargs='?', type=argparse.FileType('r'), help=\
    'File containing the Average properties of grains analysed from the Big Data file')
parser.add_argument('Neighbour_Grains',default=[],nargs='?',type=argparse.FileType('r'),help='File containing the Nearest Neighbouring Grains')
#The User input are read into the 'sys.argv' list in their original format
args = parser.parse_args()
for l in range(1,7):
    sys.argv[l]=int(sys.argv[l])
print(sys.argv)
##==============================================
max_struct_typ = sys.argv[3]  #struct max is 2 if there is only other defects, FCC and HCP
data = np.loadtxt(sys.argv[8],comments=['#'])
if sys.argv[7]=='A':
    counter = 0
elif sys.argv[7]=='B':
    counter = 1
elif sys.argv[7]=='C':
    counter = 2
elif sys.argv[7]=='D':
    counter = max_struct_typ+1
elif sys.argv[7]=='E':
    counter = max_struct_typ+2
elif sys.argv[7]=='F':
    counter = max_struct_typ+3  
avg_prop_data = data[sys.argv[2] * (counter):(sys.argv[2]*(counter+1)),:]    #Extracting interested Data from the average properties  
vol = avg_prop_data[:,10]    #Extracting the average volume for a structure type
stress = avg_prop_data[:,4:10]   #Extracting average stresses for a structure type
vm_stress = avg_prop_data[:,12]
#Reading the Neighbour_Grains data from the file
words=list()
t = list()
fname = open(sys.argv[9],"r")
for line in fname:
    t = line.strip().split('\t')
    words.append(t)
nng = [[int(j) for j in i] for i in words]  #converting everything into integers
##-------------------------------------------------------
##STARSTS THE ACTUAL CALCULATION FOR DATA REQUIRED FOR VISUALIZATION:
#GRAIN VOL FRACTION, AVG. STRESS, INDIVIDUAL STRESS, TRI AXIALITY
temp = list()
tot_cong_vol = list()  #total volume of each conglomerate
tot_cong_stress = list()  #list of stresses in each conglomerate
stress_in_cong = list()
hydstress_cong = list()
gr_vol_frac = list()
vol_list= list()        
for a in range(0,len(nng)):
    temp = nng[a]
    t_vol = 0.0
    temp_vol = list()
    t_stress= list()
    temp_stress = [0.0,0.0,0.0,0.0,0.0,0.0]
    for b in range(0,len(temp)):
        c = temp[b]-1
        temp_vol.append(vol[c])
        t_stress.append(vm_stress[c])
        t_vol += vol[c]
        temp_stress += (stress[c,:]*vol[c])
    temp_stress = temp_stress/t_vol
    tot_cong_vol.append(t_vol)
    vol_list.append(temp_vol)           #list of lists of volumes of each conglomerate
    gr_vol_frac.append(temp_vol/t_vol)      #list of lists of Grain Volume Fraction for each conglomerate
    tot_cong_stress.append(t_stress)        # list of lists of vm stresses directly from data for each conglomerate 
    #Calculating the Von Mises Stress of a conglomerate individually by below equation and appending a tensor of stresses for each conglomerate
    stress_in_cong.append(np.sqrt((0.5)*((temp_stress[0]-temp_stress[1])**2+(temp_stress[1]-temp_stress[2])**2+(temp_stress[2]-temp_stress[0])**2)+3*((temp_stress[3])**2+(temp_stress[4])**2+(temp_stress[5])**2)))
    hydstress_cong.append((temp_stress[0]+temp_stress[1]+temp_stress[2])/3)   #Calculating the Hydrostatic Stresses from the stress tensor for each conglomerate
    #the only difference between stress_in_cong and tot_cong_stress is that the former is calculated via equation whereas the latter is from available data 
vol_mean = np.mean(vol)     #useful for plotting in function 1 for average volume plots
##==================================================#Calculations for 3rd Function Distribution Plots!!Nothing to do with Neighbouring Grains
gr_size = list() 
hyd_stress = list()
triaxiality=list()
int_stress = [0.0,0.0,0.0,0.0,0.0,0.0]
for i in range(0,len(vol)):
    int_stress += stress[i,:]*vol[i]
    hyd_stress.append((stress[i,0]+stress[i,1]+stress[i,2])/3)      #list of hydrostatic stresses
vol_sum = np.sum(vol)
int_stress = int_stress/vol_sum
#Single values below
Macro_VM = np.sqrt((0.5)*((int_stress[0]-int_stress[1])**2+(int_stress[1]-int_stress[2])**2+(int_stress[2]-int_stress[0])**2)+3*((int_stress[3])**2+(int_stress[4])**2+(int_stress[5])**2))
Hydro_stress = np.mean(hyd_stress)
#list of triaxiality factor below 
triaxiality = hyd_stress/vm_stress
TDF = Hydro_stress/Macro_VM
#Below are the Polarized Values
polarized_stress_vm = vm_stress-Macro_VM
polarized_stress_hyd = hyd_stress-Hydro_stress
polarized_triaxiality = triaxiality-TDF
#Calculate the Equivalent Grain Size using Volume of a sphere equation
gr_size = np.cbrt(vol*3/(4*np.pi))
D = np.sum(gr_size)
grain_frac = gr_size/D
#====================================================
#Looping Condition over Parsed Data !!!
#Set 1
if sys.argv[1]==1:
    if sys.argv[5]==0:
        f = sys.argv[4]
        Volume_Plots(f,gr_vol_frac,vol_list,tot_cong_stress,nng,vol_mean)
    else:
        for f in range(sys.argv[4],sys.argv[6],sys.argv[5]):
            Volume_Plots(f,gr_vol_frac,vol_list,tot_cong_stress,nng,vol_mean)
#Set 2
if sys.argv[1]==2:
    if sys.argv[5]==0:
        q = sys.argv[4]
        Stress_Plots(q,gr_vol_frac,vol_list,tot_cong_stress,nng,Macro_VM,stress_in_cong,hydstress_cong)
    else:
        for q in range(sys.argv[4],sys.argv[6],sys.argv[5]):
            Stress_Plots(q,gr_vol_frac,vol_list,tot_cong_stress,nng,Macro_VM,stress_in_cong,hydstress_cong)
#Set 3
if sys.argv[1]==3:
    for r in range(1,4):
        Distribution_Plots(r,gr_size,polarized_stress_hyd,hyd_stress,triaxiality,polarized_stress_vm,vm_stress)
#Set 4
#call all functions
if sys.argv[1]==4:
    for r in range(1,4):
        Distribution_Plots(r,gr_size,polarized_stress_hyd,hyd_stress,triaxiality,polarized_stress_vm,vm_stress)
    if sys.argv[5]==0:
        f = sys.argv[4]
        q = sys.argv[4]
        Volume_Plots(f,gr_vol_frac,vol_list,tot_cong_stress,nng,vol_mean)
        Stress_Plots(q,gr_vol_frac,vol_list,tot_cong_stress,nng,Macro_VM,stress_in_cong,hydstress_cong)
    else:
        for f,q in range(sys.argv[4],sys.argv[6],sys.argv[5]):
            Volume_Plots(f,gr_vol_frac,vol_list,tot_cong_stress,nng,vol_mean) 
            Stress_Plots(q,gr_vol_frac,vol_list,tot_cong_stress,nng,Macro_VM,stress_in_cong,hydstress_cong)
#End of the Program!!!!!!
#==================================================!!!

