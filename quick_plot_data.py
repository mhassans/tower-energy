import pandas as pd
import matplotlib.pyplot as plt
import math

sensorCell = True #sensor cell or trigger cell

u = 2
v = 2
layer = 26

eta_min = 32 
eta_max = 38  

phi_min = 8 
phi_max = 16  

eta_ticks = [a*(2*math.pi/72) for a in range(eta_min, eta_max)] 
phi_ticks = [b*(2*math.pi/72) for b in range(phi_min, phi_max)] 

if (sensorCell):
    cells = pd.read_csv('TCPositions/sensorCell_positions.txt', sep=' ')
    cells.columns= ["layer","waferu","waferv","triggercellu","triggercellv","SC_eta","SC_phi"]
else:
    cells = pd.read_csv('TCPositions/TCPositions_Zminus_siliconOnly.csv', sep=' ')

cell_u_v_layer = cells[(cells['waferu'] == u) & (cells['waferv'] == v) & (cells['layer'] == layer)]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

if(sensorCell):
    plt.scatter(cell_u_v_layer['SC_eta']*-1.0, cell_u_v_layer['SC_phi'])
else:
    plt.scatter(cell_u_v_layer['triggercelleta']*-1.0, cell_u_v_layer['triggercellphi'])
    
ax.set_xticks(eta_ticks)
ax.set_yticks(phi_ticks)
ax.grid(which='both')
plt.show()

