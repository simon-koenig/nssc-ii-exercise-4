import matplotlib.pyplot as plt
import numpy as np




def plot_graph(elements, nodes, T, dT, qi, nx, ny, variation):
    element_coords = []
    for el in elements:
        _coords = []
        for node_ind in el:
            _coords.append(nodes[int(node_ind)])
        element_coords.append(np.mean(np.array(_coords), axis = 0 ))

    el_x, el_y = zip(*element_coords)
    dT_x, dT_y = zip(*dT)
    qi_x, qi_y = zip(*qi)

    T = T.round(1)
    x,y = np.array(nodes).T
    T_rs = T.reshape((nx,ny))
    x_rs = x.reshape((nx,ny))
    y_rs = y.reshape((nx,ny))

    def trunc(values, decs=0):
        return np.trunc(values*10**decs)/(10**decs)
    levels = list(set(trunc(T, decs=4).flatten()))
    levels.sort()
    clev = np.arange(T.min(),T.max() + 0.5,0.5) #Adjust the .001 to get finer gradient


    #PLOTS
    fig,ax=plt.subplots(2,1,figsize = (16,30))
    cp0 = ax[0].contourf(
        x_rs, y_rs, T_rs,
        levels = clev, #old version: levels
        cmap = 'coolwarm'
        )
    cp1 = ax[1].contourf(
        x_rs, y_rs, T_rs,
        levels = clev, #old version: levels
        cmap = 'coolwarm'
        )
    fig.colorbar(cp0,
        ax=ax[0], shrink=0.7,
        anchor = (-0.4,0.3),
        use_gridspec=False)

    fig.colorbar(cp1,
        ax=ax[1], shrink=0.7,
        anchor = (-0.4,0.3),
        use_gridspec=False
        )

    #ax.clabel(cp, cp.levels, inline=True, fontsize=10)
    ax[0].triplot(*zip(*nodes), triangles=elements, color='black', lw=0.8)
    ax[1].triplot(*zip(*nodes), triangles=elements, color='black', lw=0.8)

    q1 = ax[0].quiver(el_x,el_y, dT_x, dT_y,
        fc = 'black', width = 0.006,
        ec = 'white', linewidth = 0.5,
        angles = 'xy',
        scale = 55000,
    #    scale_units = 'xy'
        )
    ax[0].quiverkey(q1, X=1.1, Y=0.85, U=6000,
                label='Temperature \n Gradient \n 6000 [K/m]', 
                labelpos='N',
                fontproperties = {'size':10})

    q2 = ax[1].quiver(el_x,el_y, qi_x, qi_y,
        fc = 'black', width = 0.006,
        ec = 'white', linewidth = 0.5,
        angles = 'xy',
        scale = 22000000
        )
    ax[1].quiverkey(q2, X=1.1, Y=0.9, U=2000000,
                label='Heat Flux \n 2000 [kW/$m^2$]',
                labelpos='N',
                fontproperties = {'size':10})


    xx,yy = zip(*nodes)
    step_x = (max(xx) - min(xx)) / len(xx)
    step_y = (max(yy) - min(yy)) / len(yy)
    #fig.set_title('Filled Contours Plot')
    fig.subplots_adjust(right = 0.7)
    ax[0].set_xlim([min(x)- 5*step_x, max(x) + 5*step_x])
    ax[0].set_ylim([min(y)- 15*step_y, max(y) + 15*step_y])
    ax[1].set_xlim([min(x)- 5*step_x, max(x) + 5*step_x])
    ax[1].set_ylim([min(y)- 15*step_y, max(y) + 15*step_y])
    ax[0].set_title('Temperature Gradient', fontsize = 20)
    ax[0].set_xlabel('x in [m]', fontsize = 15)
    ax[0].set_ylabel('y in [m]', fontsize = 15)
    ax[1].set_title('Heat Flux', fontsize = 20)
    ax[1].set_xlabel('x (m)', fontsize = 15)
    ax[1].set_ylabel('y (m)', fontsize = 15)

    #plt.show()
    plt.savefig(f"Temperature_Gradient_Heat_Flux_{str(variation).replace('.','_')}")
