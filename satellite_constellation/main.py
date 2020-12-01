import matplotlib.pyplot as plt
import numpy as np
from satellite_constellation.Constellation import *
from satellite_constellation.SceneCreator import *
from satellite_constellation.utils import *




if __name__ == '__main__':

    walker_constellation = WalkerConstellation(3, 1, 1, 30, 1000, 0, 20)
    r = walker_constellation.altitude + heavenly_body_radius[walker_constellation.focus]
    # print(walker_constellation.perigee_positions[2])
    print(walker_constellation.perigee_positions)
    RAAN = (walker_constellation.raan[1]) * math.pi / 180
    roty = 0
    rotx = walker_constellation.inclination * math.pi / 180

    t = np.linspace(0, 2 * math.pi, 100)
    plt.plot(r * np.cos(t), r * np.sin(t), 0)
    x, y, z = r * np.cos(t), r * np.sin(t), 0 * t

    for idx in range(100):
        temp_mat = np.array([x[idx], y[idx], z[idx]])
        rot_mat = rotate(temp_mat,rotx,'x')
        x[idx] = rot_mat[0]
        y[idx] = rot_mat[1]
        z[idx] = rot_mat[2]

    x1, y1, z1 = r * np.cos(t), r * np.sin(t), 0 * t

    fig = plt.figure()

    sat_x, sat_y, sat_z = polar2cart(r, 70 * math.pi / 180, 45 * math.pi / 180)
    print(sat_x, sat_y, sat_z)
    perspectives = [[0,0],[90,0],[45,45]]

    ax = [plt.subplot(2, 2, 1, projection='3d'), plt.subplot(2, 2, 2, projection='3d'), plt.subplot(2, 2, 3, projection='3d')]
    for idx in range(3):
        ax[idx].view_init(elev=perspectives[idx][0], azim=perspectives[idx][1])
        ax[idx].set_xlim(-r, r)
        ax[idx].set_ylim(-r, r)
        ax[idx].set_zlim(-r, r)
        ax[idx].plot(x, y, z, '--', linewidth=0.5)
        ax[idx].plot(x1, y1, z1, '--', linewidth=0.5, color='r')
        ax[idx].zaxis.set_tick_params(labelsize=3)
        ax[idx].xaxis.set_tick_params(labelsize=3)
        ax[idx].yaxis.set_tick_params(labelsize=3)
        ax[idx].set_xlabel("X", fontsize=3)
        ax[idx].set_ylabel("Y", fontsize=3)
        ax[idx].set_zlabel("Z", fontsize=3)
        for idy in range(walker_constellation.num_sats):
            x_i, y_i, z_i = sat_x, sat_y, sat_z = polar2cart(r, (90) * math.pi / 180,
                                                             (walker_constellation.perigee_positions[
                                                                 idy] +
                                                              + walker_constellation.ta[idy]) * math.pi / 180)
            coords = np.array([x_i, y_i, z_i])
            coord_rot = rotate(coords, 90 * math.pi / 180, 'z')
            coord_rot = rotate(coord_rot, walker_constellation.inclination * math.pi / 180, 'x')
            # coord_rot = rotate(coords, walker_constellation.raan[idy]* math.pi / 180, 'z')
            ax[idx].scatter(coord_rot[0], coord_rot[1], coord_rot[2])



    plt.savefig('../../plot.png', dpi=300, bbox_inches='tight')
