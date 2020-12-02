import matplotlib.pyplot as plt
import numpy as np
from satellite_constellation.Constellation import *
from satellite_constellation.SceneCreator import *
from satellite_constellation.utils import *


def draw_walker(walker_constellation):
    r = walker_constellation.altitude + heavenly_body_radius[walker_constellation.focus]

    t = np.linspace(0, 2 * math.pi, 100)

    x1, y1, z1 = r * np.cos(t), r * np.sin(t), 0 * t

    fig = plt.figure()

    perspectives = [[0, 0], [90, 0], [45, 45]]

    ax = [plt.subplot(2, 2, 1, projection='3d'), plt.subplot(2, 2, 2, projection='3d'),
          plt.subplot(2, 2, 3, projection='3d')]
    for idx in range(3):
        ax[idx].view_init(elev=perspectives[idx][0], azim=perspectives[idx][1])
        ax[idx].set_xlim(-r, r)
        ax[idx].set_ylim(-r, r)
        ax[idx].set_zlim(-r, r)
        ax[idx].plot(x1, y1, z1, '--', linewidth=0.1, color='r')  # Plot equatorial circle
        ax[idx].zaxis.set_tick_params(labelsize=3)
        ax[idx].xaxis.set_tick_params(labelsize=3)
        ax[idx].yaxis.set_tick_params(labelsize=3)
        ax[idx].set_xlabel("X", fontsize=3)
        ax[idx].set_ylabel("Y", fontsize=3)
        ax[idx].set_zlabel("Z", fontsize=3)

        for idy in range(walker_constellation.num_planes):  # Plot orbital planes
            ang = idy * 360 / walker_constellation.num_planes
            t = np.linspace(0, 2 * math.pi, 100)
            plt.plot(r * np.cos(t), r * np.sin(t), 0)
            x, y, z = r * np.cos(t), r * np.sin(t), 0 * t
            for idz in range(100):
                coords = np.array([x[idz], y[idz], z[idz]])
                rot_coords = rotate(coords, walker_constellation.inclination * math.pi / 180, 'x')
                rot_coords = rotate(rot_coords, ang * math.pi / 180, 'z')
                x[idz] = rot_coords[0]
                y[idz] = rot_coords[1]
                z[idz] = rot_coords[2]
            ax[idx].plot(x, y, z, '--', linewidth=0.5)

        for idy in range(walker_constellation.num_planes):  # Plot satellites
            for idz in range(walker_constellation.sats_per_plane):
                ctr = idz + idy * (walker_constellation.sats_per_plane)

                x_i, y_i, z_i = sat_x, sat_y, sat_z = polar2cart(r, (90) * math.pi / 180,
                                                                 (walker_constellation.perigee_positions[
                                                                      ctr] +
                                                                  + walker_constellation.ta[ctr]
                                                                  + walker_constellation.raan[ctr]) * math.pi / 180)
                coords = np.array([x_i, y_i, z_i])
                coords = rotate(coords, 90 * math.pi / 180, 'z')
                coords = rotate(coords, walker_constellation.inclination * math.pi / 180, 'x')
                coords = rotate(coords, (walker_constellation.raan[ctr]) * math.pi / 180, 'z')
                ax[idx].scatter(coords[0], coords[1], coords[2])

    plt.savefig('../../walker_plot.png', dpi=300, bbox_inches='tight')


def draw_flower(flower_constellation):
    a = flower_constellation.semi_major
    b = a * math.sqrt(1 - math.pow(flower_constellation.eccentricity, 2))
    f = (flower_constellation.altitude + heavenly_body_radius[flower_constellation.focus]) * 10 ** 3
    disp = a - f

    t = np.linspace(0, 2 * math.pi, 100)

    r = heavenly_body_radius[flower_constellation.focus] * 10 ** 3
    x1, y1, z1 = r * np.cos(t), r * np.sin(t), 0 * t
    x2, y2, z2 = r * np.cos(t), 0 * t, r * np.sin(t)

    fig = plt.figure()

    r = a

    perspectives = [[0, 0], [90, 0], [45, 45]]

    ax = [plt.subplot(2, 2, 1, projection='3d'), plt.subplot(2, 2, 2, projection='3d'),
          plt.subplot(2, 2, 3, projection='3d')]
    for idx in range(3):
        ax[idx].view_init(elev=perspectives[idx][0], azim=perspectives[idx][1])
        ax[idx].set_xlim(-3 / 2 * a, 3 / 2 * a)
        ax[idx].set_ylim(-3 / 2 * a, 3 / 2 * a)
        ax[idx].set_zlim(-3 / 2 * a, 3 / 2 * a)
        ax[idx].plot(x1, y1, z1, '--', linewidth=1.5, color='r')  # Plot equatorial circle
        ax[idx].plot(x2, y2, z2, '--', linewidth=1.5, color='r')  # Plot equatorial circle

        ax[idx].zaxis.set_tick_params(labelsize=3)
        ax[idx].xaxis.set_tick_params(labelsize=3)
        ax[idx].yaxis.set_tick_params(labelsize=3)
        ax[idx].set_xlabel("X", fontsize=3)
        ax[idx].set_ylabel("Y", fontsize=3)
        ax[idx].set_zlabel("Z", fontsize=3)
        for idy in range(flower_constellation.num_orbits):  # Plot orbital planes
            x, y, z = disp + a * np.cos(t), b * np.sin(t), 0 * t
            for idz in range(100):
                coords = np.array([x[idz], y[idz], z[idz]])
                coords = rotate(coords, flower_constellation.raan[idy] * math.pi / 180, 'z')
                coords = rotate(coords, flower_constellation.inclination * math.pi / 180, 'x')
                x[idz] = coords[0]
                y[idz] = coords[1]
                z[idz] = coords[2]
            ax[idx].plot(x, y, z, '--', linewidth=0.5)

        for idy in range(flower_constellation.num_satellites):  # Plot satellites
            ang = (flower_constellation.mean_anomaly[idy] + 180) * math.pi / 180
            x_i, y_i, z_i = disp + a * np.cos(ang), b * np.sin(ang), 0
            coords = np.array([x_i, y_i, z_i])
            coords = rotate(coords, flower_constellation.raan[idy] * math.pi / 180, 'z')
            coords = rotate(coords, flower_constellation.inclination * math.pi / 180, 'x')
            ax[idx].scatter(coords[0], coords[1], coords[2], s=2)

    plt.savefig('../../flower_plot.png', dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    draw_walker(WalkerConstellation(8, 2, 1, 45, 1000, 0, 20))
    draw_flower(FlowerConstellation(31, 11, 30, 7, 10, 0, 45, 25000, 0))
