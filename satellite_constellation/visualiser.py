import matplotlib.pyplot as plt
from satellite_constellation.Constellation import *
from satellite_constellation.utils import *


def draw_walker(walker_constellation):
    r = walker_constellation.altitude + heavenly_body_radius[walker_constellation.focus]

    if walker_constellation.inclination % 90 == 0:
        plane_range = 180
    else:
        plane_range = 360

    t = np.linspace(0, 2 * math.pi, 100)

    x1, y1, z1 = r * np.cos(t), r * np.sin(t), 0 * t

    fig = plt.figure()

    perspectives = [[0, 0], [0, 45], [90, 0], [60,60]]

    ax = [plt.subplot(2, 2, 1, projection='3d'), plt.subplot(2, 2, 2, projection='3d'),
          plt.subplot(2, 2, 3, projection='3d'), plt.subplot(2, 2, 4, projection='3d')]
    for idx in range(4):

        ax[idx].view_init(elev=perspectives[idx][0], azim=perspectives[idx][1])
        ax[idx].set_xlim(-r, r)
        ax[idx].set_ylim(-r, r)
        ax[idx].set_zlim(-r, r)
        ax[idx].plot(x1, y1, z1, '--', linewidth=1, color='r')  # Plot equatorial circle
        ax[idx].zaxis.set_tick_params(labelsize=3)
        ax[idx].xaxis.set_tick_params(labelsize=3)
        ax[idx].yaxis.set_tick_params(labelsize=3)
        ax[idx].set_xlabel("X", fontsize=3)
        ax[idx].set_ylabel("Y", fontsize=3)
        ax[idx].set_zlabel("Z", fontsize=3)

        # Plot target at revisit time
        xp, yp, zp = (r, 0, 0)
        target_ang = (walker_constellation.revisit_time / (24*60*60)) * 2 * math.pi
        target_coords = np.array([xp, yp, zp])
        ax[idx].scatter(xp, yp, zp, color='black', marker='x')
        target_coords = rotate(target_coords, target_ang, 'z')
        ax[idx].scatter(target_coords[0], target_coords[1], target_coords[2], color='green', marker='x')

        for idy in range(walker_constellation.num_planes):  # Plot orbital planes
            ang = idy * plane_range / walker_constellation.num_planes
            t = np.linspace(0, 2 * math.pi, 100)
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
                ctr = idz + idy * walker_constellation.sats_per_plane

                x_i, y_i, z_i = polar2cart(r, 90 * math.pi / 180,
                                           (walker_constellation.perigee_positions[
                                                ctr] +
                                            + walker_constellation.ta[ctr]
                                            ) * math.pi / 180)
                coords = np.array([x_i, y_i, z_i])
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
    x3, y3, z3 = 0 * t, r * np.cos(t), r * np.sin(t)

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
        ax[idx].plot(x1, y1, z1, '--', linewidth=0.1, color='r')  # Plot equatorial circle
        ax[idx].plot(x2, y2, z2, '--', linewidth=0.1, color='r')  # Plot equatorial circle
        ax[idx].plot(x3, y3, z3, '--', linewidth=0.1, color='r')  # Plot equatorial circle

        # Plot target at revisit time
        xp, yp, zp = (6371*10**3, 0, 0)
        target_ang = (flower_constellation.revisit_time)* 2 * math.pi
        target_coords = np.array([xp, yp, zp])
        ax[idx].scatter(xp, yp, zp, color='black', marker='x')
        target_coords = rotate(target_coords, target_ang, 'z')
        ax[idx].scatter(target_coords[0], target_coords[1], target_coords[2], color='green', marker='x')


        ax[idx].zaxis.set_tick_params(labelsize=3)
        ax[idx].xaxis.set_tick_params(labelsize=3)
        ax[idx].yaxis.set_tick_params(labelsize=3)
        ax[idx].set_xlabel("X", fontsize=3)
        ax[idx].set_ylabel("Y", fontsize=3)
        ax[idx].set_zlabel("Z", fontsize=3)
        for idy in range(
                min(flower_constellation.num_orbits, flower_constellation.num_satellites)):  # Plot orbital planes
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
            ang = (flower_constellation.true_anomaly[idy] + 180) * math.pi / 180
            x_i, y_i, z_i = disp + a * np.cos(ang), b * np.sin(ang), 0
            coords = np.array([x_i, y_i, z_i])
            coords = rotate(coords, flower_constellation.raan[idy] * math.pi / 180, 'z')
            coords = rotate(coords, flower_constellation.inclination * math.pi / 180, 'x')
            ax[idx].scatter(coords[0], coords[1], coords[2], s=2)

    plt.savefig('../../flower_plot.png', dpi=300, bbox_inches='tight')
