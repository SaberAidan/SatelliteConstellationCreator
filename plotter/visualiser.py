from satellite_constellation.utils import *
import plotly.graph_objects as go


def draw_walker_plotly(walker_constellation, satellites=True, orbits=True, links=False, sensor_regions=False):
    sat_coords = np.empty([walker_constellation.num_satellites, 3])
    r = walker_constellation.altitude + heavenly_body_radius[walker_constellation.focus]
    fig = go.Figure()

    max_dist = np.array([])

    if orbits:
        for idy in range(walker_constellation.num_planes):  # Plot orbital planes
            ang = idy * 360 / walker_constellation.num_planes
            t = np.linspace(0, 2 * math.pi, 100)
            x, y, z = r * np.cos(t), r * np.sin(t), 0 * t
            for idz in range(100):
                coords = np.array([x[idz], y[idz], z[idz]])
                rot_coords = rotate(coords, walker_constellation.inclination * math.pi / 180, 'x')
                rot_coords = rotate(rot_coords, ang * math.pi / 180, 'z')
                x[idz] = rot_coords[0]
                y[idz] = rot_coords[1]
                z[idz] = rot_coords[2]

            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='Orbit ' + str(idy + 1)))
            max_dist = np.append(math.sqrt(np.sum(coords ** 2)), np.max(y))

    if satellites:
        d_sat = 0

        for sats in walker_constellation.satellites:
            coord = sat_to_xyz(sats)
            sat_coords[d_sat] = coord
            fig.add_trace(
                go.Scatter3d(x=[coord[0]], y=[coord[1]], z=[coord[2]], mode='markers', name='real_sat '
                                                                                            + str(
                    d_sat), showlegend=False))
            d_sat += 1

    if links:
        for idy in range(0, sat_coords.shape[0]):  # Draw line of sight between satellites
            ctr = 0
            for idz in range(0, sat_coords.shape[0]):
                if idz != idy:
                    temp_coords = np.append([sat_coords[idy, :]], [sat_coords[idz, :]], axis=0)
                    if not sphere_intercept(temp_coords[0], temp_coords[1],
                                            heavenly_body_radius[walker_constellation.focus]):
                        ctr += 1
                        fig.add_trace(
                            go.Scatter3d(x=temp_coords[:, 0], y=temp_coords[:, 1], z=temp_coords[:, 2], mode='lines',
                                         name='link {0} to {1}'.format(idy, idz), showlegend=False,
                                         line=dict(color='green', width=0.5, dash='dash')))

    if sensor_regions:

        for idy in range(walker_constellation.num_satellites):
            rE = heavenly_body_radius[walker_constellation.focus]
            v_r = rE * sat_coords[idy] / math.sqrt(np.sum(sat_coords[idy] ** 2))
            t = np.linspace(0, 2 * math.pi, 100)
            projection_coords = np.array([[0, 0, 0]])

            basis = v_r / math.sqrt(np.sum(v_r ** 2))

            basis = np.array(basis)
            x = np.array([0, 1, 0])
            basis2 = np.cross(basis, x)

            half_width = (walker_constellation.beam_width/2 * math.pi/180)
            if half_width > np.arcsin(rE / r):
                half_width = np.arcsin(rE / r)
            theta_2 = np.arcsin(np.sin(half_width) / (rE / r)) - half_width

            walker_constellation.earth_coverage_angle = theta_2

            start = rotate(v_r, walker_constellation.earth_coverage_angle, 'custom', basis=basis2)

            for idk in range(0, 100):
                coords = rotate(start, t[idk], 'custom', basis=basis)
                coords = rE * coords/math.sqrt(np.sum(coords**2))

                projection_coords = np.append(projection_coords, [coords], axis=0)

            fig.add_trace(
                go.Scatter3d(x=projection_coords[:, 0], y=projection_coords[:, 1], z=projection_coords[:, 2],
                             mode='lines',
                             name='Projection ' + str(idy), line=dict(color='red', width=5), showlegend=False))

    r = heavenly_body_radius[walker_constellation.focus]
    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    fig.add_trace(go.Surface(x=x, y=y, z=z, name='Earth', showlegend=False, colorscale='blues'))

    max_dist = np.append(max_dist, np.max(y))
    max_val = np.max(max_dist)

    try:
        spherical_earth_map = np.load('map_sphere.npy')
    except FileNotFoundError:
        pass
    else:
        xm, ym, zm = spherical_earth_map.T * 6371
        fig.add_trace(
            go.Scatter3d(x=xm, y=ym, z=zm, name='Earth', mode='lines', showlegend=False,
                         line=dict(color='white', width=1)))

    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            yaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            zaxis=dict(range=[-max_val - 5000, max_val + 5000]), aspectmode='cube')
    )

    fig.show()


def draw_flower_plotly(flower_constellation, satellites=True, orbits=True, links=False):
    sat_coords = np.array([[0, 0, 0]])

    fig = go.Figure()

    a = flower_constellation.semi_major
    b = a * math.sqrt(1 - math.pow(flower_constellation.eccentricity, 2))
    f = (flower_constellation.altitude + heavenly_body_radius[flower_constellation.focus]) * 10 ** 3
    disp = a - f

    t = np.linspace(0, 2 * math.pi, 100)

    max_dist = np.array([])

    if orbits:
        for idy in range(
                min(flower_constellation.num_orbits, flower_constellation.num_satellites)):  # Plot orbital planes
            x, y, z = disp + a * np.cos(t), b * np.sin(t), 0 * t
            for idz in range(100):
                coords = np.array([x[idz], y[idz], z[idz]]) * 10 ** -3
                coords = rotate(coords, flower_constellation.raan[idy] * math.pi / 180, 'z')
                coords = rotate(coords, flower_constellation.inclination * math.pi / 180, 'x')
                x[idz] = coords[0]
                y[idz] = coords[1]
                z[idz] = coords[2]

            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='Orbit ' + str(idy + 1), showlegend=False))
            max_dist = np.append(math.sqrt(np.sum(coords ** 2)), np.max(y))

    if satellites:

        for sats in flower_constellation.satellites:
            coord = sat_to_xyz(sats)
            fig.add_trace(
                go.Scatter3d(x=[coord[0]], y=[coord[1]], z=[coord[2]], mode='markers', name='real_sat '
                                                                                            + "flower",
                             showlegend=False))

        flower_constellation.as_cartesian()

        if satellites:
            d_sat = 0

            for coords in flower_constellation.as_cartesian():
                sat_coords = np.append(sat_coords, [coords], axis=0)

                fig.add_trace(
                    go.Scatter3d(x=[coords[0]], y=[coords[1]], z=[coords[2]], mode='markers', name='real_sat '
                                                                                                   + str(
                        idy), showlegend=False))
                d_sat += 1

    if satellites and links:
        for idy in range(1, sat_coords.shape[0]):  # Draw line of sight between satellites
            for idz in range(1, sat_coords.shape[0]):
                if idz != idy:
                    temp_coords = np.append([sat_coords[idy, :]], [sat_coords[idz, :]], axis=0)
                    if not sphere_intercept(temp_coords[0], temp_coords[1],
                                            heavenly_body_radius[flower_constellation.focus]):
                        fig.add_trace(
                            go.Scatter3d(x=temp_coords[:, 0], y=temp_coords[:, 1], z=temp_coords[:, 2], mode='lines',
                                         name='link {0} to {1}'.format(idy, idz), showlegend=False,
                                         line=dict(color='green', width=0.5, dash='dash')))

    r = heavenly_body_radius[flower_constellation.focus]

    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    fig.add_trace(go.Surface(x=x, y=y, z=z, name='Earth', showlegend=False, colorscale='blues'))

    max_dist = np.append(max_dist, np.max(y))
    max_val = np.max(max_dist)

    try:
        spherical_earth_map = np.load('map_sphere.npy')
    except FileNotFoundError:
        pass
    else:
        xm, ym, zm = spherical_earth_map.T * 6371
        fig.add_trace(
            go.Scatter3d(x=xm, y=ym, z=zm, name='Earth', mode='lines', showlegend=False,
                         line=dict(color='white', width=1)))

    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            yaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            zaxis=dict(range=[-max_val - 5000, max_val + 5000]), aspectmode='cube')
    )

    fig.show()


def draw_soc_plotly(SOC_Constellation, satellites=True, orbits=True, links=False, sensor_regions=False):
    sat_coords = np.array([[0, 0, 0]])
    r = SOC_Constellation.altitude + heavenly_body_radius[SOC_Constellation.focus]
    rE = heavenly_body_radius[SOC_Constellation.focus]
    fig = go.Figure()

    max_dist = np.array([])

    t = np.linspace(0, 2 * math.pi, 100)

    if orbits:
        for idy in range(SOC_Constellation.num_streets):  # Plot orbital planes
            # ang = idy * 360 / SOC_Constellation.num_streets
            ang = SOC_Constellation.raan[idy]
            t = np.linspace(0, 2 * math.pi, 100)
            x, y, z = r * np.cos(t), r * np.sin(t), 0 * t
            for idz in range(100):
                coords = np.array([x[idz], y[idz], z[idz]])
                rot_coords = rotate(coords, SOC_Constellation.inclination * math.pi / 180, 'x')
                rot_coords = rotate(rot_coords, ang * math.pi / 180, 'z')
                x[idz] = rot_coords[0]
                y[idz] = rot_coords[1]
                z[idz] = rot_coords[2]

            fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='Orbit ' + str(idy + 1)))
            max_dist = np.append(max_dist, math.sqrt(np.sum(coords ** 2)))

    if satellites:
        d_sat = 0

        for sats in SOC_Constellation.satellites:
            coord = sat_to_xyz(sats)
            fig.add_trace(
                go.Scatter3d(x=[coord[0]], y=[coord[1]], z=[coord[2]], mode='markers', name='real_sat '
                                                                                            + str(
                    d_sat), showlegend=False))
            d_sat += 1

    if sensor_regions:
        for idy in range(SOC_Constellation.num_streets):  # Plot sensor regions
            for idz in range(SOC_Constellation.sats_per_street):
                projection_coords = np.array([[0, 0, 0]])
                ctr = idz + idy * SOC_Constellation.sats_per_street
                for idk in range(0, 100):
                    coords = np.array(polar2cart(rE, SOC_Constellation.earth_coverage_angle, t[idk]))
                    coords = rotate(coords, (SOC_Constellation.perigee_positions[ctr] + SOC_Constellation.ta[
                        ctr] + 90) * math.pi / 180, 'x')
                    coords = rotate(coords, (SOC_Constellation.raan[idy] + 90) * math.pi / 180, 'z')

                    projection_coords = np.append(projection_coords, [coords], axis=0)

                fig.add_trace(
                    go.Scatter3d(x=projection_coords[:, 0], y=projection_coords[:, 1], z=projection_coords[:, 2],
                                 mode='lines',
                                 name='Projection ' + str(idz), line=dict(color='red', width=5), showlegend=False))
        street_width = heavenly_body_radius[SOC_Constellation.focus] * SOC_Constellation.street_width * math.pi / 180

        t = np.linspace(0, 2 * math.pi, 100)
        x, y, z = rE * np.cos(t), rE * np.sin(t), 0 * t
        for idz in range(100):
            coords = np.array([x[idz], y[idz], z[idz]])
            coords = rotate(coords, SOC_Constellation.inclination * math.pi / 180, 'x')
            x[idz] = coords[0]
            y[idz] = coords[1] + street_width / 2
            z[idz] = coords[2]
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='Street', line=dict(color='green', width=5)))

        x, y, z = rE * np.cos(t), rE * np.sin(t), 0 * t
        for idz in range(100):
            coords = np.array([x[idz], y[idz], z[idz]])
            coords = rotate(coords, SOC_Constellation.inclination * math.pi / 180, 'x')
            x[idz] = coords[0]
            y[idz] = coords[1] - street_width / 2
            z[idz] = coords[2]
        fig.add_trace(go.Scatter3d(x=x, y=y, z=z, mode='lines', name='Street', line=dict(color='green', width=5)))

    if satellites and links:  # Draw line of sight between satellites
        for idy in range(1, sat_coords.shape[0]):
            for idz in range(1, sat_coords.shape[0]):
                if idz != idy:
                    temp_coords = np.append([sat_coords[idy, :]], [sat_coords[idz, :]], axis=0)
                    if not sphere_intercept(temp_coords[0], temp_coords[1],
                                            heavenly_body_radius[SOC_Constellation.focus]):
                        fig.add_trace(
                            go.Scatter3d(x=temp_coords[:, 0], y=temp_coords[:, 1], z=temp_coords[:, 2], mode='lines',
                                         name='link {0} to {1}'.format(idy, idz), showlegend=False,
                                         line=dict(color='green', width=0.5, dash='dash')))

    r = heavenly_body_radius[SOC_Constellation.focus]

    pi = np.pi
    cos = np.cos
    sin = np.sin
    phi, theta = np.mgrid[0:pi:101j, 0:2 * pi:101j]
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)

    fig.add_trace(go.Surface(x=x, y=y, z=z, name='Earth', showlegend=False, colorscale='blues'))

    try:
        spherical_earth_map = np.load('map_sphere.npy')
    except FileNotFoundError:
        pass
    else:
        xm, ym, zm = spherical_earth_map.T * 6371
        fig.add_trace(
            go.Scatter3d(x=xm, y=ym, z=zm, name='Earth', mode='lines', showlegend=False,
                         line=dict(color='white', width=1)))

    max_dist = np.append(max_dist, np.max(y))
    max_val = np.max(max_dist)

    fig.update_layout(
        scene=dict(
            xaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            yaxis=dict(range=[-max_val - 5000, max_val + 5000], ),
            zaxis=dict(range=[-max_val - 5000, max_val + 5000]), aspectmode='cube')
    )

    fig.show()
