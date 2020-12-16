import math
import numpy as np
import time

heavenly_body_radius = {  # [km]
    "earth": 6371,
    "luna": 1737,
    "mars": 3390,
    "venus": 6052,
    "mercury": 2440,
    "sol": 695700,
    "jupiter": 69911,
    "saturn": 58232,
    "uranus": 25362,
    "neptune": 24622,
    "pluto": 1188,
}

heavenly_body_mass = {  # [kg]
    "earth": 5.972 * 10 ** 24,
    "luna": 73.46 * 10 ** 21,
    "mars": 641.71 * 10 ** 21,
    "venus": 4867.5 * 10 ** 21,
    "mercury": 330.11 * 10 ** 21,
    "sol": 1.9885 * 10 ** 30,
    "jupiter": 1.8982 * 10 ** 27,
    "saturn": 5.6834 * 10 ** 26,
    "uranus": 8.6810 * 10 ** 25,
    "neptune": 1.02413 * 10 ** 26,
    "pluto": 13.03 * 10 ** 21,
}

heavenly_body_period = {  # [days]
    "earth": 1,
    "luna": 27.321661,
    "mars": 1.02595675,
    "venus": 243.0187,
    "mercury": 58.6462,
    "sol": 25.379995,
    "jupiter": 0.41007,
    "saturn": 0.426,
    "uranus": 0.71833,
    "neptune": 0.67125,
    "pluto": 6.38718,
}

constants = {
    "G": 6.67408 * 10 ** (-11),  # Gravitational constant [m^3 kg^-1 s^-2]
    "wE": 7.2921159 * 10 ** (-5),  # Earth angular velocity [rad/s]
    "J2E": 10826269 * 10 ** (-3),  # Earth J2 constant
}


def mod(x, y):
    """
    Mappable modulo function
    :param x: First number
    :param y: Second number
    :return: x % y
    """
    return [a % b for a, b in zip(x, y)]


def proper_round(num, dec=0):  # Add exception check for no decimal point found

    num = str(num)[:str(num).index('.') + dec + 2]
    if num[-1] >= '5':
        return float(num[:-2 - (not dec)] + str(int(num[-2 - (not dec)]) + 1))
    return float(num[:-1])


def rotate_np(vecs, angs, ax='x', rpy=None, basis=None):
    if ax == 'x':
        rotation = np.array([[np.ones(len(angs)), np.zeros(len(angs)), np.zeros(len(angs))],
                             [np.zeros(len(angs)), np.cos(angs), -1 * np.sin(angs)],
                             [np.zeros(len(angs)), np.sin(angs), np.cos(angs)]])
        rotation_mats = np.einsum('ijk -> kij', rotation)
        rotated_vecs = np.einsum('ijk,ik->ij', rotation_mats, vecs)
        return rotated_vecs
    if ax == 'y':
        rotation = np.array([[np.cos(angs), np.zeros(len(angs)), np.sin(angs)],
                             [np.zeros(len(angs)), np.ones(len(angs)), np.zeros(len(angs))],
                             [-np.sin(angs), np.zeros(len(angs)), np.cos(angs)]])
        rotation_mats = np.einsum('ijk -> kij', rotation)
        rotated_vecs = np.einsum('ijk,ik->ij', rotation_mats, vecs)
        return rotated_vecs
    if ax == 'z':
        rotation = np.array([[np.cos(angs), -np.sin(angs), np.zeros(len(angs))],
                             [np.sin(angs), np.cos(angs), np.zeros(len(angs))],
                             [np.zeros(len(angs)), np.zeros(len(angs)), np.ones(len(angs))]])
        rotation_mats = np.einsum('ijk -> kij', rotation)
        rotated_vecs = np.einsum('ijk,ik->ij', rotation_mats, vecs)
        return rotated_vecs

    elif ax == 'c':

        zeros = np.zeros(len(rpy))
        ones = np.ones(len(rpy))

        ang_yaw, ang_pitch, ang_roll = rpy[:, 2], rpy[:, 1], rpy[:, 0]
        ang_yaw *= math.pi / 180
        ang_pitch *= math.pi / 180
        ang_roll *= math.pi / 180

        r_yaw = np.array([[ones, zeros, zeros],
                          [zeros, np.cos(ang_yaw), -1 * np.sin(ang_yaw)],
                          [zeros, np.sin(ang_yaw), np.cos(ang_yaw)]])
        r_pitch = np.array([[np.cos(ang_pitch), zeros, np.sin(ang_pitch)],
                            [zeros, ones, zeros],
                            [-np.sin(ang_pitch), zeros, np.cos(ang_pitch)]])
        r_roll = np.array([[np.cos(ang_roll), zeros, np.sin(ang_roll)],
                           [zeros, ones, zeros],
                           [-np.sin(ang_roll), zeros, np.cos(ang_roll)]])
        r_c = np.einsum('ijk,ijk -> ijk', r_yaw, r_pitch)
        rotation = np.einsum('ijk,ijk -> ijk', r_c, r_roll)

    elif ax == "custom":
        ux = basis[:, 0]
        uy = basis[:, 1]
        uz = basis[:, 2]
        a = (1 - np.cos(angs))
        rotation = np.array([
            [np.cos(angs) + np.power(ux, 2) * a, ux * uy * a - uz * np.sin(angs), ux * uz * a + uy * np.sin(angs)],
            [ux * uy * a + uz * np.sin(angs), np.cos(angs) + np.power(uy, 2) * a, uy * uz * a - ux * np.sin(angs)],
            [uz * ux * a - uy * np.sin(angs), uz * uy * a + ux * np.sin(angs), np.cos(angs) + np.power(uz, 2) * a]])

        rotation_mats_foo = np.einsum('ijk -> kij', rotation)
        rotated_vecs_foo = np.einsum('ijk,ik->ij', rotation_mats_foo, vecs)

        # print(np.shape(rotation_mats_foo), np.shape(vecs))
        #
        # print(np.shape(rotation))

    rotation_mats = np.einsum('ijk -> kij', rotation)
    rotated_vecs = np.einsum('ijk,ik->ij', rotation_mats, vecs)
    return rotated_vecs


def rotate(vec, ang, ax='x', rpy=[0, 0, 0], basis=None):
    if ax == 'x':
        r_x = np.array([[1, 0, 0],
                        [0, math.cos(ang), -1 * math.sin(ang)],
                        [0, math.sin(ang), math.cos(ang)]])
        return np.matmul(r_x, vec)
    elif ax == 'y':
        r_y = np.array([[math.cos(ang), 0, math.sin(ang)],
                        [0, 1, 0],
                        [-math.sin(ang), 0, math.cos(ang)]])
        return np.matmul(r_y, vec)
    elif ax == 'z':
        r_z = np.array([[math.cos(ang), -math.sin(ang), 0],
                        [math.sin(ang), math.cos(ang), 0],
                        [0, 0, 1]])
        return np.matmul(r_z, vec)
    elif ax == 'c':
        ang_yaw, ang_pitch, ang_roll = rpy[2], rpy[1], rpy[0]
        ang_yaw *= math.pi / 180
        ang_pitch *= math.pi / 180
        ang_roll *= math.pi / 180
        r_yaw = np.array([[1, 0, 0],
                          [0, math.cos(ang_yaw), -1 * math.sin(ang_yaw)],
                          [0, math.sin(ang_yaw), math.cos(ang_yaw)]])
        r_pitch = np.array([[math.cos(ang_pitch), 0, math.sin(ang_pitch)],
                            [0, 1, 0],
                            [-math.sin(ang_pitch), 0, math.cos(ang_pitch)]])
        r_roll = np.array([[math.cos(ang_roll), 0, math.sin(ang_roll)],
                           [0, 1, 0],
                           [-math.sin(ang_roll), 0, math.cos(ang_roll)]])
        r_c = np.matmul(r_yaw, r_pitch)
        r_c = np.matmul(r_c, r_roll)
        r_c = np.matmul(r_c, vec)
        return r_c
    elif ax == "custom":
        ux = basis[0]
        uy = basis[1]
        uz = basis[2]
        a = (1 - math.cos(ang))
        R = np.array([
            [math.cos(ang) + math.pow(ux, 2) * a, ux * uy * a - uz * math.sin(ang), ux * uz * a + uy * math.sin(ang)],
            [ux * uy * a + uz * math.sin(ang), math.cos(ang) + math.pow(uy, 2) * a, uy * uz * a - ux * math.sin(ang)],
            [uz * ux * a - uy * math.sin(ang), uz * uy * a + ux * math.sin(ang), math.cos(ang) + math.pow(uz, 2) * a]])
        return np.matmul(R, vec)


def sphere_intercept(P1, P2, R):
    x1 = P1[0]
    x2 = P2[0]
    y1 = P1[1]
    y2 = P2[1]
    z1 = P1[2]
    z2 = P2[2]

    a = math.pow(x2 - x1, 2) + math.pow(y2 - y1, 2) + math.pow(z2 - z1, 2)
    b = 2 * (x1 * (x2 - x1) + y1 * (y2 - y1) + z1 * (z2 - z1))
    c = math.pow(x1, 2) + math.pow(y1, 2) + math.pow(z1, 2) - math.pow(R, 2)

    determinant = math.pow(b, 2) - 4 * a * c
    if determinant < 0:
        return False
    elif determinant == 0:
        return True
    else:
        return True


def sphere_intercept_np(P1, P2, R):
    x1 = P1[0]
    x2 = P2[:, 0]
    y1 = P1[1]
    y2 = P2[:, 1]
    z1 = P1[2]
    z2 = P2[:, 2]

    a = np.power(x2 - x1, 2) + np.power(y2 - y1, 2) + np.power(z2 - z1, 2)
    b = 2 * (x1 * (x2 - x1) + y1 * (y2 - y1) + z1 * (z2 - z1))
    c = np.power(x1, 2) + np.power(y1, 2) + np.power(z1, 2) - np.power(R, 2)

    determinant = np.power(b, 2) - 4 * a * c

    results = np.full((len(P2)), True, dtype=bool)
    results[determinant < 0] = False
    results[determinant >= 0] = True
    #
    return results
    # if determinant < 0:
    #     return False
    # elif determinant == 0:
    #     return True
    # else:
    #     return True


def geographic_distance_np(target_lat, target_lon, lats, lons, radius, radians=False):
    if not radians:
        target_lat = target_lat * math.pi / 180
        lats = lats * math.pi / 180
        target_lon = target_lon * math.pi / 180
        lons = lons * math.pi / 180

    a = np.power(np.sin((lats - target_lat) / 2), 2) + np.cos(target_lat) * np.cos(lats) * np.power(
        np.sin((lons - target_lon) / 2), 2)

    return radius * 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))


def geographic_distance(lat1, lon1, lat2, lon2, radius, radians=False):
    if not radians:
        lat1 = lat1 * math.pi / 180
        lat2 = lat2 * math.pi / 180
        lon2 = lon2 * math.pi / 180
        lon1 = lon1 * math.pi / 180

    a = math.pow(math.sin((lat2 - lat1) / 2), 2) + math.cos(lat1) * math.cos(lat2) * math.pow(
        math.sin((lon2 - lon1) / 2), 2)

    return radius * 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))


def geographic_area(lat1, lon1, lat2, lon2, radius, radians=False):
    if not radians:
        lat1 = lat1 * math.pi / 180
        lat2 = lat2 * math.pi / 180
        lon2 = lon2 * math.pi / 180
        lon1 = lon1 * math.pi / 180

    area = math.pow(radius, 2) * abs(math.sin(lat1) - math.sin(lat2)) * abs(lon1 - lon2)

    return area


def sat_to_xyz(satellite):
    r = satellite.true_alt

    if satellite.eccentricity > 0:

        # print('ra',satellite.right_ascension_r)
        # print('i',satellite.inclination_r)

        a = satellite.semi_major
        b = a * math.sqrt(1 - math.pow(satellite.eccentricity, 2))
        f = (satellite.altitude + heavenly_body_radius[satellite._focus]) * 10 ** 3
        disp = a - f

        ang = deg_2_rad(satellite.ta + 180)

        x_i, y_i, z_i = disp + a * np.cos(ang), b * np.sin(ang), 0
        coords = np.array([x_i, y_i, z_i]) * 10 ** -3
        coords = rotate(coords, deg_2_rad(satellite.right_ascension), 'z')
        coords = rotate(coords, deg_2_rad(satellite.inclination), 'x')

    else:

        ax1 = np.array([r, 0, 0])
        ax1 = rotate(ax1, satellite.right_ascension_r, 'z')
        ax2 = rotate(ax1, math.pi / 2, 'z')
        ax2 = rotate(ax2, satellite.inclination_r, 'custom',
                     basis=ax1 / math.sqrt(np.sum(ax1 ** 2)))

        basis = np.cross(ax1 / math.sqrt(np.sum(ax1 ** 2)), ax2 / math.sqrt(np.sum(ax2 ** 2)))

        coords = np.array([r, 0, 0])
        coords = rotate(coords, satellite.right_ascension_r, 'z')
        coords = rotate(coords, (satellite.perigee_r + satellite.ta_r), 'custom', basis=basis)
    #
    return coords


def sat_to_xyz_np(satellites):
    r = satellites[:, 1]
    eccentricity = satellites[:, 2]
    semi_major = satellites[:, 8]

    right_ascension = satellites[:, 4] * math.pi / 180
    inclination = satellites[:, 3] * math.pi / 180
    ta = satellites[:, 6] * math.pi / 180
    perigee = satellites[:, 5] * math.pi / 180

    if satellites[0, 2] == 0:
        ax1 = np.column_stack((r, np.zeros(len(r)), np.zeros(len(r))))
        ax1 = rotate_np(ax1, right_ascension, 'z')
        ax2 = rotate_np(ax1, np.full(len(ax1), math.pi / 2), 'z')

        basis = np.sqrt(np.sum(np.square(ax1), axis=1))

        ax2 = rotate_np(ax2, inclination, 'custom', basis=ax1 / basis[:, None])

        b1 = np.sqrt(np.sum(np.square(ax1), axis=1))
        b1 = ax1 / basis[:, None]
        b2 = np.sqrt(np.sum(np.square(ax2), axis=1))
        b2 = ax2 / basis[:, None]

        basis = np.cross(b1, b2)

        coords = np.column_stack((r, np.zeros(len(r)), np.zeros(len(r))))
        coords = rotate_np(coords, right_ascension, 'z')
        coords = rotate_np(coords, (perigee + ta), 'custom', basis=basis)

    else:

        b = semi_major * np.sqrt(1 - np.power(eccentricity, 2))
        f = r * 10 ** 3
        disp = semi_major - f

        ang = ta + math.pi

        coords = np.column_stack((disp + semi_major * np.cos(ang), b * np.sin(ang), np.zeros(len(b)))) * 10 ** -3
        coords = rotate_np(coords, right_ascension, 'z')
        coords = rotate_np(coords, inclination, 'x')

    return coords


def polar2cart(r, phi, theta):
    return [
        r * math.sin(phi) * math.cos(theta),
        r * math.sin(theta) * math.sin(phi),
        r * math.cos(phi)
    ]


def polar2cart_np(vs):
    r = vs[:, 0]
    phi = vs[:, 1]
    theta = vs[:, 2]
    coords = np.column_stack(
        (r * math.sin(phi) * math.cos(theta), r * math.sin(theta) * math.sin(phi), r * math.cos(phi)))

    return coords


def cart2polar(x, y, z):
    v = np.array([x, y, z])
    r = math.sqrt(np.sum(v ** 2))
    azimuth = math.atan2(y, x)
    inclination = math.acos(z / r)

    return r, inclination, azimuth


def cart2polar_np(vs):
    r = np.sqrt(np.sum(vs ** 2, axis=1))
    azimuth = np.arctan2(vs[:, 1], vs[:, 0])
    inclination = np.arccos(vs[:, 2] / r)

    coords = np.column_stack((r, inclination, azimuth))

    return coords


def spherical2geographic(polar, azimuth, radians):
    if radians:
        polar = polar * 180 / math.pi
        azimuth = azimuth * 180 / math.pi

    azimuth = (azimuth + 360) % 360

    if azimuth > 180:
        azimuth -= 360

    if polar > 180:
        polar = 360 - polar

    latitude = 90 - polar
    longitude = azimuth

    return latitude, longitude


def spherical2geographic_np(coordinates, radians):
    if radians:
        coordinates = coordinates * 180 / math.pi

    coordinates = (coordinates + 360) % 360

    coordinates[:, 1][coordinates[:, 1] > 180] -= 360
    coordinates[:, 0][coordinates[:, 0] > 180] = 360 - coordinates[:, 0][coordinates[:, 0] > 180]

    # latitude = 90 - coordinates[:, 0]
    # longitude = coordinates[:, 1]

    return np.column_stack((90 - coordinates[:, 0], coordinates[:, 1]))


def geographic2spherical(latitude, longitude, altitude):
    polar = deg_2_rad(90 - latitude)
    azimuth = deg_2_rad(longitude)
    return polar, azimuth


def rad_2_deg(angle):
    return angle * 180 / math.pi


def deg_2_rad(angle):
    return angle * math.pi / 180
