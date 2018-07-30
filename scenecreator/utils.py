"""

"""


def mod(x, y):
    """
    Mappable modulo function
    :param x: First number
    :param y: Second number
    :return: x % y
    """
    return [a % b for a, b in zip(x, y)]


heavenly_body_radius = {
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
