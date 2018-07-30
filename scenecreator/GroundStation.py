import warnings


class GroundStation(object):
    def __init__(self, name, lat, long, elevation, beam_width):
        self.__name = name
        self.__lat = lat
        self.__long = long
        self.__elevation = elevation
        self.__beam = beam_width

    @property
    def name(self):
        return self.__name

    @name.setter
    def name(self, new_name):
        self.__name = new_name

    @property
    def lat(self):
        return self.__lat

    @lat.setter
    def lat(self, new_lat):
        if (new_lat < -90) or (new_lat > 90):
            return ValueError("Latitude must be between -90 and 90")
        else:
            self.__lat = new_lat

    @property
    def long(self):
        return self.__long

    @long.setter
    def long(self, new_long):
        if (new_long < -180) or (new_long > 180):
            return ValueError("Longitude must be between -180 and 180")
        else:
            self.__long = new_long

    @property
    def elevation(self):
        return self.__elevation

    @elevation.setter
    def elevation(self, new_elev):
        if new_elev < 0:
            return ValueError("Elevation must be above 0")
        if new_elev > 8900:
            return ValueError("Elevation must be on the ground")
        else:
            self.__elevation = new_elev

    @property
    def beam(self):
        return self.__beam

    @beam.setter
    def beam(self, new_beam):
        if (new_beam < 0) or (new_beam > 180):
            return ValueError("Beam width must be between 0 and 180 degrees")
        self.__beam = new_beam

    def as_xml(self):
        warnings.warn("XML support is depreciated and not supported from PIGI 0.8.5 onward", DeprecationWarning)
        return '\t\t<Entity Type="GroundStation" Name="{0}">\n' \
               '\t\t\t<PropertySection Name="UserProperties">\n' \
               '\t\t\t\t<FloatPropertyValue name="Latitude" value="{1}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Longitude" value="{2}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Elevation" value="{3}"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="BeamWidth" value="{4}"/>\n' \
               '\t\t\t\t<StringPropertyValue name="PlanetName" value="Earth"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Animation">\n' \
               '\t\t\t\t<ArrayPropertyValue name="Position" value="[6339.69, -699.193, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Orientation" value="[1, 0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Scale" value="[1, 1, 1]"/>\n' \
               '\t\t\t\t<EnumPropertyValue name="Debug" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Time Input">\n' \
               '\t\t\t\t<TimestampPropertyValue name="Timepoint" value="2016-May-07 08:32:21.059611"/>\n' \
               '\t\t\t\t<DurationPropertyValue name="Duration" value="2"/>\n' \
               '\t\t\t\t<DurationPropertyValue name="StartOffset" value="-1"/>\n' \
               '\t\t\t\t<DurationPropertyValue name="Timestep" value="0.000694444"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Favourite">\n' \
               '\t\t\t\t<EnumPropertyValue name="favourite" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Mesh">\n' \
               '\t\t\t\t<ArrayPropertyValue name="Position" value="[0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Orientation" value="[1, 0, 0, 0]"/>\n' \
               '\t\t\t\t<ArrayPropertyValue name="Scale" value="[500, 500, 500]"/>\n' \
               '\t\t\t\t<StringPropertyValue name="name" value="GroundStation.mesh"/>\n' \
               '\t\t\t\t<EnumPropertyValue name="Debug" value="0"/>\n' \
               '\t\t\t\t<StringPropertyValue name="group" value="SolarSystem"/>\n' \
               '\t\t\t\t<EnumPropertyValue name="visibility" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t\t<PropertySection Name="Look Angles">\n' \
               '\t\t\t\t<ArrayPropertyValue name="LatLon" value="[0, 0]"/>\n' \
               '\t\t\t\t<FloatPropertyValue name="Elevation" value="0"/>\n' \
               '\t\t\t</PropertySection>\n' \
               '\t\t</Entity>\n'.format(self.name, self.lat, self.long, self.elevation, self.beam)

    def __repr__(self):
        return "{0}, {1}, {2}, {3}, {4}".format(self.name, self.lat, self.long, self.elevation, self.beam)

    def __str__(self):
        return "Ground Station: {0}\n" \
               "Latitude: {1}, Longitude: {2}" \
               "Elevation: {3}, Beam Width: [4}".format(self.name, self.lat, self.long, self.elevation, self.beam)

    def as_dict(self):
        return {"Name": self.name,
                "Latitude": self.lat,
                "Longitude": self.long,
                "Elevation": self.elevation,
                "Beam Width": self.beam,
                "Type": 'station'}
