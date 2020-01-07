# -*- coding: utf-8 -*-

import arcpy
import datetime
import numpy as np
from typing import NamedTuple

class Bessel1841(NamedTuple):
    # parameters from https://www.geoportal.sk/files/gz/etrs89_s-jtsk_tech_sprava_2014_ver3_0.pdf
    # or https://www.geoportal.sk/sk/geodeticke-zaklady/geodeticke-systemy-transformacie/
    a = 6377397.155 # semi-major elipsoid axis
    f = 1/(299.1528154) # reciprocal flattening
    b = a*(1 - f) # secondary elipsoid axis
    e = np.sqrt((a**2 - b**2) / a**2) # first excentricity

class Grs80(NamedTuple):
    # parameters from https://en.wikipedia.org/wiki/Geodetic_datum#From_ECEF_to_geodetic
    # or 
    a = 6378137.000
    f = 1/(298.257222101)
    b = a*(1 - f)
    e = np.sqrt((a**2 - b**2) / a**2)

class Wgs84(NamedTuple):
    # parameters from https://en.wikipedia.org/wiki/World_Geodetic_System
    a = 6378137.000
    f = 1/(298.257223563)
    b = a*(1 - f)
    e = np.sqrt((a**2 - b**2) / a**2)

class Transform(object):

    @staticmethod
    def plh_xyz(plh, elipsoid):
        """ See https://en.wikipedia.org/wiki/Reference_ellipsoid for reference.
            Elipsoid parameter defines the elipsoid with 2 properties - a and f (semi-major axis and reciprocal flattening)."""
        phi, lam, hei = plh
        N = elipsoid.a/np.sqrt(1 - elipsoid.e**2 * (np.sin(phi)**2))
        X = (N+hei)*np.cos(phi)*np.cos(lam)
        Y = (N+hei)*np.cos(phi)*np.sin(lam)
        Z = (N*(1-elipsoid.e**2)+hei)*np.sin(phi)
        return X, Y, Z

    @staticmethod
    def quadrant_lon(lam, x, y, timeoffset):
        # returns -pi; pi in the correct quadrant
        a = [0, 3/2*np.pi, 1/2*np.pi, np.pi]
        i = 2*(x < 0) + (y < 0)
        while lam < a[i]:
            lam += 1/2*np.pi
        lam += 2*np.pi/86400 * timeoffset
        if lam > np.pi:
            lam = -2*np.pi + lam
        return lam

    @staticmethod
    def xyz_plh(xyz, elipsoid, timeoffset=0):
        X, Y, Z = xyz
        lam = np.arctan(Y/X)          # longitude
        lam = Transform.quadrant_lon(lam, X, Y, timeoffset)
        phi = np.arctan(Z/(np.sqrt(X**2+Y**2)*(1-elipsoid.e**2)))   # latitude
        diff = 1
        while diff > 10E-12:
            N = elipsoid.a / np.sqrt(1 - elipsoid.e**2 * (np.sin(phi)**2))
            h = np.sqrt(X**2 + Y**2)/np.cos(phi) - N
            phi1 = np.arctan(Z/np.sqrt(X**2+Y**2) * ((N+h)/(N+h - elipsoid.e**2*N)))
            diff = abs(phi-phi1)
            phi = phi1

        return phi, lam, h

class Project(object):

    @staticmethod
    def igs_wgs(xyz, timeoffset=0):
        plh = Transform.xyz_plh(xyz, Grs80, timeoffset)
        plh = np.append(np.degrees(plh[:2]), plh[2])
        return plh

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "SpaceVehicle"
        self.alias = "Projects the IGS/ITRF coordinates to WGS84"

        # List of tool classes associated with this toolbox
        self.tools = [ConvertSP3Tool]


class ConvertSP3Tool(object):
    def __init__(self):
        """Define the tool (tool name is the name of the class)."""
        self.label = "Convert SP3 file"
        self.description = "Converts SP3 file with satellite coordinates to new feature class. SP3 is obtained from ftp://cddis.nasa.gov/gnss/products/"
        self.canRunInBackground = True

    def getParameterInfo(self):
        """Define parameter definitions"""

        param0 = arcpy.Parameter(
            displayName="IGS SP3 file",
            name="sp3_file",
            datatype="DEFile",
            parameterType="Required",
            direction="Input")

        param1 = arcpy.Parameter(
            displayName="Output workspace",
            name="output_workspace",
            datatype="DEWorkspace",
            parameterType="Required",
            direction="Input")

        param2 = arcpy.Parameter(
            displayName="Output Featureclass Name",
            name="output_fc_name",
            datatype="GPString",
            parameterType="Required",
            direction="Input")

        param3 = arcpy.Parameter(
            displayName="Truncate if exists",
            name="truncate",
            datatype="GPBoolean",
            parameterType="Optional",
            direction="Input")            
        param3.value = True

        return [param0, param1, param2, param3]

    def isLicensed(self):
        """Set whether tool is licensed to execute."""
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return

    def create_featureclass(self, workspace, fc_name):
        fc = arcpy.management.CreateFeatureclass(workspace, fc_name, "POINT", None, "DISABLED", "ENABLED", "GEOGCS['GCS_WGS_1984',DATUM['D_WGS_1984',SPHEROID['WGS_1984',6378137.0,298.257223563]],PRIMEM['Greenwich',0.0],UNIT['Degree',0.0174532925199433]];-400 -400 1000000000;-100000 10000;-100000 10000;8.98315284119521E-09;0.001;0.001;IsHighPrecision", '', 0, 0, 0, '')
        arcpy.management.AddField(fc, "Name", "TEXT", 38, 8, 10, "Name", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "Time", "DATE", 38, 8, None, "Time", "NON_NULLABLE", "NON_REQUIRED", '')        
        arcpy.management.AddField(fc, "X", "DOUBLE", 38, 8, None, "IGS X [km]", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "Y", "DOUBLE", 38, 8, None, "IGS Y [km]", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "Z", "DOUBLE", 38, 8, None, "IGS Z [km]", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "T", "DOUBLE", 38, 8, None, "T [us]", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "a", "SHORT", 3, 8, 10, "a", "NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "b", "SHORT", 3, 8, 10, "b", "NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "c", "SHORT", 3, 8, 10, "c", "NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "d", "SHORT", 3, 8, 10, "d", "NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "lon", "DOUBLE", 38, 8, None, "lon", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "lat", "DOUBLE", 38, 8, None, "lat", "NON_NULLABLE", "NON_REQUIRED", '')
        arcpy.management.AddField(fc, "alt", "DOUBLE", 38, 8, None, "alt [m]", "NON_NULLABLE", "NON_REQUIRED", '')

        return fc

    def _splitline(self, line, delimiters):
        idx = np.cumsum([0] + list(delimiters))
        slices = [slice(i, j) for (i, j) in zip(idx[:-1], idx[1:])]
        data = [line[s] for s in slices]
        return data

    def _intornone(self, value):
        return int(value) if value.strip() else None

    def readsp3(self, sp3_file):
        timestamp = None
        with open(sp3_file) as f:
            for line in f:
                if line.startswith('*'):
                    # timestamp line
                    delimiter = (2, 5, 3, 3, 3, 3, 3)
                    data = self._splitline(line, delimiter)
                    ts = [int(d) for d in data[1:]]
                    timestamp = datetime.datetime(ts[0], ts[1], ts[2], ts[3], ts[4], ts[5])
                elif line.startswith('EOF'):
                    return
                elif timestamp:
                    delimiters = (4, 14, 14, 14, 14, 3, 3, 3, 4)
                    data = self._splitline(line, delimiters)            
                    row = [timestamp, data[0], float(data[1]), float(data[2]), float(data[3]), float(data[4]), 
                        self._intornone(data[5]), self._intornone(data[6]), self._intornone(data[7]), self._intornone(data[8])]
                    yield row

    def append(self, sp3_file, featureclass):
        with arcpy.da.InsertCursor(featureclass, ['Shape@XYZ', 'Time', 'Name', 'X', 'Y', 'Z', 'T', 'a', 'b', 'c', 'd', 'lon', 'lat', 'alt']) as icur:
            for row in self.readsp3(sp3_file):
                plh = Project.igs_wgs((row[2]*1000, row[3]*1000, row[4]*1000), row[0].hour*3600 + row[0].minute*60 + row[0].second)
                data = list(row)
                data.insert(0, (plh[1], plh[0], plh[2]))
                data.append(plh[1])
                data.append(plh[0])
                data.append(plh[2])
                icur.insertRow(tuple(data))
            
    def execute(self, parameters, messages):
        sp3 = parameters[0].valueAsText
        ws = parameters[1].valueAsText
        fc_name = parameters[2].valueAsText
        truncate = parameters[3].valueAsText
        fc = ws + '/' + fc_name
        if truncate == 'true':
            arcpy.AddMessage("Truncating " + fc)
            arcpy.management.TruncateTable(fc)
        else:
            fc = self.create_featureclass(ws, fc_name)
        self.append(sp3, fc)
