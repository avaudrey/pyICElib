#!/usr/bin/python
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------
#   Copyright (C) 2017 <Alexandre Vaudrey>                                  |
#                                                                           |
#   This program is free software: you can redistribute it and/or modify    |
#   it under the terms of the GNU General Public License as published by    |
#   the Free Software Foundation, either version 3 of the License, or       |
#   (at your option) any later version.                                     |
#                                                                           |
#   This program is distributed in the hope that it will be useful,         |
#   but WITHOUT ANY WARRANTY; without even the implied warranty of          |
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           |
#   GNU General Public License for more details.                            |
#                                                                           |
#   You should have received a copy of the GNU General Public License       |
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.   |
#---------------------------------------------------------------------------|

import numpy as np
import scipy.optimize as sp

__docformat__ = "restructuredtext en"
__author__ = "Alexandre Vaudrey <alexandre.vaudrey@gmail.com>"
__date__ = "18/09/2017"

class EngineCylinder:
    """ 
    Main geometrical and operating parameters of the concerned engine
    cylinders. As attributes:
    - bore, stroke, connecting rod length (both in [mm]) and compression ratio.
    - angles at which intake and exhaust valves are open and closed.
    - actual crank angle.
    And as methods:
    - displaced, clearance, maximum and actual volumes (both in [l]).
    - states of intake and exhaust valves, i.e. open or closed.
    """
    # TODO: the cross sectional areas of valves will probably be required in the
    # future.
    # The default values of these attributes are coming from an existing
    # 4-strokes and 4-cylinders spark ignited engine from NISSAN.
    def __init__(self):
        # If it's necessary to give a name to the engine geometry 
        self.name = 'cylinder1'
        # Diameter of each cylinder, in [mm].
        self.bore = 75.
        # Stroke, in [mm].
        self.stroke = 90.
        # Compression ratio
        self.compression_ratio = 10.
        # Connecting rod length, in [mm] 
        self._connecting_rod_length = 130.
        # Angles at which the intake valve is open and closed, respectively.
        # Theses values correspond to ideal Beau de Rochas and Diesel cycles.
        self.IVO_angle = -360.
        self.IVC_angle = -180.
        # Angles at which the exhaust valve is open and closed, respectively
        self.EVO_angle = 180.
        self.EVC_angle = 360.
    # Connecting rod length is defined as a property to avoid any mechanical
    # interference with the engine shaft.
    @property
    def connecting_rod_length(self):
        """ Length of the connecting rod, in mm. """
        return self._connecting_rod_length
    @connecting_rod_length.setter
    def connecting_rod_length(self, l):
        """ New length of the connecting rod, in mm. """
        # Check if the entered value is actually greater than the stroke
        if (l <= self.stroke):
            raise ValueError("The connecting rod length cannot be lower than"
                             "the stroke.")
        self._connecting_rod_length = l
        pass
    # And methods
    def displaced_volume(self):
        """ Total displaced volume, or swept volume of the engine, in liter. """
        return 1e-6*0.25*np.pi*pow(self.bore,2)*self.stroke
    def clearance_volume(self):
        """ Total clearance volume of the engine, in liter. """
        return self.displaced_volume()/(self.compression_ratio-1.0)
    def maximum_volume(self):
        """The maximum volume inside the cylinders, equal to the sum of the
        displaced volume and of the clearance one, is used as an horizontal
        limit within the indicated diagram."""
        return self.displaced_volume()+self.clearance_volume()
    def actual_volume(self, crank_angle):
        """Actual volume within which the working gas is enclosed, as a function
        of the crank angle, in °. (0° corresponds to top-center and 180° to
        bottom-center)"""
        # Angle in radians
        radangle = np.radians(crank_angle)
        # Crank radius
        r = 0.5*self.stroke
        # Ratio of the crank radius on the connecting rod length
        ratio = r/self.connecting_rod_length
        # Fraction of the stroke corresponding to the actual volume
        x = r*((1-np.cos(radangle))\
               +(1-np.sqrt(1-pow(ratio*np.sin(radangle),2)))/ratio)
        # Actual volume
        V = self.clearance_volume()+self.displaced_volume()/self.stroke*x
        return V
    def volume_angle_variation(self, crank_angle):
        """ Variation of the actual volume V (in [l]) with the crank angle (in
        [°]) for a given value of the latter.""" 
        # Angle in radians
        radangle = np.radians(crank_angle)
        # Ratio of the crank radius on the connecting rod length
        ratio = 0.5*self.stroke/self.connecting_rod_length
        # Seeked variation
        dVdtheta = 0.5*self.displaced_volume()*(np.sin(radangle)\
                                              +0.5*ratio*np.sin(2*radangle))
        return dVdtheta*np.pi/180.
    # FIXME: these two methods must be done taking into account the periodic
    # motion of the crank
    def is_intake_valve_open(self, crank_angle):
        """ Is the intake valve open at this crank angle. """
        if (crank_angle < self.IVC_angle):
            IV_open = True
        else:
            IV_open = False
        return IV_open
    def is_exhaust_valve_open(self, crank_angle):
        """ Is the intake valve open at this crank angle. """ 
        if (crank_angle < self.EVC_angle):
            EV_open = True
        else:
            EV_open = False
        return EV_open

# === Old version, to cancel once used ========================================
class EngineGeometry:
    """
    Geometric model of a reciprocating internal combustion engine.
    """
    # The default values of these attributes are coming from an existing
    # 4-strokes and 4-cylinders spark ignited engine from NISSAN.
    def __init__(self):
        # If it's necessary to give a name to the engine geometry 
        self.name = 'geometry1'
        # Number of cylinders
        self.number_of_cylinders = 4
        # Diameter of each cylinder, in mm.
        self.bore = 73.6
        # Stroke
        self.stroke = 88.0
        # Compression ratio
        self.compression_ratio = 9.5
        # Ratio of connecting rod length to the crank radius, sometimes a
        # starting parameter
        self._connecting_rod_ratio = 3.0
        # Sometimes, the connecting rod length is known
        self._connecting_rod_length = 132.0
    # Attributes defined as properties
    @property
    def connecting_rod_ratio(self):
        """ Ratio of the connecting rod length to the crank radius,
        dimensionless. """
        return self._connecting_rod_ratio
    @connecting_rod_ratio.setter
    def connecting_rod_ratio(self, R):
        """ New value of the Ratio of the connecting rod length to the crank
        radius, dimensionless. """
        # A wrong value of this parameter, because of a possible confusion with
        # its inverse --- which is as much often used as this one --- would
        # create NaN values in the rest.
        if R<=1.0:
            raise ValueError("The 'connecting rod ratio' is the ratio of the"
                             "connecting rod length to the crank radius, it"
                             "must be greater than one !")
        self._connecting_rod_ratio = R
        # And the related connecting rod length
        self._connecting_rod_length = 0.5*self._connecting_rod_ratio*self.stroke
        pass
    @property
    def connecting_rod_length(self):
        """ Length of the connecting rod, in mm. """
        return self._connecting_rod_length
    @connecting_rod_length.setter
    def connecting_rod_length(self, l):
        """ New length of the connecting rod, in mm. """
        # Check if the entered value is actually greater than the stroke
        if (l <= self.stroke):
            raise ValueError("The connecting rod length cannot be lower than"
                             "the stroke.")
        self._connecting_rod_length = l
        # Calculation of the resulting connecting rod ratio
        self._connecting_rod_ratio = l/(0.5*self.stroke)
        pass
    # And other methods
    def displaced_volume(self):
        """ Total displaced volume, or swept volume of the engine, in liter. """
        return 1e-6*self.number_of_cylinders*0.25*np.pi*pow(self.bore,2)*\
                self.stroke
    def clearance_volume(self):
        """ Total clearance volume of the engine, in liter. """
        return self.displaced_volume()/(self.compression_ratio-1.0)
    def maximum_volume(self):
        """The maximum volume inside the cylinders, equal to the sum of the
        displaced volume and of the clearance one, is used as an horizontal
        limit within the indicated diagram."""
        return self.displaced_volume()+self.clearance_volume()
    def piston_position(self, angle):
        """ Relative position of the piston, =1 at TDC and =0 at BDC, regarding
        to the crank angle in degres. """
        # Angle in radians
        radangle = np.radians(angle)
        # Ratio of the crank radius on the connecting rod length
        ratio = 1/self.connecting_rod_ratio
        return 1-0.5*((1-np.cos(radangle))+\
                      (1-np.sqrt(1-pow(ratio*np.sin(radangle),2)))/ratio)

class EngineCycle(EngineGeometry):
    """Internal cylinder pressure vs. crank angle data series."""
    # TODO : Computation of the mean indicated pressure using the (p,V) datas
    def __init__(self):
        super(EngineCycle, self).__init__()
        # Name of the cycle if required
        self.cycle_name = 'cycle1'
        # List of angles values
        self.angle = [-360,-180,0,180,360]
        # List of corresponding pressure values inside the cylinder
        self.pressure = [1.0,1.0,2.46,1.0,1.0]
        # Ratio of connecting rod length to the crank radius, sometimes a
        # starting parameter
        connecting_rod_ratio =\
                property(EngineGeometry.connecting_rod_ratio.__get__)
        # Sometimes, the connecting rod length is known
        connecting_rod_length =\
                property(EngineGeometry.connecting_rod_length.__get__)
    def piston_position(self):
        """ Relative position of the piston, =1 at TDC and =0 at BDC, regarding
        to the crank angle in degres. """
        # Angle in radians
        radangle = np.radians(self.angle)
        # Ratio of the crank radius on the connecting rod length
        ratio = 1/self.connecting_rod_ratio
        return 1-0.5*((1-np.cos(radangle))+\
                      (1-np.sqrt(1-pow(ratio*np.sin(radangle),2)))/ratio)
    def actual_volume(self):
        """Actual volume in the cylinders, in liters."""
        return self.clearance_volume()+(1-self.piston_position())*\
                self.displaced_volume()

if __name__ == '__main__':
    cylindre1 = EngineCylinder()
    print(cylindre1.__dict__)
    # Technical parameters of the engine
#    cummins = EngineCycle()
#    cummins.name = 'Cummins QSB6.7'
#    cummins.stroke = 124.0
#    cummins.bore = 107.0
#    cummins.compression_ratio = 17.2
#    cummins.number_of_cylinders = 6
#    cummins._connecting_rod_length = 261.5
