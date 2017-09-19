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

if __name__ == '__main__':
    pass
