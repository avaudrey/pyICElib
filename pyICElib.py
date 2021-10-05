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

class SimpleChemicalReaction:
    """
    Chemical reaction transforming a set of reactants in a set of products. Each
    set is represented as a dictionary of type 'set = {'Compound1':(nu1, n1(0)),
    ...}' with:
        - 'Compound1' the chemical formula of compound 1, as e.g.  'CH4', 'H2',
        'O2', 'H2Ov' (for water vapour), 'H2Ol' (for liquid water), etc.
        - nu1 the (always) positive stoichiometric coefficient of compound 1.
        - n1(0)>0 the initial molar amount of compound 1, in mol.
    """
    # Attributes
    def __init__(self):
        # If it's necessary to give a name to the engine geometry 
        self.name = 'reaction1'
        # Set of reactants, here for example hydrogen and oxygen in
        # stoichiometric proportions
        self.reactants = {'H2':(1.0,1.0), 'O2':(0.5,0.5)}
        # Set of products, here only water vapour
        self.products = {'H2Ov':(1.0,0.0)}
    # Methods
    def maximum_extent_of_reaction(self):
        """
        Calculate the maximum extent of reaction that can be reached according
        to the initial amounts of reactants available.
        """
        # Calculation of the maximum reaction progress corresponding to each
        # reactant
        list_of_xi_max = []
        for values in self.reactants.values():
            if (values[0] != 0):
                list_of_xi_max.append(values[1]/values[0])
        # And return of the minimum
        return min(list_of_xi_max)
    def molar_amounts(self, relative_xi):
        """
        Calculation of the molar amounts of chemical compounds at an
        intermediate extent of reaction represented by 'relative_xi' with:
            relative_xi = extent/maximum_extent
        """ 
        if (relative_xi < 0) or (relative_xi > 1.0):
            raise ValueError('Relative extent of reaction is from 0% to 100%')
        else:
            # Maximum extent of reaction
            xi_max = self.maximum_extent_of_reaction()
            # Dictionary containing the final molar amounts of reactants and
            # products
            amounts = {}
            # We start with the calculation of final amounts of reactants
            for react, values in zip(self.reactants.keys(),\
                                     self.reactants.values()):
                amounts[react] = values[1]-relative_xi*xi_max*values[0]
                # And for the products
            for prod, values in zip(self.products.keys(),\
                                    self.products.values()):
                amounts[prod] = values[1]+relative_xi*xi_max*values[0]
            return amounts
    def final_molar_amounts(self):
        """
        Final molar amounts of reactants and products.
        """
        # The final composition is just the one of reaction but with a 100%
        # relative extent of reaction.
        return self.molar_amounts(1.0)
    def molar_concentrations(self, relative_xi):
        """ 
        Molar concentration of each compound in the reactive mixture at an
        intermediate extent of reaction represented by 'relative_xi' with:
            relative_xi = extent/maximum_extent
        """ 
        if (relative_xi < 0) or (relative_xi > 1.0):
            raise ValueError('Relative extent of reaction is from 0% to 100%')
        else:
            # Molar quantities
            n = self.molar_amounts(relative_xi)
            # Total amount of moles at the end of the reaction
            ntot = sum(list(n.values()))
            # Final concentrations
            x = {}
            for elem, ni in zip(n.keys(), n.values()):
                x[elem] = ni/ntot
            return x
    def final_molar_concentrations(self):
        """
        Molar concentration of each components in the final mixture.
        """
        return self.molar_concentrations(1.0)

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
        # Angles at which the intake valve opening and closing process start,
        # respectively. Theses values correspond to ideal Beau de Rochas and
        # Diesel cycles.
        self.IVO_angle = -360.
        self.IVC_angle = -180.
        # Angles at which the exhaust valve starts to open and closed,
        # respectively.
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

class NISTShomateEquation:
    """
    Shomate equation used to manipulate temperature depending thermochemical
    properties of various compounds, see : https://webbook.nist.gov/. Concerned
    properties are:
        - Specific molar heat capacity at constant pressure, in kJ/(mol.K).
        - Standard molar enthalpy, in kJ/mol.
        - Standard molar entropy at standard pressure of 1 bar, in kJ/(mol.K).
    """
    def __init__(self, compound):
        # We load the database of NIST coefficients
        self.__NIST_datas = np.load("NIST_data.npy", allow_pickle=True).all()
        # Specific compound under concern
        self._name = compound
        # NIST coefficients related to the compound
        self._coefficients = self.__NIST_datas[self._name]
        # Temperature at which the physical properties are considered
        self._temperature = [298.15]
        # Universal constant of ideal gases
        self.__R = 8.314462618
    # Attributes as properties
    @property
    def name(self):
        """ Specific compound under concern. """
        return self._name
    @name.setter
    def name(self, compound):
        """ Change the compound under concern. """
        # Check if the chose compound is in the list
        clist = self.compounds_list()
        if (compound not in clist):
            raise ValueError("The compound you choose is not in the current list!")
        # And import the NIST coefficients related to this new compound
        self._coefficients = self.__NIST_datas[compound]
    @property
    def temperature(self):
        """ Temperature at which the asked properties are considered. """
        return self._temperature
    @temperature.setter
    def temperature(self, newT):
        """ Change the temperature value. """
        if (type(newT) == float):
            self._temperature = [newT]
        elif (type(newT) == np.ndarray) or (type(newT) == tuple):
            self._temperature = list(newT)
        elif (type(newT) == list):
            self._temperature = newT
        else:
            raise ValueError("Temperature must be a numerical value, as float,\
                             list, tuple or array")
    # Methods
    def compounds_list(self):
        """ List the chemical compounds that can be manipulated. """
        return list(self.__NIST_datas.keys())
    def molar_heat_capacity_at_cst_pressure(self):
        """
        Specific molar heat capacity at constant pressure Cp, in J/(mol.K).
        """
        specific_heat_capacity = []
        # For each temperature value
        for temp in self._temperature:
            # We look for the set of coefficient corresponding to the range of
            # temperature that concerns us
            for rangeT in self._coefficients:
                if (rangeT[0][0] <= temp <= rangeT[0][1]):
                    # And we extract the coefficients for the calculation of the
                    # cp
                    coeffs = rangeT[1:]
                    # Reduced temperature
                    t = 1e-3*temp
                    specific_heat_capacity.append(coeffs[0]+coeffs[1]*t\
                                                  +coeffs[2]*pow(t,2)\
                                                  +coeffs[3]*pow(t,3)\
                                                  +coeffs[4]/pow(t,2))
        return specific_heat_capacity
    def molar_heat_capacity_at_cst_volume(self):
        """
        Specific molar heat capacity at constant volume CV, in J/(mol.K).
        """
        # TODO : this formula is supposed to work only for ideal gas, how to do
        # otherwise?
        # Use of the Mayer's relation
        cpressure = self.molar_heat_capacity_at_cst_pressure()
        return list(np.array(cpressure)-self.__R)
    def heat_capacity_ratio(self):
        """
        Heat capacity ratio 'gamma'.
        """
        cpressure = self.molar_heat_capacity_at_cst_pressure()
        cvolume = self.molar_heat_capacity_at_cst_volume()
        return list(np.array(cpressure)/np.array(cvolume))
    def molar_enthalpy(self):
        """ Molar enthalpy, in kJ/mol."""
        enthalpy = []
        # For each temperature value
        for temp in self._temperature:
            # We look for the set of coefficient corresponding to the range of
            # temperature that concerns us
            for rangeT in self._coefficients:
                if (rangeT[0][0] <= temp <= rangeT[0][1]):
                    # And we extract the coefficients for the calculation of the
                    # cp
                    coeffs = rangeT[1:]
                    # Reduced temperature
                    t = 1e-3*temp
                    enthalpy.append(coeffs[0]*t+0.5*coeffs[1]*pow(t,2)\
                                    +0.333*coeffs[2]*pow(t,3)\
                                    +0.25*coeffs[3]*pow(t,4)\
                                    -coeffs[4]/t+coeffs[5]+coeffs[7])
        return enthalpy
    def molar_entropy(self):
        pass

if __name__ == '__main__':
    pass
#    cylindre1 = EngineCylinder()
    # Technical parameters of the engine
#    cummins = EngineCycle()
#    cummins.name = 'Cummins QSB6.7'
#    cummins.stroke = 124.0
#    cummins.bore = 107.0
#    cummins.compression_ratio = 17.2
#    cummins.number_of_cylinders = 6
#    cummins._connecting_rod_length = 261.5
