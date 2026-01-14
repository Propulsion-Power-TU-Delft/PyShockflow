import numpy as np
from scipy.interpolate import RectBivariateSpline


class Grids:

    def __init__(self, dim0, dim1, variables):

        self.Variables = variables
        self.Dim0 = dim0
        self.Dim1 = dim1
        if not isinstance(variables, list):
            raise Exception('variables should be a list')
        self.PropsValues = []

    def set_values(self, x0_array, x1_array, inputspec, abstractstate, import_node_values, node_values_file_name, kx = 1, ky = 1):
        if import_node_values:
            self.PropsValues = np.load(node_values_file_name, allow_pickle=True)
        else:
            for variable in self.Variables:
                self.PropsValues.append(np.ndarray([self.Dim0, self.Dim1]))

            for i in range(self.Dim0):
                for j in range(self.Dim1):
                    try:
                        abstractstate.update(inputspec, x0_array[i], x1_array[j])
                    except:
                        print(inputspec, x0_array[i], x1_array[j])
                        raise
                    for variable in self.Variables:
                        loc = self.Variables.index(variable)
                        if variable == 'drhomassdPcT':
                            # turbosim compatibility
                            self.PropsValues[loc][i, j] = abstractstate.first_partial_deriv('rhomass', 'P', 'T')
                        else:
                            self.PropsValues[loc][i, j] = getattr(abstractstate, variable)()

        self.MaxX0 = np.max(x0_array)
        self.MinX0 = np.min(x0_array)
        self.MaxX1 = np.max(x1_array)
        self.MinX1 = np.min(x1_array)
        for variable in self.Variables:
            loc = self.Variables.index(variable)
            grid = RectBivariateSpline(x0_array, x1_array, self.PropsValues[loc], kx=kx, ky=ky)
            setattr(self, variable, grid)

    def export_node_values(self, filename):
        np.save(filename, self.PropsValues)
