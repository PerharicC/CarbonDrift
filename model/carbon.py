from opendrift.models.oceandrift import Lagrangian3DArray
import numpy as np

class Carbon(Lagrangian3DArray):

    variables = Lagrangian3DArray.add_variables([
    ('mass', {'dtype': np.float32,
                      'units': 'g',
                      'default': 1})])