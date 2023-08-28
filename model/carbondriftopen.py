import logging; logger = logging.getLogger(__name__)

def open(filename, distribution, plot_distribution = False, times=None, elements=None, load_history=True):
    '''Import netCDF output file as OpenDrift object of correct class'''

    import os
    import pydoc
    from netCDF4 import Dataset
    if not os.path.exists(filename):
        logger.info('File does not exist, trying to retrieve from URL')
        import urllib
        try:
            urllib.urlretrieve(filename, 'opendrift_tmp.nc')
            filename = 'opendrift_tmp.nc'
        except:
            raise ValueError('%s does not exist' % filename)
    n = Dataset(filename)
    try:
        module_name = n.opendrift_module
        class_name = n.opendrift_class
    except:
        logger.warning(filename + ' does not contain global attributes '
                    'opendrift_module and opendrift_class, defaulting to OceanDrift')
        module_name = 'oceandrift'
        class_name = 'OceanDrift'
    n.close()

    if class_name != "CarbonDrift":
        print("The imported dataset is not a CarbonDrift module.")
        from time import sleep
        sleep(3)

        if class_name == 'OpenOil3D':
            class_name = 'OpenOil'
            module_name = 'opendrift.models.openoil'
        if class_name == 'OceanDrift3D':
            class_name = 'OceanDrift'
            module_name = 'opendrift.models.oceandrift'
    
    cls = pydoc.locate(module_name + '.' + class_name)
    
    if cls is None:
        from opendrift.models import oceandrift
        cls = oceandrift.OceanDrift
    
    if class_name == "CarbonDrift":
        o = cls(distribution= distribution, plot_distribution = plot_distribution)
    else:
        o = cls()
    
    o.io_import_file(filename, times=times, elements=elements, load_history=load_history)
    logger.info('Returning ' + str(type(o)) + ' object')
    return o

