import configparser
import os

class Config:
    def __init__(self, config_file='input.ini'):
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"Config file '{config_file}' not found.")
        
        self.config_parser = configparser.ConfigParser()
        self.config_parser.read(config_file)
        
    def get_config_value(self, section, option, default=None):
        """
        Helper method to retrieve a configuration value with a default fallback.
        """
        try:
            return self.config_parser.get(section, option)
        except (configparser.NoSectionError, configparser.NoOptionError):
            return default
    
    def getNumberOfPoints(self):
        return int(self.config_parser.get('SIMULATION', 'NUMBER_POINTS')) 
    
    def getLength(self):
        return float(self.config_parser.get('GEOMETRY', 'LENGTH')) 
    
    def getPressureLeft(self):
        return float(self.config_parser.get('SIMULATION', 'PRESSURE_LEFT')) 
    
    def getPressureRight(self):
        return float(self.config_parser.get('SIMULATION', 'PRESSURE_RIGHT')) 
    
    def getDensityLeft(self):
        return float(self.config_parser.get('SIMULATION', 'DENSITY_LEFT')) 
    
    def getDensityRight(self):
        return float(self.config_parser.get('SIMULATION', 'DENSITY_RIGHT')) 
    
    def getTemperatureLeft(self):
        return float(self.config_parser.get('SIMULATION', 'TEMPERATURE_LEFT')) 
    
    def getTemperatureRight(self):
        return float(self.config_parser.get('SIMULATION', 'TEMPERATURE_RIGHT')) 
    
    def getVelocityLeft(self):
        return float(self.config_parser.get('SIMULATION', 'VELOCITY_LEFT')) 
    
    def getVelocityRight(self):
        return float(self.config_parser.get('SIMULATION', 'VELOCITY_RIGHT')) 
    
    def getCFLMax(self):
        return float(self.config_parser.get('SIMULATION', 'CFL_MAX')) 
    
    def getTimeMax(self):
        return float(self.config_parser.get('SIMULATION', 'TIME_MAX')) 
    
    def getFluidName(self):
        return str(self.config_parser.get('FLUID', 'FLUID_NAME'))
    
    def getFluidModel(self):
        return str(self.config_parser.get('FLUID', 'FLUID_MODEL')).lower() 
    
    def getFluidGamma(self):
        return float(self.config_parser.get('FLUID', 'FLUID_GAMMA'))
    
    def getGasRConstant(self):
        return float(self.config_parser.get('FLUID', 'GAS_R_CONSTANT'))
    
    def getInterfaceLocation(self):
        return float(self.config_parser.get('GEOMETRY', 'INTERFACE_LOCATION'))
    
    def getBoundaryConditions(self):
        left = str(self.config_parser.get('SIMULATION', 'BOUNDARY_CONDITION_LEFT')).lower() 
        right = str(self.config_parser.get('SIMULATION', 'BOUNDARY_CONDITION_RIGHT')).lower() 
        return (left, right)
    
    def getNumericalScheme(self):
        return str(self.config_parser.get('SIMULATION', 'NUMERICAL_SCHEME')).lower() 
    
    # def getTimeStepMethod(self):
    #     try:
    #         return str(self.config_parser.get('SIMULATION', 'TIME_STEP_METHOD')).lower()
    #     except:
    #         return 'constant'
    
    def getOutputFolder(self):
        return str(self.config_parser.get('OUTPUT', 'FOLDER_NAME')) 
    
    def getOutputFileName(self):
        return str(self.config_parser.get('OUTPUT', 'FILE_NAME')) 
    
    def showAnimation(self):
        res = str(self.config_parser.get('OUTPUT', 'SHOW_ANIMATION')).lower() 
        if res=='yes' or res=='true':
            return True
        else:
            return False
    
    def getMUSCLReconstruction(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'MUSCL_RECONSTRUCTION')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False # false by default
    
    
    def isMeshRefined(self):
        try:
            res = str(self.config_parser.get('SIMULATION', 'MESH_REFINEMENT')).lower() 
            if res=='yes' or res=='true':
                return True
            else:
                return False
        except:
            return False # false by default
    
    def getRefinementBoundaries(self):
        start = float(self.config_parser.get('SIMULATION', 'X_START_REFINEMENT')) 
        end = float(self.config_parser.get('SIMULATION', 'X_END_REFINEMENT')) 
        return (start, end)
    
    def getNumberPointsRefinement(self):
        return int(self.config_parser.get('SIMULATION', 'NUMBER_POINTS_REFINEMENT')) 
        
    
    def getTopology(self):
        try:
            return str(self.config_parser.get('GEOMETRY', 'TOPOLOGY')).lower()
        except:
            return 'default' 
    
    def getNozzleFilePath(self):
        return str(self.config_parser.get('GEOMETRY', 'NOZZLE_FILEPATH'))
    
    
    def getAreaReference(self):
        try:
            return float(self.config_parser.get('GEOMETRY', 'REFERENCE_AREA')) 
        except:
            return 1.0 # default
    
    