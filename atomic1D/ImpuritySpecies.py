class ImpuritySpecies(object):
	# For storing OpenADAS data related to a particular impurity species
	# Loosely based on cfe316/atomic/atomic_data.py/AtomicData class (although with much less code since
	# all of the F77 importing is done in the seperate <<make json_update>> code since BOUT++ protocol 
	# requires fortran code be isolated from main operation)

	def __init__(self,symbol,adas_files_dict={},rate_coefficients={},impurity_fraction=None):
		# Searches for year, atomic_number, has_charge_exchange from user_input.json
		# 
		# Default initialiser for class
		# symbol              : (str)                    | element symbol (e.g. 'C')
		# name                : (str)                    | full name of element (for printing only)
		# year                : (int)                    | year for which OpenADAS data was searched (1996)
		# has_charge_exchange : (bool)                   | whether cx_power (prc) was found for this element-year combination (True)
		# atomic_number       : (int)                    | number of protons for impurity species (6)
		# adas_files_dict     : (str -> str)             | dictionary of OpenADAS files, indexed by file-type ('ionisation': 'scd96_c', ...)
		# rate_coefficients   : (str -> RateCoefficient) | dictionary of RateCoefficient objects corresponding to adas files ('ionisation': <RateCoefficientObject>, ...)
				
		import json

		self = ImpuritySpecies

		with open('user_input.json','r') as fp:
			data_dict = json.load(fp)

		element_dict           = data_dict[symbol]

		assert symbol            == element_dict['symbol']
		self.symbol              = symbol
		self.name                = element_dict['name']
		self.year                = element_dict['year']
		self.has_charge_exchange = element_dict['has_charge_exchange']
		self.atomic_number       = element_dict['atomic_number']
		self.adas_files_dict     = adas_files_dict
		self.rate_coefficients   = rate_coefficients

	def __str__(self):
		# Printing method, for easier inspection of object data
		
		_print_adas_dict = ''
		if len(self.adas_files_dict) == 0:
			_print_adas_check = 'Not initialised'
		else:
			_print_adas_check = 'Initialised'
			for key, value in self.adas_files_dict.items():
				_print_adas_dict = _print_adas_dict + '{:>25} -> {}\n'.format(key,value)
		if len(self.rate_coefficients) == 0:
			_print_rate_check = 'Not initialised'
		else:
			_print_rate_check = 'Initialised'
		
		_printing_string = 'ImpuritySpecies object with attributes'+\
		'\n{:>25} = {}'.format('symbol',			 		self.symbol)+\
		'\n{:>25} = {}'.format('year',			 		self.year)+\
		'\n{:>25} = {}'.format('has_charge_exchange',	self.has_charge_exchange)+\
		'\n{:>25} = {}'.format('atomic_number',	 		self.atomic_number)+\
		'\n{:>25} = {}'.format('adas_files_dict',		_print_adas_check)+\
		'\n{:>25} = {}'.format('rate_coefficients',		_print_rate_check)

		if len(self.adas_files_dict) != 0:
			_printing_string += '\n--------------------------------------------------\n'+_print_adas_dict

		return _printing_string

	def addJSONFiles(self,physics_process,filetype_code,JSON_database_path):
		# 1. Make the filename string expected for the json adas file
		# 2. Check that this file exists in the JSON_database_path/json_data directory
		# 3. Add this file to the atomic data .adas_files_dict attribute
		import os.path

		filename = '{}{}_{}.json'.format(filetype_code,str(self.year)[-2:],self.symbol)
		full_path = '{}/json_data/{}'.format(JSON_database_path,filename)

		if not(os.path.isfile(full_path)):
			raise FileNotFoundError('File {} not found in {}/json_data'.format(filename,JSON_database_path))

		self.adas_files_dict[physics_process] = filename

	def makeRateCoefficients(self,JSON_database_path):
		# Calls the RateCoefficient.__init__ method for each entry in the .adas_files_dict
		# Generates a dictionary of RateCoefficient objects as .rate_coefficients
		from atomic1D import RateCoefficient

		for physics_process, filename in self.adas_files_dict.items():
			full_path = '{}/json_data/{}'.format(JSON_database_path,filename)
			self.rate_coefficients[physics_process] = RateCoefficient(self,full_path)
