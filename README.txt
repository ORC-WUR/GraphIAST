GraphIAST is a user-friendly graphical user interface to predict mixed gas adsorption isotherms by Ideal Adsorption Solution Theory (IAST). GraphIAST uses pyIAST as programming basis and makes three-dimentional selectivity predictions easily accessible within just a few clicks. 
GraphIAST can be used on Windows, and linux.
The respective folders contain the executable and an image with a short help guide. Keep these files in the same folder while runnin the executable or script for optimal usage.
If the files are not in the same folder the "help" function in the software will not be operational.

When using the script make sure that the following python packages are installed: Math, Matplotlib, Numpy, Pandas, PIL, pyIAST, Tkinter (checked with python 3.9)

This package includes:
- A folder "script", which includes the python files and code of all versions so far.
        - the executables for the different platforms
        - the python script

- GraphIAST manual

- MIT-Licence

- A folder "Case studies" contains various cases e.g. selectivity vs mole fraction by constant pressure using test input data.
	-the folder "test data" contains the following test data:
		IRMOF-1_methane_isotherm_298K.csv (isotherm data obtained from the pyIAST package)
		IRMOF-1_ethane_isotherm_298K.csv (isotherm data obtained from the pyIAST package)
	-the folder "test run" contains a small tutorial with expected output:

