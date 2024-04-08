from setuptools import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'Algebra Stuff'
LONG_DESCRIPTION = ''

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="algebra_stuff", 
        version=VERSION,
        author="Matthias Schuller",
        author_email="<matthiasshuller92@gmail.com>",
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        packages=find_packages(),
        install_requires=['sympy'], # add any additional packages that 
        
        keywords=['python'],
        # classifiers= [
        #     "Programming Language :: Python :: 3",
        #     "Operating System :: Microsoft :: Windows",
        # ]
)