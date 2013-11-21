from distutils.core import setup

setup(
    name='stata_dta',
    version='0.1.0',
    author='James Fiedler',
    author_email='jrfiedler@gmail.com',
    packages=['stata_dta', 'stata_dta.stata_missing'],
    license='MIT',
    description='Work with Stata .dta files in Python'
)
