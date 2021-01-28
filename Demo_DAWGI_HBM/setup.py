from setuptools import setup

setup(
    name='Demo_DAWGI_HBM',
    version='0.1.6',
    description='Demo of hierarchical Bayesian modelling using Gaia eDR3 data for Pleiades',
    url='https://github.com/Link4138/Demo_DAWGI_HBM',
    author='Jairo Andrés Alzate Trujillo',
    author_email='andresalzate3146@gmail.com',
    license='BSD 3-clause',
    packages=['Pleiades'],
    install_requires=['matplotlib',
                      'numpy',
                      'scipy',
                      'pystan',
                      'pickle-mixin',
                      ],
    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
