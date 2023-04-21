import os
from setuptools import setup, find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='coso',
    version='1.0.3',    
    url='https://github.com/PietroTotis/CoSo',
    author='Pietro Totis',
    author_email='pietro.totis@kuleuven.be',
    license='GNU General Public License (GPL)',

    # packages=find_packages(),
    packages=['coso'],
    package_dir={'coso':'src'},
    package_data={"src": ["VisCoSo/*"]},
    install_requires=["portion",
                      "clingo==5.5.1",
                      "psutil",
                      "yattag",
                      "ply"],    
    description = ("A solver for combinatorics math word problems."),
    long_description=read('README.md'),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Artificial Intelligence",
        "License :: OSI Approved :: GNU General Public License (GPL)",
    ],
    keywords="combinatorics automated reasoning",
)
