import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name='coso',
    version='1.0.0',    
    url='https://github.com/PietroTotis/CoSo',
    author='Pietro Totis',
    author_email='pietro.totis@kuleuven.be',
    license='GPL',
    packages=['src'],
    install_requires=["portion==2.1.1",
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
        "License :: OSI Approved :: GPL",
    ],
    keywords="combinatorics automated reasoning",
)