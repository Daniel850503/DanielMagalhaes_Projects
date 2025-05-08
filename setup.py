# -*- coding: utf-8 -*-
"""
Created on Wed May  7 11:23:14 2025

@author: daniel.magalhaes_xel
"""

from setuptools import setup, find_packages

setup(
    name='power_converter_tools',
    version='0.7.0',
    packages=find_packages(include=['Buck_library', 'Buck_library.*',
                                    'Signals_library', 'Signals_library.*'
                                    'Plot_library', 'Plot_library.*']),
    description='tool to support the design of a DC/DC buck converter',
    author='Daniel Magalhães',
    author_email='seu@email.com',
    license='MIT'
)