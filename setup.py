#!/usr/bin/env python
#-*- coding:utf-8 -*-

#############################################
# File Name: setup.py
# Author: Songyan Zhu
# Mail: zhusy93@gmail.com
# Created Time:  2025-04-18
#############################################


from setuptools import setup, find_packages

setup(
	name = "bigleaf",
	version = "0.0.1",
	keywords = ("Geosciences, geospatial, data science"),
	description = "R bigleaf module in Python",
	long_description = "Processing academic data and drawing figures for geosciences researchers.",
	license = "MIT Licence",

	url="https://github.com/soonyenju/bigleaf",
	author = "Songyan Zhu",
	author_email = "zhusy93@gmail.com",

	packages = find_packages(),
	include_package_data = True,
	platforms = "any",
	install_requires=[

	]
)