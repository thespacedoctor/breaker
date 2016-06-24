#!/usr/local/bin/python
# encoding: utf-8
"""
*Get common file and folder paths for the host package*

:Author:
    David Young

:Date Created:
    October 29, 2015
"""
import os


def getpackagepath():
    """
    *Get the local path to the pacakge*
    """
    moduleDirectory = os.path.dirname(__file__)
    packagePath = os.path.dirname(__file__) + "/../"

    return packagePath
