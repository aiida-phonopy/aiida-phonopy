# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

import json

if __name__ == "__main__":
    with open("setup.json", "r") as info:
        kwargs = json.load(info)

    setup(reentry_register=True, include_package_data=True, packages=find_packages(exclude=["aiida"]), **kwargs)
