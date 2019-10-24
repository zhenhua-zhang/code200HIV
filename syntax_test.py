#!/usr/bin/env python3
# -*- UTF-8 -*-

import os
import sys
import subprocess

import numpy as np
import scipy as sp
import pandas as pd
import seaborn as sbn
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch

# phenotype::Phenotypes::remove_outliers
dtfm = pd.DataFrame(
    dict(
        normal_a=np.random.standard_normal(200),
        normal_b=np.random.standard_normal(200),
        normal_c=np.linspace(13, 10, 200) + np.random.standard_normal(200),
        normal_d=np.linspace(15, 20, 200) + np.random.standard_normal(200),
        normal_e=np.linspace(-5, 10, 200) + np.random.standard_normal(200),
        normal_f=np.linspace(10, 20, 200) + np.random.standard_normal(200),
        normal_g=np.linspace(1, 10, 200) + np.random.standard_normal(200),
        normal_h=np.linspace(10, 20, 200) + np.random.standard_normal(200),
    ),
    index=["col_{}".format(x) for x in range(200)]
)

class Animals(object):
    def __init__(self, name, gender):
        self.name = name
        self.gender = gender

    def get_name(self):
        return self.name
    
    def get_gender(self):
        return self.gender
    
class Dog(Animals):
    def __init__(self, name, gender, age):
        super().__init__(name, gender)
        self.age = age
    
    def bark(self):
        return self.get_name() + " is barking!"
    
    def get_age(self):
        return self.age
    
    def summarize(self):
        print("Here is the info of the dog")
        print("Name   : {}".format(self.get_name()))
        print("Age    : {}".format(self.get_age()))
        print("Gender : {}".format(self.get_gender()))

dog = Dog("Jack", "male", 10)
dog.summarize()