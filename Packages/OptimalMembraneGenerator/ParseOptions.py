#!/usr/bin/env python
import configparser

def get_config(filename):
    conf = configparser.SafeConfigParser()
    conf.read(filename)

    # add validation routines here -- for now, we run fast and loose

    return conf
