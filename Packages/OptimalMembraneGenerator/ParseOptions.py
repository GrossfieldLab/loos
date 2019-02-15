#!/usr/bin/env python3
import configparser

def get_config(filename):
    conf = configparser.SafeConfigParser()
    conf.read(filename)

    # add validation routines here -- for now, we run fast and loose

    return conf
