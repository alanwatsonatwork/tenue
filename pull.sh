#!/bin/sh

rsync -ahv --exclude=".*" --exclude="__pycache__" colibri-data:tenue/src/ . 
