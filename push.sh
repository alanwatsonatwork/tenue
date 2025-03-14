#!/bin/sh

rsync -ahv --exclude=".*" --exclude="__pycache__" . oan-data:tenue/
