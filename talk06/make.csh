#!/bin/bash
g++ -o execute_name `root-config --cflags --glibs` loopPlot.C
g++ -o execute_name16 `root-config --cflags --glibs` loopPlot16.C
g++ -o execute_name_all `root-config --cflags --glibs` loopPlot_all.C
