#!/bin/sh
cd v2
python3 project_v2.py True 0.2 100
python3 project_v2.py False 0.2 100
cd ../v1
python3 project_v1.py
