#!/bin/sh
cd v2
python3 project_v2.py True 1 50
python3 project_v2.py False 1 50
cd ../v1
python3 project_v1.py