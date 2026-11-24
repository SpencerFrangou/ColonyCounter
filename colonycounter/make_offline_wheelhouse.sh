#!/bin/bash
# Run once on a computer with internet access to fill colonycounter/wheels with every package ColonyCounter
# needs, so that the launcher can install them on a computer without internet (same operating system and
# Python version). ensure_env.py uses the cache automatically when the folder exists.
cd "$(dirname "$0")"; mkdir -p wheels
python3 -m pip download -r ../requirements.txt -d wheels && echo "wheel cache written to colonycounter/wheels ($(ls wheels | wc -l) files)"
