#!/bin/bash

grep -i "mux" log_ipi_1 | awk '{print $3 " " $6 " " $9}' > dipole_traj
