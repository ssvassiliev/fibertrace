#!/bin/bash

../tfield < ../tfield.conf
../streamline < ../streamline.conf
../load_streamline_mol2.sh -f streamline_XYZ.mol2
