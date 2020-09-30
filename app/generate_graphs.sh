#!/bin/bash

set -ex

legacy/legacy_rocket > legacy_rocket.out
refurbished/refurbished_rocket > refurbished_rocket.out

gnuplot plot.gp
