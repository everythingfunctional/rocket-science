#!/bin/bash

set -ex

fpm run --target legacy_rocket > legacy_rocket.out
fpm run --target refurbished_rocket > refurbished_rocket.out

gnuplot app/plot.gp
