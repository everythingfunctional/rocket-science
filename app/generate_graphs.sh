#!/bin/bash

set -ex

legacy/legacy_rocket > legacy_rocket.out

gnuplot plot.gp
