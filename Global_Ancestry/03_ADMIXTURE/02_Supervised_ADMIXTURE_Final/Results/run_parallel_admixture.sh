#!/bin/bash

module load admixture

parallel --jobs 30 < admixture_supervised_jobs.sh
