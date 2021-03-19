#!/usr/bin/env bash
bsub.py 10 "backTrans" "snakemake -j --latency-wait 200 --cluster-config cluster_config.json --cluster \"bsub -q {cluster.queue} -n {cluster.nCPUs} -M {cluster.memory} -R {cluster.resources} -J {cluster.name} -o {cluster.output} -e {cluster.error}\""
