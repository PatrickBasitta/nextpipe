#!/bin/bash -ue
FastRemap -f bed -c hg38ToHg19.over.chain -i output.bed -u unmapped.bed -o remapped_to_hg19
