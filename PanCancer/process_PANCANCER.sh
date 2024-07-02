#!/bin/bash
# execution of different scripts to process PANCANCER
./IDcatcher_filetransfer.sh > /PAT-Sequenzer/PANCANCER/ARCHIVE_AND_LOGS/RENAME_LOGS/$(date +"%Y%m%d%H%M".log) && \
./RUN_PAN.sh && \
./archive_transfer.sh > /PAT-Sequenzer/PANCANCER/ARCHIVE_AND_LOGS/ARCHIVE_LOGS/$(date +"%Y%m%d%H%M".log)
