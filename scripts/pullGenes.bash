#!/usr/bin/env bash

grep -P "\tgene\t" /zfs/tillers/Reference_Genomes/BTx623/v3.1.1/annotation/Sbicolor_454_v3.1.1.gene.gff3 | grep -f GenesAndNeighbors_700_Sobic.txt | awk '{print $1":"$4"-"$5}'
