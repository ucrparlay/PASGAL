#!python3
import os

graphs = [
    # Social
    "soc-LiveJournal1.bin",
    "twitter.bin",

    # Web
    "enwiki-2023.bin",
    "eu-2015-host.bin",
    "sd_arc.bin",
    "clueweb.bin",
    "hyperlink2014.bin",
    "hyperlink2012.bin",

    # Road
    "africa.bin",
    "north-america.bin",
    "asia.bin",
    "europe.bin",

    # k-NN
    "CHEM_5.bin",
    "GeoLifeNoScale_5.bin",
    "GeoLifeNoScale_10.bin",
    "Cosmo50_5.bin",

    # Synthetic
    "grid_1000_100000.bin",
    "grid_1000_100000_03.bin",
]

numactl = "numactl -i all"
graph_dir = "/data0/graphs/links/"


algorithms = ["tarjan"]

os.system("ulimit -s unlimited")
for algo in algorithms:
    os.system("echo \"### " + algo + "\" >> logs.txt")
    os.system("make " + algo)
    for graph in graphs:
        flags = "-i " + graph_dir + graph

        os.system("printf  \"%s\\t\" " + graph + " >> " + algo + ".tsv")
        os.system("drop_caches")

        cmd = numactl + " ./" + algo + " " + flags + " >> logs.txt"
        print(cmd)
        if os.system(cmd):
            exit(1)
