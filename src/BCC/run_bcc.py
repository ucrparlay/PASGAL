#!python3
import os

graphs = [
    # Social
    "soc-LiveJournal1_sym.bin",
    "socfb-konect_sym.bin",
    "com-orkut_sym.bin",
    "sinaweibo_sym.bin",
    "twitter_sym.bin",
    "friendster_sym.bin",

    # Web
    "enwiki-2023_sym.bin",
    "eu-2015-host_sym.bin",
    "sd_arc_sym.bin",
    "clueweb_sym.bin",
    "hyperlink2014_sym.bin",
    "hyperlink2012_sym.bin",

    # Road
    "africa_sym.bin",
    "north-america_sym.bin",
    "asia_sym.bin",
    "europe_sym.bin",

    # k-NN
    "CHEM_5_sym.bin",
    "GeoLifeNoScale_5_sym.bin",
    "GeoLifeNoScale_10_sym.bin",
    "Cosmo50_5_sym.bin",

    # Synthetic
    "grid_1000_100000_sym.bin",
    "grid_1000_100000_03_sym.bin",
    "hugetrace-00020_sym.bin",
    "hugebubbles-00020_sym.bin"
]

numactl = "numactl -i all"
graph_dir = "/data0/graphs/links/"


algorithms = ["fast-bcc", "tarjan-vishkin", "hopcroft-tarjan"]
# algorithms = ["hopcroft-tarjan"]

os.system("ulimit -s unlimited")
for algo in algorithms:
    os.system("echo \"### " + algo + "\" >> logs.txt")
    os.system("make " + algo)
    for graph in graphs:
        flags = "-s -i " + graph_dir + graph

        os.system("printf  \"%s\\t\" " + graph + " >> " + algo + ".tsv")
        os.system("drop_caches")

        cmd = numactl + " ./" + algo + " " + flags + " >> logs.txt"
        print(cmd)
        if os.system(cmd):
            exit(1)
