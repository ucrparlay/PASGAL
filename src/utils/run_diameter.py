#!python3
import os

graphs = [
    ############# undirected ##############
    # # Social
    # "soc-LiveJournal1_sym.bin",
    # "socfb-konect_sym.bin",
    # "com-orkut_sym.bin",
    # "sinaweibo_sym.bin",
    # "twitter_sym.bin",
    # "friendster_sym.bin",

    # # Web
    # "enwiki-2023_sym.bin",
    # "eu-2015-host_sym.bin",
    # "sd_arc_sym.bin",
    # "clueweb_sym.bin",
    # "hyperlink2014_sym.bin",
    # "hyperlink2012_sym.bin",

    # # Road
    # "africa_sym.bin",
    # "north-america_sym.bin",
    # "asia_sym.bin",
    # "europe_sym.bin",

    # # k-NN
    # "CHEM_5_sym.bin",
    # "GeoLifeNoScale_5_sym.bin",
    # "GeoLifeNoScale_10_sym.bin",
    # "Cosmo50_5_sym.bin",

    # # Synthetic
    # "grid_1000_100000_sym.bin",
    # "grid_1000_100000_03_sym.bin",
    # "hugetrace-00020_sym.bin",
    # "hugebubbles-00020_sym.bin",

    ############### directed ###############
    # # Social
    # "soc-LiveJournal1.bin",
    # "twitter.bin",

    # # Web
    # "enwiki-2023.bin",
    # "eu-2015-host.bin",
    # "sd_arc.bin",
    # "clueweb.bin",
    # "hyperlink2014.bin",
    # "hyperlink2012.bin",

    # # Road
    # "africa.bin",
    # "north-america.bin",
    # "asia.bin",
    # "europe.bin",

    # # k-NN
    # "CHEM_5.bin",
    # "GeoLifeNoScale_5.bin",
    # "GeoLifeNoScale_10.bin",
    "Cosmo50_5.bin",

    # Synthetic
    "grid_1000_100000.bin",
    "grid_1000_100000_03.bin",
]



numactl = "numactl -i all"
graph_dir = "/localdata/xdong038/"


os.system("make get_diameter")
for graph in graphs:
    flags = "-i " + graph_dir + graph
    if "sym" in graph:
        flags += " -s"
    os.system("printf  \"" + graph + "\\t\" >> diameter.tsv")
    cmd = numactl + " ./get_diameter" + " " + flags
    print(cmd)
    if os.system(cmd):
        exit(1)
