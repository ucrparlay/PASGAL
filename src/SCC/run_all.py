# This script will run our SCC algorithms on the directed. 
# Test our SCC running time under 96 cores 192 hyperthreads 
# and also the sequential running time with only one core.
# The script also test the baseline algorithms: GBBS(both parallel and sequential)
# iSpan, Multi-Step and Tarjan's SCC algorithm
# The directed output of executing the code will be output to folder '../log/exp1'
import subprocess
import os
import re

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
 

global par_rounds
global seq_rounds
GRAPH_DIR = "/data0/graphs/links"
LOG_DIR= "/data0/lwang323/PASGAL/SCC"
os.makedirs(LOG_DIR, exist_ok=True) 
dir_graphs = {
    # # Social Graphs
    # "LJ": ('soc-LiveJournal1', 'https://www.dropbox.com/scl/fi/836oq1mpruk6y0pul4d4t/soc-LiveJournal1.bin?rlkey=a97bmdoi31a0v87x2j6lmrldg&dl=0'),
    # "TW": ('twitter', 'https://www.dropbox.com/scl/fi/033q9zfjon0jylnevejm5/twitter.bin?rlkey=2j9jcrzp55dm6gg0oypxw0r85&dl=0'),

    # # # Web Graphs
    # "enWiki":("enwiki-2023", ),
    # "EU":("eu-2015-host",),
    # "SD": ('sd_arc', 'https://www.dropbox.com/scl/fi/3lf3eoff3nafyz5cxd4m8/sd_arc.bin?rlkey=m1cdvqe1obnbkd68xksuaaa1u&dl=0'),
    # "CW": ('clueweb', ),
    "HL14": ('hyperlink2014',),
    "HL12": ('hyperlink2012', ),

    # # Road graphs
    "Africa":("africa",),
    "NorthAmerica": ("north-america",),
    "Asia": ("asia",),
    "Europe": ("europe",),

    # KNN Graphs
    # "HH5": ('Household.lines_5', 'https://www.dropbox.com/scl/fi/ndpp5pq2jjlhuqhpx91r7/Household.lines_5.bin?rlkey=n8s0s9wqyzqfxdl82mmqui1z7&dl=0'),
    "CH5": ('CHEM_5', 'https://www.dropbox.com/scl/fi/f0n4o9mrbcfqhw54oa95i/CHEM_5.bin?rlkey=cji03kmi4e5tzmrdj36jsk062&dl=0'),
    # # "GL2": ('GeoLifeNoScale_2', 'https://www.dropbox.com/scl/fi/kckxpp3wcmpbkphikqer4/GeoLifeNoScale_2.bin?rlkey=ybmgdkymbq4i06fvsumk7w9ak&dl=0'),
    "GL5": ('GeoLifeNoScale_5', 'https://www.dropbox.com/scl/fi/x3rrl7imfok5742c8g22v/GeoLifeNoScale_5.bin?rlkey=pltytc4h08oyb93hgf3aho0lo&dl=0'),
    "GL10": ('GeoLifeNoScale_10', 'https://www.dropbox.com/scl/fi/5fw3fhmaqekhihevxvj5e/GeoLifeNoScale_10.bin?rlkey=ed5qn87vdc9m0fjgey7bcopzx&dl=0'),
    # # "GL15": ('GeoLifeNoScale_15', 'https://www.dropbox.com/scl/fi/roc5tryp828wj3sdlzceu/GeoLifeNoScale_15.bin?rlkey=20tfzbv69xsh841sobh5tny05&dl=0'),
    # # "GL20": ('GeoLifeNoScale_20', 'https://www.dropbox.com/scl/fi/2ufpj603of52ygv8ppqg3/GeoLifeNoScale_20.bin?rlkey=qo9csmyg9bryl3kxejwfhiiwr&dl=0'),
    "COS5": ('Cosmo50_5', 'https://www.dropbox.com/scl/fi/z6xpgzyb4vhhw42fy2xqt/Cosmo50_5.bin?rlkey=mcmlrx2ga8f6y5tfhu1wyduqi&dl=0'),

    # # Lattice Graphs
    # "SQR": ('grid_4000_4000', 'https://www.dropbox.com/scl/fi/7jt1isj9oxejjroeailvk/grid_4000_4000.bin?rlkey=boftsrs25u9gngbmg7mfuaao0&dl=0'),
    "REC": ('grid_1000_100000', 'https://www.dropbox.com/scl/fi/wxc3cbdg3i3i6kuqh8ydf/grid_1000_10000.bin?rlkey=ouj3ulurgeomxoruc7g9f6tur&dl=0'),
    # # "SQR_s": ('grid_4000_4000_03', 'https://www.dropbox.com/scl/fi/jsyw9lqsfawpk7642sfcg/grid_4000_4000_03.bin?rlkey=rdjgqc8roqfh3spos858h4aon&dl=0'),
    "REC_s": ('grid_1000_100000_03', 'https://www.dropbox.com/scl/fi/cctseke6v0d5lgugmg89z/grid_1000_10000_03.bin?rlkey=hicjbj4z9y6sj9uogvo9em4s5&dl=0'),
}

# Parallel Running Time
def Our_scc():
  print("Testing Our Parallel SCC Running Time")
  OUT_DIR=f"{LOG_DIR}/OURS"
  os.makedirs(OUT_DIR, exist_ok=True)
  for key, val in dir_graphs.items():
    graph = val[0]
    graph_in = f"{GRAPH_DIR}/{graph}.bin"
    log_out = f"{OUT_DIR}/{key}.txt"
    print(f"Running on {key}")
    # dlarge = "-large" if (key=="HL12") else ""
    cmd = f'''numactl -i all {CURRENT_DIR}/scc -i {graph_in} | tee -a {log_out}'''
    print(cmd)
    subprocess.call(cmd, shell=True)

def collect_data_inner(file_in, key_words):
    f = open(file_in,'r')
    res = f.read()
    f.close()
    data_lines = re.findall(f"{key_words}.*", res)
    data = list(map(lambda x: eval(x.split(" ")[-1]), data_lines))
    return data

def collect_data(folder, key_words):
  print("collect data")
  data=[]
  for key, val in dir_graphs.items():
    graph = val[0]
    log_out = f"{folder}/{key}.txt"
    try:
      data.append(collect_data_inner(log_out, key_words))
    except:
      print(f"run {graph} fail")
  print(data)


if __name__ == '__main__':
  global par_rounds, seq_rounds
  # par_rounds = 5
  # seq_rounds = 3
  Our_scc()
  # collect_data(f"{LOG_DIR}/OURS", "number of scc:")
  # collect_data(f"{LOG_DIR}/OURS", "max scc:")