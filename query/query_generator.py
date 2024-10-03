import random
import argparse
import os

percent = 100

def generator(input_file, output_file, query_num = 1000000):
    fin = open(input_file, "r")
    fout = open(output_file, "w")

    min_start = 1e11
    max_end = 0

    line = fin.readline()
    tmp = line.split()
    upper_num = int(tmp[1])
    lower_num = int(tmp[2])

    line = fin.readline()

    while line:
        tmp = line.split()
        start = int(tmp[2])
        end = int(tmp[3])
        min_start = min(min_start, start)
        max_end = max(max_end, end)
        line = fin.readline()
        
    print(f"min_start: {min_start}, max_end: {max_end}")

    fin.close()

    for i in range(query_num):
        src = random.randint(0, upper_num - 1)
        dst = random.randint(0, lower_num - 1)
        start = random.randint(min_start, max_end)
        end = random.randint(start, max_end)
        fout.write(str(src) + " " + str(dst) + " " + str(start) + " " + str(end) + "\n")
        
    fout.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=str, default="../data/wikiquote_pl.txt")
    parser.add_argument("--query_num", type=int, default=1000000)
    args = parser.parse_args()
    
    filename = os.path.basename(args.input).split(".")[0]  
    output_file = filename + f".query"

    generator(args.input, output_file, args.query_num)