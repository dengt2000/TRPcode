# A Contraction Framework for Answering Reachability Queries on Temporal Bipartite Graphs

## 1. Data preparation
The initial data source is the konect project, which provides a collection of network datasets. The datasets were available at [konect.cc](http://konect.cc/). However, early this year, the website was down, and the datasets were not available. We have downloaded the datasets before the website went down. 
THerefore, we uploaded the used datasets to the following link: [TRP_datasets](https://www.kaggle.com/datasets/wilfriedsun/datasets-for-trp)
> Note: The initial datasets contains only the endpoints and a single timestamp for each edge. We have extended the datasets by adding another timestamp for each edge as the end time. The end time is the start time plus a random number  following the  power low distribution. 

## 2. Running the program
### 2.1 Clone the repository
Clone the repository and run the following command to build the execution file:
```bash
madir -p build
cmake ..&& make -j
```
### 2.2 Create directories
The directory structure should be as follows:
```
|--TRPcode
   |--build
   |--data
   |--query
   |--results
   |--logs
   |--include
   |--src
   |--third
   |--tools
```
move the date file to the date directory under TRPcode repository.

### 2.3 Prepare the query file
To run the program, you need to prepare a query file. A simple query file generator is provided in the `./query` directory. Run the following command to generate a query file:
```bash
python3 query_generator.py --input ../data/wikiquote_pl.txt  --query_num 1000000
```

### 2.4 Run the program
After building the execution file, you can run the following command to execute the program in the build directory:
```bash
./Computation --inputPath=../data/wikiquote_pl.txt --queryFilePath=../query/wikiquote_pl.query --minBlock=64 --algorithm=PathTree --outputPath=../results/json-wikiquote_pl_64.json > ../logs/wikiquote_pl_64.log
```

### 2.5 Parallel computation
To run the program in parallel, please uncomment the following line(`Line 42`) in the `./src/partition.cpp`:
```C++
// #define MTHREAD
```
Then, rebuild the execution file and run the program as mentioned in the previous step.
```bash
./Computation --inputPath=../data/wikiquote_pl.txt --queryFilePath=../query/wikiquote_pl.query --minBlock=64 --numThread 8 --algorithm=PathTree --outputPath=../results/json-wikiquote_pl_64.json > ../logs/wikiquote_pl_64.log
```

### 2.6 $\delta$-constrained reachability
A parameter $\delta$ is used as input, by default, $\delta$ is set to ifinity, meaning that there is no constraint on the time interval. To set the $\delta$ value, please use the following command:
```bash
./Computation --inputPath=../data/wikiquote_pl.txt --queryFilePath=../query/wikiquote_pl.query --minBlock=64 --delta 86400 --algorithm=PathTree --outputPath=../results/json-wikiquote_pl_64.json > ../logs/wikiquote_pl_64.log
```
The unit of $\delta$ is in seconds.