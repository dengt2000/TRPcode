#!/bin/bash

if [ -z $1 ];then
    echo "1.build : generate cmakefiles and executable file"
    echo "2.test : run auto-test"
    echo "    a.without paras, test all datas"
    echo "    b.test with given file and corresponding mode"
    echo "3.clean : clear build files";
    echo "4.cleanAll : cleaning build and test files"
elif [ "$1" = "build" ];then
    # generating cmakefiles and executable file
    if [ ! -d "build" ]; then
        mkdir build
        cd build
        cmake ..
        make
        cd ..
    else 
        cd build
        rm -rf **
        cmake ..
        make
        cd ..
    fi
    if [ ! -d "exec" ]; then
        mkdir exec
    fi
    mv build/Computation exec/

    # test
    chmod +x exec/Computation
    
elif [ "$1" = "test" ];then
    if [ ! -d "result" ]; then
        mkdir result
    fi
    if [ -z $2 ];then
        echo "testing all datas"
        for file in data/rfile/*;do
            echo -n "running with $file;"
            filename=$(basename "$file" | cut -d. -f1)
            if [ ! -d "result/$filename" ]; then
                mkdir result/$filename
            fi
            python3 query_gen.py $file result/$filename/queryinfo.txt 1000
            ./exec/Computation $file r result/$filename/queryinfo.txt &> result/$filename/$filename.log
            echo "finished"
        done
        for file in data/gfile/*;do
            echo -n "running with $file;"
            filename=$(basename "$file" | cut -d. -f1)
            if [ ! -d "result/$filename" ]; then
                mkdir result/$filename
            fi
            python3 query_gen.py $file result/$filename/queryinfo.txt 1000
            ./exec/Computation $file g result/$filename/queryinfo.txt &> result/$filename/$filename.log
            echo "finished"
        done
    else
        if [ -z $3 ];then
            echo "read-mode or gen-mode ?"
            exits
        fi
        if [ -f "$2" ];then
            echo -n "runing with $2;"
            filename=$(basename "$2" | cut -d. -f1)
            if [ ! -d "result/$filename" ]; then
                mkdir result/$filename
            fi
            python3 query_gen.py $2 result/$filename/queryinfo.txt 1000
            ./exec/Computation $2 $3 result/$filename/queryinfo.txt &> result/$filename/$filename.log
            echo "finished"
        else
            echo "file not exits"
        fi
    fi

elif [ "$1" = "clean" ];then
    echo "cleaning build files"
    rm -rf build/
    rm -rf exec/
    
elif [ "$1" = "cleanAll" ];then
    echo "cleaning build and test files"
    rm -rf build/
    rm -rf exec/
    rm -rf result/
fi