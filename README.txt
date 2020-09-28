CSCI-654: Project 1

Alberto Serrano-Calva (axs4986)

Sequential program:

    compile:
        gcc gaussian_elimination_sequential.c

    run:
        ./a.out n
            where n is number of dimensions

Parallel program:

    compile:
        gcc gaussian_elimination_parallel.c -lpthread -lm

    run:
        ./a.out d n
            where d is the number of dimensions and n is the number of threads


Additional Materials:

    An additional script, written using Python, is included, which is used to plot
    graph seen in the project report. The program will display the graph during execution
    and save a copy as a PNG file.

    requirements:
        matplotlib, pandas, seaborn

    run:
        python3 plot.py

