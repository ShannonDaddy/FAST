# FAST: Fast Annotation with Synchronized Threads

Dongjae Kim, Aria S. Hahn, Niels W. Hanson, Kishori M. Konwar, and Steven J. Hallam

Comparative genomic analysis relies heavily on sequence alignment, searching reference databases to infer the metabolic potential of individual organisms or entire communities. Rapid advances in next-generation sequencing technologies have generated a torrent of sequence information that typically requires some form of database search for taxonomic or functional annotation. This paper introduces FAST, a multi-threaded, I/O optimized \se algorithm with thread synchronization enabling efficient use of multiple CPU processes while requiring only 4--8 GB of memory. FAST advances database formatting through parallel and incremental construction of Suffix Array-based database indexes, significantly improving sequence alignment workflows requiring frequent database updates. For example, the RefSeq database (9.8GB) can be formatted in under 5 minutes using 40 threads. FAST scales linearly with the number of cores used and can align 1.28 GB of sequence data to the RefSeq database in under 9 hours using 16 threads. FAST is written in C++ 03 as a portable software and is available on GitHub along with tutorials and examples: http://www.github.com/hallamlab/FAST.

More information can be found on the [Wiki](https://github.com/hallamlab/FAST/wiki).

 Test datasets used to produce the claimed results can be found on Dropbox https://www.dropbox.com/sh/xvvavweuzgqybc4/AAAWCrRnol67ZsXi2qLQWUuOa?.
