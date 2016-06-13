#include "externalsort.hh"

const std::size_t BATCHSIZE = 200;

/* Sort the input sequences and divide them into blocks */
void disk_sort_file(const string &outputdir,
                    const string &outputFile,
                    const std::vector<std::string> &mergelist) {

    // Nothing to merge, rename from tmp dir over to the intended location and exit
    if (mergelist.size() <= 1) {
        std::rename(mergelist.back().c_str(), outputFile.c_str());
        return;
    }

    std::size_t num_files = mergelist.size();
    if (num_files > BATCHSIZE) {
        // It may be that num/BATCHSIZE > BATCHSIZE so we need to iterate until we
        // have merged them to roughly BATCHSIZE files
        std::vector<std::string> files;
        while (num_files > BATCHSIZE) {
            files = mergeFilesInBatches(mergelist, outputdir);
            num_files = files.size();
        }
        merge_sorted_files(files, outputFile, outputdir);
    } else {
        merge_sorted_files(mergelist, outputFile, outputdir);
    }
}

// Generate random string for output in order to allow mutiple LAST+ binaries to run simultaneously on a single machine.
// Check if the directory structure already exists. If it does we need to generate a new randstr
std::vector<std::string> mergeFilesInBatches(const std::vector<std::string> &mergelist,
                                             const string &tmpdir) {

    // Recursively merge in batches of BATCHSIZE until all the remaining files can be fit into one batch
    TempFiles fileptr(tmpdir, generate_directory_name(tmpdir) + "LASTtemp0");

    std::size_t rounds = mergelist.size() / BATCHSIZE + 1;
    for (std::size_t i = 0; i <= rounds - 1; i++) {
        std::vector<std::string>::const_iterator it = mergelist.begin() + i * BATCHSIZE;
        if (i == rounds - 1) {
            std::vector<std::string> batch(it, mergelist.end());
            merge_sorted_files(batch, fileptr.nextFileName(), tmpdir);
        } else {
            std::vector<std::string> batch(it, mergelist.begin() + (i + 1) * BATCHSIZE);
            merge_sorted_files(batch, fileptr.nextFileName(), tmpdir);
        }
    }
    return fileptr.getFileNames();
}

void merge_sorted_files(const vector<string> &filenames,
                        const string &outputFile,
                        const string &tmpdir) {

    Line *curr_lines = new Line[filenames.size()]();
    vector<istream_iterator<Line> > f_its;
    vector<ifstream *> ifstreams;
    vector<pair<int, Line *> > values;

    openFileHandlers(filenames, f_its,
                     ifstreams, values,
                     curr_lines);

    heapSort(outputFile, values,
             f_its, curr_lines);

    delete[] curr_lines;
    vector<ifstream *>::iterator it = ifstreams.begin();
    for (; it != ifstreams.end(); ++it) {
        delete *it;
    }
    removeFiles(filenames);
}

void heapSort(const string &sortedFile,
              vector<pair<int, Line *> > &values,
              vector<istream_iterator<Line> > &f_its,
              Line *curr_lines) {

    build_heap(f_its.size(), values);

    // Open the output file
    ofstream outputfile(sortedFile.c_str());
    if (!outputfile.is_open()) {
        cerr << "Failed to open output file " << sortedFile << " ... exiting";
        exit(0);
    }

    istream_iterator<Line> empty_it;
    std::size_t S = f_its.size();
    while (S > 0) {
        int iter_id = values[0].first;
        // Get the top item off the heap (the longest remaining sequence in any file)
        Line line = *(f_its[iter_id]);
        line.print(outputfile);
        f_its[iter_id]++;

        // Add the next sequence to the top of the heap
        if (f_its[iter_id] != empty_it) {
            curr_lines[iter_id] = *(f_its[iter_id]);
            values[0] = make_pair(iter_id, curr_lines + iter_id);
        } else {
            values[0] = values[S - 1];
            S = S - 1;
        }

        // Re-heapify to make the top value percolate down if it isn't the longest
        if (S > 0) {
            heapify(values, 0, S);
        }
    }
}

// "Constructor" sort of function for readability
void openFileHandlers(const vector<string> &filenames,
                      vector<istream_iterator<Line> > &f_its,
                      vector<ifstream *> &ifstream_for_filenames,
                      vector<pair<int, Line *> > &values,
                      Line *curr_lines) {
    // Open an istream_iterator for each fasta file
    std::vector<std::string>::const_iterator it = filenames.begin();
    for (; it != filenames.end(); ++it) {
        // Make sure to keep the file stream "alive" outside
        // this loop. This is to ensure the iterator will work
        ifstream *f_str = new ifstream(it->c_str());
        if (f_str->peek() != std::ifstream::traits_type::eof()) {
            f_its.push_back(istream_iterator<Line>(*f_str));
            ifstream_for_filenames.push_back(f_str);
        } else {
            delete f_str;
        }
    }

    istream_iterator<Line> empty_it;
    for (int i = 0; i < f_its.size(); i++) {
        if (f_its[i] != empty_it) {
            curr_lines[i] = *(f_its[i]);
            values.push_back(make_pair(i, curr_lines + i));
        } else {
            cerr << "Iterator is empty";
        }
    }
}

void removeFiles(const std::vector<std::string> &batch) {
    std::vector<std::string>::const_iterator it = batch.begin();
    for (; it != batch.end(); ++it) {
        if (remove(it->c_str()) != 0) {
            cerr << "Error deleting file " << *it << "\n";
        }
    }
}