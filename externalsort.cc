#include "externalsort.hh"

/* Sort the input sequences and divide them into blocks */
void disk_sort_file(const string &outputdir,
                    const string &tobe_sorted_file_name,
                    const std::vector<std::string> &mergelist) {

	// Nothing to merge, rename from tmp dir over to the intended location and exit
	if(mergelist.size() <= 1) {
		std::rename( mergelist.back().c_str(), tobe_sorted_file_name.c_str() );
		return;
	}

	std::size_t num_files = mergelist.size();
	const std::size_t BATCHSIZE = 200;
	if ( num_files > BATCHSIZE ) {
		// There may be so many files that num/200 > 200 so we need to iterate until we
		// have merged them to roughly 200 files
		std::vector<std::string> files;
		while ( num_files > BATCHSIZE ) {
			files = mergeFilesInBatches(mergelist, outputdir);
			num_files = files.size();
		}
		merge_sorted_files(files, tobe_sorted_file_name, outputdir);
	} else {
		merge_sorted_files( mergelist, tobe_sorted_file_name, outputdir );
	}
}

// Generate random string for output in order to allow mutiple LAST+ binaries to run simultaneously on a single machine.
// Check if the directory structure already exists. If it does we need to generate a new randstr
std::vector<std::string> mergeFilesInBatches(const std::vector<std::string> &mergelist,
                                             const string &tmpdir){

	// Recursively merge in batches of 200 until all the remaining files can be fit into one batch
	TempFiles fileptr(tmpdir, generate_directory_name(tmpdir) + "LASTtemp0");
	fileptr.clear();

	std::size_t rounds = mergelist.size() / 200 + 1;
	for(std::size_t i=0; i<=rounds-1; i++){
		std::vector<std::string>::const_iterator it = mergelist.begin() + i*200;
		if(i == rounds-1){
			std::vector<std::string> batch(it, mergelist.end());
			merge_sorted_files(batch, fileptr.nextFileName(), tmpdir);
			removeFiles(batch);
		} else {
			std::vector<std::string> batch(it, mergelist.begin() + (i+1)*200);
			merge_sorted_files(batch, fileptr.nextFileName(), tmpdir);
			removeFiles(batch);
		}
	}
	return fileptr.getFileNames();
}

void merge_sorted_files(const vector<string> &filenames,
                        const string &sorted_file_name,
                        const string &tmpdir) {

	vector<istream_iterator<Line> > f_its;
	vector<ifstream *> ifstream_for_filenames;
	vector<pair<int, Line *> > values;
	Line *curr_lines = new Line[f_its.size()];

	openFileHandlers(filenames, f_its,
	                 ifstream_for_filenames, values,
	                 curr_lines);

	heapSort(sorted_file_name, values, f_its, curr_lines);

	delete[] curr_lines;
	vector<ifstream *>::iterator it = ifstream_for_filenames.begin();
	for(; it != ifstream_for_filenames.end(); ++it){
		delete *it;
	}
}

void heapSort(const string &sorted_file_name,
              vector<pair<int, Line *> > &values,
              vector<istream_iterator<Line> > &f_its,
              Line *curr_lines){

	std::size_t S = f_its.size();
	build_heap(S, values);

	// Open the output file
	ofstream outputfile( (sorted_file_name).c_str() );
	if (!outputfile.is_open()) {
		cout << "Failed to  open output  file " << sorted_file_name << " ... exiting.\n";
		exit(0);
	}

	istream_iterator<Line> empty_it;
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
                      Line *curr_lines)
{
	// Open an istream_iterator for each fasta file
	for (int i = 0; i < filenames.size(); i++) {
		ifstream *f_str = new ifstream(filenames[i].c_str()); // Make sure to keep the file stream "alive" outside this loop
		if ( f_str->peek() != std::ifstream::traits_type::eof() ){
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
			pair<int, Line *> mod_line(i, curr_lines + i);
			values.push_back(mod_line);
		}else {
			exit(0);
		}
	}
}

void removeFiles(const std::vector<std::string> &batch){
	std::vector<std::string>::const_iterator it = batch.begin();
	for( ; it != batch.end(); ++it) {
		if ( remove(it->c_str()) != 0 ) {
			cerr << "Error deleting file " << *it << "\n";
		}
	}
}