#include "externalsort.hh"

// Generate random string for output in order to allow mutiple LAST+ binaries to run simultaneously on a single machine.
// Check if the directory structure already exists. If it does we need to generate a new randstr

std::vector<std::string> merge_some_files(const std::vector<std::string> &mergelist,
                                          std::vector<TempFiles> &directories,
                                          const string &tmpdir){

	// Recursively merge in batches of 200 until all the remaining files can be fit into one batch
	std::size_t rounds = mergelist.size() / 200 + 1;

	TempFiles fileptr(tmpdir, generate_directory_name(tmpdir) + "LASTtemp0");
	directories.push_back(fileptr);
	fileptr.clear();

	for(int i=0; i<=rounds-1; i++){
		std::vector<std::string>::const_iterator it = mergelist.begin() + i*200;
		std::string next_name = fileptr.nextFileName();
		if(i == rounds-1){
			std::vector<std::string> batch(it, mergelist.end());
			//print_vector(batch);
			merge_sorted_files(batch, next_name, tmpdir);
		}else{
			std::vector<std::string> batch(it, mergelist.begin() + (i+1)*200);
			//print_vector(batch);
			merge_sorted_files(batch, next_name, tmpdir);
		}
	}
	return fileptr.getFileNames();
}

/* Sort the input sequences and divide them into blocks */
void disk_sort_file(const string &outputdir,
                    const string &tobe_sorted_file_name,
                    const string &sorted_file_name,
                    const std::vector<std::string> &mergelist) {

	std::size_t num_files = mergelist.size();
	std::vector<std::string> files;
	std::vector<TempFiles> directories;

	/*
     std::cout << "NUM FILES : " << mergelist.size() << std::endl;
     std::cout << "NUM ROUNDS : " << mergelist.size()/200+1 << std::endl;
     */

	if(num_files > 200){
		while(num_files > 200){
			files = merge_some_files(mergelist, directories, outputdir);
			num_files = files.size();
		}
		merge_sorted_files( files, tobe_sorted_file_name, outputdir );
	}else{
		merge_sorted_files( mergelist, tobe_sorted_file_name, outputdir );
	}

	for(int i=0; i<directories.size(); i++){
		directories[i].clear();
	}
}

void heapSort(const string &sorted_file_name,
              vector<pair<int, Line *> > &values,
              vector<istream_iterator<Line> > &f_its,
              Line *curr_lines){

	int S = f_its.size();
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
			pair<int, Line *> mod_line;
			mod_line.first = iter_id;
			mod_line.second = curr_lines + iter_id;

			values[0] = mod_line;
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
