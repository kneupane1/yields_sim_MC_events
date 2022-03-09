#include "clas12_yields.hpp"
#include <future>
#include <thread>

// #include <glob.h>    // glob(), globfree()
// #include <string.h>  // memset()
// #include <sstream>
// #include <stdexcept>
// #include <string>
// #include <vector>

// std::vector<std::string> glob(const std::string& pattern) {
//   using namespace std;

//   // glob struct resides on the stack
//   glob_t glob_result;
//   memset(&glob_result, 0, sizeof(glob_result));

//   // do the glob operation
//   int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
//   if (return_value != 0) {
//     globfree(&glob_result);
//     stringstream ss;
//     ss << "glob() failed with return_value " << return_value << endl;
//     throw std::runtime_error(ss.str());
//   }

//   // collect all the filenames into a std::list<std::string>
//   vector<string> filenames;
//   for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
//     filenames.push_back(string(glob_result.gl_pathv[i]));
//   }

//   // cleanup
//   globfree(&glob_result);

//   // done
//   return filenames;
// }

int main(int argc, char** argv) {
        // Need this to make sure root doesn't break
        ROOT::EnableThreadSafety();
        // std::ios::sync_with_stdio(false);

        // Make sure we don't create more threads than files
        int NUM_THREADS = 4;
        if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
        if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

        // Make a vector of vectors of strings the size of the number of threads
        std::vector<std::vector<std::string> > infilenames(NUM_THREADS);
        // Get the output file name
        std::string outfilename;

        if (argc >= 2) {
                // First argument is the output file
                outfilename = argv[1];
                // All other files are split evently by the under of threads
                for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
        } else {
                return 1;
        }

        // Make your histograms object as a shared pointer that all the threads will have
        auto csv_output_file = std::make_shared<SyncFile>(outfilename);
        csv_output_file->write(csv_data::header());
        // auto run_files = [&csv_output_file](auto&& inputs, auto&& thread_id) mutable {
        auto run_files = [&csv_output_file](std::vector<std::string> inputs, auto&& thread_id) mutable {

                                 // Called once for each thread
                                 // Make a new chain to process for this thread
                                 auto chain = std::make_shared<TChain>("clas12");

                                 // Add every file to the chain
                                 for (auto in : inputs) chain->Add(in.c_str());

                                 // Run the function over each thread
                                 // return run(chain, csv_output_file, thread_id);
                                 return run<uconn_Cuts>(std::move(chain), csv_output_file, thread_id);


                         };

        // Make a set of threads (Futures are special threads which return a value)
        std::future<size_t> threads[NUM_THREADS];

        // Define events to be used to get Hz later
        size_t events = 0;

        // Start timer
        auto start = std::chrono::high_resolution_clock::now();
        // For each thread
        for (size_t i = 0; i < NUM_THREADS; i++) {
                // Set the thread to run a task A-Syncroisly
                // The function we run is the first argument (run_files)
                // The functions areruments are all the remaining arguments
                threads[i] = std::async(run_files, infilenames.at(i), i);
        }

        // For each thread
        for (size_t i = 0; i < NUM_THREADS; i++) {
                // Get the information from the thread in this case how many events each thread actually computed
                events += threads[i].get();
        }

        // Timer and Hz calculator functions that print at the end
        std::cout.imbue(std::locale("")); // Puts commas in
        std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
        std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
        std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;

        csv_output_file->writeToFile();

        std::chrono::duration<double> elapsed_full_write = (std::chrono::high_resolution_clock::now() - start);
        std::cout << RED << elapsed_full_write.count() << " sec" << DEF << std::endl;
        std::cout << BOLDYELLOW << events / elapsed_full_write.count() << " Hz" << DEF << std::endl;

        return 0;
}
