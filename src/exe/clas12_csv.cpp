#include "clas12_csv.hpp"
#include <fstream>
#include <future>
#include <iostream>
#include <thread>

int main(int argc, char** argv) {
  // Need this to make sure root doesn't break
  ROOT::EnableThreadSafety();
  bool mc = false;

  int NUM_THREADS = 4;
  if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
  if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

  // Make a vector of vectors of strings the size of the number of threads
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
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

  // Make a set of threads (Futures are special threads which return a value)
  std::future<std::string> threads[NUM_THREADS];

  auto run_files = [](std::vector<std::string> inputs, int thread_id, bool mc) {
    // Called once for each thread
    // Make a new chain to process for this thread
    auto chain = std::make_shared<TChain>("clas12");
    // Add every file to the chain
    for (auto in : inputs) chain->Add(in.c_str());

    // Run the function over each thread
    return run<Cuts>(chain, thread_id, mc);
  };

  // Define events to be used to get Hz later
  std::string events;

  // Start timer
  auto start = std::chrono::high_resolution_clock::now();
  // For each thread
  for (size_t i = 0; i < NUM_THREADS; i++) {
    // Set the thread to run a task A-Syncroisly
    // The function we run is the first argument (run_files)
    // The functions areruments are all the remaining arguments
    threads[i] = std::async(run_files, infilenames.at(i), i, mc);
  }

  // For each thread
  for (size_t i = 0; i < NUM_THREADS; i++) {
    // Get the information from the thread in this case how many events each thread actually computed
    events += threads[i].get();
  }

  std::ofstream outfile;
  outfile.open(outfilename, std::ios::out | std::ios::trunc);
  if (mc == true)
    outfile << "e_rec_p,e_rec_theta,e_rec_phi,e_sec,e_thrown_p,e_thrown_theta,e_thrown_phi" << std::endl;
  else
    outfile << "e_rec_p,e_rec_theta,e_rec_phi,e_sec" << std::endl;
  outfile << events << std::endl;
  outfile.close();

  // Timer and Hz calculator functions that print at the end
  std::cout.imbue(std::locale(""));  // Puts commas in
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  return 0;
}
