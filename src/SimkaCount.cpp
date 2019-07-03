/*****************************************************************************
 *   Simka: Fast kmer-based method for estimating the similarity between numerous metagenomic datasets
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2015  INRIA
 *   Authors: G.Benoit, C.Lemaitre, P.Peterlongo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "SimkaPotara.hpp"
#include "minikc/MiniKC.hpp"
//#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;

//#define NB_COUNT_CACHE 1
//#define TRACK_DISK_USAGE



class SimkaCount : public Tool
{
public:

	SimkaCount () : Tool ("SimkaCount")
    {
        //getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "output file",           true));
        //getParser()->push_back (new OptionOneParam (STR_ID,   "dataset id", true));
        //getParser()->push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size", true));
        getParser()->push_back (new OptionOneParam ("-out-tmp-simka",   "tmp output", true));
        getParser()->push_back (new OptionOneParam ("-bank-name",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-bank-index",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MAX_READS,   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-datasets",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-partitions",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));

        getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);
        if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    }

    void execute ()
    {


    	//size_t datasetId =  getInput()->getInt(STR_ID);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	string bankName =  getInput()->getStr("-bank-name");
    	size_t bankIndex =  getInput()->getInt("-bank-index");
    	size_t minReadSize =  getInput()->getInt(STR_SIMKA_MIN_READ_SIZE);
    	double minReadShannonIndex =  getInput()->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
    	u_int64_t maxReads =  getInput()->getInt(STR_SIMKA_MAX_READS);
    	size_t nbDatasets =   getInput()->getInt("-nb-datasets");
    	size_t nbPartitions =   getInput()->getInt("-nb-partitions");
    	CountNumber abundanceMin =   getInput()->getInt(STR_KMER_ABUNDANCE_MIN);
    	CountNumber abundanceMax =   getInput()->getInt(STR_KMER_ABUNDANCE_MAX);

    	Parameter params(*this, kmerSize, outputDir, bankName, minReadSize, minReadShannonIndex, maxReads, nbDatasets, nbPartitions, abundanceMin, abundanceMax, bankIndex);

        Integer::apply<Functor,Parameter> (kmerSize, params);

    }


    struct Parameter
    {
        Parameter (SimkaCount& tool, size_t kmerSize, string outputDir, string bankName, size_t minReadSize, double minReadShannonIndex, u_int64_t maxReads, size_t nbDatasets, size_t nbPartitions, CountNumber abundanceMin, CountNumber abundanceMax, size_t bankIndex) :
        	tool(tool), kmerSize(kmerSize), outputDir(outputDir), bankName(bankName), minReadSize(minReadSize), minReadShannonIndex(minReadShannonIndex), maxReads(maxReads), nbDatasets(nbDatasets), nbPartitions(nbPartitions), abundanceMin(abundanceMin), abundanceMax(abundanceMax), bankIndex(bankIndex)  {}
        SimkaCount& tool;
        size_t kmerSize;
        string outputDir;
        string bankName;
        size_t minReadSize;
        double minReadShannonIndex;
        u_int64_t maxReads;
        size_t nbDatasets;
        size_t nbPartitions;
        CountNumber abundanceMin;
        CountNumber abundanceMax;
        size_t bankIndex;
    };

    template<size_t span> struct Functor  {

        typedef typename Kmer<span>::Type  Type;
        typedef typename Kmer<span>::Count Count;
        typedef typename SimkaCompressedProcessor<span>::Kmer_BankId_Count Kmer_BankId_Count;

    	void operator ()  (Parameter p){


			IProperties* props = p.tool.getInput();
			vector<string> outInfo;



			IBank* bank = Bank::open(p.outputDir + "/input/" + p.bankName);
			LOCAL(bank);


			vector<u_int64_t> nbKmerPerParts(p.nbPartitions, 0);
			vector<u_int64_t> nbDistinctKmerPerParts(p.nbPartitions, 0);
			vector<u_int64_t> chordNiPerParts(p.nbPartitions, 0);


			Configuration config;
			{
				Repartitor* repartitor = new Repartitor();
				LOCAL(repartitor);

				{
					Storage* storage = StorageFactory(STORAGE_HDF5).load (p.outputDir + "/" + "config.h5");
					LOCAL (storage);
					config.load(storage->getGroup(""));
					repartitor->load(storage->getGroup(""));
				}

				vector<Bag<Kmer_BankId_Count>* > bags;
				vector<Bag<Kmer_BankId_Count>* > cachedBags;
		    	for(size_t i=0; i<p.nbPartitions; i++){
					string outputFilename = p.outputDir + "/solid/part_" + Stringify::format("%i", i) + "/__p__" + Stringify::format("%i", p.bankIndex) + ".gz";
					Bag<Kmer_BankId_Count>* bag = new BagGzFile<Kmer_BankId_Count>(outputFilename);
					Bag<Kmer_BankId_Count>* cachedBag = new BagCache<Kmer_BankId_Count>(bag, 10000);
					cachedBags.push_back(cachedBag);
					//BagCache bagCache(*bag, 10000);
		        	bags.push_back(bag);
		    	}


				string tempDir = p.outputDir + "/temp/" + p.bankName;
				System::file().mkdir(tempDir, -1);

				SimkaSequenceFilter sequenceFilter(p.minReadSize, p.minReadShannonIndex);
				IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, p.maxReads, p.nbDatasets);
				LOCAL(filteredBank);
				//LOCAL(bank);

				SimkaCompressedProcessor<span>* proc = new SimkaCompressedProcessor<span>(cachedBags, nbKmerPerParts, nbDistinctKmerPerParts, chordNiPerParts, p.abundanceMin, p.abundanceMax, p.bankIndex);

				u_int64_t nbReads = 0;

				if(p.kmerSize <= 15){
					MiniKC<span> miniKc(p.tool.getInput(), p.kmerSize, filteredBank, *repartitor, proc);
					miniKc.execute();

					nbReads = miniKc._nbReads;
				}
				else{
					std::vector<ICountProcessor<span>* > procs;
					procs.push_back(proc);
					SortingCountAlgorithm<span> algo (filteredBank, config, repartitor,
							procs,
							props);

					algo.execute();

					nbReads = algo.getInfo()->getInt("seq_number");
				}


				u_int64_t nbDistinctKmers = 0;
				u_int64_t nbKmers = 0;
				u_int64_t chord_N2 = 0;
				for(size_t i=0; i<p.nbPartitions; i++){
					nbDistinctKmers += nbDistinctKmerPerParts[i];
					nbKmers += nbKmerPerParts[i];
					chord_N2 += chordNiPerParts[i];
				}
				outInfo.push_back(Stringify::format("%llu", nbReads));
				outInfo.push_back(Stringify::format("%llu", nbDistinctKmers));
				outInfo.push_back(Stringify::format("%llu", nbKmers));
				outInfo.push_back(Stringify::format("%llu", chord_N2));



#ifdef TRACK_DISK_USAGE
				string command = "du -sh " +  p.outputDir;
				system(command.c_str());
#endif

				System::file().rmdir(tempDir);

		    	for(size_t i=0; i<p.nbPartitions; i++){
		    		delete cachedBags[i];
		    	}

			}

			string contents = "";
			for(size_t i=0; i<nbDistinctKmerPerParts.size(); i++){
				contents += Stringify::format("%llu", nbDistinctKmerPerParts[i]) + "\n";
			}
			IFile* nbKmerPerPartFile = System::file().newFile(p.outputDir + "/kmercount_per_partition/" + p.bankName + ".txt", "w");
			nbKmerPerPartFile->fwrite(contents.c_str(), contents.size(), 1);
			nbKmerPerPartFile->flush();
			delete nbKmerPerPartFile;


			writeFinishSignal(p, outInfo);
		}

		void writeFinishSignal(Parameter& p, const vector<string>& outInfo){

			string finishFilename = p.outputDir + "/count_synchro/" +  p.bankName + ".ok";
			IFile* file = System::file().newFile(finishFilename, "w");
			string contents = "";

			for(size_t i=0; i<outInfo.size(); i++){
				contents += outInfo[i] + "\n";
			}
			file->fwrite(contents.c_str(), contents.size(), 1);
			file->flush();

			delete file;
		}

    };

};

/********************************************************************************/
/*                       Dump solid kmers in ASCII format                       */
/********************************************************************************/
int main (int argc, char* argv[])
{
    try
    {
    	SimkaCount().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}
//! [snippet1]
