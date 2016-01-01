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
//#include <gatb/gatb_core.hpp>

// We use the required packages
using namespace std;













template<typename Filter> class SimkaPotaraBankFiltered : public BankDelegate
{
public:

	SimkaPotaraBankFiltered (IBank* ref, const Filter& filter, u_int64_t maxReads, size_t nbDatasets) : BankDelegate (ref), _filter(filter)  {
		//_nbReadsPerDataset = nbReadsPerDataset;
		_maxReads = maxReads;
		_nbDatasets = nbDatasets;
	}


    Iterator<Sequence>* iterator ()
    {

        Iterator<Sequence>* it = _ref->iterator ();
        //std::vector<Iterator<Sequence>*> iterators = it->getComposition();
        return new SimkaInputIterator<Sequence, Filter> (it, _nbDatasets, _maxReads, _filter);
    	//return filterIt;

    }

private:

	//vector<u_int64_t> _nbReadsPerDataset;
    u_int64_t _maxReads;
    Filter _filter;
    u_int64_t _nbReadToProcess;
    size_t _datasetId;
    size_t _nbDatasets;
};


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
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SIZE,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_READ_SHANNON_INDEX,   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MAX_READS,   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-datasets",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        //getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));

        getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);
        if (Option* p = dynamic_cast<Option*> (getParser()->getParser(STR_KMER_ABUNDANCE_MIN)))  {  p->setDefaultValue ("0"); }
    }

    void execute ()
    {


    	//size_t datasetId =  getInput()->getInt(STR_ID);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	//cout << kmerSize << endl;

    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	string bankName =  getInput()->getStr("-bank-name");
    	size_t minReadSize =  getInput()->getInt(STR_SIMKA_MIN_READ_SIZE);
    	double minReadShannonIndex =  getInput()->getDouble(STR_SIMKA_MIN_READ_SHANNON_INDEX);
    	u_int64_t maxReads =  getInput()->getInt(STR_SIMKA_MAX_READS);
    	size_t nbDatasets =   getInput()->getInt("-nb-datasets");
    	size_t nbPartitions =   getInput()->getInt("-nb-partitions");

    	Parameter params(*this, kmerSize, outputDir, bankName, minReadSize, minReadShannonIndex, maxReads, nbDatasets, nbPartitions);

        Integer::apply<Functor,Parameter> (kmerSize, params);



		//SimkaBankId* bank = new SimkaBankId(_banks, i);
		//cout << config._nb_partitions << endl;
		//KmerCountCompressor<span>* kmerCountCompressor = new KmerCountCompressor<span>(outputDir, config._nb_partitions, 1);

		//SimkaCompProcessor<span>* processor = new SimkaCompProcessor<span>(kmerCountCompressor);
		//vector<ICountProcessor<span>*> procs;
		//procs.push_back(processor);

		//algo.addProcessor(processor);

		//algo.execute();

		//delete kmerCountCompressor;
		//itBanks[i]->


        // We get a handle on the HDF5 storage object.
        // Note that we use an auto pointer since the StorageFactory dynamically allocates an instance
        //Storage* storage = StorageFactory(DSK::getStorageMode()).load (getInput()->getStr(STR_URI_FILE));
        //LOCAL (storage);

        //string kmerSizeStr = storage->getGroup("params").getProperty ("kmer_size");

        //if (kmerSizeStr.empty())  { throw Exception ("unable to get the kmer size"); }

        //size_t kmerSize = atoi (kmerSizeStr.c_str());

    }


    struct Parameter
    {
        Parameter (SimkaCount& tool, size_t kmerSize, string outputDir, string bankName, size_t minReadSize, double minReadShannonIndex, u_int64_t maxReads, size_t nbDatasets, size_t nbPartitions) :
        	tool(tool), kmerSize(kmerSize), outputDir(outputDir), bankName(bankName), minReadSize(minReadSize), minReadShannonIndex(minReadShannonIndex), maxReads(maxReads), nbDatasets(nbDatasets), nbPartitions(nbPartitions)  {}
        SimkaCount& tool;
        //size_t datasetId;
        size_t kmerSize;
        string outputDir;
        string bankName;
        size_t minReadSize;
        double minReadShannonIndex;
        u_int64_t maxReads;
        size_t nbDatasets;
        size_t nbPartitions;
    };

    template<size_t span> struct Functor  {

        typedef typename Kmer<span>::Type  Type;

    	void operator ()  (Parameter p){


			IProperties* props = p.tool.getInput();





			IBank* bank = Bank::open(p.outputDir + "/input/" + p.bankName);
			LOCAL(bank);
			u_int64_t nbSeqs = 1;
	        IBank* sampleBank = new SimkaBankSample(bank, nbSeqs);
			SortingCountAlgorithm<span> sortingCount (sampleBank, props);
			SimkaNullProcessor<span>* proc = new SimkaNullProcessor<span>();
			sortingCount.addProcessor (proc);
			sortingCount.execute();
			Configuration config = sortingCount.getConfig();
			//_nbPartitions = _maxJobMerge;
			config._nb_partitions = p.nbPartitions;

			uint64_t memoryUsageCachedItems;
			config._nb_cached_items_per_core_per_part = 1 << 8; // cache at least 256 items (128 here, then * 2 in the next while loop)
			do
			{
				config._nb_cached_items_per_core_per_part *= 2;
				memoryUsageCachedItems = 1LL * config._nb_cached_items_per_core_per_part *config._nb_partitions * config._nbCores * sizeof(Type);
			}
			while (memoryUsageCachedItems < config._max_memory * MBYTE / 10);








			//Configuration* config = new Configuration();
			{
				Repartitor* repartitor = new Repartitor();
				LOCAL(repartitor);

				{
					Storage* storage = StorageFactory(STORAGE_HDF5).load (p.outputDir + "/" + "config.h5");
					LOCAL (storage);
					//config->load(storage->getGroup(""));
					repartitor->load(storage->getGroup(""));
				}

				//config._abundanceUserNb = 1;
				//config._abundance.clear();
				//CountRange range(props->getInt(STR_KMER_ABUNDANCE_MIN), 100000);
				//config._abundance.push_back(range);

				string tempDir = p.outputDir + "/temp/" + p.bankName;
				System::file().mkdir(tempDir, -1);
				//cout << i << endl;
				//string outputDir = p.outputDir + "/comp_part" + to_string(p.datasetId) + "/";

				//cout << "\tinput: " << p.outputDir + "/input/" + p.bankName << endl;

				SimkaSequenceFilter sequenceFilter(p.minReadSize, p.minReadShannonIndex);
				IBank* filteredBank = new SimkaPotaraBankFiltered<SimkaSequenceFilter>(bank, sequenceFilter, p.maxReads, p.nbDatasets);
				// = new SimkaPotaraBankFiltered(bank)
				LOCAL(filteredBank);
				//LOCAL(bank);

				Storage* solidStorage = 0;

				string solidsName = p.outputDir + "/solid/" +  p.bankName + ".h5";
				bool autoDelete = false; // (solidsName == "none") || (solidsName == "null");
				solidStorage = StorageFactory(STORAGE_HDF5).create (solidsName, true, autoDelete);
				LOCAL(solidStorage);

				solidStorage->root().setCompressLevel (1);
				//props->add(1, STR_HISTOGRAM_MAX, "0");
				//props->add(1, STR_KMER_ABUNDANCE_MIN_THRESHOLD, "0");
				//props->add(1, STR_SOLIDITY_KIND, "sum");
				//props->add(1, STR_URI_OUTPUT_TMP, tempDir);
				//cout << tempDir << endl;

				SortingCountAlgorithm<span> algo (filteredBank, config, repartitor,
						SortingCountAlgorithm<span>::getDefaultProcessorVector (config, props, solidStorage),
						props);

				algo.execute();

				System::file().rmdir(tempDir);
			}

			//cout << "heo" << endl;
			//delete config;
			//cout << "heo" << endl;
			writeFinishSignal(p);
			//cout << "heo" << endl;
		}

		void writeFinishSignal(Parameter& p){

			string finishFilename = p.outputDir + "/count_synchro/" +  p.bankName + ".ok";
			IFile* file = System::file().newFile(finishFilename, "w");
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
