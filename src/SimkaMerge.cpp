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


#include <gatb/gatb_core.hpp>
#include <SimkaAlgorithm.hpp>
#include <SimkaDistance.hpp>
#include <fstream>
#include <gzstream.h>
#include "json.hpp"
// We use the required packages
using namespace std;



using namespace gatb::core::system;
using namespace gatb::core::system::impl;
using json = nlohmann::json;

#define MERGE_BUFFER_SIZE 1000
#define SIMKA_MERGE_MAX_FILE_USED 200

struct sortItem_Size_Filename_ID{

	u_int64_t _size;
	size_t _datasetID;

	sortItem_Size_Filename_ID(){}

	sortItem_Size_Filename_ID(u_int64_t size, size_t datasetID){
		_size = size;
		_datasetID = datasetID;
	}
};

bool sortFileBySize (sortItem_Size_Filename_ID i, sortItem_Size_Filename_ID j){
	return ( i._size < j._size );
}

u_int64_t getFileSize(const string& filename){
	std::ifstream in(filename.c_str(), std::ifstream::ate | std::ifstream::binary);
	u_int64_t size = in.tellg();
	in.close();
	return size;
}




struct Parameter
{
    Parameter (IProperties* props, string inputFilename, string outputDir, size_t partitionId, size_t kmerSize, double minShannonIndex, bool computeSimpleDistances, bool computeComplexDistances, size_t nbCores, string f_matrix, string d_matrix, bool is_pipe, string json_path) : props(props), inputFilename(inputFilename), outputDir(outputDir), partitionId(partitionId), kmerSize(kmerSize), minShannonIndex(minShannonIndex), computeSimpleDistances(computeSimpleDistances), computeComplexDistances(computeComplexDistances), nbCores(nbCores), f_matrix(f_matrix), d_matrix(d_matrix), is_pipe(is_pipe), json_path(json_path) {}
    IProperties* props;
    string inputFilename;
    string outputDir;
    size_t partitionId;
    size_t kmerSize;
    double minShannonIndex;
    bool computeSimpleDistances;
    bool computeComplexDistances;
    size_t nbCores;
    string f_matrix;
    string d_matrix;
    bool is_pipe;
    string json_path;
};


template<size_t span=KMER_DEFAULT_SPAN>
class StorageIt
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;

    struct Kmer_BankId_Count{
		Type _type;
		u_int32_t _bankId;
		u_int64_t _count;

		Kmer_BankId_Count(){

		}

		Kmer_BankId_Count(Type type, u_int64_t bankId, u_int64_t count){
			_type = type;
			_bankId = bankId;
			_count = count;
		}
	};

    StorageIt(Iterator<Kmer_BankId_Count>* it, size_t bankId, size_t partitionId);

    ~StorageIt(){
    	delete _it;
    }


	bool next(){
		_it->next();
		return !_it->isDone();
	}

	Type& value(){
		return _it->item()._type;
	}

	u_int16_t getBankId(){
		return _it->item()._bankId;
	}

	u_int64_t& abundance(){
		return _it->item()._count;
	}


	u_int16_t _bankId;
	u_int16_t _partitionId;
    Iterator<Kmer_BankId_Count>* _it;
};

template<size_t span>
StorageIt<span>::StorageIt (Iterator<StorageIt::Kmer_BankId_Count> *it, size_t bankId, size_t partitionId)
{
    _it = it;
    _bankId = bankId;
    _partitionId = partitionId;
}


class SimkaCounterBuilderMerge
{
public:

    /** Constructor.
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
	SimkaCounterBuilderMerge (CountVector& abundancePerBank)  :  _abundancePerBank(abundancePerBank)  {}

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank, CountNumber abundance)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= abundance;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank, CountNumber abundance)  {  _abundancePerBank [idxBank] += abundance;  }

    /** Set the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    //void set (CountNumber val, size_t idxBank=0)  {  _abundancePerBank [idxBank] = val;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    //CountNumber operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** */
    //const CountVector& get () const { return _abundancePerBank; }

    void print(const string& kmer){
		cout << kmer << ": ";
    	for(size_t i=0; i<size(); i++){
    		cout << _abundancePerBank[i] << " ";
    	}
    	cout << endl;
    }

private:
    CountVector& _abundancePerBank;
};


template<size_t span>
class DiskBasedMergeSort
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
	typedef typename StorageIt<span>::Kmer_BankId_Count Kmer_BankId_Count;

	struct kxp{
		Type _type;
		u_int32_t _bankId;
		u_int64_t _count;
		StorageIt<span>* _it;

		kxp(){

		}

		kxp(Type type, u_int64_t bankId, u_int64_t count, StorageIt<span>* it){
			_type = type;
			_bankId = bankId;
			_count = count;
			_it = it;
		}
	};

	struct kxpcomp { bool operator() (kxp& l, kxp& r) { return (r._type < l._type); } } ;

	string _outputDir;
	string _outputFilename;
	vector<size_t>& _datasetIds;
	size_t _partitionId;
	Bag<Kmer_BankId_Count>* _outputGzFile;
	Bag<Kmer_BankId_Count>* _cachedBag;



    DiskBasedMergeSort(size_t mergeId, const string& outputDir, vector<size_t>& datasetIds, size_t partitionId):
    	_datasetIds(datasetIds)
    {
    	_outputDir = outputDir;
    	_partitionId = partitionId;

    	_outputFilename = _outputDir + "/solid/part_" + Stringify::format("%i", partitionId) + "/__p__" + Stringify::format("%i", mergeId) + ".gz.temp";
    	_outputGzFile = new BagGzFile<Kmer_BankId_Count>(_outputFilename);
    	_cachedBag = new BagCache<Kmer_BankId_Count>(_outputGzFile, 10000);

    }

    ~DiskBasedMergeSort(){
    }

    void execute(){

		vector<IterableGzFile<Kmer_BankId_Count>* > partitions;
		vector<StorageIt<span>*> its;

		size_t _nbBanks = _datasetIds.size();

		for(size_t i=0; i<_nbBanks; i++){
			string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);
			partitions.push_back(partition);
			its.push_back(new StorageIt<span>(partition->iterator(), i, _partitionId));
		}

		Type previous_kmer;

		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;
		StorageIt<span>* bestIt;


		for(size_t i=0; i<_nbBanks; i++){
			StorageIt<span>* it = its[i];
			it->_it->first();
		}

		//fill the  priority queue with the first elems
		for (size_t ii=0; ii<_nbBanks; ii++)
		{
			pq.push(kxp(its[ii]->value(), its[ii]->getBankId(), its[ii]->abundance(), its[ii]));
		}

		if (pq.size() != 0) // everything empty, no kmer at all
		{
			//get first pointer
			bestIt = pq.top()._it; pq.pop();
			_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
			//best_p = get<1>(pq.top()) ; pq.pop();
			//previous_kmer = bestIt->value();
			//solidCounter->init (bestIt->getBankId(), bestIt->abundance());
			//nbBankThatHaveKmer = 1;

			while(1){

				if (! bestIt->next())
				{
					//reaches end of one array
					if(pq.size() == 0){
						break;
					}

					//otherwise get new best
					//best_p = get<1>(pq.top()) ; pq.pop();
					bestIt = pq.top()._it; pq.pop();
				}

				pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

		    	bestIt = pq.top()._it; pq.pop();
		    	_cachedBag->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
		    	//cout << bestIt->value().toString(31) << " " << bestIt->getBankId() <<  " "<< bestIt->abundance() << endl;
				//bestIt = get<3>(pq.top()); pq.pop();


				//pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt));

			}
		}

		for(size_t i=0; i<partitions.size(); i++){
			delete partitions[i];
		}

		for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}


		_cachedBag->flush();
    	delete _cachedBag;

		for(size_t i=0; i<_nbBanks; i++){
			string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			System::file().remove(filename);
		}

		string newOutputFilename = _outputFilename;
		newOutputFilename.erase(_outputFilename.size()-5, 5);
    	System::file().rename(_outputFilename, newOutputFilename); //remove .temp at the end of new merged file
    }

};


template<size_t span>
class SimkaMergeAlgorithm : public Algorithm
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
	typedef typename DiskBasedMergeSort<span>::Kmer_BankId_Count Kmer_BankId_Count;
	typedef typename DiskBasedMergeSort<span>::kxp kxp;

    struct kxpcomp { bool operator() (kxp& l,kxp& r) { return (r._type < l._type); } } ;

	Parameter& p;

	SimkaMergeAlgorithm(Parameter& p) :
		Algorithm("SimkaMergeAlgorithm", p.nbCores, p.props), p(p)
	{
		_abundanceThreshold.first = 0;
		_abundanceThreshold.second = 999999999;

		_computeSimpleDistances = p.computeSimpleDistances;
		_computeComplexDistances = p.computeComplexDistances;
		_kmerSize = p.kmerSize;
		_minShannonIndex = p.minShannonIndex;
	}

	~SimkaMergeAlgorithm(){
		//delete _progress;
	}

	void execute(){
	    //Alexandre
        _output_matrix = p.f_matrix;
        _is_pipe = p.is_pipe;
        _output_dir_m = p.d_matrix;
        _nbCores = p.nbCores;

        _json_file = p.json_path;
        bool _groups = (_json_file != "None");
        json _j_groups;
        if (_groups)
        {
            std::ifstream ifs(_json_file);
            ifs >> _j_groups;
        }

		//removeStorage(p);

		_partitionId = p.partitionId;

		ofstream matrix_pipe;
        ogzstream matrix_file;
        char buffer[2048];

        if ( _is_pipe )
        {
            matrix_pipe.rdbuf()->pubsetbuf(buffer, sizeof(buffer));
		    matrix_pipe.open(_output_matrix, std::ios::app);
        }

        else
        {
            const std::string matrix_part = _output_dir_m + "/" + Stringify::format("%i", _partitionId) + ".gz";
            matrix_file.open(matrix_part.c_str());
        }
        //Alexandre
		createDatasetIdList(p);
		_nbBanks = _datasetIds.size();

		string partDir = p.outputDir + "/solid/part_" + Stringify::format("%i", _partitionId) + "/";
		vector<string> filenames = System::file().listdir(partDir);
		vector<string> partFilenames;
		vector<sortItem_Size_Filename_ID> filenameSizes;

		for(size_t i=0; i<filenames.size(); i++){
			if(filenames[i].find("__p__") != std::string::npos){


				string id = string(filenames[i]);
				id.erase(0, 5);
				std::string::size_type pos = id.find(".gz");
				id.erase(pos, 3);

				size_t datasetId = atoll(id.c_str());

				filenameSizes.push_back(sortItem_Size_Filename_ID(getFileSize(partDir+filenames[i]), datasetId));
			}
		}

		while(filenameSizes.size() > SIMKA_MERGE_MAX_FILE_USED){

			sort(filenameSizes.begin(),filenameSizes.end(),sortFileBySize);

			vector<size_t> mergeDatasetIds;
			vector<size_t> toRemoveItem;


			for(size_t i=0; i<SIMKA_MERGE_MAX_FILE_USED; i++){
				sortItem_Size_Filename_ID sfi = filenameSizes[i];
				mergeDatasetIds.push_back(sfi._datasetID);
			}

			for(size_t i=0; i<mergeDatasetIds.size(); i++){
				filenameSizes.erase(filenameSizes.begin());
			}

			size_t mergedId = mergeDatasetIds[0];
			DiskBasedMergeSort<span> diskBasedMergeSort(mergedId, p.outputDir, mergeDatasetIds, _partitionId);
			diskBasedMergeSort.execute();

			filenameSizes.push_back(sortItem_Size_Filename_ID(getFileSize(diskBasedMergeSort._outputFilename), mergedId));
		}

		_stats = new SimkaStatistics(_nbBanks, p.computeSimpleDistances, p.computeComplexDistances, p.outputDir, _datasetIds);

		string line;
		vector<IterableGzFile<Kmer_BankId_Count>* > partitions;
		vector<StorageIt<span>*> its;
		u_int64_t nbKmers = 0;

    	for(size_t i=0; i<filenameSizes.size(); i++){
    		size_t datasetId = filenameSizes[i]._datasetID;
    		string filename = p.outputDir + "/solid/part_" + Stringify::format("%i", p.partitionId) + "/__p__" + Stringify::format("%i", datasetId) + ".gz";
    		IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);

    		partitions.push_back(partition);
    		its.push_back(new StorageIt<span>(partition->iterator(), i, _partitionId));

    		size_t currentPart = 0;
	    	ifstream file((p.outputDir + "/kmercount_per_partition/" +  _datasetIds[i] + ".txt").c_str());
			while(getline(file, line)){
				if(line == "") continue;
				if(currentPart == _partitionId){
					nbKmers += strtoull(line.c_str(), NULL, 10);
					break;
				}
				currentPart += 1;
			}
			file.close();
    	}

		
		_nbDistinctKmers = 0;
		_nbSharedDistinctKmers = 0;
		u_int64_t nbKmersProcessed = 0;
		size_t nbBankThatHaveKmer = 0;
		u_int16_t best_p = 0;
		Type previous_kmer;
	    CountVector abundancePerBank;
		abundancePerBank.resize(_nbBanks, 0);
		SimkaCounterBuilderMerge* solidCounter = new SimkaCounterBuilderMerge(abundancePerBank);;
		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;

    	StorageIt<span>* bestIt;

		for(size_t i=0; i<its.size(); i++){
			StorageIt<span>* it = its[i];
			it->_it->first();
		}

	    //fill the  priority queue with the first elems
	    for (size_t ii=0; ii<its.size(); ii++)
	    {
	    	pq.push(kxp(its[ii]->value(), its[ii]->getBankId(), its[ii]->abundance(), its[ii]));
	    }

	    if (pq.size() != 0) // everything empty, no kmer at all
	    {
	        //get first pointer
	    	bestIt = pq.top()._it; pq.pop();
	        //best_p = get<1>(pq.top()) ; pq.pop();
	        previous_kmer = bestIt->value();
	        solidCounter->init (bestIt->getBankId(), bestIt->abundance());
	        nbBankThatHaveKmer = 1;

	        unsigned int counter = 0;
			while(1){

				if (! bestIt->next())
				{
					//reaches end of one array
					if(pq.size() == 0){
						break;
					}

					//otherwise get new best
					//best_p = get<1>(pq.top()) ; pq.pop();
			    	bestIt = pq.top()._it; pq.pop();
				}


				if (bestIt->value() != previous_kmer )
				{
					//if diff, changes to new array, get new min pointer
					pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

			    	bestIt = pq.top()._it; pq.pop();
					//best_p = get<1>(pq.top()) ; pq.pop();

					//if new best is diff, this is the end of this kmer
					if(bestIt->value()!=previous_kmer )
					{
						insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
						//alexandre
						if ( _is_pipe )
                        {
                            if (_groups) matrix_pipe << toMatrix(previous_kmer, abundancePerBank, _j_groups);
                            else matrix_pipe << toMatrix (previous_kmer, abundancePerBank);
                        }
                        else
                        {
                            if (_groups) matrix_file << toMatrix(previous_kmer, abundancePerBank, _j_groups);
                            else matrix_file << toMatrix (previous_kmer, abundancePerBank);
                        }

						solidCounter->init (bestIt->getBankId(), bestIt->abundance());
						nbBankThatHaveKmer = 1;
						previous_kmer = bestIt->value();
					}
					else
					{
						solidCounter->increase (bestIt->getBankId(), bestIt->abundance());
						nbBankThatHaveKmer += 1;
					}
				}
				else
				{
					solidCounter->increase (bestIt->getBankId(), bestIt->abundance());
					nbBankThatHaveKmer += 1;
				}
			}

			insert(previous_kmer, abundancePerBank, nbBankThatHaveKmer);
			//Alexandre
            if ( _is_pipe )
            {
                if (_groups) matrix_pipe << toMatrix(previous_kmer, abundancePerBank, _j_groups);
                else matrix_pipe << toMatrix (previous_kmer, abundancePerBank);
            }
            else
            {
                if (_groups) matrix_file << toMatrix(previous_kmer, abundancePerBank, _j_groups);
                else matrix_file << toMatrix (previous_kmer, abundancePerBank);
            }
        }


	    matrix_file.close();
		matrix_pipe.close();

	    for ( size_t i = 0; i < partitions.size(); i++ )
		{
			delete partitions[i];
		}

		saveStats(p);

		delete _stats;
		delete solidCounter;
		
        for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}

		writeFinishSignal(p);
	}
	
    void insert(const Type& kmer, const CountVector& counts, size_t nbBankThatHaveKmer)
    {
		//_stats->_nbDistinctKmers += 1;
        if ( nbBankThatHaveKmer > 1 ) { _stats->_nbSharedKmers += 1; }
	}

	std::string toMatrix (const Type& kmer, const CountVector& counts) {
        std::string new_line;

        int sumLine = 0;
        for ( auto& n : counts )
        {
            sumLine += n;
            if ( sumLine > 1) goto keep;
        }
        return new_line;

        keep:
            _stats->_nbDistinctKmers += 1;
            new_line += kmer.toString(_kmerSize);
            new_line += " ";
            for ( auto& i : counts )
            {
                if (i) new_line += "1";
                else { new_line += "0";}
            }
            new_line += "\n";
            return new_line;
    }

    std::string toMatrix (const Type& kmer, const CountVector& counts, const json& groups)
    {
	    std::string new_line(kmer.toString(_kmerSize));
	    new_line += " ";
        bool keep_kmers = false;
	    for (int i=0; i<counts.size(); i++)
        {
	        if (counts[i] == 0) new_line += "0";
	        else if (counts[i] > 1)
            {
                new_line += "1";
                keep_kmers = true;
            }
	        else if (counts[i] == 1)
            {
	            std::cout << "enter" << std::endl;
	            bool in_grp = check_group(counts, groups, i);
	            if (in_grp)
                {
                    keep_kmers = true;
	                new_line += "1";
                }
	            else new_line += "0";
            }
        }

        if (keep_kmers) _stats->_nbDistinctKmers += 1;

	    new_line += "\n";
	    return new_line;
    }

    bool check_group(const CountVector& counts, const json& groups, const int& exp)
    {
	    auto l_groups = groups[std::to_string(exp)];
	    int sum_in_group = 0;
	    for ( auto& pos : l_groups )
        {
	        sum_in_group += counts[pos.get<int>()];
	        if ( sum_in_group > 1 ) return true;
        }
	    return false;
    }
    //Alexandre
    void createDatasetIdList(Parameter& p)
    {

    	string datasetIdFilename = p.outputDir + "/" + "datasetIds";
        IFile* inputFile = System::file().newFile(datasetIdFilename, "rb");

        inputFile->seeko(0, SEEK_END);
		u_int64_t size = inputFile->tell();
		inputFile->seeko(0, SEEK_SET);
		char buffer2[size];
		inputFile->fread(buffer2, size, size);
		string fileContents(buffer2, size);

		string line;
		string linePart;
		vector<string> linePartList;
		stringstream fileContentsStream(fileContents);

		while(getline(fileContentsStream, line)){

			if(line == "") continue;

			_datasetIds.push_back(line);
		}

		delete inputFile;
	}


	//void removeStorage(Parameter& p){
	//	//Storage* storage = 0;
	//	//storage = StorageFactory(STORAGE_HDF5).create (p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".stats", true, true);
	//	//LOCAL (storage);
	//}

	void saveStats(Parameter& p){

		string filename = p.outputDir + "/stats/part_" + SimkaAlgorithm<>::toString(p.partitionId) + ".gz";

		_stats->save(filename); //storage->getGroup(""));

	}

	void writeFinishSignal(Parameter& p){
		string finishFilename = p.outputDir + "/merge_synchro/" +  SimkaAlgorithm<>::toString(p.partitionId) + ".ok";
		IFile* file = System::file().newFile(finishFilename, "w");
		delete file;
	}

private:
	size_t _nbBanks;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	size_t _kmerSize;
	float _minShannonIndex;
    string _output_matrix;
    string _output_dir_m;
    bool _is_pipe;
    string _json_file;
	pair<size_t, size_t> _abundanceThreshold;
	vector<string> _datasetIds;
	size_t _partitionId;

	IteratorListener* _progress;

    //vector<ICommand*> _cmds;
	//ICommand* _mergeCommand;
	size_t _nbCores;


	SimkaStatistics* _stats;
	//SimkaCountProcessorSimple<span>* _processor;
	u_int64_t _nbDistinctKmers;
	u_int64_t _nbSharedDistinctKmers;
};


class SimkaMerge : public Tool
{
public:

	SimkaMerge () : Tool ("SimkaMerge")
    {
		//Original input filename given to simka. Used to recreate dataset id list
        getParser()->push_back (new OptionOneParam (STR_NB_CORES,   "nb cores", true));
        getParser()->push_back (new OptionOneParam (STR_KMER_SIZE,   "kmer size", true));
        getParser()->push_back (new OptionOneParam (STR_URI_INPUT,   "input filename", true));
        getParser()->push_back (new OptionOneParam ("-out-tmp-simka",   "tmp output", true));
        getParser()->push_back (new OptionOneParam ("-partition-id",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-nb-cores",   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-max-memory",   "bank name", true));
        getParser()->push_back (new OptionOneParam (STR_SIMKA_MIN_KMER_SHANNON_INDEX,   "bank name", true));
        getParser()->push_back (new OptionOneParam ("-matrix", "output matrix", true));
        getParser()->push_back (new OptionOneParam ("-dir-matrix", "dir output matrix", false, "./simka_results"));
        getParser()->push_back (new OptionOneParam ("-pipe", "if pipe", false, "false"));
        getParser()->push_back (new OptionOneParam ("-groups", "json file", false, "None"));

        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES.c_str(), "compute simple distances"));
        getParser()->push_back (new OptionNoParam (STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES.c_str(), "compute complex distances"));
    }

    void execute ()
    {


    	size_t nbCores =  getInput()->getInt(STR_NB_CORES);
    	size_t kmerSize =  getInput()->getInt(STR_KMER_SIZE);
    	size_t partitionId =  getInput()->getInt("-partition-id");
    	string inputFilename =  getInput()->getStr(STR_URI_INPUT);
    	string outputDir =  getInput()->getStr("-out-tmp-simka");
    	double minShannonIndex =   getInput()->getDouble(STR_SIMKA_MIN_KMER_SHANNON_INDEX);
    	bool computeSimpleDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_SIMPLE_DISTANCES);
    	bool computeComplexDistances =   getInput()->get(STR_SIMKA_COMPUTE_ALL_COMPLEX_DISTANCES);
        string f_matrix = getInput()->getStr("-matrix");
        string d_matrix = getInput()->getStr("-dir-matrix");
        string pipe = getInput()->getStr("-pipe");
        string json_path = getInput()->getStr("-groups");

        bool is_pipe;
        if (pipe == "true") is_pipe = true;
        else is_pipe = false;

        Parameter params(getInput(), inputFilename, outputDir, partitionId, kmerSize, minShannonIndex, computeSimpleDistances, computeComplexDistances, nbCores, f_matrix, d_matrix, is_pipe, json_path);

        Integer::apply<Functor,Parameter> (kmerSize, params);

    }


    template<size_t span>
    struct Functor  {

    	void operator ()  (Parameter& p)
		{
    		SimkaMergeAlgorithm<span>(p).execute();
		}

    };
};


int main (int argc, char* argv[])
{
    try
    {
    	SimkaMerge().run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}


//! [snippet1]






