/*
 * Simka2Utils.hpp
 *
 *  Created on: 7 nov. 2016
 *      Author: gbenoit
 */

#ifndef GATB_SIMKA_SRC_SIMKA2_SIMKA2UTILS_HPP_
#define GATB_SIMKA_SRC_SIMKA2_SIMKA2UTILS_HPP_

#include "../core/SimkaDistance.hpp"
//#define SIMKA2_NB_PARTITIONS 200

const string STR_SIMKA2_DATASET_ID = "-id";
const string STR_SIMKA2_DATABASE_DIR = "-database-dir";
const string STR_SIMKA2_NB_PARTITION = "-nb-partitions";
const string STR_SIMKA2_PARTITION_ID = "-partition-id";
const string STR_SIMKA2_INPUT_IDS = "-in-ids";

//const string STR_SIMKA2_DISTANCE_INPUT_1 = "-in-already-computed";
//const string STR_SIMKA2_DISTANCE_INPUT_2 = "-in-to-compute";

string getDatasetID(const string& kmerSpectrumDir){
	string datasetID = System::file().getBaseName(kmerSpectrumDir);
	//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
	return datasetID;
}



u_int64_t simka2_getFileSize(const string& filename){
	std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
	u_int64_t size = in.tellg();
	in.close();
	return size;
}
/*
u_int64_t simka2_getDatasetSize(const string& kmerSpectrumDir){
	u_int64_t datasetSize = 0;
	for(size_t partitionID=0; partitionID<SIMKA2_NB_PARTITIONS; partitionID++){
		string partFilename = kmerSpectrumDir + "/" + Stringify::format("%i", partitionID) + ".gz";
		datasetSize += simka2_getFileSize(partFilename);
	}
	return datasetSize;
}
*/



void simka2_writeDatasetInfo(ofstream& file, const string& datasetID, u_int64_t nbReads, u_int64_t nbDistinctKmers, u_int64_t nbKmers, u_int64_t chord_N2){
    simka2_writeString(datasetID, file);
    file.write((char const*)(&nbReads), sizeof(nbReads));
    file.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
    file.write((char const*)(&nbKmers), sizeof(nbKmers));
    file.write((char const*)(&chord_N2), sizeof(chord_N2));
}

void simka2_readDatasetInfo(ifstream& file, string& datasetID, u_int64_t& nbReads, u_int64_t& nbDistinctKmers, u_int64_t& nbKmers, u_int64_t& chord_N2){
	simka2_readString(datasetID, file);
	file.read((char *)(&nbReads), sizeof(nbReads));
	file.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	file.read((char *)(&nbKmers), sizeof(nbKmers));
	file.read((char *)(&chord_N2), sizeof(chord_N2));
}

void simka2_transferDatasetInfo(ifstream& source, ofstream& dest){
	string datasetID;
	u_int64_t nbReads;
	u_int64_t nbDistinctKmers;
	u_int64_t nbKmers;
	u_int64_t chord_N2;

	simka2_readString(datasetID, source);
	source.read((char *)(&nbReads), sizeof(nbReads));
	source.read((char *)(&nbDistinctKmers), sizeof(nbDistinctKmers));
	source.read((char *)(&nbKmers), sizeof(nbKmers));
	source.read((char *)(&chord_N2), sizeof(chord_N2));

    simka2_writeString(datasetID, dest);
    dest.write((char const*)(&nbReads), sizeof(nbReads));
    dest.write((char const*)(&nbDistinctKmers), sizeof(nbDistinctKmers));
    dest.write((char const*)(&nbKmers), sizeof(nbKmers));
    dest.write((char const*)(&chord_N2), sizeof(chord_N2));
}











void simka2_loadStatInfos(const string& databaseDir, const set<string>& uniqDirs, const vector<string>& ids, vector<string>& kmerSpectrumDirs, SimkaStatistics* stats, map<string, string>& datasetFilenames){
	map<string, vector<u_int64_t> > datasetInfos;
	map<string, string> test;

	//for(size_t i=0; i<_database._uniqKmerSpectrumDirs.size(); i++){
	for (set<string>::iterator i = uniqDirs.begin(); i != uniqDirs.end(); i++) {
		string dir = *i;

		u_int64_t nbMergedBanks = 0;
		string mergeInfoFilename = databaseDir + "/" + dir + "/merge_info.bin";

		ifstream mergedLinkFile(mergeInfoFilename.c_str(), std::ios::binary);
		mergedLinkFile.read((char*)(&nbMergedBanks), sizeof(nbMergedBanks));

		//cout << dir << "  " << nbMergedBanks << endl;

		for(size_t i=0; i<nbMergedBanks; i++){

			string datasetID;
			u_int64_t nbReads;
			u_int64_t nbDistinctKmers;
			u_int64_t nbKmers;
			u_int64_t chord_N2;
			simka2_readDatasetInfo(mergedLinkFile, datasetID, nbReads, nbDistinctKmers, nbKmers, chord_N2);

			vector<u_int64_t> infos;
			infos.push_back(nbReads);
			infos.push_back(nbDistinctKmers);
			infos.push_back(nbKmers);
			infos.push_back(chord_N2);

			datasetInfos[datasetID] = infos;
		}
		mergedLinkFile.close();

		//cout << nbMergedBanks << endl;
	}

	set<string> dirPresent;

	for(size_t i=0; i<ids.size(); i++){
		string id = ids[i];
		vector<u_int64_t> infos = datasetInfos[id];

		stats->_datasetNbReads[i] = infos[0];
		stats->_nbSolidDistinctKmersPerBank[i] = infos[1];
		stats->_nbSolidKmersPerBank[i] = infos[2];


		if(stats->_computeSimpleDistances){
			stats->_chord_sqrt_N2[i] = sqrt(infos[3]);
		}

		string dir = databaseDir + "/" + datasetFilenames[id];
		if(dirPresent.find(dir) != dirPresent.end()) continue;

		kmerSpectrumDirs.push_back(dir);
		dirPresent.insert(dir);

		//cout << kmerSpectrumDirs.size() << endl;
	}


}







template<size_t span=KMER_DEFAULT_SPAN>
class StorageIt
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    u_int64_t _bankIdOffset;
    //typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    //typedef typename ModelCanonical::Kmer                                   KmerType;

    StorageIt(Iterator<Kmer_BankId_Count>* it, size_t bankId, size_t partitionId, u_int64_t bankIdOffset){
    	_it = it;
    	//cout << h5filename << endl;
    	_bankId = bankId;
    	_partitionId = partitionId;
    	_bankIdOffset = bankIdOffset;


		//Iterator<Count>* it2 = partition1.iterator();
		//Collection<Count>& kmers1 = (*partition1)[_partitionId];
		//collections.push_back(&kmers1);

		//_it = kmers1.iterator();

		//_nbKmers = it->estimateNbItems();
		//it2->first();
		//while(!it2->isDone()){
		//	cout << it2->item().value.toString(31) << endl;
		//	it2->next();
		//}
    }

    ~StorageIt(){
    	delete _it;
    }

    //void setPartitionId(size_t partitionId){
    //	_partitionId = partitionId;
    //}

	bool next(){
		_it->next();

		//cout << "is done?" <<  _it->isDone() << endl;
		return !_it->isDone();
	}

	Type& value(){
		return get<0>(_it->item());
	}

	u_int64_t getBankId(){
		return (((u_int64_t)(get<1>(_it->item()))) + _bankIdOffset);
	}

	u_int64_t& abundance(){
		return get<2>(_it->item());
	}



	//u_int64_t getNbKmers(){
	//	return _nbKmers;
	//}

	u_int16_t _bankId;
	u_int16_t _partitionId;
    Iterator<Kmer_BankId_Count>* _it;
    //u_int64_t _nbKmers;
};






template<size_t span>
class DiskBasedMergeSort
{

public:

	typedef typename Kmer<span>::Type                                       Type;
	typedef typename Kmer<span>::Count                                      Count;
    typedef tuple<Type, u_int64_t, u_int64_t> Kmer_BankId_Count;
    typedef tuple<Type, u_int64_t, u_int64_t, StorageIt<span>*> kxp;
	struct kxpcomp { bool operator() (kxp l,kxp r) { return (get<0>(r) < get<0>(l)); } } ;



	//string _outputDir;
	vector<string>& _datasetToMergeDirs;
	size_t _partitionId;
	size_t _nbBanks;
	vector<string> _currentDatasetIds;

    DiskBasedMergeSort(size_t partitionId, vector<string>& datasetToMergeDirs):
    	_datasetToMergeDirs(datasetToMergeDirs)
    {
    	//_outputDir = outputDir;
    	_partitionId = partitionId;

    	_nbBanks = 0;
		//_nbBanks = _datasetIds.size();

    	//cout << _outputFilename << endl;
    }

    ~DiskBasedMergeSort(){
    }

    void execute(){

    	//cout << endl << "start merging" << endl;
		vector<IterableGzFile<Kmer_BankId_Count>* > partitions;
		vector<StorageIt<span>*> its;

		u_int64_t bankIdOffset = 0;

		for(size_t i=0; i<_datasetToMergeDirs.size(); i++){
			string kmerSpectrumDir = _datasetToMergeDirs[i];

			string filename = kmerSpectrumDir + "/" + Stringify::format("%i", _partitionId) + ".gz";

			//cout << filename << endl;

			u_int64_t nbMergedBanks = 0;
			//size_t bankOffset = 0;
	    	string mergeInfoFilename = kmerSpectrumDir + "/merge_info.bin";

			//cout << "\tmerge info:    "<< mergeInfoFilename << "    " << System::file().doesExist(mergeInfoFilename) << endl;

	    	//if(System::file().doesExist(mergeInfoFilename)){
				ifstream mergedLinkFile(mergeInfoFilename.c_str(), std::ios::binary);
				mergedLinkFile.read((char*)(&nbMergedBanks), sizeof(nbMergedBanks));
				//cout << nbMergedBanks << endl;
				for(size_t i=0; i<nbMergedBanks; i++){

					string datasetID;
					u_int64_t nbReads;
					u_int64_t nbDistinctKmers;
					u_int64_t nbKmers;
					u_int64_t chord_N2;
					simka2_readDatasetInfo(mergedLinkFile, datasetID, nbReads, nbDistinctKmers, nbKmers, chord_N2);
					//string datasetId = _allDatasetIds[_datasetIds[i]];
					//u_int64_t size = datasetId.size();
					//mergedLinkFile.write((char const*)(&size), sizeof(size));
					//mergedLinkFile.write(datasetId.c_str(), size);
					//cout << datasetId << endl;
					//string linkedDatasetID;
					//simka2_readString(linkedDatasetID, mergedLinkFile);
					//_currentDatasetIds.push_back(linkedDatasetID);
					//cout << "\tmerge info:" << linkedDatasetID << endl;
				}
				mergedLinkFile.close();

				//_nbBanks += nbMergedBanks;
				//}
				//else{

				//string datasetID = getDatasetID(System::file().getDirectory(filename));
				//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
				//string datasetID = System::file().getBaseName(filename);
				//datasetID.erase(datasetID.end()-string("_kmerSpectrum").size(), datasetID.end());
				//_currentDatasetIds.push_back(datasetID);

				//nbMergedBanks = 1;
				//}

	    	_nbBanks += nbMergedBanks;

			//cout << _datasetIds[i] << endl;
			//string filename = _outputDir + "/solid/part_" +  Stringify::format("%i", _partitionId) + "/__p__" + Stringify::format("%i", _datasetIds[i]) + ".gz";
			//cout << "\t\t" << filename << endl;
			IterableGzFile<Kmer_BankId_Count>* partition = new IterableGzFile<Kmer_BankId_Count>(filename, 10000);
			partitions.push_back(partition);
			//cout << "\tbank offset: " << bankIdOffset << endl;
			its.push_back(new StorageIt<span>(partition->iterator(), i, _partitionId, bankIdOffset));
			bankIdOffset += nbMergedBanks;
			//nbKmers += partition->estimateNbItems();

			//size_t currentPart = 0;
			//ifstream file((_outputDir + "/kmercount_per_partition/" +  _datasetIds[i] + ".txt").c_str());
			//while(getline(file, line)){
			//	if(line == "") continue;
			//	if(currentPart == _partitionId){
			//		//cout << stoull(line) << endl;
			//		nbKmers += strtoull(line.c_str(), NULL, 10);
			//		break;
			//	}
			//	currentPart += 1;
			//}
			//file.close();
		}

		//cout << _nbBanks << endl;
		//u_int64_t progressStep = nbKmers / 1000;
		//_progress = new ProgressSynchro (
		//	createIteratorListener (nbKmers, "Merging kmers"),
		//	System::thread().newSynchronizer());
		//_progress->init ();



		//_nbDistinctKmers = 0;
		//_nbSharedDistinctKmers = 0;
		//u_int64_t nbKmersProcessed = 0;
		//size_t nbBankThatHaveKmer = 0;
		//u_int16_t best_p = 0;
		Type previous_kmer;
		//CountVector abundancePerBank;
		//abundancePerBank.resize(_nbBanks, 0);
		//SimkaCounterBuilderMerge* solidCounter = new SimkaCounterBuilderMerge(abundancePerBank);;
		std::priority_queue< kxp, vector<kxp>,kxpcomp > pq;
		StorageIt<span>* bestIt;


		for(size_t i=0; i<its.size(); i++){
			StorageIt<span>* it = its[i];
				it->_it->first();
			//}
		}

		//fill the  priority queue with the first elems
		for (size_t ii=0; ii<its.size(); ii++)
		{
			//pq.push(Kmer_BankId_Count(ii,its[ii]->value()));
			if (!its[ii]->_it->isDone()){
				pq.push(kxp(its[ii]->value(), its[ii]->getBankId(), its[ii]->abundance(), its[ii]));
			}
		}

		if (pq.size() != 0) // everything empty, no kmer at all
		{
			//get first pointer
			bestIt = get<3>(pq.top()); pq.pop();
	    	process(bestIt->value(), bestIt->getBankId(), bestIt->abundance());

	    	//previous_kmer = bestIt->value();
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
					bestIt = get<3>(pq.top()); pq.pop();
				}

				//if(bestIt->value() == previous_kmer){
				//	process(bestIt->value(), bestIt->getBankId(), bestIt->abundance());
				//}
				//else{
					pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt)); //push new val of this pointer in pq, will be counted later

			    	bestIt = get<3>(pq.top()); pq.pop();
			    	//previous_kmer = bestIt->value();
			    	process(bestIt->value(), bestIt->getBankId(), bestIt->abundance());
			    	//}

		    	//cout << bestIt->value().toString(31) << " " << bestIt->getBankId() <<  " "<< bestIt->abundance() << endl;
		    	//cout << "\t" << bestIt->value().toString(31) << " " << bestIt->getBankId() << endl;
				//bestIt = get<3>(pq.top()); pq.pop();


				//pq.push(kxp(bestIt->value(), bestIt->getBankId(), bestIt->abundance(), bestIt));

			}

			//process(bestIt->value(), bestIt->getBankId(), bestIt->abundance());
	    	//_outputGzFile->insert(Kmer_BankId_Count(bestIt->value(), bestIt->getBankId(), bestIt->abundance()));
	    	//cout << "lol  " << bestIt->value().toString(31) << " " << previous_kmer.toString(31) << endl;
		}

		for(size_t i=0; i<partitions.size(); i++){
			delete partitions[i];
		}

		for(size_t i=0; i<its.size(); i++){
			delete its[i];
		}

		end();
    }

	virtual void process(Type& kmer, u_int64_t bankId, u_int64_t abundance) = 0;

	virtual void end() = 0;

};


#endif /* GATB_SIMKA_SRC_SIMKA2_SIMKA2UTILS_HPP_ */