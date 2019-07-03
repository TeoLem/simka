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

#ifndef TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_
#define TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_

#include <gatb/gatb_core.hpp>
#include "SimkaCommons.hpp"
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include<stdio.h>
#include <iostream>
#include <fstream>
//#define PRINT_STATS
//#define CHI2_TEST
//#define SIMKA_POTARA
//#define BOOTSTRAP
#define MAX_BOOTSTRAP 50
#define NB_BOOTSTRAP 45
//#define SIMKA_FUSION
//#define MULTI_PROCESSUS
//#define MULTI_DISK
//#define SIMKA_MIN
#include "SimkaDistance.hpp"



enum SIMKA_SOLID_KIND{
	RANGE,
	SUM,
};



typedef u_int16_t bankIdType;





class SimkaCounterBuilder
{
public:

    /** Constructor.
     * \param[in] nbBanks : number of banks parsed during kmer counting.
     */
    SimkaCounterBuilder (size_t nbBanks=1)  :  _abundancePerBank(nbBanks)  {}

    /** Get the number of banks.
     * \return the number of banks. */
    size_t size() const  { return _abundancePerBank.size(); }

    /** Initialization of the counting for the current kmer. This method should be called
     * when a kmer is seen for the first time.
     * \param[in] idxBank : bank index where the new current kmer has been found. */
    void init (size_t idxBank=0)
    {
        for (size_t k=0; k<_abundancePerBank.size(); k++)  { _abundancePerBank[k]=0; }
        _abundancePerBank [idxBank]= 1;
    }

    /** Increase the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void increase (size_t idxBank=0)  {  _abundancePerBank [idxBank] ++;  }

    /** Set the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank */
    void set (CountNumber val, size_t idxBank=0)  {  _abundancePerBank [idxBank] = val;  }

    /** Get the abundance of the current kmer for the provided bank index.
     * \param[in] idxBank : index of the bank
     * \return the abundance of the current kmer for the given bank. */
    CountNumber operator[] (size_t idxBank) const  { return _abundancePerBank[idxBank]; }

    /** */
    const CountVector& get () const { return _abundancePerBank; }

private:
    CountVector _abundancePerBank;
};


template<typename Filter> class SimkaBankFiltered : public BankDelegate
{
public:

	u_int64_t _refNbReads;
	u_int64_t _refTotalSeqSize;
	u_int64_t _refMaxReadSize;

	/** Constructor.
     * \param[in] ref : referred bank.
     * \param[in] filter : functor that filters sequence.
     */
	SimkaBankFiltered (IBank* ref, const Filter& filter, const vector<size_t>& nbPaireds, u_int64_t maxReads) : BankDelegate (ref), _filter(filter)  {

		_nbPaireds = nbPaireds;
		_maxReads = maxReads;
		_nbBanks = ref->getCompositionNb();
		ref->estimate(_refNbReads, _refTotalSeqSize, _refMaxReadSize);
	}


    void estimate (u_int64_t& number, u_int64_t& totalSize, u_int64_t& maxSize){


    	if(_maxReads == 0){
    		number = _refNbReads;
    		totalSize = _refTotalSeqSize;
    		maxSize = _refMaxReadSize;
    	}
    	else{

    		u_int64_t maxReads = 0;
    		for(size_t i=0; i<_nbBanks; i++){
    			maxReads += _maxReads * _nbPaireds[i];
    		}
    		maxReads = min (maxReads, _refNbReads);

			if(maxReads == _refNbReads){
	    		number = _refNbReads;
	    		totalSize = _refTotalSeqSize;
	    		maxSize = _refMaxReadSize;
			}
			else{
				number = maxReads;
				double factor =  (double)maxReads / (double)_refNbReads;
				totalSize = _refTotalSeqSize * factor;
				maxSize = _refMaxReadSize;
			}
    	}

    	cout << number << endl;
    }

    /** \copydoc tools::collections::Iterable::iterator */
    Iterator<Sequence>* iterator ()
    {

        // We create one iterator from the reference
        Iterator<Sequence>* it = _ref->iterator ();

        // We get the composition for this iterator
        std::vector<Iterator<Sequence>*> iterators = it->getComposition();

        //if (iterators.size() == 1)  { return new FilterIterator<Sequence,Filter> (it, _filter); }
        //else
        //{
            // We are going to create a new CompositeIterator, we won't need the one we just got from the reference
		LOCAL(it);

		// We may have to encapsulate each sub iterator with the filter.
		for (size_t i=0; i<iterators.size(); i++)  {
			iterators[i] = new SimkaInputIterator<Sequence, Filter> (iterators[i], _nbPaireds[i], _maxReads, _filter);
            }

		return new CompositeIterator<Sequence> (iterators);
    }

private:

	vector<size_t> _nbPaireds;
    Filter _filter;
    u_int64_t _maxReads;
    size_t _nbBanks;
};


/*********************************************************************
* ** SimkaAlgorithm
*********************************************************************/
template<size_t span=KMER_DEFAULT_SPAN>
class SimkaAlgorithm : public Algorithm
{

public:


    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::Count                                      Count;
    typedef typename Kmer<span>::ModelCanonical                             ModelCanonical;
    typedef typename ModelCanonical::Kmer                                   KmerType;

	SimkaAlgorithm(IProperties* options);
	~SimkaAlgorithm();
	void execute();
	void print();

    static string toString(u_int64_t value){
    	char buffer[40];
    	snprintf(buffer, 30, "%llu", value);
    	return string(buffer);
    }

protected:


    bool setup();
    bool isInputValid();
    void parseArgs();
    bool createDirs();
    void computeMaxReads();
	void layoutInputFilename();
	void createBank();
	void count();

	void outputMatrix();

	void clear();
    
	u_int64_t _maxMemory;
	size_t _nbCores;
	string _outputDir;
	string _outputDirTemp;
	size_t _nbBanks;
	string _inputFilename;
	size_t _kmerSize;
	pair<CountNumber, CountNumber> _abundanceThreshold;
	SIMKA_SOLID_KIND _solidKind;
	bool _soliditySingle;
	int64_t _maxNbReads;
	size_t _minReadSize;
	double _minReadShannonIndex;
	double _minKmerShannonIndex;
	size_t _nbMinimizers;

    std::string _output_m;
    bool _pipe;
    std::string _json_path;

	SimkaStatistics* _stats;
	//SimkaDistance* _simkaDistance;

	string _banksInputFilename;
	vector<string> _tempFilenamesToDelete;
	IBank* _banks;
	IProperties* _options;

	vector<string> _bankNames;
	//vector<u_int64_t> _nbReadsPerDataset;

	string _outputFilenameSuffix;

	u_int64_t _totalKmers;
    vector<size_t> _nbBankPerDataset;

	string _largerBankId;
	bool _computeSimpleDistances;
	bool _computeComplexDistances;
	bool _keepTmpFiles;

};


#endif /* TOOLS_SIMKA_SRC_SIMKAALGORITHM_HPP_ */
