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

/*
TODO:
	- Faire la config dans un job a part (job_count.bash) pour avoir la même config pour les job de comptage et la config
	- Verifier les paramètre passer au jobs généré (nbcores, maxmemory...)

*/
SimkaPotara::SimkaPotara(const string& execFilename)  : Tool ("Simka")
{

	_execFilename = execFilename;

	Simka::createOptionsParser(getParser());

	//Kmer parser
    IOptionsParser* coreParser = getParser()->getParser("core");

    coreParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_COUNT, "maximum number of simultaneous counting jobs (a higher value improve execution time but increase temporary disk usage)", false));
    coreParser->push_back (new OptionOneParam (STR_SIMKA_NB_JOB_MERGE, "maximum number of simultaneous merging jobs (1 job = 1 core)", false));


    IOptionsParser* clusterParser = new OptionsParser ("cluster");
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_COMMAND, "command to submit counting job", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_COMMAND, "command to submit merging job", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_COUNT_FILENAME, "filename to the couting job template", false ));
    clusterParser->push_back (new OptionOneParam (STR_SIMKA_JOB_MERGE_FILENAME, "filename to the merging job template", false ));


	//getParser()->push_back(coreParser);
	getParser()->push_back(clusterParser);
}

struct Parameter
{
    Parameter (IProperties* props, const string& execFilename) : _props(props), _execFilename(execFilename) {}
    IProperties* _props;
    string _execFilename;
};

template<size_t span> struct Functor  {  void operator ()  (Parameter p)
{
	SimkaPotaraAlgorithm<span> simkaAlgorithm (p._props, p._execFilename);
	simkaAlgorithm.execute();
}};

void SimkaPotara::execute ()
{
	IProperties* input = getInput();
	Parameter params(input, _execFilename);

	size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    Integer::apply<Functor,Parameter> (kmerSize, params);
}


int main (int argc, char* argv[])
{
    try
    {
        // We run the tool with the provided command line arguments.

        SimkaPotara(string(argv[0])).run (argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
