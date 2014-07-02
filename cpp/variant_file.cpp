/*
 * variant_file.cpp
 *
 *  Created on: Dec 11, 2012
 *      Author: amarcketta
 */

#include "variant_file.h"

variant_file::~variant_file() {}

// Return the number of individuals that have not been filtered out
int variant_file::N_kept_individuals() const
{
	int N_kept = 0;
	for (unsigned int ui=0; ui<include_indv.size(); ui++)
		if (include_indv[ui] == true)
			N_kept++;
	return N_kept;
}

// Return the number of sites that have not been filtered out
int variant_file::N_kept_sites() const
{
	return N_kept_entries;
}

// Return the total number of sites in the file
int variant_file::N_total_sites() const
{
	return N_entries;
}

void variant_file::ByteSwap(unsigned char *b, int n) const
{
   register int i = 0;
   register int j = n-1;
   while (i<j)
   {
      std::swap(b[i], b[j]);
      i++, j--;
   }
}

void variant_file::get_contigs(const string &contigs_file, vector<string> &contig_vector)
{
	if (contigs_file == "")
		LOG.error("Contig declarations in header are necessary for BCF conversion. Use --contigs <filename> to add contigs to the header.");

	ifstream contigs(contigs_file.c_str());
	if (!contigs.is_open())
		LOG.error("Could not open contigs file: " + contigs_file);

	string line;
	int contig_lines = 0;
	contig_vector.resize(0);

	while (getline(contigs, line))
	{
		if (line.find("##contig=")==string::npos)
			LOG.error("Contigs file must contain only contig header lines.");

		contig_vector.push_back(line);
		contig_lines++;
	}

	contigs.close();
	LOG.printLOG("Including "+header::int2str(contig_lines)+" header lines from the contig file.\n");
}
