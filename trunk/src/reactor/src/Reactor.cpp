#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <getopt.h>
#include "holders.h"

using namespace std;

int main (int argc, char **argv) {
	int c;
	string smile_file = "", smirk_file = "";
	string output_prefix = "result";
	string usage = "Usage: Reactor -m <SMILES_MOLECUES> -r <SMIRKS_REACTIONS> [-o <OUTPUT_PREFIX>] [-h]";

    if (2 > argc) {
        cout << usage << endl;
        return 1;
    }

	while((c =  getopt(argc, argv, ":m:r:o:h")) != -1) {
		switch (c) {
			case 'm':
				smile_file = optarg;
				break;
			case 'r':
				smirk_file = optarg;
				break;
			case 'o':
				output_prefix = optarg;
			 	break;
			case 'h':
				cout << usage << endl;
				return 0;
			default:
				cout << "Unrecognized operation " << c << ". Ignoring" << endl;
				break;
		}
	}
    if (smile_file == "" || smirk_file == "") {
    	cout << usage << endl;
    	return 1;
    }

    ReactionHolder rxn;
    SmilesHolder total_smiles;

    string line;
    stringstream csv, output_list;
    csv << "Smiles,";

    // Parse the smirks file
    ifstream smirks(smirk_file);

    while (getline(smirks, line)) {
    	rxn.addReaction((char*) line.c_str(), line.size());
    	csv << line << ',';
    }
    csv << endl;

    smirks.close();

    // Create the reactions for the smiles
    ifstream smiles(smile_file);

    while (getline(smiles, line)) {
    	auto product_list = rxn.runReactionsOverReactant((char*) line.c_str(), line.size());
    	csv << line << ',';
    	for (auto it = product_list.begin(); it != product_list.end(); ++it) {
    		csv << it->getListAsString('.') << ',';
    		total_smiles.addSmiles(*it);
    	}
    	csv << endl;
    }

    smiles.close();

    // Save the outputs
    ofstream out_csv(output_prefix + "_table.csv");
    out_csv << csv.str();
    out_csv.close();

    cout << "Table saved to file " << output_prefix << "_table.csv" << endl;

    ofstream out_list(output_prefix + "_list.txt");
    out_list << total_smiles.getListAsString('\n');
    out_list.close();

    cout << "Smiles list saved to file " << output_prefix << "_list.txt" << endl;

    ofstream out_rxn(output_prefix + "_reactions.txt");
    out_rxn << rxn.collectAllReactions().getListAsString('\n');
    out_rxn.close();

	cout << "Smiles list saved to file " << output_prefix << "_reactions.txt" << endl;


    return 0;
}

